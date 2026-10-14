"""公开合成基准和案例的运行入口；不读取私人输入或预存解。"""
from pathlib import Path
import argparse,json,hashlib,sys,subprocess,shutil,csv
ROOT=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify():
 m=json.loads((ROOT/'benchmark_manifest.json').read_text(encoding='utf8'));assert len(m['tasks'])==21
 for key,value in m['input_sha256'].items():assert sha(ROOT/key)==value,key
 assert sum(x['denominator'] for x in m['tasks'])==264
 print('21任务、264请求与全部输入SHA校验通过。')
def run_checked(runtime,tag,domain,N,B,budget,out):
 command=[sys.executable,'-B',str(runtime/'solve.py'),'--tag',tag,'--domain',domain,'--N',str(N),'--B',str(B),'--budget',str(budget),'--out',str(out)]
 if domain=='input_aware_48':
  # 既有函数支持此域；保留原命令行的两个域选项，不改模型。
  command=[sys.executable,'-B','-c',"import sys;sys.path.insert(0,sys.argv[1]);from solve import run;run(sys.argv[2],sys.argv[3],int(sys.argv[4]),int(sys.argv[5]),float(sys.argv[6]),sys.argv[7])",str(runtime),tag,domain,str(N),str(B),str(budget),str(out)]
 subprocess.run(command,check=True);subprocess.run([sys.executable,'-B',str(runtime/'audit.py'),str(out)],check=True)
def benchmark(budget):
 verify();runtime=ROOT/'benchmark';subprocess.run([sys.executable,'-B',str(runtime/'prepare_domains.py'),'--check'],check=True);tags=[x['tag'] for x in json.loads((ROOT/'benchmark_manifest.json').read_text())['tasks']];records=[]
 for tag in tags:
  for domain in ['baseline','full_mother']:
   for N,B in [(1,2),(1,3),(2,2),(2,3)]:
    out=runtime/'results/finite_pool'/tag/domain/f'N{N}_B{B}.json';out.parent.mkdir(parents=True,exist_ok=True);run_checked(runtime,tag,domain,N,B,budget,out);records.append(json.loads(out.read_text()))
  for N in [1,2]:
   out=runtime/'results/equal48'/tag/'input_aware_48'/f'N{N}_B3.json';out.parent.mkdir(parents=True,exist_ok=True);run_checked(runtime,tag,'input_aware_48',N,3,budget,out);records.append(json.loads(out.read_text()))
  ready=tag+'__ready';dest=runtime/'inputs'/ready
  if not dest.exists():shutil.copytree(runtime/'inputs'/tag,dest)
  p=dest/'contract.json';j=json.loads(p.read_text());j['tag']=ready;j['充满时长']=0.;j['就绪假设']='瞬时就绪控制，最多四架次分别提供全充模块；不是实测充电性能';p.write_text(json.dumps(j,ensure_ascii=False,indent=2),encoding='utf8')
  for domain in ['baseline','full_mother']:
   out=runtime/'results/instant_ready'/ready/domain/'N2_B4.json';out.parent.mkdir(parents=True,exist_ok=True);run_checked(runtime,ready,domain,2,4,budget,out);records.append(json.loads(out.read_text()))
 keys=['tag','domain','N','B','denominator','on','delivered','transport_necessary_on_upper','transport_necessary_delivery_upper','search_seconds','stop'];report=ROOT/'reports';report.mkdir(exist_ok=True)
 with (report/'benchmark_results.csv').open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows({k:v[k] for k in keys} for v in records)
 assert len(records)==252;print('252配置完成并逐方案独立验收；结果为当前机器和预算的可行下界。')
def cases(data_root,budget):
 assert data_root is not None,'案例需指定合法取得的公开数据目录'
 subprocess.run([sys.executable,'-B',str(ROOT/'case_processing/process_public_data.py'),'--data-root',str(data_root)],check=True)
 runtime=ROOT/'case_processing';subprocess.run([sys.executable,'-B',str(runtime/'generate_cases.py')],check=True);subprocess.run([sys.executable,'-B',str(runtime/'prepare_domains.py')],check=True)
 for tag in ['dingri_communication','rasuwa_dhunche','rasuwa_betrawati']:
  for domain in ['baseline','input_aware_48','full_mother']:
   for N,B in [(1,2),(1,3),(2,2),(2,3)]:
    out=runtime/'results'/tag/domain/f'N{N}_B{B}.json';out.parent.mkdir(parents=True,exist_ok=True);run_checked(runtime,tag,domain,N,B,budget,out)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--mode',choices=['verify','benchmark','cases','all'],default='verify');p.add_argument('--budget',type=float,default=60);p.add_argument('--data-root',type=Path);a=p.parse_args()
 if a.mode=='verify':verify()
 if a.mode in ['benchmark','all']:benchmark(a.budget)
 if a.mode in ['cases','all']:cases(a.data_root,a.budget)
