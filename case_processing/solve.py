"""便携合成输入的新轨迹求解；连续反例反馈，不使用历史热启动。"""
from pathlib import Path
import json,time,copy,argparse,warnings
from model import Task,ROOT,sha
from service_candidate_solver import solve_service_candidate
warnings.filterwarnings('ignore',message='Unrecognized options detected.*')
def value(c):return c['on']+.001*c['delivered']
def run(tag,domain,N,B,budget=30,out=None):
 task=Task(tag,domain);start=time.perf_counter();deadline=start+budget;best=copy.deepcopy(task.direct);best.update({'N':N,'B':B})
 trace=[{'seconds':0.,'on':best['on'],'delivered':best['delivered']}];calls=[];cycle=0;stop='TIME_LIMIT'
 skeletons=[(1,[1]),(1,[2])]+([(2,[1,1]),(2,[2,1]),(2,[2,2])] if N==2 else [])
 while time.perf_counter()<deadline:
  if best['on']==task.upper and best['delivered']==task.delivery_upper:stop='TRANSPORT_NECESSARY_UPPER_REACHED';break
  advanced=False
  for used,counts in skeletons:
   left=deadline-time.perf_counter()
   if left<=.01:break
   if not len(task.pool):stop='NO_SUPPORT_POOL';break
   rows,info=solve_service_candidate(task.p,used,B,1,counts,task.pool,task.queries,task.cover,task.weights,seconds=min(5.,left),margin=.02)
   entry={'cycle':cycle,'seconds':time.perf_counter()-start,'N_used':used,'counts':counts,'solver':info};calls.append(entry)
   if rows is None:continue
   before=len(task.feedback);state=[(r['point_index'],r['link_done'],r['end']) for r in rows];cert=task.certify(state,N,B)
   entry['certification']={'status':cert['status'],'on':cert['on'],'delivered':cert['delivered'],'feedback_added':len(task.feedback)-before};entry['end_seconds']=time.perf_counter()-start
   if cert['status']=='CERTIFIED' and value(cert)>value(best) and time.perf_counter()<=deadline:
    best=cert;trace.append({'seconds':time.perf_counter()-start,'on':cert['on'],'delivered':cert['delivered']});advanced=True
   if len(task.feedback)>before:advanced=True
   if best['on']==task.upper and best['delivered']==task.delivery_upper:break
  cycle+=1
  if not advanced:stop='NO_NEW_CERTIFICATE_OR_COUNTEREXAMPLE';break
 obj={'status':'COMPLETE_CERTIFIED_SEARCH','tag':tag,'domain':domain,'N':N,'B':B,'full_charge_seconds':task.p.full,'denominator':len(task.ledger),'on':best['on'],'delivered':best['delivered'],'transport_necessary_on_upper':task.upper,'transport_necessary_delivery_upper':task.delivery_upper,'certificate':best,'budget_seconds':budget,'preprocessing_seconds':task.preprocess_seconds,'search_seconds':time.perf_counter()-start,'stop':stop,'trace':trace,'calls':calls,'counterexamples':task.feedback,'input_hashes':task.bindings,'code_hashes':{str(f.relative_to(ROOT)):sha(f) for f in ROOT.glob('*.py')},'scope':'事件公开点位和地形上的预设服务；既有MILP，每轮48点/每体最多2槽，连续验收；上界仅适用于给定任务剖面。'}
 path=Path(out) if out else ROOT/'results'/tag/domain/('N'+str(N)+'_B'+str(B)+'.json');path.parent.mkdir(parents=True,exist_ok=True)
 path.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf-8');print(tag,domain,N,B,best['on'],best['delivered'],len(task.ledger),stop,flush=True);return obj
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--tag',required=True);p.add_argument('--domain',choices=['baseline','input_aware_48','full_mother'],required=True);p.add_argument('--N',type=int,required=True);p.add_argument('--B',type=int,required=True);p.add_argument('--budget',type=float,default=30);p.add_argument('--out')
 a=p.parse_args();run(a.tag,a.domain,a.N,a.B,a.budget,a.out)

