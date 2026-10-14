"""由必要通信查询确定性选择48点，不读取已有解或服务结果。"""
from pathlib import Path
import argparse,json,numpy as np
from model import Task
ROOT=Path(__file__).resolve().parent
def main(check):
 tags=[r['tag'] for r in json.loads((ROOT.parent/'benchmark_manifest.json').read_text())['tasks']]
 for tag in tags:
  q=Task(tag,'full_mother');ids=list(map(int,q.pool));ids+=[int(i) for i in np.flatnonzero(q.p.cap>1) if int(i) not in ids][:48-len(ids)];ids+=[int(i) for i in range(len(q.p.palette)) if int(i) not in ids][:48-len(ids)];ids=ids[:48];assert len(ids)==48
  path=ROOT/'inputs'/tag/'input_aware_48.npz'
  if check:
   with np.load(path) as z:assert np.array_equal(z['palette'],q.p.palette[ids]) and np.array_equal(z['agl'],q.p.agl[ids]),tag
  else:np.savez_compressed(path,palette=q.p.palette[ids],agl=q.p.agl[ids])
 print('21任务的等基数48点规则与输入一致。' if check else '21任务的48点输入重建完成。')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');main(p.parse_args().check)
