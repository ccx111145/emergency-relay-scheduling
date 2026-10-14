"""只读便携合成输入与完整请求账本；不导入旧工作区。"""
from pathlib import Path
from fractions import Fraction as F
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import json,copy,time,math,hashlib
import numpy as np
import independent_radio_geometry as G
from general_solver_core import choose_pool,verify
from frontier_engine import resource_options,coverage_sortie
from channel_support_vector_v7 import support_vector
ROOT=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def make_rows(p,state):
 return [{'sid':'中继架次'+str(j),'point_index':int(i),'link_done':float(S),'end':float(E),'lon':float(p.palette[i,0]),'lat':float(p.palette[i,1]),'hover_alt':float(p.palette[i,2]),'prepare':float(S-p.out[i]-p.link-p.prep),'depart':float(S-p.out[i]-p.link),'back':float(E+p.back[i])} for j,(i,S,E) in enumerate(state)]
class Problem:
 def __init__(self,folder,domain):
  self.folder=Path(folder);self.contract=json.loads((self.folder/'contract.json').read_text(encoding='utf-8'));self.task=json.loads((self.folder/'task.json').read_text(encoding='utf-8'))
  self.tag=self.contract['tag'];self.rt=self.contract['中继参数'];self.cp=self.contract['无线参数'];self.node=self.contract['基地'];self.radius=self.contract['地球半径'];self.full=self.contract['充满时长'];self.limit=self.contract['滞空上限']
  with np.load(self.folder/'dem.npz') as z:self.dem=G.Raster(z['z'],z['lon'],z['lat'],float(z['nodata']))
  with np.load(self.folder/(domain+'.npz')) as z:self.palette=z['palette'];self.agl=z['agl']
  flights=np.array([G.relay_flight(v,float(h),self.dem,self.node,self.rt,self.radius,self.contract['重力加速度'],self.contract['每kWh焦耳数'])[:3] for v,h in zip(self.palette,self.agl)])
  self.out,self.back,self.fly=flights.T;self.power=self.rt['p_hover']+self.rt['p_comm'];self.prep=self.rt['t_prep'];self.link=self.rt['t_link'];self.turn=self.rt['t_turn'];self.emax=self.rt['e_use']*(1-self.rt['rho'])
  self.cap=np.minimum(self.limit-self.out-self.back-self.link,(self.emax-self.fly)*3600/self.power-self.link)
  self.gw=np.array([self.node['lon'],self.node['lat'],self.node['elev']+self.cp['gw_ant_h']])
  self.budgets={name:G.budget_pair(self.cp[a],self.cp[b],self.cp) for name,a,b in [('direct','transport','gateway'),('access','transport','relay_access'),('backhaul','relay_backhaul','gateway')]}
  self.phases=self.task['phases'];self.bh_cache={};self.geometry_cache={}
 def position(self,ph,t):return np.array(ph['P0'])+(np.array(ph['P1'])-ph['P0'])*((t-ph['t0'])/(ph['t1']-ph['t0']))
 def link_ok(self,a,b,budget):return G.point_link(a,b,self.dem,self.cp,budget,self.radius)
class Task:
 def __init__(self,tag,domain):
  start=time.perf_counter();self.p=Problem(ROOT/'inputs'/tag,domain);self.tag=tag;self.domain=domain;p=self.p
  self.ledger=p.task['账本'];self.sidK={ph['sid']:ph['K'] for ph in p.phases}
  self.eligible={b['计划架次'] for b in self.ledger if b['计划架次'] is not None and b['计划硬截止满足']}
  self.ks=sorted({self.sidK[s] for s in self.eligible});self.sid_for_K={self.sidK[s]:s for s in self.eligible}
  self.counts={K:(sum(b['计划按期'] for b in self.ledger if b['计划架次']==self.sid_for_K[K]),sum(b['计划架次']==self.sid_for_K[K] for b in self.ledger)) for K in self.ks}
  self.upper=sum(a for a,b in self.counts.values());self.delivery_upper=sum(b for a,b in self.counts.values());self.weights={K:a+.001*b for K,(a,b) in self.counts.items()}
  self.queries=[];self.feedback=[];self.bindings={str(f.relative_to(ROOT)):sha(f) for f in p.folder.iterdir() if f.is_file()}
  for j,ph in enumerate(p.phases):
   if ph['K'] not in self.ks:continue
   for t in np.linspace(ph['t0'],ph['t1'],17):
    P=p.position(ph,float(t))
    if G.certified_point_link(p.gw,P,p.dem,p.cp,p.budgets['direct'],p.radius) is False:self.queries.append({'K':ph['K'],'sid':ph['sid'],'t':float(t),'P':P.tolist(),'phase_index':j})
  self.update_cover();self.direct=self.certify([],0,0,feedback=False);assert self.direct['status']=='CERTIFIED'
  self.preprocess_seconds=time.perf_counter()-start
 def update_cover(self):
  p=self.p
  self.cover=np.column_stack([support_vector(p,q,p.geometry_cache) for q in self.queries]) if self.queries else np.zeros((len(p.palette),0),bool)
  self.cover[p.cap<=1,:]=False
  self.pool=choose_pool(p,self.queries,self.cover,48) if self.queries else np.flatnonzero(p.cap>1.)[:48]
 def certify(self,state,N,B,feedback=True):
  p=self.p;options=[r for r in resource_options(p,make_rows(p,state),p.full,1) if r['N']<=N and r['B']<=B]
  if not options:return {'status':'FAIL_RESOURCE','on':0,'delivered':0}
  rows=options[0]['rows'];old=p.phases;p.phases=[]
  resources=verify(p,rows,N,B,1);p.phases=old
  if resources['状态']!='PASS':return {'status':resources['状态'],'on':0,'delivered':0}
  checks={K:coverage_sortie(p,rows,K) for K in self.ks};served={K for K,c in checks.items() if c['状态']=='PASS'}
  if feedback:
   existing={(q['K'],q['t']) for q in self.queries};added=[]
   for K,c in checks.items():
    if c['状态']!='FAIL':continue
    w=c['反例'];ph=old[w['phase_index']];t=w['t']
    if (K,t) in existing:continue
    P=p.position(ph,t)
    if G.certified_point_link(p.gw,P,p.dem,p.cp,p.budgets['direct'],p.radius) is not False:continue
    existing.add((K,t));added.append({'K':K,'sid':ph['sid'],'t':t,'P':P.tolist(),'phase_index':w['phase_index']})
   if added:self.queries.extend(added);self.feedback.extend(added);self.update_cover()
  boxes=[b for b in self.ledger if b['计划架次'] in self.eligible and self.sidK[b['计划架次']] in served]
  p.phases=[ph for ph in old if ph['K'] in served];final=verify(p,rows,N,B,1);p.phases=old
  if final['状态']!='PASS':return {'status':final['状态'],'on':0,'delivered':0}
  return {'status':'CERTIFIED','N':N,'B':B,'on':sum(b['计划按期'] for b in boxes),'delivered':len(boxes),'denominator':len(self.ledger),'on_ids':sorted(b['箱ID'] for b in boxes if b['计划按期']),'delivered_ids':sorted(b['箱ID'] for b in boxes),'served_K':sorted(served),'rows':rows,'checks':checks,'final_check':final}

