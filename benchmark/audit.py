"""便携独立验收：重算箱截止、运输能量、中继物理、组件占用与连续通信。"""
from pathlib import Path
from fractions import Fraction as F
import json,math,argparse,time
import numpy as np
from model import Problem,ROOT,sha
from frontier_engine import coverage_sortie
import independent_radio_geometry as G
ff=lambda x:F.from_float(float(x))
def audit(path):
 started=time.perf_counter();path=Path(path);obj=json.loads(path.read_text(encoding='utf-8'))
 for file,h in obj['input_hashes'].items():assert sha(ROOT/file)==h
 for file,h in obj['code_hashes'].items():assert sha(ROOT/file)==h
 p=Problem(ROOT/'inputs'/obj['tag'],obj['domain']);cert=obj['certificate'];ledger=p.task['账本'];requests={b['箱ID']:b for b in ledger}
 assert len(requests)==len(ledger)==obj['denominator']
 served=set(cert['served_K']);chosen=[r for r in p.task['transport'] if r['K'] in served];assert {r['K'] for r in chosen}==served
 ids=[];on=[];c=p.contract['货运参数'];cargo=[]
 for r in chosen:
  phases=[ph for ph in p.phases if ph['K']==r['K']];assert len(phases)==7
  assert phases[0]['t0']==r['start']+c['prep'] and phases[-1]['t1']==r['end']
  boxes=[requests[key] for key in r['deliver']];payload=len(boxes)*c['payload_per_box']
  h1=phases[0]['P1'][2]-phases[0]['P0'][2];h2=phases[4]['P1'][2]-phases[4]['P0'][2];travel=phases[1]['t1']-phases[1]['t0'];d=G.distance(phases[1]['P0'],phases[1]['P1'],p.radius)
  assert abs(travel-d/c['v_cruise'])<1e-8
  e=2*c['p_cruise']*travel/3600+c['p_hover']*c['dwell']/3600+p.contract['重力加速度']*((c['mtow']+payload)*h1+c['mtow']*h2)/(c['eta_up']*p.contract['每kWh焦耳数'])
  assert abs(e-r['energy'])<1e-8 and e<=c['e_use']*(1-c['rho'])+1e-10
  for ph in phases:assert G.ray_closed_clear(ph['P0'],ph['P1'],p.dem)
  ratio=ff(r['energy'])/F(str(c['e_use']));charge=F(str(c['full_charge']))*F(7,2)*ratio if ratio<=F(1,10) else F(str(c['full_charge']))*(F(5,18)+F(13,18)*ratio)
  cargo.append((r['drone'],r['battery'],ff(r['start']),ff(r['end']),ff(r['end'])+charge))
  for box,t in r['deliver'].items():
   b=requests[box];assert b['计划架次']==r['sid'] and abs(t-b['计划送达秒'])<1e-8
   assert b['硬截止秒'] is None or ff(t)<=ff(b['硬截止秒'])
   ids.append(box)
   if ff(t)<=ff(b['期望送达秒']):on.append(box)
 for i,(body,module,A,R,C) in enumerate(cargo):
  assert A>=0
  for body2,module2,A2,R2,C2 in cargo[i+1:]:
   assert body!=body2 or R<=A2 or R2<=A
   assert module!=module2 or C<=A2 or C2<=A
 assert len(ids)==len(set(ids)) and sorted(ids)==cert['delivered_ids'] and sorted(on)==cert['on_ids']
 assert len(ids)==obj['delivered'] and len(on)==obj['on']
 bodies={};modules={};physics=[]
 for r in cert['rows']:
  i=r['point_index'];point=p.palette[i];assert np.max(abs(point-np.array([r['lon'],r['lat'],r['hover_alt']])))<1e-8
  base=np.array([p.node['lon'],p.node['lat'],p.node['elev']]);d=G.distance(base,[point[0],point[1],base[2]],p.radius);f=np.linspace(0,1,max(2,int(math.ceil(d/25))+1));xy=base[:2]+f[:,None]*(point[:2]-base[:2])
  cc=np.rint((xy[:,0]-p.dem.lon[0])/p.dem.dx).astype(int);rr=np.rint((xy[:,1]-p.dem.lat[0])/p.dem.dy).astype(int)
  cruise=max(float(p.dem.z[rr,cc].max())+50,base[2]+50,point[2]);h1=cruise-base[2];h2=cruise-point[2];travel=d/p.rt['v_cruise']
  out=h1/p.rt['v_up']+travel+h2/p.rt['v_dn'];back=h2/p.rt['v_up']+travel+h1/p.rt['v_dn'];fly=2*p.rt['p_cruise']*travel/3600+p.rt['mtow']*p.contract['重力加速度']*(h1+h2)/(p.rt['eta_up']*p.contract['每kWh焦耳数'])
  assert max(abs(out-p.out[i]),abs(back-p.back[i]),abs(fly-p.fly[i]))<1e-8
  S=ff(r['link_done']);E=ff(r['end']);A=S-ff(p.out[i])-F(str(p.link))-F(str(p.prep));R=E+ff(p.back[i]);energy=ff(p.fly[i])+F(str(p.power))*(E-S+F(str(p.link)))/3600
  assert A>=0 and E>S and R-A-F(str(p.prep))<=F(str(p.limit)) and energy<=F(str(p.rt['e_use']))*(1-F(str(p.rt['rho'])))
  used=energy/F(str(p.rt['e_use']));T=F(str(p.full));charge=T*F(7,2)*used if used<=F(1,10) else T*(F(5,18)+F(13,18)*used)
  record={'A':A,'depart':A+F(str(p.prep)),'S':S,'E':E,'R':R,'C':R+charge}
  bodies.setdefault(r['drone'],[]).append(record);modules.setdefault(r['comp'],[]).append(record)
  assert G.certified_point_link(p.gw,point,p.dem,p.cp,p.budgets['backhaul'],p.radius) is True
  physics.append({'site':i,'out':out,'back':back,'energy':fly})
 assert len(bodies)<=obj['N'] and len(modules)<=obj['B']
 windows=[]
 for chain in bodies.values():
  chain.sort(key=lambda r:r['S']);windows +=[(chain[0]['depart'],chain[0]['S']),(chain[-1]['E'],chain[-1]['R'])]
  for a,b in zip(chain,chain[1:]):assert b['A']>=a['R']+F(str(p.turn));windows.append((a['E'],b['S']))
 for chain in modules.values():
  chain.sort(key=lambda r:r['A'])
  for a,b in zip(chain,chain[1:]):assert b['A']>=a['C']
 events={}
 for a,b in windows:events[a]=events.get(a,0)+1;events[b]=events.get(b,0)-1
 peak=used=0
 for t,change in sorted(events.items()):used+=change;peak=max(peak,used)
 assert peak<=1
 checks={K:coverage_sortie(p,cert['rows'],K) for K in served};assert all(c['状态']=='PASS' for c in checks.values())
 receipt={'status':'PASS','result_sha256':sha(path),'on':len(on),'delivered':len(ids),'denominator':len(ledger),'separate_charge_formula':True,'separate_cargo_energy':True,'separate_relay_physics':physics,'continuous_coverage':checks,'blackout_peak':peak,'seconds':time.perf_counter()-started,'shared_radio_kernel':True,'scope':'重新读输入、另一计分/能量/充电与库存实现；数值传播内核共享，不代表实测或机器区间证明。'}
 path.with_name(path.stem+'_独立验收.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2),encoding='utf-8')
 print(obj['tag'],obj['domain'],obj['N'],obj['B'],'独立验收PASS',flush=True);return receipt
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('file');audit(parser.parse_args().file)

