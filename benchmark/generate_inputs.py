"""生成有发布权的合成地形和请求；不读取受限第三方输入，不筛选有利随机种子。"""
from pathlib import Path
import json,hashlib
import numpy as np
import independent_radio_geometry as G
ROOT=Path(__file__).resolve().parent
REGIMES=['flat_reference','transport_energy','deadline','ridge_communication','long_duty','battery_readiness','battery_outage']
SEEDS=[2026101301,2026101302,2026101303]
RADIUS=6371000.;DEG=RADIUS*np.pi/180
def save(p,o):p.write_text(json.dumps(o,ensure_ascii=False,indent=2),encoding='utf-8')
def charge(e,E,T):
 r=e/E
 return 3.5*T*r if r<=.1 else T*(5/18+13/18*r)
def build(regime,seed):
 rng=np.random.default_rng(seed);folder=ROOT/'inputs'/(regime+'_'+str(seed));folder.mkdir(parents=True,exist_ok=True)
 x=np.linspace(-2000,14000,65);y=np.linspace(-6000,6000,49);X,Y=np.meshgrid(x,y)
 ridge=regime in ['ridge_communication','long_duty','battery_readiness','battery_outage']
 Z=80+8*np.sin(X/1700)*np.sin(Y/2100)+(600*np.exp(-((X-4000)/550)**2) if ridge else 0)
 dem=G.Raster(Z,x/DEG,y/DEG,-9999.)
 ground=lambda xx,yy:float(dem.z[int(round((yy-y[0])/(y[1]-y[0]))),int(round((xx-x[0])/(x[1]-x[0])))])
 node={'lon':0.,'lat':0.,'elev':ground(0,0)}
 rt={'gid':'示例中继','mtow':23.5,'v_cruise':15.,'p_cruise':1.15,'e_use':.9 if regime in ['battery_readiness','battery_outage'] else 3.2,'rho':.2,'t_prep':180.,'t_link':30.,'t_turn':300.,'v_up':4.,'v_dn':3.,'eta_up':.72,'p_hover':1.05,'p_comm':.05,'h_max_agl':300.}
 cp={'f_mhz':2400.,'l_sys':3.,'l_obs':10.,'p_sens':-98.,'margin':8.,'transport':[20.,3.,3.],'relay_access':[20.,6.,6.],'relay_backhaul':[19.,8.,8.],'gateway':[27.,12.,12.],'gw_ant_h':20.}
 contract={'tag':folder.name,'基地':node,'中继参数':rt,'无线参数':cp,'地球半径':RADIUS,'重力加速度':9.80665,'每kWh焦耳数':3600000.,'充满时长':43200. if regime=='battery_outage' else (5400. if regime=='battery_readiness' else 1800.),'滞空上限':7233.,'模型来源':'全部参数是明确选择的合成示例假设，未经设备校准，不代表实测数据。','请求来源':'固定随机种子程序生成，非真实居民需求','地形来源':'解析函数生成的合成高程栅格，非真实地区DSM','货运参数':{'e_use':.35 if regime=='transport_energy' else 4.,'rho':.2,'p_cruise':1.4,'p_hover':1.1,'mtow':18.,'payload_per_box':1.5,'v_cruise':20.,'v_up':3.,'v_down':2.,'eta_up':.72,'prep':90.,'dwell':90.,'full_charge':2400.}}
 palettes={};heights={}
 for name,xs,ys in [('baseline',np.array([0.,2000.,10000.,12000.]),np.linspace(-3000,3000,4)),('full_mother',np.linspace(0,12000,13),np.linspace(-3000,3000,9))]:
  coords=[];agl=[]
  for h in [50.,150.,300.]:
   for xx in xs:
    for yy in ys:
     coords.append([xx/DEG,yy/DEG,ground(xx,yy)+h]);agl.append(h)
  palettes[name]=np.array(coords);heights[name]=np.array(agl)
 n=8 if regime=='long_duty' else 6
 interval=1600. if regime=='long_duty' else (2500. if regime in ['battery_readiness','battery_outage'] else 1800.)
 phases=[];ledger=[];transports=[];exclusions=[]
 c=contract['货运参数']
 for j in range(n):
  xx=([4000.,8000.,12000.][j%3] if regime=='transport_energy' else 8000.)+rng.uniform(-350,350)
  yy=rng.uniform(-700,700);dest=np.array([xx/DEG,yy/DEG,ground(xx,yy)+60.]);base=np.array([0.,0.,node['elev']+20.])
  d=G.distance(base,[dest[0],dest[1],base[2]],RADIUS);nseg=max(2,int(np.ceil(d/25))+1)
  f=np.linspace(0,1,nseg);xy=base[:2]+f[:,None]*(dest[:2]-base[:2])
  cols=np.rint((xy[:,0]-dem.lon[0])/dem.dx).astype(int);rows=np.rint((xy[:,1]-dem.lat[0])/dem.dy).astype(int)
  cruise=max(float(dem.z[rows,cols].max())+120,base[2]+120,dest[2]+120)
  up=(cruise-base[2])/c['v_up'];down=(cruise-dest[2])/c['v_down'];travel=d/c['v_cruise']
  back_up=(cruise-dest[2])/c['v_up'];back_down=(cruise-base[2])/c['v_down']
  loading=2*c['payload_per_box'];energy=2*c['p_cruise']*travel/3600+c['p_hover']*c['dwell']/3600+contract['重力加速度']*((c['mtow']+loading)*(cruise-base[2])+c['mtow']*(cruise-dest[2]))/(c['eta_up']*contract['每kWh焦耳数'])
  start=900.+j*interval;depart=start+c['prep'];arrival=depart+up+travel+down;end=arrival+c['dwell']+back_up+travel+back_down
  expected=arrival+(120 if j%2==0 else -60)
  hard=arrival-60 if regime=='deadline' and j%2==1 else arrival+600
  usable=energy<=c['e_use']*(1-c['rho'])
  sid='航班'+str(j+1);K='运输'+str(j+1)
  if usable:
   P0=base;P1=np.array([base[0],base[1],cruise]);P2=np.array([dest[0],dest[1],cruise]);P3=dest
   segments=[(depart,depart+up,P0,P1),(depart+up,depart+up+travel,P1,P2),(depart+up+travel,arrival,P2,P3),(arrival,arrival+c['dwell'],P3,P3),(arrival+c['dwell'],arrival+c['dwell']+back_up,P3,P2),(arrival+c['dwell']+back_up,arrival+c['dwell']+back_up+travel,P2,P1),(arrival+c['dwell']+back_up+travel,end,P1,P0)]
   for a,b,A,B in segments:
    phases.append({'K':K,'sid':sid,'t0':float(a),'t1':float(b),'P0':A.tolist(),'P1':B.tolist()})
   transports.append({'sid':sid,'K':K,'type':'示例货运','drone':'货运机'+str(j+1),'battery':'货运电池'+str(j+1),'start':start,'end':end,'energy':energy,'chg':charge(energy,c['e_use'],c['full_charge']),'deliver':{sid+'-箱1':arrival,sid+'-箱2':arrival}})
  else:exclusions.append({'sid':sid,'reason':'固定同基地直飞能量超预算','required_kwh':energy,'available_kwh':c['e_use']*(1-c['rho'])})
  for b in [1,2]:
   ledger.append({'箱ID':sid+'-箱'+str(b),'请求坐标':dest.tolist(),'期望送达秒':expected,'硬截止秒':hard,'计划架次':sid if usable else None,'计划送达秒':arrival if usable else None,'计划按期':bool(arrival<=expected) if usable else False,'计划硬截止满足':bool(arrival<=hard) if usable else False,'排除原因':None if usable else '固定同基地能量超限'})
 save(folder/'contract.json',contract);save(folder/'task.json',{'账本':ledger,'transport':transports,'phases':phases,'energy_exclusions':exclusions,'seed':seed,'regime_design':regime,'regime_label_status':'预设机制方向，不代表实验已证实该瓶颈'})
 np.savez_compressed(folder/'dem.npz',z=Z,lon=x/DEG,lat=y/DEG,nodata=-9999.)
 for name in palettes:np.savez_compressed(folder/(name+'.npz'),palette=palettes[name],agl=heights[name])
 return {'tag':folder.name,'regime_design':regime,'seed':seed,'denominator':len(ledger),'energy_excluded_boxes':2*len(exclusions),'hard_deadline_excluded_boxes':sum(b['计划架次'] is not None and not b['计划硬截止满足'] for b in ledger),'source':'完全合成，自行生成；全部固定种子保留'}
def main():
 rows=[build(r,s) for r in REGIMES for s in SEEDS]
 files={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'inputs').rglob('*') if p.is_file()}
 save(ROOT/'输入清单.json',{'tasks':rows,'hashes':files,'scope':'21个固定种子的合成基准任务，独立分母；不是历史救援数据。'})
 print('生成',len(rows),'个完整合成任务，不代表实验已经完成',flush=True)
if __name__=='__main__':main()

