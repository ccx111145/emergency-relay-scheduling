from pathlib import Path
import json,math,shutil,sys
import numpy as np,pandas as pd
from pyproj import Transformer
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT));import independent_radio_geometry as G
O=ROOT.parent;D=O/'source_data';RADIUS=6371000.
def save(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
def charge(e,E,T):
 r=e/E;return T*(3.5*r if r<=.1 else 5/18+13/18*r)
def build(tag,points,node,terrain,mode):
 f=ROOT/'inputs'/tag;f.mkdir(parents=True,exist_ok=True);shutil.copy2(D/'地形裁切'/(terrain+'_dem.npz'),f/'dem.npz');z=np.load(f/'dem.npz');dem=G.Raster(z['z'],z['lon'],z['lat'],float(z['nodata']))
 def ground(lon,lat):return float(dem.z[int(round((lat-dem.lat[0])/dem.dy)),int(round((lon-dem.lon[0])/dem.dx))])
 node=dict(node);node['elev']=ground(node['lon'],node['lat'])
 rt={'gid':'R','mtow':23.5,'v_cruise':15.,'p_cruise':1.15,'e_use':3.2,'rho':.2,'t_prep':180.,'t_link':30.,'t_turn':300.,'v_up':4.,'v_dn':3.,'eta_up':.72,'p_hover':1.05,'p_comm':.05,'h_max_agl':300.}
 cp={'f_mhz':2400.,'l_sys':3.,'l_obs':10.,'p_sens':-98.,'margin':8.,'transport':[20.,3.,3.],'relay_access':[20.,6.,6.],'relay_backhaul':[19.,8.,8.],'gateway':[27.,12.,12.],'gw_ant_h':20.}
 c={'e_use':4.,'rho':.2,'p_hover':8.003571428571428,'mtow':65.,'payload_per_box':1.5,'v_cruise':15.,'v_up':3.,'v_down':2.5,'eta_up':.72,'prep':300.,'dwell':90.,'full_charge':2400.}
 # 沿用参考Type B的航程能量标定，载荷1.5kg：在空载与满载功率间线性插值。
 c['p_cruise']=4*15*3600/28000+(4*15*3600/16000-4*15*3600/28000)*(1.5/30)
 contract={'tag':tag,'基地':node,'中继参数':rt,'无线参数':cp,'地球半径':RADIUS,'重力加速度':9.80665,'每kWh焦耳数':3600000.,'充满时长':1800.,'滞空上限':7233.,'货运参数':c,'案例模式':mode,'地形来源':'Copernicus GLO-30公开DSM，真实裁切','模型解释':'中继和无线参数沿用主稿参考；货物量/服务窗口/交付期限为协议设定，不是观测救援任务；未校准高海拔性能。'}
 lo=min([x['lon'] for x in points]+[node['lon']])-.015;hi=max([x['lon'] for x in points]+[node['lon']])+.015;bottom=min([x['lat'] for x in points]+[node['lat']])-.015;top=max([x['lat'] for x in points]+[node['lat']])+.015
 xs=np.linspace(lo,hi,13);ys=np.linspace(bottom,top,9)
 for domain,ix,iy in [('baseline',[0,4,8,12],[0,3,5,8]),('full_mother',range(13),range(9))]:
  pal=[];agl=[]
  for h in [50.,150.,300.]:
   for i in ix:
    for j in iy:pal.append([float(xs[i]),float(ys[j]),ground(xs[i],ys[j])+h]);agl.append(h)
  np.savez_compressed(f/(domain+'.npz'),palette=np.array(pal),agl=np.array(agl))
 phases=[];ledger=[];transports=[];excluded=[];next_ready=[0.]*4
 for j,q in enumerate(points):
  sid='request'+str(j+1);K='service'+str(j+1);dest=np.array([q['lon'],q['lat'],ground(q['lon'],q['lat'])+(2 if mode=='communication' else 60)]);base=np.array([node['lon'],node['lat'],node['elev']+20.]);usable=True
  if mode=='communication':
   start=900.+j*120;arrival=start;end=start+300;energy=0.;segments=[(start,end,dest,dest)];body='ground-terminal-'+str(j);battery=body
  else:
   dist=G.distance(base,[dest[0],dest[1],base[2]],RADIUS);parts=max(2,int(np.ceil(dist/25))+1);xy=base[:2]+np.linspace(0,1,parts)[:,None]*(dest[:2]-base[:2]);cc=np.rint((xy[:,0]-dem.lon[0])/dem.dx).astype(int);rr=np.rint((xy[:,1]-dem.lat[0])/dem.dy).astype(int);cruise=max(float(dem.z[rr,cc].max())+120,base[2]+120,dest[2]+120)
   up=(cruise-base[2])/c['v_up'];down=(cruise-dest[2])/c['v_down'];travel=dist/c['v_cruise'];bu=(cruise-dest[2])/c['v_up'];bd=(cruise-base[2])/c['v_down'];start=max(900.+j*600,next_ready[j%4]);depart=start+c['prep'];arrival=depart+up+travel+down;end=arrival+c['dwell']+bu+travel+bd
   energy=2*c['p_cruise']*travel/3600+c['p_hover']*c['dwell']/3600+9.80665*((65+1.5)*(cruise-base[2])+65*(cruise-dest[2]))/(c['eta_up']*3600000.)
   usable=bool(energy<=c['e_use']*(1-c['rho']));P0=base;P1=np.array([base[0],base[1],cruise]);P2=np.array([dest[0],dest[1],cruise]);P3=dest
   segments=[(depart,depart+up,P0,P1),(depart+up,depart+up+travel,P1,P2),(depart+up+travel,arrival,P2,P3),(arrival,arrival+c['dwell'],P3,P3),(arrival+c['dwell'],arrival+c['dwell']+bu,P3,P2),(arrival+c['dwell']+bu,arrival+c['dwell']+bu+travel,P2,P1),(arrival+c['dwell']+bu+travel,end,P1,P0)];body='TypeB-'+str(j%4);battery=body+'-battery'
   if usable:next_ready[j%4]=end+charge(energy,c['e_use'],c['full_charge'])
  if usable:
   for a,b,P,Q in segments:phases.append({'K':K,'sid':sid,'t0':float(a),'t1':float(b),'P0':P.tolist(),'P1':Q.tolist()})
   transports.append({'sid':sid,'K':K,'type':mode,'drone':body,'battery':battery,'start':start,'end':end,'energy':energy,'chg':0 if mode=='communication' else charge(energy,c['e_use'],c['full_charge']),'deliver':{sid:arrival}})
  else:excluded.append({'sid':sid,'reason':'给定往返飞行剖面能量超限','required_kwh':energy,'available_kwh':3.2,'horizontal_necessary_kwh':2*c['p_cruise']*travel/3600})
  ledger.append({'箱ID':sid,'来源点位':q,'请求坐标':dest.tolist(),'期望送达秒':arrival+120,'硬截止秒':arrival+600,'计划架次':sid if usable else None,'计划送达秒':arrival if usable else None,'计划按期':usable,'计划硬截止满足':usable,'排除原因':None if usable else '给定往返剖面能量超限'})
 save(f/'contract.json',contract);save(f/'task.json',{'账本':ledger,'transport':transports,'phases':phases,'energy_exclusions':excluded,'scope':'点位和地形来自公开数据，服务单位和时间为预设模型条件；不代表历史行动结果'});print(tag,len(ledger),len(excluded),'能量排除',flush=True)
 return tag

def main():
 a=pd.read_csv(D/'定日真实居民点.csv');nodes=pd.read_csv(D/'CN_案例居民点.csv');n=nodes[nodes.geonameid.eq(9089058)].iloc[0];p=[{'lon':r.lon,'lat':r.lat,'geonameid':int(r.geonameid),'name':r.asciiname} for r in a.itertuples()]
 build('dingri_communication',p,{'lon':n.lon,'lat':n.lat,'name':n.asciiname,'geonameid':int(n.geonameid)},'dingri','communication')
 b=pd.read_csv(D/'Rasuwa确认受损建筑点.csv');tr=Transformer.from_crs(4326,32645,always_xy=True);x,y=tr.transform(b['经度'].values,b['纬度'].values);b['grid_x']=np.floor(x/2000).astype(int);b['grid_y']=np.floor(y/2000).astype(int);b['x']=x;b['y']=y;points=[]
 for (ix,iy),g in b.groupby(['grid_x','grid_y'],sort=True):
  r=g.assign(d2=(g.x-g.x.mean())**2+(g.y-g.y.mean())**2).sort_values(['d2','来源FID']).iloc[0];points.append({'lon':float(r['经度']),'lat':float(r['纬度']),'source_file':r['来源文件'],'source_layer':r['图层'],'source_fid':int(r['来源FID']),'grid_x':int(ix),'grid_y':int(iy),'confirmed_buildings_in_cell':len(g)})
 save(D/'Rasuwa标准服务点.json',points)
 nodes=pd.read_csv(D/'NP_案例居民点.csv')
 for key,id in [('dhunche',1283447),('betrawati',1283631)]:
  n=nodes[nodes.geonameid.eq(id)].iloc[0];build('rasuwa_'+key,points,{'lon':n.lon,'lat':n.lat,'name':n.asciiname,'geonameid':int(id)},'rasuwa','delivery')
if __name__=='__main__':main()
