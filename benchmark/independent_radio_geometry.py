"""独立数值内核：双事件流逐像元射线、重心坐标二维线性规划连续扇面。"""
from dataclasses import dataclass
from itertools import combinations
import math
import numpy as np

@dataclass
class Raster:
 z:np.ndarray
 lon:np.ndarray
 lat:np.ndarray
 nodata:float
 def __post_init__(self):
  self.dx=float(self.lon[1]-self.lon[0]);self.dy=float(self.lat[1]-self.lat[0])
 def xy(self,p):return np.array([(p[0]-self.lon[0])/self.dx,(p[1]-self.lat[0])/self.dy])

def distance(a,b,radius):
 # 独立写出球面中心角及高差，单位均为米。
 la,lb=map(math.radians,(a[1],b[1]));dl=math.radians(b[0]-a[0])
 h=math.sin((lb-la)/2)**2+math.cos(la)*math.cos(lb)*math.sin(dl/2)**2
 horizontal=2*radius*math.atan2(math.sqrt(max(0.,h)),math.sqrt(max(0.,1-h)))
 return math.hypot(horizontal,b[2]-a[2])

def loss(distance_m,frequency):
 return 32.45+20*math.log10(frequency)+20*math.log10(max(distance_m/1000,1e-9))

def axis_events(start,finish):
 if finish==start:return iter(())
 lo,hi=min(start,finish),max(start,finish)
 events=((i+.5-start)/(finish-start) for i in range(math.ceil(lo-.5),math.floor(hi-.5)+1))
 values=[v for v in events if 0<v<1]
 if finish<start:values.reverse()
 return iter(values)

def ray_status(a,b,dem):
 """开放射线，每个像元内检查线性高度下确界；未知优先于遮挡。"""
 a,b=np.asarray(a,float),np.asarray(b,float)
 if not np.isfinite(np.r_[a,b]).all():return -1
 start,end=dem.xy(a),dem.xy(b);streams=[axis_events(start[k],end[k]) for k in range(2)]
 nexts=[next(it,1.) for it in streams];left=0.;blocked=False
 while left<1:
  right=min(nexts);mid=(left+right)/2
  c,row=map(lambda v:int(round(v)),start+(end-start)*mid)
  if not (0<=row<dem.z.shape[0] and 0<=c<dem.z.shape[1]):return -1
  terrain=float(dem.z[row,c])
  if not math.isfinite(terrain) or terrain==dem.nodata:return -1
  z0=a[2]+(b[2]-a[2])*left;z1=a[2]+(b[2]-a[2])*right
  if min(z0,z1)<terrain-1e-9:blocked=True
  left=right
  if right>=1:break
  for k in range(2):
   if nexts[k]<=right:nexts[k]=next(streams[k],1.)
 return int(blocked)

def point_link(a,b,dem,cp,budget,radius):
 l=loss(distance(a,b,radius),cp['f_mhz'])
 if l>budget+1e-9:return False
 if l+cp['l_obs']<=budget+1e-9:return True
 flag=ray_status(a,b,dem)
 if flag<0:return None
 return l+flag*cp['l_obs']<=budget+1e-9

def ray_closed_clear(a,b,dem,margin=1e-5):
 """保守点证书包含像元边界和孤立角点接触，不仅检查正长度射线段。"""
 a,b=np.asarray(a,float),np.asarray(b,float);u,v=dem.xy(a),dem.xy(b)
 streams=[axis_events(u[k],v[k]) for k in range(2)];nexts=[next(it,1.) for it in streams];left=0.
 def terrain_safe(cols,rows,height):
  for row in rows:
   for col in cols:
    if not(0<=row<dem.z.shape[0] and 0<=col<dem.z.shape[1]):return False
    z=float(dem.z[row,col])
    if not math.isfinite(z) or z==dem.nodata or height<z+margin:return False
  return True
 def touching(t):
  point=u+(v-u)*t;cells=[]
  for x in point:
   nearest=int(round(x));boundary=math.floor(x)+.5
   cells.append([math.floor(x),math.floor(x)+1] if abs(x-boundary)<=1e-9 else [nearest])
  return terrain_safe(cells[0],cells[1],float(a[2]+(b[2]-a[2])*t))
 if not touching(0.):return False
 while left<1:
  right=min(nexts);mid=(left+right)/2;col,row=map(lambda x:int(round(x)),u+(v-u)*mid)
  height=min(a[2]+(b[2]-a[2])*left,a[2]+(b[2]-a[2])*right)
  if not terrain_safe([col],[row],float(height)) or not touching(right):return False
  left=right
  if right>=1:break
  for k in range(2):
   if nexts[k]<=right:nexts[k]=next(streams[k],1.)
 return True

def certified_point_link(a,b,dem,cp,budget,radius):
 """上界单独使用保守点证书；下界可能支持集合仍用point_link，不强删候选。"""
 possible=point_link(a,b,dem,cp,budget,radius)
 if possible is not True:return possible
 if loss(distance(a,b,radius),cp['f_mhz'])+cp['l_obs']<=budget-1e-6:return True
 return True if ray_closed_clear(a,b,dem) else None

def interval_distance_bound(r,a,b,radius):
 """中点到固定端点距离，加上路径半程长度上界，采用乘积空间三角不等式。"""
 mid=(np.asarray(a)+b)/2
 near=0. if min(a[1],b[1])<=0<=max(a[1],b[1]) else min(abs(a[1]),abs(b[1]))
 lat=radius*math.radians(b[1]-a[1]);lon=radius*math.cos(math.radians(near))*math.radians(b[0]-a[0])
 half=math.sqrt(lat*lat+lon*lon+(b[2]-a[2])**2)/2
 return distance(r,mid,radius)+half

def fan_lp_clear(r,a,b,dem):
 """在重心参数(u,v)中对每个像元枚举二维LP顶点，未复用原多边形裁剪代码。"""
 r,a,b=np.asarray(r),np.asarray(a),np.asarray(b)
 xy=dem.xy(r);delta=np.column_stack([dem.xy(a)-xy,dem.xy(b)-xy]);zdelta=np.array([a[2]-r[2],b[2]-r[2]])
 points=np.array([dem.xy(p) for p in (r,a,b)])
 lower=np.ceil(points.min(axis=0)-.5-1e-10).astype(int);upper=np.floor(points.max(axis=0)+.5+1e-10).astype(int)
 if lower[0]<0 or lower[1]<0 or upper[0]>=dem.z.shape[1] or upper[1]>=dem.z.shape[0]:return False,{'原因':'扇面矩形含越界像元'}
 cols,rows=np.meshgrid(np.arange(lower[0],upper[0]+1),np.arange(lower[1],upper[1]+1))
 cols,rows=cols.ravel(),rows.ravel();count=0;minimum_margin=math.inf
 # u>=0,v>=0,u+v<=1，加上像元的四个矩形边界。
 A=np.array([[-1.,0.],[0.,-1.],[1.,1.],delta[0],-delta[0],delta[1],-delta[1]])
 pairs=[]
 for i,j in combinations(range(7),2):
  det=A[i,0]*A[j,1]-A[i,1]*A[j,0]
  if det!=0:
   if abs(det)<1e-14:return False,{'原因':'近退化顶点，需细分或有理数复核'}
   pairs.append((i,j,det))
 for start in range(0,len(rows),4096):
  rr,cc=rows[start:start+4096],cols[start:start+4096];n=len(rr)
  B=np.column_stack([np.zeros(n),np.zeros(n),np.ones(n),cc+.5-xy[0]+1e-10,xy[0]-cc+.5+1e-10,rr+.5-xy[1]+1e-10,xy[1]-rr+.5+1e-10])
  low=np.full(n,np.inf)
  for i,j,det in pairs:
   u=(B[:,i]*A[j,1]-A[i,1]*B[:,j])/det;v=(A[i,0]*B[:,j]-B[:,i]*A[j,0])/det
   # 允许更宽可行域并留高度安全余量，作为保守数值检查。
   valid=np.ones(n,bool)
   for k in range(7):valid &= A[k,0]*u+A[k,1]*v<=B[:,k]+1e-9
   height=r[2]+zdelta[0]*u+zdelta[1]*v
   low=np.minimum(low,np.where(valid,height,np.inf))
  mask=np.isfinite(low);count+=int(mask.sum())
  ground=dem.z[rr,cc]
  if np.any(mask & (~np.isfinite(ground)|(ground==dem.nodata))):return False,{'原因':'扇面含nodata'}
  margin=low-ground
  if mask.any():minimum_margin=min(minimum_margin,float(margin[mask].min()))
  bad=np.flatnonzero(mask & (margin<1e-5))
  if len(bad):
   j=int(bad[0]);return False,{'原因':'未能证明整个扇面净空','像元':[int(rr[j]),int(cc[j])],'高度下界':float(low[j]),'地面':float(ground[j])}
 return count>0,{'算法':'重心坐标LP顶点枚举','像元数':count,'最小净空':minimum_margin}

def interval_link(r,a,b,dem,cp,budget,radius):
 bound=interval_distance_bound(r,a,b,radius);l=loss(bound,cp['f_mhz'])
 if l+cp['l_obs']<=budget-1e-6:return True,{'依据':'所有二值遮挡状态均满足预算','距离上界':bound}
 if l>budget-1e-6:return False,{'原因':'距离上界未闭合'}
 okay,proof=fan_lp_clear(r,a,b,dem)
 return okay,{'距离上界':bound,**proof}

def budget_pair(a,b,cp):
 threshold=cp['p_sens']+cp['margin']
 return min(a[0]+a[1]+b[2],b[0]+b[1]+a[2])-cp['l_sys']-threshold

def relay_flight(point,agl,dem,node,rt,radius,gravity,joule):
 start=[node['lon'],node['lat'],node['elev']];d=distance(start,[point[0],point[1],node['elev']],radius)
 n=max(2,math.ceil(d/25)+1);f=np.arange(n,dtype=float)/(n-1)
 lon=start[0]+(point[0]-start[0])*f;lat=start[1]+(point[1]-start[1])*f
 col=np.rint((lon-dem.lon[0])/dem.dx).astype(int);row=np.rint((lat-dem.lat[0])/dem.dy).astype(int)
 if np.any((col<0)|(col>=dem.z.shape[1])|(row<0)|(row>=dem.z.shape[0])):raise ValueError('飞行剖面越界')
 gmax=float(dem.z[row,col].max());c0,r0=map(lambda v:int(round(v)),dem.xy(point));alt=float(dem.z[r0,c0])+agl
 cruise=max(gmax+50,node['elev']+50,alt);h1=cruise-node['elev'];h2=cruise-alt;horizontal=d/rt['v_cruise']
 tout=h1/rt['v_up']+horizontal+h2/rt['v_dn'];tback=h2/rt['v_up']+horizontal+h1/rt['v_dn']
 efly=2*rt['p_cruise']*horizontal/3600+rt['mtow']*gravity*(h1+h2)/(rt['eta_up']*joule)
 return tout,tback,efly,alt
