"""逐像元开放射线的数组实现，只用于候选支持；最终验收保留原内核。"""
import numpy as np
import independent_radio_geometry as G

def ray_status_vector(a,b,dem):
 a,b=np.asarray(a,float),np.asarray(b,float)
 if not np.isfinite(np.r_[a,b]).all():return -1
 start,end=dem.xy(a),dem.xy(b)
 events=np.unique(np.array([0.,1.]+list(G.axis_events(start[0],end[0]))+list(G.axis_events(start[1],end[1])),float))
 left,right=events[:-1],events[1:];mid=(left+right)/2
 cells=np.rint(start[None,:]+(end-start)[None,:]*mid[:,None]).astype(np.int64);c,row=cells[:,0],cells[:,1]
 if np.any(row<0) or np.any(c<0) or np.any(row>=dem.z.shape[0]) or np.any(c>=dem.z.shape[1]):return -1
 terrain=dem.z[row,c].astype(float)
 if not np.isfinite(terrain).all() or np.any(terrain==dem.nodata):return -1
 z0=a[2]+(b[2]-a[2])*left;z1=a[2]+(b[2]-a[2])*right
 return int(np.any(np.minimum(z0,z1)<terrain-1e-9))

def support_vector(p,q,cache):
 def link(a,b,budget):
  loss=G.loss(G.distance(a,b,p.radius),p.cp['f_mhz'])
  if loss>budget+1e-9:return False
  if loss+p.cp['l_obs']<=budget+1e-9:return True
  key=(tuple(np.asarray(a,float)),tuple(np.asarray(b,float)))
  if key not in cache:cache[key]=ray_status_vector(a,b,p.dem)
  flag=cache[key]
  if flag<0:return None
  return loss+flag*p.cp['l_obs']<=budget+1e-9
 P=np.array(q['P']);column=[]
 for i,point in enumerate(p.palette):
  access=link(point,P,p.budgets['access'])
  if access is False:column.append(False);continue
  if i not in p.bh_cache:p.bh_cache[i]=link(p.gw,point,p.budgets['backhaul'])
  column.append(p.bh_cache[i] is not False)
 return np.array(column,bool)
