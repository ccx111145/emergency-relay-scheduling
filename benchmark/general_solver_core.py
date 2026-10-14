"""便携的既有候选MILP及连续验收函数；不依赖私人工作簿或旧实验路径。"""
from itertools import combinations
from fractions import Fraction as F
import math,time
import numpy as np
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix
import independent_radio_geometry as G

def charge(e,rt,full):
 used=np.asarray(e,dtype=float)/rt['e_use']
 return np.minimum(full*3.5*used,full*(.2777777777777778+.65/.9*used))

def dominance_pool(problem,cover):
 """只对当前必要点模型作支配消元，每轮新点加入后从全部场址重新构建。"""
 masks=[int.from_bytes(np.packbits(row,bitorder='little').tobytes(),'little') for row in cover]
 score=problem.out+problem.back+1000*problem.fly
 order=np.argsort(score,kind='stable');kept=[]
 for p in order:
  if problem.cap[p]<0:continue
  dominated=False
  for q in kept:
   if problem.out[q]<=problem.out[p] and problem.back[q]<=problem.back[p] and problem.fly[q]<=problem.fly[p] and (masks[p]&masks[q])==masks[p]:dominated=True;break
  if not dominated:kept.append(int(p))
 return np.array(kept,int)

def choose_pool(problem,queries,cover,budget):
 """纯输入驱动选点：不给已知见证场址硬编码，不将选择失败用于下界。"""
 if budget==0:return dominance_pool(problem,cover)
 score=cover.sum(axis=1)/(1.+(problem.out+problem.back)/2000)
 ranked=np.argsort(-score,kind='stable');chosen=[]
 # 先覆盖尚未支持的必要请求，再按需求分块选择局部强点。
 missing=np.ones(cover.shape[1],bool)
 while missing.any() and len(chosen)<budget:
  gain=cover[:,missing].sum(axis=1)/(1.+(problem.out+problem.back)/2000)
  if chosen:gain[chosen]=-1
  p=int(np.argmax(gain))
  if gain[p]<=0:break
  chosen.append(p);missing &= ~cover[p]
 chunks=np.array_split(np.arange(len(queries)),max(1,budget//2))
 for chunk in chunks:
  if not len(chunk):continue
  gain=cover[:,chunk].sum(axis=1)/(1.+(problem.out+problem.back)/2000)
  for p in np.argsort(-gain,kind='stable')[:3]:
   if gain[p]>0 and int(p) not in chosen:chosen.append(int(p));break
  if len(chosen)>=budget:break
 for p in ranked:
  if len(chosen)>=budget:break
  if score[p]>0 and int(p) not in chosen:chosen.append(int(p))
 return np.array(chosen,int)

def solve_candidate(problem,N,B,k,counts,pool,queries,cover,seconds=30,margin=.01):
 """连续S/E与悬停选择MILP，精确线性化充满模型和组件身份接续。

 架次数量、机体归属顺序和候选池受限，故负答不能证明原模型不可行。
 """
 rt=problem.rt;M=max(p['t1'] for p in problem.phases)+20000.;M*=4;horizon=M/2;slots=sum(counts)
 lower=[];upper=[];cost=[];integer=[];constraints=[];rowlo=[];rowhi=[]
 def var(lo=0.,hi=horizon,c=0.,binary=False):
  i=len(lower);lower.append(lo);upper.append(hi);cost.append(c);integer.append(int(binary));return i
 def row(c,lo=-math.inf,hi=math.inf):constraints.append(c);rowlo.append(lo);rowhi.append(hi)
 def add(a,b,scale=1.):
  d=dict(a)
  for i,v in b.items():d[i]=d.get(i,0.)+scale*v
  return d
 # 原来只奖励E-S，整条轨迹整体后移仍同价，会贴着覆盖断点反复平移。
 # 保留驻留奖励并对绝对开航时刻加正成本；仅为候选搜索目标，不改变可行域。
 S=[var(c=1e-3) for _ in range(slots)];E=[var(c=-9e-4) for _ in range(slots)]
 X=[[var(hi=1,binary=True) for p in pool] for _ in range(slots)]
 energy=[var(hi=problem.emax) for _ in range(slots)];chg=[var(hi=problem.full) for _ in range(slots)]
 pack=[[var(hi=1,binary=True) for _ in range(B)] for _ in range(slots)]
 out=[];back=[];prepare=[];returned=[];released=[]
 for j in range(slots):
  row({x:1 for x in X[j]},lo=1,hi=1);row({x:1 for x in pack[j]},lo=1,hi=1)
  o={x:float(problem.out[p]) for x,p in zip(X[j],pool)};b={x:float(problem.back[p]) for x,p in zip(X[j],pool)};out.append(o);back.append(b)
  a=add({S[j]:1},o,-1);prepare.append(a)
  r=add({E[j]:1},b);returned.append(r);released.append(add(r,{chg[j]:1}))
  row(a,lo=problem.prep+problem.link+margin)
  row({E[j]:1,S[j]:-1},lo=.1)
  row(add({E[j]:1,S[j]:-1},{x:float(problem.cap[p]) for x,p in zip(X[j],pool)},-1),hi=-margin)
  coefficients=add({energy[j]:1,E[j]:-problem.power/3600,S[j]:problem.power/3600},{x:float(problem.fly[p]) for x,p in zip(X[j],pool)},-1)
  row(coefficients,lo=problem.power*problem.link/3600,hi=problem.power*problem.link/3600)
  fast=var(hi=1,binary=True);threshold=rt['e_use']*.1
  row({energy[j]:1,fast:-M},lo=threshold-M);row({energy[j]:1,fast:-M},hi=threshold)
  slow=problem.full*3.5/rt['e_use'];rapid=problem.full*.65/.9/rt['e_use'];intercept=problem.full*.2777777777777778
  row({chg[j]:1,energy[j]:-slow,fast:-M},hi=0);row({chg[j]:1,energy[j]:-slow,fast:M},lo=0)
  row({chg[j]:1,energy[j]:-rapid,fast:M},hi=intercept+M);row({chg[j]:1,energy[j]:-rapid,fast:-M},lo=intercept-M)
 for i,q in enumerate(queries):
  options=[]
  for j in range(slots):
   y=var(hi=1,binary=True);options.append(y)
   row(add({y:1},{x:float(cover[p,i]) for x,p in zip(X[j],pool)},-1),hi=0)
   row({S[j]:1,y:M},hi=q['t']+M-margin);row({E[j]:1,y:-M},lo=q['t']+margin-M)
  row({y:1 for y in options},lo=1)
 bodyid=[];windows=[];startslot=0
 for body,num in enumerate(counts):
  chain=list(range(startslot,startslot+num));startslot+=num;bodyid.extend([body]*num)
  for a,b in zip(chain,chain[1:]):row(add(add({S[b]:1,E[a]:-1},out[b],-1),back[a],-1),lo=problem.turn+problem.prep+problem.link+margin)
  first,last=chain[0],chain[-1]
  windows.append((add({S[first]:1},out[first],-1),-problem.link,{S[first]:1},0.))
  for a,b in zip(chain,chain[1:]):windows.append(({E[a]:1},0.,{S[b]:1},0.))
  windows.append(({E[last]:1},0.,returned[last],0.))
 # 容量k：任意k+1个区间不能有共同正长度交集；一维区间等价于其中存在一对有序不交。
 if k<N:
  for subset in combinations(range(len(windows)),k+1):
   selectors=[]
   for a,b in ((a,b) for a in subset for b in subset if a!=b):
    y=var(hi=1,binary=True);selectors.append(y);wa,wb=windows[a],windows[b]
    c=add(wa[2],wb[0],-1);c[y]=M;row(c,hi=M+wb[1]-wa[3]-margin)
   row({y:1 for y in selectors},lo=1)
 # 同一个组件从准备至返航后充满均被占用；精确组件身份而非0/1路径占用裁剪。
 for a,b in combinations(range(slots),2):
  y=var(hi=1,binary=True)
  for component in range(B):
   c=add(released[a],prepare[b],-1);c[y]=M;c[pack[a][component]]=M;c[pack[b][component]]=M
   row(c,hi=3*M-problem.prep-problem.link-margin)
   c=add(released[b],prepare[a],-1);c[y]=-M;c[pack[a][component]]=M;c[pack[b][component]]=M
   row(c,hi=2*M-problem.prep-problem.link-margin)
 # 组件标签的对称性：第一架次使用组件0。
 row({pack[0][0]:1},lo=1,hi=1)
 rr=[];cc=[];vv=[]
 for i,c in enumerate(constraints):
  for j,v in c.items():
   if v:rr.append(i);cc.append(j);vv.append(v)
 matrix=coo_matrix((vv,(rr,cc)),shape=(len(constraints),len(lower))).tocsc()
 start=time.perf_counter();answer=milp(np.array(cost),integrality=np.array(integer),bounds=Bounds(lower,upper),constraints=LinearConstraint(matrix,np.array(rowlo),np.array(rowhi)),options={'time_limit':seconds,'mip_rel_gap':.05})
 info={'求解器状态':int(answer.status),'消息':answer.message,'耗时秒':time.perf_counter()-start,'变量数':len(lower),'约束数':len(constraints),'候选数':len(pool),'架次骨架':counts,'N':N,'B':B,'k':k,'范围':'候选搜索；时间连续，但候选池及架次骨架受限；失败一律UNKNOWN'}
 if answer.x is None:return None,info
 rows=[]
 for j in range(slots):
  p=int(pool[np.argmax(answer.x[X[j]])]);start,end=answer.x[S[j]],answer.x[E[j]];e=problem.fly[p]+problem.power*(end-start+problem.link)/3600
  rows.append({'sid':f'{problem.tag}-G{j}','drone':f'R{bodyid[j]+1:02}','point_index':p,'lon':float(problem.palette[p,0]),'lat':float(problem.palette[p,1]),'hover_alt':float(problem.palette[p,2]),'link_done':float(start),'end':float(end),'depart':float(start-problem.out[p]-problem.link),'prepare':float(start-problem.out[p]-problem.link-problem.prep),'back':float(end+problem.back[p]),'e_total':float(e),'chg':float(charge(e,rt,problem.full)),'comp':f'组件{np.argmax(answer.x[pack[j]])}'})
 return rows,info

def verify(problem,rows,N,B,k):
 """独立全阶段验收：资源分数化重算，连续区间LP与真实端点。"""
 issues=[];windows=[];bodies={};packs={};exact=[];rt=problem.rt
 for r in rows:
  p=r['point_index'];S,E=F.from_float(r['link_done']),F.from_float(r['end']);a=S-F.from_float(problem.out[p])-F(str(problem.link))-F(str(problem.prep));b=E+F.from_float(problem.back[p]);energy=F.from_float(problem.fly[p])+F(str(problem.power))*(E-S+F(str(problem.link)))/3600
  # 充电曲线独立按SOC定义重算；Fraction保证时间与库存边界不被EPS放过。
  used=energy/F(str(rt['e_use']));slow=F(str(problem.full))*F('3.5')*used;fast=F(str(problem.full))*(F(5,18)+F(13,18)*used);cg=min(slow,fast)
  exact_budget=F(str(rt['e_use']))*(1-F(str(rt['rho'])))
  if a<0 or E<=S or b-(a+F(str(problem.prep)))>problem.limit or energy>exact_budget:issues.append({'类别':'准备/滞空/能量不满足','sid':r['sid']})
  record={'sid':r['sid'],'prepare':a,'depart':a+F(str(problem.prep)),'S':S,'E':E,'return':b,'pack_release':b+cg};exact.append(record);bodies.setdefault(r['drone'],[]).append(record);packs.setdefault(r['comp'],[]).append(record)
 if len(bodies)>N or len(packs)>B:issues.append({'类别':'机体或组件数超限'})
 for chain in bodies.values():
  chain.sort(key=lambda r:r['S']);windows.extend([(chain[0]['depart'],chain[0]['S']),(chain[-1]['E'],chain[-1]['return'])])
  for a,b in zip(chain,chain[1:]):
   if b['prepare']<a['return']+F(str(problem.turn)):issues.append({'类别':'机体接续违反'})
   windows.append((a['E'],b['S']))
 for chain in packs.values():
  chain.sort(key=lambda r:r['prepare'])
  for a,b in zip(chain,chain[1:]):
   if b['prepare']<a['pack_release']:issues.append({'类别':'组件未充满即复用'})
 events={}
 for a,b in windows:
  events[a]=events.get(a,0)+1;events[b]=events.get(b,0)-1
 count=0;peak=0
 for t,delta in sorted(events.items()):count+=delta;peak=max(peak,count)
 if peak>k:issues.append({'类别':'同时黑障窗口数超k','峰值':peak})
 if issues:return {'状态':'FAIL_RESOURCE','问题':issues}
 available=[]
 for r in rows:
  p=np.array([r['lon'],r['lat'],r['hover_alt']])
  bh=G.certified_point_link(problem.gw,p,problem.dem,problem.cp,problem.budgets['backhaul'],problem.radius)
  if bh is not True:return {'状态':'UNKNOWN_BACKHAUL' if bh is None else 'FAIL_BACKHAUL','sid':r['sid']}
  available.append((r,p))
 leaves=0;endpoints=0;unknown=[]
 for j,p in enumerate(problem.phases):
  cuts=sorted({p['t0'],p['t1']}|{r[key] for r in rows for key in ['link_done','end'] if p['t0']<r[key]<p['t1']})
  for t in cuts:
   endpoint=problem.position(p,t);carriers=[(problem.gw,problem.budgets['direct'])]+[(loc,problem.budgets['access']) for r,loc in available if r['link_done']<=t<r['end']];endpoints+=1
   checks=[G.certified_point_link(a,endpoint,problem.dem,problem.cp,budget,problem.radius) for a,budget in carriers]
   if not any(v is True for v in checks):
    if all(v is False for v in checks):return {'状态':'FAIL_COVERAGE','witness':{'t':t,'P':endpoint.tolist(),'K':p['K'],'phase_index':j}}
    unknown.append([j,t,t])
  for a,b in zip(cuts,cuts[1:]):
   mid=(a+b)/2;carriers=[(problem.gw,problem.budgets['direct'])]+[(loc,problem.budgets['access']) for r,loc in available if r['link_done']<=mid<r['end']];stack=[(a,b,0)]
   while stack:
    x,y,depth=stack.pop();p0,p1=problem.position(p,x),problem.position(p,y)
    if any(G.interval_link(loc,p0,p1,problem.dem,problem.cp,budget,problem.radius)[0] for loc,budget in carriers):leaves+=1;continue
    t=(x+y)/2;pt=problem.position(p,t);checks=[problem.link_ok(loc,pt,budget) for loc,budget in carriers]
    if all(v is False for v in checks):return {'状态':'FAIL_COVERAGE','witness':{'t':t,'P':pt.tolist(),'K':p['K'],'phase_index':j}}
    if y-x<.002 or depth>=24:unknown.append([j,x,y]);continue
    stack.extend([(t,y,depth+1),(x,t,depth+1)])
 return {'状态':'UNKNOWN' if unknown else 'PASS','区间叶数':leaves,'端点数':endpoints,'UNKNOWN':unknown,'组件数':len(packs),'机体数':len(bodies),'黑障峰值':peak,'范围':'独立数值几何和精确分数资源验收；非机器区间或野外实测。'}
