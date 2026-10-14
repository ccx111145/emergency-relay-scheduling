"""TRE统一验收：精确机体/电池资源与完整飞行阶段通信。"""
from fractions import Fraction as F
from itertools import product
import numpy as np
import independent_radio_geometry as G

ff=lambda x:F.from_float(float(x))
def resource_options(p,rows,full,k=1):
    """同一固定时序的全部机体划分；共享电池区间最优着色。"""
    rt=p.rt;exact=[]
    for r in rows:
        site=r['point_index'];S=ff(r['link_done']);E=ff(r['end'])
        a=S-ff(p.out[site])-F(str(p.link))-F(str(p.prep));b=E+ff(p.back[site])
        e=ff(p.fly[site])+F(str(p.power))*(E-S+F(str(p.link)))/3600
        if a<0 or E<=S or b-(a+F(str(p.prep)))>p.limit or e>F(str(rt['e_use']))*(1-F(str(rt['rho']))):return []
        cg=min(F(full)*F(7,2)*e/F(str(rt['e_use'])),F(full)*(F(5,18)+F(13,18)*e/F(str(rt['e_use']))))
        exact.append({'sid':r['sid'],'prepare':a,'depart':a+F(str(p.prep)),'S':S,'E':E,'return':b,'release':b+cg,'energy':e})
    order=sorted(range(len(exact)),key=lambda i:exact[i]['prepare'])
    releases=[];pack={}
    for i in order:
        found=next((j for j,v in enumerate(releases) if v<=exact[i]['prepare']),None)
        if found is None:found=len(releases);releases.append(F(0))
        releases[found]=exact[i]['release'];pack[i]=found
    if not exact:return [{'N':0,'B':0,'k_peak':0,'rows':[],'energy':0.,'资源分数':[]}]
    answers=[]
    # 限制增长序列枚举无标签划分，每种机体划分只计算一次。
    def partitions(labels):
        if len(labels)==len(exact):yield labels;return
        for label in range(max(labels,default=-1)+2):yield from partitions(labels+[label])
    for labels in partitions([]):
        chains={};valid=True
        for i,label in enumerate(labels):chains.setdefault(label,[]).append(i)
        windows=[]
        for chain in chains.values():
            chain.sort(key=lambda i:exact[i]['S'])
            for a,b in zip(chain,chain[1:]):
                if exact[b]['prepare']<exact[a]['return']+F(str(p.turn)):valid=False;break
            if not valid:break
            windows += [(exact[chain[0]]['depart'],exact[chain[0]]['S']),(exact[chain[-1]]['E'],exact[chain[-1]]['return'])]
            windows += [(exact[a]['E'],exact[b]['S']) for a,b in zip(chain,chain[1:])]
        if not valid:continue
        events={}
        for a,b in windows:
            events[a]=events.get(a,0)+1;events[b]=events.get(b,0)-1
        occupancy=0;peak=0
        for _,delta in sorted(events.items()):occupancy+=delta;peak=max(peak,occupancy)
        if peak>k:continue
        allocated=[{**r,'drone':f'R{labels[i]+1:02}','comp':f'电池{pack[i]}'} for i,r in enumerate(rows)]
        answers.append({'N':len(chains),'B':len(releases),'k_peak':peak,'rows':allocated,
            'energy':float(sum(r['energy'] for r in exact)),
            '资源分数':[{key:str(v) if isinstance(v,F) else v for key,v in r.items()} for r in exact]})
    # 机体数最小的合法划分足以用于容量网格，同计划电池着色与划分无关。
    return sorted(answers,key=lambda r:r['N'])[:1]

def coverage_sortie(p,rows,K):
    """逐运输架次查全连续阶段，返回PASS、FAIL或UNKNOWN。"""
    available=[]
    for r in rows:
        loc=p.palette[r['point_index']]
        bh=G.certified_point_link(p.gw,loc,p.dem,p.cp,p.budgets['backhaul'],p.radius)
        if bh is not True:continue
        available.append((r,loc))
    leaves=0;endpoints=0;unknown=[]
    phases=[(j,ph) for j,ph in enumerate(p.phases) if ph['K']==K]
    if not phases:return {'状态':'UNKNOWN','原因':'无飞行阶段'}
    for j,ph in phases:
        cuts=sorted({ph['t0'],ph['t1']}|{r[key] for r in rows for key in ['link_done','end'] if ph['t0']<r[key]<ph['t1']})
        for t in cuts:
            pos=p.position(ph,t)
            carriers=[(p.gw,p.budgets['direct'])]+[(loc,p.budgets['access']) for r,loc in available if r['link_done']<=t<r['end']]
            checks=[G.certified_point_link(a,pos,p.dem,p.cp,budget,p.radius) for a,budget in carriers];endpoints+=1
            if not any(v is True for v in checks):
                if all(v is False for v in checks):return {'状态':'FAIL','反例':{'K':K,'phase_index':j,'t':t,'P':pos.tolist()}}
                unknown.append([j,t,t])
        for a,b in zip(cuts,cuts[1:]):
            mid=(a+b)/2
            carriers=[(p.gw,p.budgets['direct'])]+[(loc,p.budgets['access']) for r,loc in available if r['link_done']<=mid<r['end']]
            stack=[(a,b,0)]
            while stack:
                x,y,depth=stack.pop();p0,p1=p.position(ph,x),p.position(ph,y)
                if any(G.interval_link(loc,p0,p1,p.dem,p.cp,budget,p.radius)[0] for loc,budget in carriers):leaves+=1;continue
                t=(x+y)/2;pos=p.position(ph,t)
                checks=[G.certified_point_link(loc,pos,p.dem,p.cp,budget,p.radius) for loc,budget in carriers]
                if all(v is False for v in checks):return {'状态':'FAIL','反例':{'K':K,'phase_index':j,'t':t,'P':pos.tolist()}}
                if y-x<.002 or depth>=24:unknown.append([j,x,y]);continue
                stack.extend([(t,y,depth+1),(x,t,depth+1)])
    return {'状态':'UNKNOWN' if unknown else 'PASS','区间叶数':leaves,'端点数':endpoints,'未知区间':unknown}
