from pathlib import Path
import sys,json,time,numpy as np
R=Path(__file__).resolve().parent;sys.path.insert(0,str(R));from model import Task
for tag in ['dingri_communication','rasuwa_dhunche','rasuwa_betrawati']:
 t=time.time();q=Task(tag,'full_mother');valid=np.flatnonzero(q.p.cap>1).tolist();pool=list(map(int,q.pool));pool+= [i for i in valid if i not in pool][:48-len(pool)];pool+= [i for i in range(len(q.p.palette)) if i not in pool][:48-len(pool)];pool=pool[:48]
 f=R/'inputs'/tag;np.savez_compressed(f/'input_aware_48.npz',palette=q.p.palette[pool],agl=q.p.agl[pool]);(f/'选址输入记录.json').write_text(json.dumps({'mother_indices':pool,'queries':len(q.queries),'direct_count':q.direct['on'],'fixed_profile_upper':q.upper,'rule':'必要通信查询的既有choose_pool，按有效点确定性补齐48；不读取最终服务成绩'},ensure_ascii=False,indent=2),encoding='utf8');print(tag,'直接',q.direct['on'],'任务剖面上界',q.upper,'必要查询',len(q.queries),'预处理秒',time.time()-t,flush=True)
