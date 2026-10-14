"""普通连续候选MILP的部分服务扩展，不宣称新算法或全域最优。"""
import inspect
from general_solver_core import solve_candidate

source=inspect.getsource(solve_candidate)
source=source.replace('def solve_candidate(problem,N,B,k,counts,pool,queries,cover,seconds=30,margin=.01):',
    'def solve_service_candidate(problem,N,B,k,counts,pool,queries,cover,weights,seconds=30,margin=.01):')
assert 'S=[var(c=1e-3)' in source and 'for i,q in enumerate(queries):' in source
source=source.replace('S=[var(c=1e-3)','S=[var(c=1e-9)').replace('E=[var(c=-9e-4)','E=[var(c=-9e-10)')
source=source.replace(' for i,q in enumerate(queries):',
    " served={K:var(hi=1 if w>0 else 0,c=-float(w),binary=True) for K,w in weights.items()}\n for i,q in enumerate(queries):")
source=source.replace('  row({y:1 for y in options},lo=1)',
    "  conditional={y:1 for y in options};conditional[served[q['K']]]=-1\n  row(conditional,lo=0)")
source=source.replace("'mip_rel_gap':.05","'mip_rel_gap':0.")
source=source.replace(' rows=[]\n for j in range(slots):',
    " info['选择执行K']=[K for K,varid in served.items() if answer.x[varid]>.5]\n rows=[]\n for j in range(slots):")
source=source.replace("'mip_rel_gap':0.","'mip_rel_gap':0.,'threads':1")
namespace=dict(solve_candidate.__globals__)
exec(compile(source,__file__,'exec'),namespace)
solve_service_candidate=namespace['solve_service_candidate']
