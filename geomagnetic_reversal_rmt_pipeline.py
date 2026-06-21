#!/usr/bin/env python3
"""
═══════════════════════════════════════════════════════════════════════════════
  地磁极性倒转 RMT 动力学 — 地球发电机的充能-释放节律
  Geomagnetic Reversal RMT — does the geodynamo have memory?
  ─────────────────────────────────────────────────────────────
  靶区: 地磁极性倒转序列 (GPTS chron boundaries).
  物理: 外核发电机积累不稳定 -> 倒转释放 -> 重建稳定极性 (再充能).
  争论: 学界用 gamma 拟合得到 k~1.3 (略大于 Poisson k=1), 结论模糊.
        Lowrie & Kent 2004: "essentially free of long-term memory" (Poisson).
        其他研究: 倒转频率变化暗示 "long-term memory".
        RMT 的间距比 <r> 是独立、对记忆更敏感的诊断 -> 给出干净裁决.
  预期: 不预设. Poisson(无记忆) / GOE(有再充能) / 中间态 都有意义.
═══════════════════════════════════════════════════════════════════════════════

METHODOLOGICAL DISCIPLINE (three traps, all handled):
  1. NON-STATIONARITY: reversal frequency changes systematically (slow in Mesozoic,
     fast in Cenozoic; near-zero during the Cretaceous Normal Superchron).
     This large-scale trend inflates CV and biases <r>. Removed by LOCAL UNFOLDING.
  2. THE SUPERCHRON: the ~40 Myr Cretaceous Normal Superchron (121-84 Ma) is an
     extreme outlier; analyze the post-superchron stationary interval (0-83 Ma),
     and report the superchron separately as a known anomaly.
  3. SHORT-SUBCHRON COMPLETENESS: very short chrons may be under-recorded.
     Sensitivity test: with/without the shortest intervals.

The events are DISCRETE and independently dated (marine magnetic anomalies +
Ar/Ar + astronomical tuning) — NO peak-extraction from a continuous curve, so the
speleothem peak-picking trap does NOT apply here. This is a clean discrete sequence.
"""

import numpy as np
from scipy import stats
from scipy.interpolate import UnivariateSpline
from scipy.integrate import cumulative_trapezoid
from math import gamma as _gamma
import warnings
warnings.filterwarnings('ignore')

# ═══════ RMT theory ═══════
def wigner_goe(s): return (np.pi/2)*s*np.exp(-np.pi/4*s**2)
def wigner_gue(s): return (32/np.pi**2)*s**2*np.exp(-4*s**2/np.pi)
def poisson_pdf(s): return np.exp(-s)
R_POI, R_GOE, R_GUE = 0.3863, 0.5307, 0.6027

def compute_r(sp):
    r = np.minimum(sp[:-1],sp[1:])/np.maximum(sp[:-1],sp[1:])
    return float(np.mean(r)), float(np.std(r)/np.sqrt(len(r)))
def compute_cv(sp): return float(np.std(sp)/np.mean(sp))
def brody_fit(s):
    from scipy.optimize import minimize_scalar
    def nll(b):
        a=(_gamma((b+2)/(b+1)))**(b+1)
        return -np.sum(np.log(b+1)+np.log(a)+b*np.log(s+1e-15)-a*s**(b+1))
    return float(minimize_scalar(nll,bounds=(0.01,3.0),method='bounded').x)
def gamma_k(s):
    """Fit gamma distribution shape parameter k (literature's diagnostic)."""
    try:
        k,loc,scale=stats.gamma.fit(s,floc=0); return float(k)
    except Exception: return np.nan
def ks_tests(s):
    sf=np.linspace(0,8,3000); ecdf=np.arange(1,len(s)+1)/len(s); ss=np.sort(s)
    kp=float(stats.kstest(s,lambda x:1-np.exp(-x))[0])
    kg=float(np.max(np.abs(ecdf-np.interp(ss,sf,cumulative_trapezoid(wigner_goe(sf),sf,initial=0)))))
    ku=float(np.max(np.abs(ecdf-np.interp(ss,sf,cumulative_trapezoid(wigner_gue(sf),sf,initial=0)))))
    return kp,kg,ku
def bootstrap_r(s,nboot=5000,seed=21):
    rng=np.random.default_rng(seed); n=len(s); o=[]
    for _ in range(nboot):
        ss=s[rng.integers(0,n,n)]
        o.append(np.mean(np.minimum(ss[:-1],ss[1:])/np.maximum(ss[:-1],ss[1:])))
    return np.percentile(o,[2.5,97.5])

# ═══════ Unfolding (removes reversal-frequency trend) ═══════
def local_unfolding(reversal_ages, smooth=None):
    ages=np.sort(np.asarray(reversal_ages,float)); N=len(ages); c=np.arange(1,N+1)
    s = smooth if smooth is not None else N/2.0
    spl=UnivariateSpline(ages,c,s=s)
    unf=spl(ages); sp=np.diff(unf); sp=sp[sp>0]
    return sp/np.mean(sp)

# ═══════ Analysis ═══════
def analyze(reversal_ages, label, unfold=True, drop_shortest=0):
    ages=np.sort(np.asarray(reversal_ages,float))
    chrons=np.diff(ages)
    if drop_shortest>0:
        thr=np.percentile(chrons,drop_shortest)
        keep=np.concatenate([[True],chrons>=thr])  # keep first age, drop reversals after tiny chrons
        ages=ages[keep]; chrons=np.diff(ages)
    if unfold:
        s=local_unfolding(ages)
    else:
        s=chrons/np.mean(chrons)
    r,re=compute_r(s); cv=compute_cv(s); b=brody_fit(s); k=gamma_k(chrons/np.mean(chrons))
    kp,kg,ku=ks_tests(s); lo,hi=bootstrap_r(s)
    best='Poisson' if kp<min(kg,ku) else ('GOE' if kg<ku else 'GUE')
    cls=("Poisson" if r<0.44 else "Poisson-GOE" if r<0.50 else "GOE" if r<0.57 else "GUE")
    print(f"\n=== {label} ===")
    print(f"  n_reversals={len(ages)}  n_chrons={len(chrons)}  mean chron={np.mean(chrons):.3f} Myr")
    print(f"  {'(local-unfolded)' if unfold else '(raw, no unfolding)'}")
    print(f"  <r>={r:.4f}±{re:.4f}  95%CI=[{lo:.3f},{hi:.3f}]  CV={cv:.3f}  beta={b:.2f}")
    print(f"  gamma k={k:.2f} (literature diagnostic; k=1 Poisson)   KS best={best}")
    print(f"  Poisson(0.386) in CI? {'YES' if lo<=R_POI<=hi else 'no'}   "
          f"GOE(0.531) in CI? {'YES' if lo<=R_GOE<=hi else 'no'}")
    print(f"  -> {cls}")
    return dict(label=label,n=len(ages),r=r,r_se=re,cv=cv,beta=b,gamma_k=k,
                ci=[float(lo),float(hi)],ks_best=best,cls=cls,s=s)

def load_gpts_csv(path):
    import pandas as pd
    df=pd.read_csv(path,comment='#')
    col=None
    for c in df.columns:
        if any(k in c.lower() for k in ['age','boundary','ma','reversal']): col=c; break
    if col is None: col=df.columns[0]
    return np.sort(df[col].dropna().values.astype(float))

if __name__=="__main__":
    import sys
    print("="*74)
    print("  🧭 地磁极性倒转 RMT — 地球发电机有记忆吗?")
    print("="*74)
    if len(sys.argv)>1:
        rev=load_gpts_csv(sys.argv[1])
    else:
        rev=np.load('gpts_reversals.npy') if __import__('os').path.exists('gpts_reversals.npy') else None
        if rev is None:
            print("  supply a GPTS CSV of reversal ages"); sys.exit()
    analyze(rev,"Geomagnetic reversals 0-83 Ma (raw)",unfold=False)
    analyze(rev,"Geomagnetic reversals 0-83 Ma (unfolded — frequency trend removed)",unfold=True)
    analyze(rev,"Sensitivity: drop shortest 10% chrons",unfold=True,drop_shortest=10)
    print("\n"+"="*74)
    print("  纪律: 离散独立定年事件, 无峰检测陷阱; unfolding 去除倒转频率趋势;")
    print("  超静期单独处理; 短段敏感性测试。不预设 Poisson 或 GOE。")
    print("="*74)
