#!/usr/bin/env python3
"""Module score (formula 10) - committed generator + stability analysis.
Gap closed 2026-09-27: results/module_score.json (commit 26843e0) had NO
committed generating script; this script reproduces it exactly from the
frozen artifacts, then runs the two stability analyses named in section
9.5: (i) bootstrap distribution of CORE6 KW p and ordinal c-index over
1,000 sample resamples (descriptive, in-sample stability - no new claim);
(ii) leave-one-out member contributions to the CORE6 c-index.
Definition (frozen, section 6b formula 10): S_j = sum_m w_m * z(x'_mj),
w_m = sign(d_severe-mild from gate a), z = per-miRNA standardization over
the 146 graded samples, x' = log2(CPM+1). Members from gate-(b) verdicts."""
import gzip, json, itertools
import numpy as np, pandas as pd
from scipy.stats import kruskal
SEED=20260926
xw=pd.read_csv('sources/GSE299582_sample_crosswalk.csv')
xw['ch']=xw['characteristics'].map(lambda s: json.loads(s.replace('""','"')) if s else {})
xw['sev']=xw['ch'].map(lambda c:c.get('chronic chagas_cardiomyopathy_severity','?'))
xw['prefix']=xw['title'].str.split(' - ').str[0].str.strip()
with gzip.open('sources/matrices/GSE299582_normalized_counts.csv.gz','rt') as f: M=pd.read_csv(f,index_col=0)
od=xw[xw.sev.isin(['mild','moderate','severe'])|(xw.label=='control')].copy()
od['y']=od.apply(lambda r:0 if r['label']=='control' else ['mild','moderate','severe'].index(r['sev'])+1,axis=1)
X=M[od['prefix']]; L=np.log2(X.div(X.sum(axis=0),axis=1)*1e6+1); y=od['y'].values
ga=pd.read_csv('results/h2_gate_a_passing.csv').set_index('mirna')
gb=pd.read_csv('results/h2_gate_b_enrichment.csv')
verd={}
for m,g in gb.groupby('mirna'):
    h=g[g.label.str.contains('hipsc')].iloc[0]; b=g[g.label.str.contains('blood')].iloc[0]
    verd[m]=(h.fdr<=0.05, b.fdr<=0.05)
core6=[m for m,(h,b) in verd.items() if h and b]
card18=[m for m,(h,b) in verd.items() if h]
def score(members):
    hm=['hsa-'+m if not m.startswith('hsa-') else m for m in members]
    w=np.array([np.sign(ga.loc[h,'d_sev_mild']) for h in hm])
    V=L.loc[hm].values
    Z=(V - V.mean(axis=1,keepdims=True)) / V.std(axis=1,keepdims=True)
    return w @ Z
def cindex(yy,s):
    c=t=0
    for i,j in itertools.combinations(range(len(yy)),2):
        if yy[i]==yy[j]: continue
        t+=1; c+=(s[i]-s[j])*(yy[i]-yy[j])>0
    return c/t
def kw(yy,s):
    return kruskal(*[s[yy==k] for k in sorted(set(yy))]).pvalue
def medians(yy,s): return {str(k):float(np.median(s[yy==k])) for k in sorted(set(yy))}
out={}
for name,members in [('CORE6',sorted(core6)),('CARDIAC18',sorted(card18))]:
    s=score(members)
    out[name]={'members':members,'kw_p':float(kw(y,s)),'ordinal_cindex':float(cindex(y,s)),'median_by_group':medians(y,s)}
print('REPRODUCTION:', json.dumps({k:{kk:v[kk] for kk in ('kw_p','ordinal_cindex')} for k,v in out.items()},indent=1))
ref=json.load(open('results/module_score.json'))
for k in out:
    assert abs(out[k]['kw_p']-ref[k]['kw_p'])<1e-15, (k,'kw mismatch')
    assert abs(out[k]['ordinal_cindex']-ref[k]['ordinal_cindex'])<1e-12, (k,'cindex mismatch')
print('EXACT REPRODUCTION of module_score.json confirmed')
# Stability (descriptive)
rng=np.random.default_rng(SEED)
s6=score(sorted(core6)); B=1000; bkw=[];bci=[]
n=len(y)
for _ in range(B):
    idx=rng.integers(0,n,n)
    if len(set(y[idx]))<4: continue
    bkw.append(kw(y[idx],s6[idx])); bci.append(cindex(y[idx],s6[idx]))
loo={}
for m in sorted(core6):
    rest=[x for x in sorted(core6) if x!=m]
    loo[m]=float(cindex(y,score(rest)))
stab={'bootstrap_B':B,'kw_p_median':float(np.median(bkw)),'kw_p_q025':float(np.quantile(bkw,0.025)),
      'kw_p_q975':float(np.quantile(bkw,0.975)),'cindex_median':float(np.median(bci)),
      'cindex_q025':float(np.quantile(bci,0.025)),'cindex_q975':float(np.quantile(bci,0.975)),
      'leave_one_out_cindex':loo,
      'note':'descriptive in-sample stability of the CORE6 module score; members/weights from same cohort, no generalization claim'}
json.dump(stab,open('results/module_score_stability.json','w'),indent=2)
print(json.dumps(stab,indent=2))
