#!/usr/bin/env python3
"""E5 per ADDENDUM_5 (+round-03 weakness-15 jackknife): influence diagnostics
on the locked H1' pipeline (severe n=30, full n=146).
(a) Cook's D per sample on the locked OOF logistic fits (max over folds x
    thresholds), 4/n_train rule, flagged samples disclosed.
(b) leave-5%-out (7 samples, 200 seeded draws) OOF c-index stability.
(c) leave-severe-out: refit on 116 (control/mild/moderate), OOF c-index +
    KW of the CORE6 module score across the 3 remaining groups.
    READING (locked): gradient claim holds only if leave-severe-out KW<0.05
    AND no Cook's D >4/n without disclosure; else downgrade.
(d) LOO jackknife: full-pipeline OOF c-index per dropped sample (146 runs)."""
import gzip, json, itertools
import numpy as np, pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from scipy.stats import ttest_ind, kruskal
SEED=20260926
xw=pd.read_csv('sources/GSE299582_sample_crosswalk.csv')
xw['ch']=xw['characteristics'].map(lambda s: json.loads(s.replace('""','"')) if s else {})
xw['sev']=xw['ch'].map(lambda c:c.get('chronic chagas_cardiomyopathy_severity','?'))
xw['prefix']=xw['title'].str.split(' - ').str[0].str.strip()
with gzip.open('sources/matrices/GSE299582_normalized_counts.csv.gz','rt') as f: M=pd.read_csv(f,index_col=0)
od=xw[xw.sev.isin(['mild','moderate','severe'])|(xw.label=='control')].copy()
od['y']=od.apply(lambda r:0 if r['label']=='control' else ['mild','moderate','severe'].index(r['sev'])+1,axis=1)
X=M[od['prefix']]; L=np.log2(X.div(X.sum(axis=0),axis=1)*1e6+1).T.values; y=od['y'].values.astype(int)
n=len(y)
def cindex(yy,s):
    c=t=0
    for i,j in itertools.combinations(range(len(yy)),2):
        if yy[i]==yy[j]: continue
        t+=1; c+=(s[i]-s[j])*(yy[i]-yy[j])>0
    return c/t if t else np.nan
def ofit(Xt,yt):
    K=sorted(set(yt)); ms=[]
    for k in K[1:]:
        m=LogisticRegression(max_iter=20000,C=1.0).fit(Xt,(yt>=k).astype(int)); ms.append(m)
    return lambda Xn: 1+np.column_stack([m.predict_proba(Xn)[:,1] for m in ms]).sum(axis=1)
def pipeline(mask):
    """Full H1' pipeline on samples mask; returns pooled OOF c-index and per-sample oof scores."""
    idx=np.where(mask)[0]; yy=y[idx]; Ll=L[idx]
    oof=np.full(len(idx),np.nan)
    for tr,te in StratifiedKFold(5,shuffle=True,random_state=SEED).split(Ll,yy):
        t,pv=ttest_ind(Ll[tr][yy[tr]==0],Ll[tr][yy[tr]!=0],axis=0,equal_var=False)
        keep=np.argsort(np.nan_to_num(pv,nan=1.0))[:100]
        sc=StandardScaler().fit(Ll[tr][:,keep]); pr=ofit(sc.transform(Ll[tr][:,keep]),yy[tr])
        oof[te]=pr(sc.transform(Ll[te][:,keep]))
    return cindex(yy,oof),oof,idx
# module score for KW (CORE6 signed, per module_score.py convention)
core6=['hsa-miR-1-3p','hsa-miR-122-5p','hsa-miR-192-5p','hsa-miR-30c-5p','hsa-miR-145-5p','hsa-miR-194-5p']
ga=pd.read_csv('results/h2_gate_a_passing.csv').set_index('mirna')
signs={m:(1 if ga.loc[m,'d_sev_mild']>0 else -1) for m in core6}
gi={g:i for i,g in enumerate(M.index)}
modscore=np.zeros(n)
for m,s in signs.items(): modscore+=s*L[:,gi[m]]
# (a) Cook's D on locked fits (max over folds/thresholds)
cooks=np.zeros(n)
for tr,te in StratifiedKFold(5,shuffle=True,random_state=SEED).split(L,y):
    t,pv=ttest_ind(L[tr][y[tr]==0],L[tr][y[tr]!=0],axis=0,equal_var=False)
    keep=np.argsort(np.nan_to_num(pv,nan=1.0))[:100]
    sc=StandardScaler().fit(L[tr][:,keep]); Xt=sc.transform(L[tr][:,keep])
    Xd=np.column_stack([np.ones(len(Xt)),Xt]); pn=Xd.shape[1]
    H=Xd@np.linalg.pinv(Xd.T@Xd)@Xd.T; h=np.diag(H)
    for k in sorted(set(y[tr]))[1:]:
        yy=(y[tr]>=k).astype(int)
        m=LogisticRegression(max_iter=20000,C=1.0).fit(Xt,yy)
        pr=m.predict_proba(Xt)[:,1]; r=(yy-pr)/np.sqrt(pr*(1-pr)+1e-12)
        d=(r**2*h)/(pn*(1-h)**2+1e-12)
        cooks[tr]=np.maximum(cooks[tr],d)
thr=4/(n*0.8)
flagged=[(od['prefix'].iloc[i],int(y[i]),float(cooks[i])) for i in np.where(cooks>thr)[0]]
# (b) leave-5%-out stability
rng=np.random.default_rng(SEED+5); l5=[]
for d in range(200):
    drop=rng.choice(n,7,replace=False); mask=np.ones(n,bool); mask[drop]=False
    c,_,_=pipeline(mask); l5.append(c)
    if (d+1)%50==0: print(f'l5 {d+1}/200',flush=True)
# (c) leave-severe-out
mask=y!=3
c_lso,oof_lso,idx_lso=pipeline(mask)
kw_lso=kruskal(*[modscore[mask&(y==g)] for g in (0,1,2)]).pvalue
c_lsom,_c,_=pipeline(mask)  # same; c-index on 3-class remainder
# (d) LOO jackknife
loo=[]
for i in range(n):
    mask=np.ones(n,bool); mask[i]=False
    c,_,_=pipeline(mask); loo.append(c)
    if (i+1)%25==0: print(f'loo {i+1}/{n}',flush=True)
out={'spec':'ADDENDUM_5 E5 + round-03 jackknife','n':n,
 'cooks':{'threshold_4_over_ntrain':float(thr),'n_flagged':len(flagged),'flagged':flagged,
   'max':float(cooks.max())},
 'leave5pct':{'draws':200,'median':float(np.median(l5)),'q025':float(np.percentile(l5,2.5)),
   'q975':float(np.percentile(l5,97.5)),'min':float(np.min(l5))},
 'leave_severe_out':{'oof_cindex_3class':float(c_lso),'module_KW_p_no_severe':float(kw_lso),
   'gradient_survives':bool(kw_lso<0.05)},
 'loo_jackknife':{'median':float(np.median(loo)),'min':float(np.min(loo)),'max':float(np.max(loo)),
   'q025':float(np.percentile(loo,2.5)),'q975':float(np.percentile(loo,97.5))}}
json.dump(out,open('results/e5_influence_diagnostics.json','w'),indent=2)
print(json.dumps(out,indent=1)[:1500])
