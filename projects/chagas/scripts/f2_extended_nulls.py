#!/usr/bin/env python3
"""F2 per ADDENDUM_6: extended null framework on the H1' pipeline.
Same folds/seed/preprocessing as h1prime_ordinal.py (5-fold stratified CV,
seed 20260926, top-100 Welch t-test selection inside train, ordinal logistic).
Nulls, 1000 draws each: (a) label shuffle; (c) random 100-feature sets.
(b) identity shuffle is degenerate here (pipeline is annotation-free; column
permutation leaves p-value-rank selection invariant) - documented, 100 draws
run for the record. Statistic: pooled OOF c-index (matches committed
h1prime_oof_predictions.csv value). Empirical p = (1+#(null>=obs))/(1+N).
READING (locked): observed must exceed the 95th pct of each null to keep the
'not random' sentence; otherwise it is cut."""
import gzip, json, itertools, sys, time
import numpy as np, pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from scipy.stats import ttest_ind
SEED=20260926
xw=pd.read_csv('sources/GSE299582_sample_crosswalk.csv')
xw['ch']=xw['characteristics'].map(lambda s: json.loads(s.replace('""','"')) if s else {})
xw['sev']=xw['ch'].map(lambda c:c.get('chronic chagas_cardiomyopathy_severity','?'))
xw['prefix']=xw['title'].str.split(' - ').str[0].str.strip()
with gzip.open('sources/matrices/GSE299582_normalized_counts.csv.gz','rt') as f: M=pd.read_csv(f,index_col=0)
od=xw[xw.sev.isin(['mild','moderate','severe'])|(xw.label=='control')].copy()
od['y']=od.apply(lambda r:0 if r['label']=='control' else ['mild','moderate','severe'].index(r['sev'])+1,axis=1)
X=M[od['prefix']]; L=np.log2(X.div(X.sum(axis=0),axis=1)*1e6+1).T.values; y=od['y'].values.astype(int)
def cindex(yy,s):
    c=t=0
    for i,j in itertools.combinations(range(len(yy)),2):
        if yy[i]==yy[j]: continue
        t+=1; c+=(s[i]-s[j])*(yy[i]-yy[j])>0
    return c/t
def ofit(Xt,yt):
    K=sorted(set(yt)); ms=[]
    for k in K[1:]:
        m=LogisticRegression(max_iter=20000,C=1.0).fit(Xt,(yt>=k).astype(int)); ms.append(m)
    return lambda Xn: 1+np.column_stack([m.predict_proba(Xn)[:,1] for m in ms]).sum(axis=1)
FOLDS=list(StratifiedKFold(5,shuffle=True,random_state=SEED).split(L,y))
def run_pipeline(yy, feat_set=None):
    """feat_set None -> per-fold top-100 t-test selection (the H1' pipeline);
    else fixed 100-column index array used in every fold. Returns pooled OOF c-index."""
    oof=np.full(len(yy),np.nan)
    for tr,te in FOLDS:
        if feat_set is None:
            t,pv=ttest_ind(L[tr][yy[tr]==0],L[tr][yy[tr]!=0],axis=0,equal_var=False)
            keep=np.argsort(np.nan_to_num(pv,nan=1.0))[:100]
        else: keep=feat_set
        sc=StandardScaler().fit(L[tr][:,keep]); pr=ofit(sc.transform(L[tr][:,keep]),yy[tr])
        oof[te]=pr(sc.transform(L[te][:,keep]))
    return cindex(yy,oof)
mode=sys.argv[1]; ndraws=int(sys.argv[2]); t0=time.time()
obs=run_pipeline(y)
print(f'observed pooled OOF c-index={obs:.4f} ({time.time()-t0:.1f}s)',flush=True)
rng=np.random.default_rng(SEED+1); nulls=[]
if mode=='label':
    for d in range(ndraws):
        nulls.append(run_pipeline(rng.permutation(y)))
        if (d+1)%50==0: print(f'{d+1}/{ndraws} ({time.time()-t0:.0f}s)',flush=True)
elif mode=='randfeat':
    G=L.shape[1]
    for d in range(ndraws):
        nulls.append(run_pipeline(y,feat_set=rng.choice(G,100,replace=False)))
        if (d+1)%50==0: print(f'{d+1}/{ndraws} ({time.time()-t0:.0f}s)',flush=True)
elif mode=='identity':
    G=L.shape[1]
    for d in range(ndraws):
        perm=rng.permutation(G)
        L[:]=L[:,perm]
        nulls.append(run_pipeline(y))
        if (d+1)%50==0: print(f'{d+1}/{ndraws} ({time.time()-t0:.0f}s)',flush=True)
else: sys.exit('unknown mode')
nulls=np.array(nulls)
p=(1+int((nulls>=obs).sum()))/(1+len(nulls))
out={'mode':mode,'draws':len(nulls),'observed_oof_cindex':float(obs),
     'null_mean':float(nulls.mean()),'null_sd':float(nulls.std()),
     'null_p95':float(np.percentile(nulls,95)),'null_max':float(nulls.max()),
     'empirical_p':p,'exceeds_p95':bool(obs>np.percentile(nulls,95)),
     'pipeline':'h1prime_ordinal.py identical folds/seed/selection','seed':SEED}
if mode=='identity':
    out['near_degenerate_note']='pipeline is annotation-free: selection is by p-value rank, so column permutation changes the selected top-100 only via floating-point tie order; null expected to concentrate at the observed value'
json.dump(out,open(f'results/f2_null_{mode}.json','w'),indent=2)
np.savetxt(f'results/f2_null_{mode}_draws.csv',nulls,fmt='%.6f')
print(json.dumps({k:out[k] for k in('mode','draws','observed_oof_cindex','null_mean','null_p95','null_max','empirical_p','exceeds_p95')}))
