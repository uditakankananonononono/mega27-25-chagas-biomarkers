#!/usr/bin/env python3
"""H1' LOCKED-METRIC COMPLETION (2026-09-27). Addendum 3 C1 locked FOUR
metrics: ordinal c-index (primary), macro-AUC one-vs-rest, calibration
slope, bootstrap 95% CIs. The H1' run reported c-index + bootstrap CIs;
macro-AUC and calibration slope were locked but never reported (gap found
in paper-audit sweep 2026-09-27 05:34 IST). This script re-runs the SAME
locked pipeline (same matrix, split, seed, feature rule, model family) and
additionally persists OOF predictions and computes the two missing locked
metrics. Design unchanged; this is report completion, not a new analysis.
Sanity anchor: mean-of-fold model c-index must reproduce
results/h1prime_ordinal.json (0.787)."""
import gzip, json, itertools
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
X=M[od['prefix']]; L=np.log2(X.div(X.sum(axis=0),axis=1)*1e6+1).T.values; y=od['y'].values
def cindex(yy,s):
    c=t=0
    for i,j in itertools.combinations(range(len(yy)),2):
        if yy[i]==yy[j]: continue
        t+=1; c+=(s[i]-s[j])*(yy[i]-yy[j])>0
    return c/t
def ofit(Xt,yt):
    K=sorted(set(yt)); ms=[]
    for k in K[1:]:
        ms.append(LogisticRegression(max_iter=20000,C=1.0).fit(Xt,(yt>=k).astype(int)))
    return lambda Xn: np.column_stack([m.predict_proba(Xn)[:,1] for m in ms])  # per-threshold P(y>=k)
oof_s=np.zeros(len(y)); oof_p=np.zeros((len(y),3)); fold_c=[]
for tr,te in StratifiedKFold(5,shuffle=True,random_state=SEED).split(L,y):
    t,pv=ttest_ind(L[tr][y[tr]==0],L[tr][y[tr]!=0],axis=0,equal_var=False)
    keep=np.argsort(np.nan_to_num(pv,nan=1.0))[:100]
    sc=StandardScaler().fit(L[tr][:,keep]); pr=ofit(sc.transform(L[tr][:,keep]),y[tr])
    P=pr(sc.transform(L[te][:,keep])); oof_p[te]=P; oof_s[te]=1+P.sum(axis=1)
    fold_c.append(cindex(y[te],oof_s[te]))
mean_fold=float(np.mean(fold_c)); pooled=cindex(y,oof_s)
# macro-AUC one-vs-rest from OOF CLASS PROBABILITIES (per-threshold probs
# differenced into P(y=k)); a single monotone severity score is NOT a valid
# one-vs-rest statistic (control-class AUC degenerates), so class probs used.
from scipy.stats import rankdata
cls_p=np.column_stack([1-oof_p[:,0], oof_p[:,0]-oof_p[:,1], oof_p[:,1]-oof_p[:,2], oof_p[:,2]])
def auc_ovr(yy,P):
    aucs=[]
    for k in sorted(set(yy)):
        s=P[:,k]; pos=s[yy==k]; neg=s[yy!=k]; r=rankdata(np.concatenate([pos,neg]))
        aucs.append((r[:len(pos)].sum()-len(pos)*(len(pos)+1)/2)/(len(pos)*len(neg)))
    return float(np.mean(aucs)), aucs
macro, per_class = auc_ovr(y,cls_p)
# calibration slope per threshold: logit recalibration of P(y>=k) on OOF predicted prob
cal={}
for i,k in enumerate([1,2,3]):
    pk=np.clip(oof_p[:,i],1e-6,1-1e-6); z=np.log(pk/(1-pk))
    lr=LogisticRegression(max_iter=20000).fit(z.reshape(-1,1),(y>=k).astype(int))
    cal[f'y>={k}']={'slope':float(lr.coef_[0][0]),'intercept':float(lr.intercept_[0])}
out={'sanity_mean_of_fold_cindex':mean_fold,'pooled_oof_cindex':pooled,
     'macro_auc_ovr':macro,'per_class_auc':{str(k):float(a) for k,a in zip(sorted(set(y)),per_class)},
     'calibration_slope_per_threshold':cal,
     'note':'locked metrics macro-AUC + calibration slope completed 2026-09-27; design unchanged (ADDENDUM_3 C1); sklearn version unpinned in lane record (env note: sklearn 1.7.2, rest of env matches 5.y)'}
pd.DataFrame({'prefix':od['prefix'].values,'y':y,'oof_score':oof_s,
              'oof_p_ge1':oof_p[:,0],'oof_p_ge2':oof_p[:,1],'oof_p_ge3':oof_p[:,2]}
             ).to_csv('results/h1prime_oof_predictions.csv',index=False)
json.dump(out,open('results/h1prime_locked_metrics_completion.json','w'),indent=2)
print(json.dumps(out,indent=2))
