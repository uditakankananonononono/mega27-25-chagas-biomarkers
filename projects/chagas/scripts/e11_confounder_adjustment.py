#!/usr/bin/env python3
"""E11 per ADDENDUM_5: age/sex-augmented ordinal logistic on the locked H1'
folds (same seed/selection; age median-imputed, sex one-hot appended AFTER
scaling of the top-100). BMI/comorbidity NOT in GSE299582 metadata ->
documented, not adjustable. Compare pooled OOF c-index vs unadjusted 0.7871."""
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
xw['age']=pd.to_numeric(xw['ch'].map(lambda c:c.get('age','nan')),errors='coerce')
xw['sex']=xw['ch'].map(lambda c:c.get('gender','?'))
xw['prefix']=xw['title'].str.split(' - ').str[0].str.strip()
with gzip.open('sources/matrices/GSE299582_normalized_counts.csv.gz','rt') as f: M=pd.read_csv(f,index_col=0)
od=xw[xw.sev.isin(['mild','moderate','severe'])|(xw.label=='control')].copy()
od['y']=od.apply(lambda r:0 if r['label']=='control' else ['mild','moderate','severe'].index(r['sev'])+1,axis=1)
X=M[od['prefix']]; L=np.log2(X.div(X.sum(axis=0),axis=1)*1e6+1).T.values; y=od['y'].values.astype(int)
A=np.column_stack([od['age'].fillna(od['age'].median()),(od['sex']=='male').astype(float)])
n_age_missing=int(od['age'].isna().sum())
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
def run(augment):
    oof=np.full(len(y),np.nan)
    for tr,te in StratifiedKFold(5,shuffle=True,random_state=SEED).split(L,y):
        t,pv=ttest_ind(L[tr][y[tr]==0],L[tr][y[tr]!=0],axis=0,equal_var=False)
        keep=np.argsort(np.nan_to_num(pv,nan=1.0))[:100]
        sc=StandardScaler().fit(L[tr][:,keep])
        Xtr,Xte=sc.transform(L[tr][:,keep]),sc.transform(L[te][:,keep])
        if augment:
            sca=StandardScaler().fit(A[tr])
            Xtr=np.hstack([Xtr,sca.transform(A[tr])]); Xte=np.hstack([Xte,sca.transform(A[te])])
        pr=ofit(Xtr,y[tr]); oof[te]=pr(Xte)
    return cindex(y,oof)
unadj=run(False); adj=run(True)
out={'spec':'ADDENDUM_5 E11','unadjusted_oof':float(unadj),'age_sex_adjusted_oof':float(adj),
 'delta':float(adj-unadj),'age_missing_imputed':n_age_missing,
 'not_adjustable':'BMI/comorbidity absent from GSE299582 series metadata (documented)'}
json.dump(out,open('results/e11_confounder_adjustment.json','w'),indent=2)
print(json.dumps(out,indent=1))
