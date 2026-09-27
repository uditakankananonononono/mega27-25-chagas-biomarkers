#!/usr/bin/env python3
"""E9 per ADDENDUM_5 (locked 2026-09-27T09:23): controlled ML baselines on
IDENTICAL outer folds (seed 20260926) and identical per-fold top-100 feature
selection as h1prime_ordinal.py; only the model swaps. Scores: expected
class (sum k*P(k)) for probabilistic models, the immediate-threshold score
for the ordinal reference. Metric: per-fold mean ordinal c-index, matching
the committed comparator. XGBoost importable (3.2.0) - no substitution
needed. READING locked: H1' advantage stated only over models it beats."""
import gzip, json, itertools, time
import numpy as np, pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier
from sklearn.svm import SVC
from xgboost import XGBClassifier
from scipy.stats import ttest_ind
SEED=20260926
t0=time.time()
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
        m=LogisticRegression(max_iter=20000,C=1.0).fit(Xt,(yt>=k).astype(int)); ms.append(m)
    return lambda Xn: 1+np.column_stack([m.predict_proba(Xn)[:,1] for m in ms]).sum(axis=1)
def expected_class(model,Xn):
    P=model.predict_proba(Xn)
    return P@np.arange(P.shape[1])
def fit_swap(name,Xt,yt):
    if name=='ordinal_logistic(H1-prime)': return lambda Xn: ofit(Xt,yt)(Xn)
    if name=='elastic_net': m=LogisticRegression(penalty='elasticnet',solver='saga',l1_ratio=0.5,C=1.0,max_iter=20000).fit(Xt,yt)
    elif name=='random_forest': m=RandomForestClassifier(n_estimators=500,random_state=SEED,n_jobs=2).fit(Xt,yt)
    elif name=='xgboost': m=XGBClassifier(n_estimators=300,max_depth=3,learning_rate=0.05,subsample=0.8,colsample_bytree=0.8,tree_method='hist',random_state=SEED).fit(Xt,yt)
    elif name=='svm_rbf': m=SVC(kernel='rbf',probability=True,random_state=SEED).fit(Xt,yt)
    return lambda Xn: expected_class(m,Xn)
models=['ordinal_logistic(H1-prime)','elastic_net','random_forest','xgboost','svm_rbf']
scores={m:[] for m in models}
for tr,te in StratifiedKFold(5,shuffle=True,random_state=SEED).split(L,y):
    t,pv=ttest_ind(L[tr][y[tr]==0],L[tr][y[tr]!=0],axis=0,equal_var=False)
    keep=np.argsort(np.nan_to_num(pv,nan=1.0))[:100]
    sc=StandardScaler().fit(L[tr][:,keep])
    Xtr=sc.transform(L[tr][:,keep]); Xte=sc.transform(L[te][:,keep])
    for mname in models:
        pr=fit_swap(mname,Xtr,y[tr])
        scores[mname].append(cindex(y[te],pr(Xte)))
    print('fold done',f'{time.time()-t0:.0f}s',flush=True)
out={m:{'per_fold_mean_cindex':float(np.mean(v)),'per_fold':[round(x,4) for x in v]} for m,v in scores.items()}
out['spec']='ADDENDUM_5 E9; identical folds+features to h1prime_ordinal.py; reading locked'
out['h1prime_committed_reference']=0.787
out['runtime_s']=round(time.time()-t0,1)
json.dump(out,open('results/e9_ml_baselines.json','w'),indent=2)
print(json.dumps(out,indent=2))
