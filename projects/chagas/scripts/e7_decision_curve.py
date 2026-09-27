#!/usr/bin/env python3
"""E7 per ADDENDUM_5 + round-03 weakness-13: decision-curve analysis on the
locked OOF predictions. Net benefit of treat-by-model vs treat-all vs
treat-none across threshold probabilities, for the severe-vs-rest contrast
(y>=3) and moderate-plus (y>=2), comparing: full H1' model, age/sex-only
model (recomputed OOF on locked folds), and combined H1'+age/sex (E11 OOF).
Descriptive; no clinical-use claim."""
import gzip, json, itertools
import numpy as np, pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from scipy.stats import ttest_ind
SEED=20260926
df=pd.read_csv('results/h1prime_oof_predictions.csv'); y=df.y.values
# age/sex-only OOF on the same folds
xw=pd.read_csv('sources/GSE299582_sample_crosswalk.csv')
xw['ch']=xw['characteristics'].map(lambda s: json.loads(s.replace('""','"')) if s else {})
xw['sev']=xw['ch'].map(lambda c:c.get('chronic chagas_cardiomyopathy_severity','?'))
xw['age']=pd.to_numeric(xw['ch'].map(lambda c:c.get('age','nan')),errors='coerce')
xw['sex']=xw['ch'].map(lambda c:c.get('gender','?'))
od=xw[xw.sev.isin(['mild','moderate','severe'])|(xw.label=='control')].copy()
od['y']=od.apply(lambda r:0 if r['label']=='control' else ['mild','moderate','severe'].index(r['sev'])+1,axis=1)
A=np.column_stack([od['age'].fillna(od['age'].median()),(od['sex']=='male').astype(float)])
def ofit(Xt,yt):
    K=sorted(set(yt)); ms=[]
    for k in K[1:]:
        m=LogisticRegression(max_iter=20000,C=1.0).fit(Xt,(yt>=k).astype(int)); ms.append(m)
    return lambda Xn: 1+np.column_stack([m.predict_proba(Xn)[:,1] for m in ms]).sum(axis=1)
oof_as=np.full(len(y),np.nan)
for tr,te in StratifiedKFold(5,shuffle=True,random_state=SEED).split(A,y):
    sc=StandardScaler().fit(A[tr]); pr=ofit(sc.transform(A[tr]),y[tr]); oof_as[te]=pr(sc.transform(A[te]))
def nb(score,yy,pt):
    pred=score>=pt; tp=(pred&(yy==1)).sum(); fp=(pred&(yy==0)).sum(); n=len(yy)
    return tp/n-fp/n*(pt/(1-pt))
pts=np.linspace(0.01,0.5,50)
out={'spec':'ADDENDUM_5 E7; descriptive, no clinical-use claim','contrasts':{}}
for lbl,yy,col in [('severe_vs_rest',(y>=3).astype(int),'oof_p_ge3'),('moderate_plus',(y>=2).astype(int),'oof_p_ge2')]:
    res={}
    for name,sc in [('full_model',df[col].values),('agesex_only',oof_as),('combined',df[col].values+oof_as)]:
        res[name]=[float(nb(sc,yy,pt)) for pt in pts]
    res['treat_all']=[float(yy.mean()-(1-yy.mean())*(pt/(1-pt))) for pt in pts]
    out['contrasts'][lbl]={'thresholds':[float(p) for p in pts],**res,
        'full_beats_agesex_frac':float(np.mean(np.array(res['full_model'])>np.array(res['agesex_only']))),
        'full_beats_all_frac':float(np.mean(np.array(res['full_model'])>np.array(res['treat_all'])))}
json.dump(out,open('results/e7_decision_curve.json','w'),indent=2)
print(json.dumps({c:{'full_beats_agesex_frac':v['full_beats_agesex_frac'],'full_beats_all_frac':v['full_beats_all_frac']} for c,v in out['contrasts'].items()},indent=1))
