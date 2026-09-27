#!/usr/bin/env python3
"""F3 per ADDENDUM_6: permutation importance of the top-100 features within
each outer fold of the H1' pipeline (1,000 permutations per feature, OOF
c-index drop on that fold's held-out samples). Descriptive; no claim change.
Reported alongside E1 selection frequencies (joined by gene where present)."""
import gzip, json
import numpy as np, pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from scipy.stats import ttest_ind
SEED=20260926; NPERM=1000
xw=pd.read_csv('sources/GSE299582_sample_crosswalk.csv')
xw['ch']=xw['characteristics'].map(lambda s: json.loads(s.replace('""','"')) if s else {})
xw['sev']=xw['ch'].map(lambda c:c.get('chronic chagas_cardiomyopathy_severity','?'))
xw['prefix']=xw['title'].str.split(' - ').str[0].str.strip()
with gzip.open('sources/matrices/GSE299582_normalized_counts.csv.gz','rt') as f: M=pd.read_csv(f,index_col=0)
od=xw[xw.sev.isin(['mild','moderate','severe'])|(xw.label=='control')].copy()
od['y']=od.apply(lambda r:0 if r['label']=='control' else ['mild','moderate','severe'].index(r['sev'])+1,axis=1)
X=M[od['prefix']]; genes=np.array(M.index)
L=np.log2(X.div(X.sum(axis=0),axis=1)*1e6+1).T.values; y=od['y'].values.astype(int)
def cindex_fast(yy,s):
    yy=np.asarray(yy); s=np.asarray(s); c=t=0
    n=len(yy)
    for i in range(n):
        d=yy[i]-yy[i+1:]; m=d!=0
        t+=m.sum(); c+=(((s[i]-s[i+1:][m])*d[m])>0).sum()
    return c/t if t else np.nan
def ofit(Xt,yt):
    K=sorted(set(yt)); ms=[]
    for k in K[1:]:
        m=LogisticRegression(max_iter=20000,C=1.0).fit(Xt,(yt>=k).astype(int)); ms.append(m)
    return lambda Xn: 1+np.column_stack([m.predict_proba(Xn)[:,1] for m in ms]).sum(axis=1)
rng=np.random.default_rng(SEED+2)
imp={}   # gene -> list of mean c-index drops across folds
sel={}   # gene -> folds selected
for fi,(tr,te) in enumerate(StratifiedKFold(5,shuffle=True,random_state=SEED).split(L,y)):
    t,pv=ttest_ind(L[tr][y[tr]==0],L[tr][y[tr]!=0],axis=0,equal_var=False)
    keep=np.argsort(np.nan_to_num(pv,nan=1.0))[:100]
    sc=StandardScaler().fit(L[tr][:,keep]); pr=ofit(sc.transform(L[tr][:,keep]),y[tr])
    Xte=sc.transform(L[te][:,keep]); base=cindex_fast(y[te],pr(Xte))
    for pos,g in enumerate(keep):
        gn=genes[g]; sel.setdefault(gn,[]).append(fi)
        drops=[]
        col=Xte[:,pos].copy()
        for _ in range(NPERM):
            Xp=Xte.copy(); Xp[:,pos]=col[rng.permutation(len(te))]
            drops.append(base-cindex_fast(y[te],pr(Xp)))
        imp.setdefault(gn,[]).append(float(np.mean(drops)))
    print(f'fold {fi} done (base={base:.3f})',flush=True)
rows=[{'gene':g,'folds_selected':len(sel[g]),'mean_oof_cindex_drop':float(np.mean(imp[g])),
       'max_fold_drop':float(np.max(imp[g]))} for g in imp]
df=pd.DataFrame(rows).sort_values('mean_oof_cindex_drop',ascending=False)
df.to_csv('results/f3_permutation_importance.csv',index=False)
e1=json.load(open('results/e1_nested_bootstrap_discovery.json'))
freq=e1.get('top20_discovery_frequency',{})
df['e1_discovery_frequency']=df['gene'].map(lambda g:freq.get(g,None))
df.to_csv('results/f3_permutation_importance.csv',index=False)
out={'spec':'ADDENDUM_6 F3; descriptive, no claim change','nperm_per_feature':NPERM,
     'n_features':len(df),'seed':SEED,
     'top10_by_mean_drop':df.head(10).to_dict('records'),
     'negative_importance_count':int((df.mean_oof_cindex_drop<0).sum())}
json.dump(out,open('results/f3_permutation_importance.json','w'),indent=2)
print(json.dumps(out['top10_by_mean_drop'],indent=1)); print('neg:',out['negative_importance_count'])
