#!/usr/bin/env python3
"""F5 per ADDENDUM_6: standard-pipeline comparison on the SAME severity screen
data (serum GSE299582, 4 groups, log2 CPM, n=146) as the frozen gate-(a)
screen (KW + BH FDR<=0.05 + |median severe-mild|>=0.5 -> 21 passing).
Alternatives: (a) plain DE (KW FDR only, no effect gate); (b) PCA clustering
(descriptive - no feature candidate set exists); (c) RF feature selection
(top-21 by importance); (d) limma-style linear-model analog (per-feature OLS
on ordinal severity, BH FDR<=0.05 - limma/R unavailable in this environment,
substitution disclosed). Overlap vs the frozen 21 reported verbatim either
way (locked READING: superiority claimed only where plain alternatives do
NOT recover the frozen candidates)."""
import gzip, json
import numpy as np, pandas as pd
from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
xw=pd.read_csv('sources/GSE299582_sample_crosswalk.csv')
xw['ch']=xw['characteristics'].map(lambda s: json.loads(s.replace('""','"')) if s else {})
xw['sev']=xw['ch'].map(lambda c:c.get('chronic chagas_cardiomyopathy_severity','?'))
xw['prefix']=xw['title'].str.split(' - ').str[0].str.strip()
with gzip.open('sources/matrices/GSE299582_normalized_counts.csv.gz','rt') as f: M=pd.read_csv(f,index_col=0)
od=xw[xw.sev.isin(['mild','moderate','severe'])|(xw.label=='control')].copy()
od['grp']=od.apply(lambda r:'control' if r['label']=='control' else r['sev'],axis=1)
order=['control','mild','moderate','severe']
X=M[od['prefix']]; L=np.log2(X.div(X.sum(axis=0),axis=1)*1e6+1)
groups=[L.loc[:,(od.grp==g).values] for g in order]
frozen=set(pd.read_csv('results/h2_gate_a_passing.csv')['mirna'])
assert len(frozen)==21, len(frozen)
genes=np.array(M.index)
# frozen screen statistics (recomputed for the alternatives' p-values)
def kw_safe(i):
    try: return kruskal(*[g.iloc[i].values for g in groups]).pvalue
    except Exception: return 1.0
pvals=np.nan_to_num(np.array([kw_safe(i) for i in range(L.shape[0])]),nan=1.0)
fdr=multipletests(pvals,method='fdr_bh')[1]
gm={g:groups[i].median(axis=1).values for i,g in enumerate(order)}
d_eff=gm['severe']-gm['mild']
# (a) plain DE: FDR only
plain=set(genes[fdr<=0.05])
# (c) RF feature selection, top-21 by importance
yord=od['grp'].map({g:i for i,g in enumerate(order)}).values
rf=RandomForestClassifier(n_estimators=500,random_state=20260926,n_jobs=4).fit(L.T.values,yord)
top21rf=set(genes[np.argsort(rf.feature_importances_)[-21:]])
# (d) limma-style analog: per-feature OLS on ordinal severity
yv=yord.astype(float); yc=yv-yv.mean()
Xc=L.values-L.values.mean(axis=1,keepdims=True)
beta=(Xc@yc)/(yc@yc)
resid=Xc-np.outer(beta,yc); s2=(resid**2).sum(axis=1)/(len(yv)-2)
se=np.sqrt(s2/(yc@yc)); tstat=beta/np.sqrt(se+1e-300)
from scipy.stats import t as tdist
p_lm=2*tdist.sf(np.abs(tstat),df=len(yv)-2)
fdr_lm=multipletests(p_lm,method='fdr_bh')[1]
lm_set=set(genes[fdr_lm<=0.05])
# (b) PCA clustering, descriptive
pc=PCA(n_components=4,random_state=20260926).fit(L.T.values)
sc=pc.transform(L.T.values)
kw_pc=[kruskal(*[sc[(od.grp==g).values,k] for g in order]).pvalue for k in range(4)]
out={'spec':'ADDENDUM_6 F5; overlap verbatim either way','frozen_n':21,
 'plain_de':{'n_fdr_only':len(plain),'overlap_with_frozen21':len(plain&frozen),
   'frozen_missed':sorted(frozen-plain)},
 'rf_top21':{'overlap_with_frozen21':len(top21rf&frozen),
   'recovered':sorted(top21rf&frozen),'frozen_missed':sorted(frozen-top21rf)},
 'limma_style_lm':{'n_fdr':len(lm_set),'overlap_with_frozen21':len(lm_set&frozen),
   'frozen_missed':sorted(frozen-lm_set),
   'substitution':'limma (R) unavailable; per-feature OLS on ordinal severity + BH FDR'},
 'pca':{'var_explained_first4':[float(v) for v in pc.explained_variance_ratio_],
   'kw_p_per_pc':[float(p) for p in kw_pc],
   'note':'descriptive; PCA yields no feature candidate set'},
 'frozen21':sorted(frozen)}
json.dump(out,open('results/f5_pipeline_comparison.json','w'),indent=2)
print(json.dumps({k:v for k,v in out.items() if k not in('frozen21','pca')},indent=1))
print('pca var:',out['pca']['var_explained_first4']); print('pca kw p:',out['pca']['kw_p_per_pc'])
