#!/usr/bin/env python3
"""E12 per ADDENDUM_5: qPCR-style reduced panel simulation. Panels = top-3
and top-6 markers by E1 discovery frequency (frozen membership - that is
the assay-transfer point). Per locked fold: sign vector estimated on TRAIN
(median severe-vs-control direction), panel score on held-out samples;
pooled OOF c-index vs the full CORE6 module evaluated identically.
Simulation-level only; no wet-lab claim."""
import gzip, json, itertools
import numpy as np, pandas as pd
from sklearn.model_selection import StratifiedKFold
SEED=20260926
xw=pd.read_csv('sources/GSE299582_sample_crosswalk.csv')
xw['ch']=xw['characteristics'].map(lambda s: json.loads(s.replace('""','"')) if s else {})
xw['sev']=xw['ch'].map(lambda c:c.get('chronic chagas_cardiomyopathy_severity','?'))
xw['prefix']=xw['title'].str.split(' - ').str[0].str.strip()
with gzip.open('sources/matrices/GSE299582_normalized_counts.csv.gz','rt') as f: M=pd.read_csv(f,index_col=0)
od=xw[xw.sev.isin(['mild','moderate','severe'])|(xw.label=='control')].copy()
od['y']=od.apply(lambda r:0 if r['label']=='control' else ['mild','moderate','severe'].index(r['sev'])+1,axis=1)
X=M[od['prefix']]; L=np.log2(X.div(X.sum(axis=0),axis=1)*1e6+1).T.values; y=od['y'].values.astype(int)
gi={g:i for i,g in enumerate(M.index)}
def cindex(yy,s):
    c=t=0
    for i,j in itertools.combinations(range(len(yy)),2):
        if yy[i]==yy[j]: continue
        t+=1; c+=(s[i]-s[j])*(yy[i]-yy[j])>0
    return c/t
panels={'top3':['hsa-miR-192-5p','hsa-miR-1-3p','hsa-miR-194-5p'],
        'top6':['hsa-miR-192-5p','hsa-miR-1-3p','hsa-miR-194-5p','hsa-miR-122-5p','hsa-miR-363-3p','hsa-miR-484'],
        'core6_full':['hsa-miR-1-3p','hsa-miR-122-5p','hsa-miR-192-5p','hsa-miR-30c-5p','hsa-miR-145-5p','hsa-miR-194-5p']}
out={'spec':'ADDENDUM_5 E12; frozen membership, train-estimated signs; simulation-level only','panels':{}}
for name,members in panels.items():
    oof=np.full(len(y),np.nan)
    for tr,te in StratifiedKFold(5,shuffle=True,random_state=SEED).split(L,y):
        sc=np.zeros(len(te))
        for m in members:
            j=gi[m]
            sgn=1 if np.median(L[tr][y[tr]!=0],axis=0)[j]-np.median(L[tr][y[tr]==0],axis=0)[j]>0 else -1
            sc+=sgn*L[te,j]
        oof[te]=sc
    out['panels'][name]={'members':members,'oof_cindex':float(cindex(y,oof))}
    print(name, round(out['panels'][name]['oof_cindex'],4),flush=True)
json.dump(out,open('results/e12_qpcr_panel_sim.json','w'),indent=2)
