#!/usr/bin/env python3
"""D3 per ADDENDUM_4 (locked 2026-09-27T08:58, NON-GATE sensitivity):
the committed gate-(b) pipeline unchanged, with TargetScan 8.0 human
predicted targets (family mapping via miR_Family_Info, species 9606)
substituted for miRTarBase validated sets, both cohorts, same 10,000
size-preserving without-replacement permutation null and BH family.
Gate-(b) verdicts stand on miRTarBase per ADDENDUM_3 C2; this reports
overlap, divergences and CORE6 membership as-is. TargetScan 8.0 pinned
+ sha256-hashed, ledgered as re-use of service 29."""
import gzip, json
import numpy as np, pandas as pd
from scipy.stats import ttest_ind
SEED=20260926
rng=np.random.default_rng(SEED)

def de_cohort_gse244827():
    with gzip.open('sources/matrices/GSE244827_CHAVArawcounts.txt.gz','rt',errors='replace') as f:
        df=pd.read_csv(f,sep='\t')
    genes=df.iloc[:,0].astype(str).str.upper()
    X=df.iloc[:,6:].apply(pd.to_numeric,errors='coerce').fillna(0).values
    with gzip.open('sources/matrices/GSE244827_CHAVArawcounts.txt.gz','rt') as f:
        mat_cols=f.readline().rstrip('\n').split('\t')[6:]
    lmap=pd.read_csv('sources/GSE244827_column_label_map.csv')
    b2l=dict(zip(lmap.b_code,lmap.label))
    labels=[b2l[c] for c in mat_cols]
    case=[i for i,l in enumerate(labels) if l=='case']; ctrl=[i for i,l in enumerate(labels) if l=='control']
    bmap=pd.read_csv('sources/services/ensembl_biomart/ensg_symbol_map.tsv',sep='\t')
    e2s=dict(zip(bmap['Gene stable ID'],bmap['Gene name'].astype(str).str.upper()))
    genes=genes.map(lambda x: e2s.get(x,x))
    return genes,X,case,ctrl

def de_cohort_gse203525():
    with gzip.open('sources/matrices/GSE203525_Counts.txt.gz','rt',errors='replace') as f:
        df=pd.read_csv(f,sep='\t')
    genes=df.iloc[:,0].astype(str).str.upper()
    X=df.iloc[:,1:].apply(pd.to_numeric,errors='coerce').fillna(0).values
    cols=pd.read_csv('sources/matrices/GSE203525_Counts.txt.gz',sep='\t',nrows=0).columns[1:]
    cc=[i for i,c in enumerate(cols) if c.startswith('CC') and '0hpi' in c]
    ind=[i for i,c in enumerate(cols) if c.startswith('IND') and '0hpi' in c]
    return genes,X,cc,ind

def de_stats(genes,X,case_idx,ctrl_idx):
    lib=X.sum(axis=0); L=np.log2(X/lib*1e6+1)
    t,p=ttest_ind(L[:,case_idx],L[:,ctrl_idx],axis=1,equal_var=False)
    fc=L[:,case_idx].mean(axis=1)-L[:,ctrl_idx].mean(axis=1)
    return genes,fc,np.nan_to_num(p,nan=1.0)

def enrichment(genes,fc,p,targets,direction,label):
    sig=direction*np.sign(fc)>0
    de=sig&(p<=0.05)
    T=set(t.upper() for t in targets)
    gene_idx={g:i for i,g in enumerate(genes)}
    tidx=[gene_idx[g] for g in T if g in gene_idx]
    if len(tidx)<5: return None
    hits=sum(de[i] for i in tidx); k=len(tidx)
    obs=hits/k
    nulls=np.empty(10000)
    pool=np.arange(len(genes))
    for b in range(10000):
        ridx=rng.choice(pool,size=k,replace=False)
        nulls[b]=de[ridx].mean()
    pval=(1+(nulls>=obs).sum())/10001
    return {'label':label,'n_targets_mapped':k,'frac_DE':obs,'null_mean':nulls.mean(),'p':pval}

fam=pd.read_csv('sources/services/targetscan/miR_Family_Info.txt.zip',sep='\t',compression='zip')
fam=fam[fam['Species ID']==9606]
mir2fam=dict(zip(fam['MiRBase ID'],fam['miR family']))
pred=pd.read_csv('sources/services/targetscan/Predicted_Targets_Info.default_predictions.txt.zip',sep='\t',compression='zip')
pred=pred[pred['Species ID']==9606]
fam2targets=pred.groupby('miR Family')['Gene Symbol'].apply(lambda s: sorted(set(s.astype(str).str.upper()))).to_dict()

cand=pd.read_csv('results/h2_gate_a_passing.csv')
cand=cand[cand['class']=='CANDIDATE']
g1,X1,c1,ct1=de_cohort_gse244827(); g1,fc1,p1=de_stats(g1,X1,c1,ct1)
g2,X2,c2,ct2=de_cohort_gse203525(); g2,fc2,p2=de_stats(g2,X2,c2,ct2)
results=[]; unmapped=[]
for _,r in cand.iterrows():
    mir=r['mirna'].replace('hsa-','')
    f_name=mir2fam.get('hsa-'+mir)
    if not f_name or f_name not in fam2targets:
        unmapped.append(mir); continue
    expected=1 if r['d_sev_mild']<0 else -1
    for lbl,gg,ff,pp in [('blood_GSE244827',g1,fc1,p1),('hipsc_GSE203525',g2,fc2,p2)]:
        res=enrichment(gg,ff,pp,fam2targets[f_name],expected,f'{mir}|{lbl}')
        if res: res['mirna']=mir; res['family']=f_name; results.append(res)
R=pd.DataFrame(results)
CORE6=['miR-1-3p','miR-122-5p','miR-192-5p','miR-30c-5p','miR-145-5p','miR-194-5p']
out={'unmapped_candidates':unmapped}
if len(R):
    o=np.argsort(R.p.values); bh=R.p.values[o]*len(R)/(np.arange(len(R))+1)
    R['fdr']=np.minimum.accumulate(bh[::-1])[::-1][np.argsort(o)]
    R.to_csv('results/d3_targetscan_sensitivity.csv',index=False)
    both=[m for m,grp in R.assign(passed=R.fdr<=0.05).groupby('mirna') if int(grp['passed'].sum())==2]
    out.update({'n_tests':len(R),'n_pass_fdr05':int((R.fdr<=0.05).sum()),
        'pass_by_cohort':R.groupby('label')['fdr'].apply(lambda s:int((s<=0.05).sum())).to_dict(),
        'both_cohort_pass':both,
        'core6_both_under_targetscan':[m for m in CORE6 if m in both],
        'sha256':{'Predicted_Targets_Info.default_predictions.txt.zip':'5f981abfc56671fab281b1b1852244472893497424edde50afbe2ac7716aa52c',
                  'miR_Family_Info.txt.zip':'d91827b9c9524f2b970d85f61050ecfa9b8cb7ca395c5f7952b9c21cb3dc3b7a'},
        'reading':'NON-GATE sensitivity (ADDENDUM_4 D3); gate-b stands on miRTarBase validated targets'})
json.dump(out,open('results/d3_targetscan_sensitivity.json','w'),indent=2)
print(json.dumps(out,indent=2))
