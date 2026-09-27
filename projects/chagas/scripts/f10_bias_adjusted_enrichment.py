#!/usr/bin/env python3
"""F10 per ADDENDUM_7 (her 11:08 directive): popularity-bias-adjusted
target-enrichment test. Per candidate, per cohort (cardiac hipsc_GSE203525,
blood GSE244827): excess = frac_DE(targets) - base rate; matched-null = 500
draws from the candidate's E3 pool (identical construction), candidate's
direction; SURVIVOR = adjusted percentile <0.05 AND raw gate-(b) FDR<=0.05.
Exact hypergeometric (E3 substitution disclosed). Both-way reading locked:
zero survivors is a reportable result. No pool changes post-results."""
import gzip, json
import numpy as np, pandas as pd
from scipy.stats import ttest_ind, hypergeom
SEED=20260926; B=500
rng=np.random.default_rng(SEED+10)
def load_counts(path, meta_cols):
    with gzip.open(path,'rt',errors='replace') as f: df=pd.read_csv(f,sep='\t')
    genes=df.iloc[:,0].astype(str).str.upper()
    X=df.iloc[:,meta_cols:].apply(pd.to_numeric,errors='coerce').fillna(0).values
    return genes,X
def de_stats(genes,X,case_idx,ctrl_idx):
    lib=X.sum(axis=0); L=np.log2(X/lib*1e6+1)
    t,p=ttest_ind(L[:,case_idx],L[:,ctrl_idx],axis=1,equal_var=False)
    fc=L[:,case_idx].mean(axis=1)-L[:,ctrl_idx].mean(axis=1)
    return genes,fc,np.nan_to_num(p,nan=1.0)
COH={}
# hipsc
g2m,X2=load_counts('sources/matrices/GSE203525_Counts.txt.gz',1)
cols=pd.read_csv('sources/matrices/GSE203525_Counts.txt.gz',sep='\t',nrows=0).columns[1:]
cc=[i for i,c in enumerate(cols) if c.startswith('CC') and '0hpi' in c]
ind=[i for i,c in enumerate(cols) if c.startswith('IND') and '0hpi' in c]
COH['hipsc_GSE203525']=de_stats(g2m,X2,cc,ind)
# blood
g,X=load_counts('sources/matrices/GSE244827_CHAVArawcounts.txt.gz',6)
lmap=pd.read_csv('sources/GSE244827_column_label_map.csv'); b2l=dict(zip(lmap.b_code,lmap.label))
with gzip.open('sources/matrices/GSE244827_CHAVArawcounts.txt.gz','rt') as f:
    mat_cols=f.readline().rstrip('\n').split('\t')[6:]
labels=[b2l[c] for c in mat_cols]
case=[i for i,l in enumerate(labels) if l=='case']; ctrl=[i for i,l in enumerate(labels) if l=='control']
bmap=pd.read_csv('sources/services/ensembl_biomart/ensg_symbol_map.tsv',sep='\t')
e2s=dict(zip(bmap['Gene stable ID'],bmap['Gene name'].astype(str).str.upper()))
COH['blood_GSE244827']=de_stats(g.map(lambda x:e2s.get(x,x)),X,case,ctrl)
tg=pd.read_csv('sources/services/mirtarbase/validated_targets.tsv',sep='\t')
tg['bare']=tg.mirna.str.upper().str.replace('HSA-','',regex=False)
mir2tg=tg.groupby('bare').target_gene.apply(lambda s:s.tolist()).to_dict()
cand=pd.read_csv('results/h2_gate_a_passing.csv'); cand=cand[cand['class']=='CANDIDATE']
gb=pd.read_csv('results/h2_gate_b_enrichment.csv')
assoc=pd.read_csv('results/h2_severity_association_all.csv')
assoc['bare']=assoc['mirna'].str.upper().str.replace('HSA-','',regex=False)
pub=json.load(open('sources/services/europepmc/mirna_pubcounts.json'))
tcount=tg.groupby('bare').size()
feat=pd.DataFrame({'tcount':tcount})
feat['expr']=assoc.set_index('bare')[['med_control','med_mild','med_moderate','med_severe']].mean(axis=1)
feat['pub']=pd.Series({k.upper():v for k,v in pub.items()})
feat=feat.dropna()
Z=(feat-feat.mean())/feat.std()
def targets_of(key):
    if key=='MIR-375-3P':
        return tg[tg.bare.isin(['MIR-375-3P','MIR-375'])].target_gene.tolist()
    return [t for m in [k for k in mir2tg if key in k] for t in mir2tg[m]]
def excess(key,targets,direction,gg,ff,pp):
    gene_idx={g:i for i,g in enumerate(gg)}
    T=[gene_idx[t.upper()] for t in targets if t.upper() in gene_idx]
    k=len(T)
    if k<5: return None,None,k
    de=(ff>0)&(pp<=0.05) if direction>0 else (ff<0)&(pp<=0.05)
    K=int(de.sum()); N=len(gg); base=K/N
    x=int(sum(de[i] for i in T))
    return x/k-base, float(hypergeom.sf(x-1,N,K,k)), k
cand_keys=[m.replace('hsa-','').upper() for m in cand['mirna']]
pools={}
for k in cand_keys:
    if k not in Z.index: continue
    d=((Z-Z.loc[k])**2).sum(axis=1)
    pools[k]=d.drop(index=[c for c in cand_keys if c in d.index],errors='ignore').nsmallest(20).index.tolist()
out={'spec':'ADDENDUM_7 F10; B=500 matched draws; exact hypergeometric','cohorts':{}}
for lbl,(gg,ff,pp) in COH.items():
    res={}
    for _,r in cand.iterrows():
        key=r['mirna'].replace('hsa-','').upper()
        direction=1 if r['d_sev_mild']<0 else -1
        tgts=targets_of(key)
        if not tgts or key not in pools: res[key]={'status':'no targets or pool'}; continue
        e_c,p_c,k_c=excess(key,tgts,direction,gg,ff,pp)
        if e_c is None: res[key]={'status':'too few mapped targets'}; continue
        nulls=[]
        for b in range(B):
            m=rng.choice(pools[key])
            mt=targets_of(m)
            if not mt: continue
            e_m,_,_=excess(m,mt,direction,gg,ff,pp)
            if e_m is not None: nulls.append(e_m)
        nulls=np.array(nulls)
        pct=float((nulls>=e_c).mean())
        raw=gb[(gb.mirna==key.lower().replace('mir','miR'))&(gb.label==f'{key.capitalize()}|{lbl}')]
        res[key]={'excess':float(e_c),'k_targets':k_c,'exact_p':p_c,
                  'null_median_excess':float(np.median(nulls)),'null_q95_excess':float(np.percentile(nulls,95)),
                  'adjusted_p':pct,'survives':bool(pct<0.05)}
        rawfdr=gb[(gb.label==f"{r['mirna'].replace('hsa-','')}|{lbl}")]
        res[key]['raw_gateb_fdr']=float(rawfdr.fdr.iloc[0]) if len(rawfdr) else None
        res[key]['survivor']=bool(pct<0.05 and len(rawfdr) and rawfdr.fdr.iloc[0]<=0.05)
    out['cohorts'][lbl]=res
    surv=[k for k,v in res.items() if v.get('survivor')]
    print(f'{lbl}: survivors {len(surv)}: {surv}',flush=True)
json.dump(out,open('results/f10_bias_adjusted_enrichment.json','w'),indent=2)
