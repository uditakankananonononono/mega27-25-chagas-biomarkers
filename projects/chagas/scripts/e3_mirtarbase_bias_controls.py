#!/usr/bin/env python3
"""E3 per ADDENDUM_5: miRTarBase-bias negative controls for the cardiac
gate-(b) claim (18/20 candidates enriched in hipsc_GSE203525 at FDR<=0.05).
Control A: matched random miRNAs (validated-target count, serum expression,
Europe PMC pubcount; z-score 20-nearest pools, excluding the 20 candidates),
full gate-(b) pipeline per draw, BH FDR across the 20, 200 draws.
Control B: size-preserving random validated targets from the SAME TargetScan
miRNA families (human 9606, pinned miR_Family_Info.txt), 200 draws/candidate.
Control C: TargetScan arm = D3 (done, results/d3_targetscan_sensitivity.json).
Computational substitution, disclosed: the size-preserving random-gene-set
null is computed EXACTLY (hypergeometric) instead of 10,000 Monte Carlo
draws - same null distribution. Observed 20-candidate cardiac p-values are
recomputed under the exact test in this script for apples-to-apples.
READING (locked): cardiac claim holds only if it exceeds all three controls
at the locked FDR; failure -> 'database-bias-consistent', disclosed."""
import gzip, json
import numpy as np, pandas as pd
from scipy.stats import ttest_ind, hypergeom
SEED=20260926; DRAWS=200
rng=np.random.default_rng(SEED+3)
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
# hipsc cohort DE (identical to gate-b script)
g2m,X2=load_counts('sources/matrices/GSE203525_Counts.txt.gz',1)
cols=pd.read_csv('sources/matrices/GSE203525_Counts.txt.gz',sep='\t',nrows=0).columns[1:]
cc=[i for i,c in enumerate(cols) if c.startswith('CC') and '0hpi' in c]
ind=[i for i,c in enumerate(cols) if c.startswith('IND') and '0hpi' in c]
g2,fc2,p2=de_stats(g2m,X2,cc,ind)
N=len(g2); gene_idx={g:i for i,g in enumerate(g2)}
de_up=(fc2>0)&(p2<=0.05); de_dn=(fc2<0)&(p2<=0.05)
K_up,K_dn=int(de_up.sum()),int(de_dn.sum())
def exact_p(targets,direction):
    T=[gene_idx[t.upper()] for t in targets if t.upper() in gene_idx]
    k=len(T)
    if k<5: return None,k,None
    de=de_up if direction>0 else de_dn; K=K_up if direction>0 else K_dn
    x=int(sum(de[i] for i in T))
    return float(hypergeom.sf(x-1,N,K,k)),k,x
tg=pd.read_csv('sources/services/mirtarbase/validated_targets.tsv',sep='\t')
tg['bare']=tg.mirna.str.upper().str.replace('HSA-','',regex=False)
cand=pd.read_csv('results/h2_gate_a_passing.csv'); cand=cand[cand['class']=='CANDIDATE']
assoc=pd.read_csv('results/h2_severity_association_all.csv')
assoc['bare']=assoc['mirna'].str.upper().str.replace('HSA-','',regex=False)
pub=json.load(open('sources/services/europepmc/mirna_pubcounts.json'))
tcount=tg.groupby('bare').size()
feat=pd.DataFrame({'tcount':tcount})
feat['expr']=assoc.set_index('bare')[['med_control','med_mild','med_moderate','med_severe']].mean(axis=1)
feat['pub']=pd.Series({k.upper():v for k,v in pub.items()})
feat=feat.dropna()
Z=(feat-feat.mean())/feat.std()
cand_keys=[m.replace('hsa-','').upper() for m in cand['mirna']]
cand_keys=[k for k in cand_keys if k in Z.index]
pools={}
for k in cand_keys:
    d=((Z-Z.loc[k])**2).sum(axis=1)
    near=d.drop(index=cand_keys,errors='ignore').nsmallest(20).index.tolist()
    pools[k]=near
# observed under exact test
obs=[]
for _,r in cand.iterrows():
    key=r['mirna'].replace('hsa-','').upper()
    sub=tg[tg.bare.str.contains(key,regex=False)]
    if key=='MIR-375-3P': sub=tg[tg.bare.isin(['MIR-375-3P','MIR-375'])]
    if len(sub)==0: obs.append(None); continue
    direction=1 if r['d_sev_mild']<0 else -1
    p,k,x=exact_p(sub.target_gene.tolist(),direction)
    obs.append({'mirna':key,'p':p,'k':k,'x':x})
obs=[o for o in obs if o]
op=np.array([o['p'] for o in obs]); o=np.argsort(op); bh=op[o]*len(op)/(np.arange(len(op))+1)
fdr_obs=np.minimum.accumulate(bh[::-1])[::-1][np.argsort(o)]
obs_pass=int((fdr_obs<=0.05).sum())
# Control A
A_passcounts=[]
cand_dirs={r['mirna'].replace('hsa-','').upper():(1 if r['d_sev_mild']<0 else -1) for _,r in cand.iterrows()}
for d in range(DRAWS):
    ps=[]
    for k in cand_keys:
        m=rng.choice(pools[k])
        sub=tg[tg.bare.str.contains(m,regex=False)]
        if len(sub)==0: continue
        p,_,_=exact_p(sub.target_gene.tolist(),cand_dirs[k])
        if p is not None: ps.append(p)
    if ps:
        pa=np.array(ps); o=np.argsort(pa); bh=pa[o]*len(pa)/(np.arange(len(pa))+1)
        fdr=np.minimum.accumulate(bh[::-1])[::-1][np.argsort(o)]
        A_passcounts.append(int((fdr<=0.05).sum()))
    if (d+1)%50==0: print(f'A {d+1}/{DRAWS}',flush=True)
# Control B
fam=pd.read_csv('sources/services/targetscan/miR_Family_Info.txt',sep='\t')
fam=fam[fam['Species ID']==9606]
m2f=dict(zip(fam['MiRBase ID'].str.upper().str.replace('HSA-','',regex=False),fam['miR family']))
f2m=fam.assign(bare=fam['MiRBase ID'].str.upper().str.replace('HSA-','',regex=False)).groupby('miR family')['bare'].apply(list).to_dict()
B={}
for _,r in cand.iterrows():
    key=r['mirna'].replace('hsa-','').upper()
    sub=tg[tg.bare.str.contains(key,regex=False)]
    if key=='MIR-375-3P': sub=tg[tg.bare.isin(['MIR-375-3P','MIR-375'])]
    if len(sub)==0: continue
    k=len(sub.target_gene.unique())
    family=m2f.get(key)
    members=[m for m in f2m.get(family,[]) if m!=key] if family else []
    ft=tg[tg.bare.isin(members)].target_gene.unique()
    if len(ft)<max(5,k): B[key]={'status':'insufficient family targets','family':family,'n_family_targets':int(len(ft))}; continue
    direction=cand_dirs[key]; ps=[]
    for d in range(DRAWS):
        samp=rng.choice(ft,size=min(k,len(ft)),replace=False)
        p,_,_=exact_p(samp.tolist(),direction)
        if p is not None: ps.append(p)
    ps=np.array(ps)
    op0=[o2['p'] for o2 in obs if o2['mirna']==key]
    B[key]={'family':family,'n_family_targets':int(len(ft)),
            'median_p_family_null':float(np.median(ps)),
            'frac_draws_p_le_obs':float((ps<=(op0[0] if op0 else 0)).mean())}
    print(f'B {key} done',flush=True)
out={'spec':'ADDENDUM_5 E3; exact-hypergeometric substitution disclosed',
 'observed_cardiac_pass_exact':obs_pass,'observed_n':len(obs),
 'control_A':{'draws':len(A_passcounts),'median_pass':float(np.median(A_passcounts)),
   'p95_pass':float(np.percentile(A_passcounts,95)),'max_pass':int(np.max(A_passcounts)),
   'frac_draws_ge_observed':float((np.array(A_passcounts)>=obs_pass).mean())},
 'control_B':B,'control_C':'D3 TargetScan arm done (results/d3_targetscan_sensitivity.json): cardiac 12/12 exact reproduction',
 'K_up':K_up,'K_dn':K_dn,'N_genes':N}
json.dump(out,open('results/e3_mirtarbase_bias_controls.json','w'),indent=2)
print(json.dumps({k:out[k] for k in ('observed_cardiac_pass_exact','observed_n','control_A')},indent=1))
