#!/usr/bin/env python3
"""E4 per ADDENDUM_5 (non-gate, descriptive): disease-pathway enrichment of
blood-cohort direction-consistent DE genes among candidate validated targets.
Union over the 20 candidates of validated targets that are direction-
consistent DE (nominal p<=0.05, candidate-predicted direction) in
GSE244827; Fisher exact vs pinned Enrichr libraries (re-use, ledgered:
KEGG_2021_Human, WikiPathways_2024_Human, GO_Biological_Process_2025),
BH FDR<=0.05; Chagas/heart-failure-relevant terms flagged. Descriptive only."""
import gzip, json
import numpy as np, pandas as pd
from scipy.stats import ttest_ind, fisher_exact
from statsmodels.stats.multitest import multipletests
def load_counts(path, meta_cols):
    with gzip.open(path,'rt',errors='replace') as f: df=pd.read_csv(f,sep='\t')
    genes=df.iloc[:,0].astype(str).str.upper()
    X=df.iloc[:,meta_cols:].apply(pd.to_numeric,errors='coerce').fillna(0).values
    return genes,X
g,X=load_counts('sources/matrices/GSE244827_CHAVArawcounts.txt.gz',6)
lmap=pd.read_csv('sources/GSE244827_column_label_map.csv'); b2l=dict(zip(lmap.b_code,lmap.label))
with gzip.open('sources/matrices/GSE244827_CHAVArawcounts.txt.gz','rt') as f:
    mat_cols=f.readline().rstrip('\n').split('\t')[6:]
labels=[b2l[c] for c in mat_cols]
case=[i for i,l in enumerate(labels) if l=='case']; ctrl=[i for i,l in enumerate(labels) if l=='control']
bmap=pd.read_csv('sources/services/ensembl_biomart/ensg_symbol_map.tsv',sep='\t')
e2s=dict(zip(bmap['Gene stable ID'],bmap['Gene name'].astype(str).str.upper()))
g=g.map(lambda x:e2s.get(x,x))
lib=X.sum(axis=0); L=np.log2(X/lib*1e6+1)
t,p=ttest_ind(L[:,case],L[:,ctrl],axis=1,equal_var=False)
fc=L[:,case].mean(axis=1)-L[:,ctrl].mean(axis=1); p=np.nan_to_num(p,nan=1.0)
de_up=set(g[(fc>0)&(p<=0.05)]); de_dn=set(g[(fc<0)&(p<=0.05)])
tg=pd.read_csv('sources/services/mirtarbase/validated_targets.tsv',sep='\t')
tg['bare']=tg.mirna.str.upper().str.replace('HSA-','',regex=False)
cand=pd.read_csv('results/h2_gate_a_passing.csv'); cand=cand[cand['class']=='CANDIDATE']
U=set()
for _,r in cand.iterrows():
    key=r['mirna'].replace('hsa-','').upper()
    sub=tg[tg.bare.str.contains(key,regex=False)]
    if key=='MIR-375-3P': sub=tg[tg.bare.isin(['MIR-375-3P','MIR-375'])]
    tt=set(t.upper() for t in sub.target_gene)
    U |= (tt&de_up) if r['d_sev_mild']<0 else (tt&de_dn)
print('blood direction-consistent DE target union:',len(U),flush=True)
BG=set(g)
rows=[]
for libf in ['KEGG_2021_Human','WikiPathways_2024_Human','GO_Biological_Process_2025']:
    libj=json.load(open(f'results/enrichr/{libf}.json'))
    data=libj.get(libf,libj)
    for row in data:
        if isinstance(row,list) and len(row)>5 and isinstance(row[5],list):
            term=row[1]; gs=set(x.upper() for x in row[5])
        elif isinstance(row,dict):
            term=row.get('term',''); gs=set(x.upper() for x in row.get('geneSymbols',[]))
        else: continue
        gs&=BG
        if len(gs)<5: continue
        a=len(U&gs)
        if a<2: continue
        orr,pv=fisher_exact([[a,len(U)-a],[len(gs)-a,len(BG)-len(U)-(len(gs)-a)]],alternative='greater')
        rows.append({'library':libf,'term':term,'overlap':a,'set_size':len(gs),'odds':float(orr),'p':float(pv)})
R=pd.DataFrame(rows)
R['fdr']=multipletests(R.p,method='fdr_bh')[1]
R=R.sort_values('fdr')
R.to_csv('results/e4_blood_pathway_enrichment.csv',index=False)
sig=R[R.fdr<=0.05]
kw=['chagas','trypanosoma','heart','cardiac','cardiomyopathy','hypertrophy','fibrosis','dilated']
flag=sig[sig.term.str.lower().str.contains('|'.join(kw))]
out={'spec':'ADDENDUM_5 E4; non-gate descriptive; pinned-library re-use',
 'union_targets_de_blood':len(U),'terms_tested':len(R),'n_fdr_pass':int(len(sig)),
 'disease_relevant_fdr_pass':flag[['library','term','overlap','set_size','fdr']].to_dict('records')[:20],
 'top10':R.head(10)[['library','term','overlap','set_size','fdr']].to_dict('records')}
json.dump(out,open('results/e4_blood_pathway.json','w'),indent=2)
print(json.dumps({'union':len(U),'tested':len(R),'fdr_pass':len(sig),'disease_relevant':len(flag)},indent=1))
print(flag[['library','term','overlap','fdr']].head(12).to_string(index=False))
