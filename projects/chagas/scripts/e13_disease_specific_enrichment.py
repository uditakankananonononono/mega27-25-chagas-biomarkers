#!/usr/bin/env python3
"""E13 per ADDENDUM_5 (non-gate): disease-specific enrichment of the CORE6
352 strong-support validated targets against Chagas/heart-failure gene sets
(pinned Enrichr libraries, re-use). Fisher exact, BH FDR<=0.05; background =
union of the pinned libraries' genes. Replaces generic-only pathway language
where supported; specificity wording per round-03 weakness 11."""
import json
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
tg=pd.read_csv('sources/services/mirtarbase/validated_targets.tsv',sep='\t')
tg['bare']=tg.mirna.str.upper().str.replace('HSA-','',regex=False)
core6=['MIR-1-3P','MIR-122-5P','MIR-192-5P','MIR-30C-5P','MIR-145-5P','MIR-194-5P']
strong=tg[tg.bare.isin(core6)&(tg.support=='Functional MTI')]
U=set(t.upper() for t in strong.target_gene.unique())
assert len(U)==352, len(U)
rows=[]
libs={}
for libf in ['KEGG_2021_Human','WikiPathways_2024_Human','GO_Biological_Process_2025']:
    libj=json.load(open(f'results/enrichr/{libf}.json'))
    data=libj.get(libf,libj); terms={}
    for row in data:
        if isinstance(row,list) and len(row)>5 and isinstance(row[5],list):
            terms[row[1]]=set(x.upper() for x in row[5])
        elif isinstance(row,dict):
            terms[row.get('term','')]=set(x.upper() for x in row.get('geneSymbols',[]))
    libs[libf]=terms
BG=set().union(*[set().union(*t.values()) for t in libs.values() if t])
for libf,terms in libs.items():
    for term,gs0 in terms.items():
        gs=gs0&BG
        if len(gs)<5: continue
        a=len(U&gs)
        if a<2: continue
        orr,pv=fisher_exact([[a,len(U)-a],[len(gs)-a,len(BG)-len(U)-(len(gs)-a)]],alternative='greater')
        rows.append({'library':libf,'term':term,'overlap':a,'set_size':len(gs),'odds':float(orr),'p':float(pv)})
R=pd.DataFrame(rows)
R['fdr']=multipletests(R.p,method='fdr_bh')[1]
R=R.sort_values('fdr')
R.to_csv('results/e13_disease_enrichment.csv',index=False)
kw=['chagas','trypanosoma','heart','cardiac','cardiomyopathy','hypertrophy','fibrosis','dilated','myocard']
sig=R[R.fdr<=0.05]
flag=sig[sig.term.str.lower().str.contains('|'.join(kw))]
out={'spec':'ADDENDUM_5 E13; non-gate; pinned-library re-use','module_targets':len(U),
 'terms_tested':len(R),'n_fdr_pass':int(len(sig)),
 'disease_relevant':flag[['library','term','overlap','set_size','odds','fdr']].to_dict('records'),
 'top10':R.head(10)[['library','term','overlap','set_size','fdr']].to_dict('records')}
json.dump(out,open('results/e13_disease_specific_enrichment.json','w'),indent=2)
print('fdr_pass',len(sig),'of',len(R)); print(flag[['library','term','overlap','set_size','fdr']].head(15).to_string(index=False))
print(R.head(8)[['library','term','overlap','fdr']].to_string(index=False))
