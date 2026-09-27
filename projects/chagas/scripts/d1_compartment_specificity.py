#!/usr/bin/env python3
"""D1 per ADDENDUM_4 (locked 2026-09-27T08:58): compartment-specificity
decisive test. Gene-level logistic model over the union of the 20
candidates' miRTarBase-validated mapped targets:
  replication ~ blood_abundance + cardiac_spec
Response: direction-consistent DE in GSE244827 early-CCC-vs-seronegative
(same machinery as h2_gate_b_enrichment.py: Welch t on log2(CPM+1),
predicted direction per the frozen rule, nominal p<=0.05). A gene targeted
by multiple candidates counts as replicated if direction-consistent for AT
LEAST ONE targeting candidate (disclosed disambiguation of the locked
spec, fixed before seeing results).
Predictors: HPA v23 PINNED files (rna_blood_cell.tsv.zip does not exist
under that name on v22-v25; HPA's blood-cell dataset is distributed as
rna_immune_cell.tsv.zip - substitution disclosed; version pin disclosed
like the miRTarBase v8.0 pin):
  blood_abundance = log1p(max nTPM across immune cell types)
  cardiac_spec    = log2((heart_muscle_nTPM+1)/(max_immune_nTPM+1))
READING locked both ways in ADDENDUM_4: blood_abundance coef positive with
Wald p<0.05 supports compartment specificity; null/negative DOWNGRADES
paper 9.3b (correction recorded, not silent)."""
import gzip, json
import numpy as np, pandas as pd
from scipy.stats import ttest_ind
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score

LEGACY_ALIASES = {'MIR-375-3P': ['MIR-375']}

# --- GSE244827 DE (frozen machinery) ---
with gzip.open('sources/matrices/GSE244827_CHAVArawcounts.txt.gz','rt',errors='replace') as f:
    df = pd.read_csv(f, sep='\t')
genes = df.iloc[:,0].astype(str).str.upper()
Xc = df.iloc[:,6:].apply(pd.to_numeric, errors='coerce').fillna(0).values
with gzip.open('sources/matrices/GSE244827_CHAVArawcounts.txt.gz','rt') as f:
    mat_cols = f.readline().rstrip('\n').split('\t')[6:]
lmap = pd.read_csv('sources/GSE244827_column_label_map.csv')
b2l = dict(zip(lmap.b_code, lmap.label))
labels = [b2l[c] for c in mat_cols]
case = [i for i,l in enumerate(labels) if l=='case']; ctrl = [i for i,l in enumerate(labels) if l=='control']
bmap = pd.read_csv('sources/services/ensembl_biomart/ensg_symbol_map.tsv', sep='\t')
e2s = dict(zip(bmap['Gene stable ID'], bmap['Gene name'].astype(str).str.upper()))
genes = genes.map(lambda x: e2s.get(x, x))
lib = Xc.sum(axis=0); Lc = np.log2(Xc/lib*1e6+1)
t,p = ttest_ind(Lc[:,case], Lc[:,ctrl], axis=1, equal_var=False)
fc = Lc[:,case].mean(axis=1) - Lc[:,ctrl].mean(axis=1)
p = np.nan_to_num(p, nan=1.0)
de_up = set(genes[(fc > 0) & (p <= 0.05)])
de_dn = set(genes[(fc < 0) & (p <= 0.05)])

# --- candidates + targets (frozen gate-b rules) ---
cand = pd.read_csv('results/h2_gate_a_passing.csv')
cand = cand[cand['class']=='CANDIDATE']
tg = pd.read_csv('sources/services/mirtarbase/validated_targets.tsv', sep='\t')
tg['key'] = tg.mirna.str.upper()
gene_expected = {}   # gene -> list of expected directions from its candidates
for _, r in cand.iterrows():
    mir = r['mirna'].replace('hsa-',''); key = mir.upper()
    sub = tg[tg.key.str.contains(key, regex=False)]
    if len(sub)==0 and key in LEGACY_ALIASES:
        sub = tg[tg.key.isin([f'HSA-{a}' for a in LEGACY_ALIASES[key]])]
    if len(sub)==0: continue
    expected = 1 if r['d_sev_mild'] < 0 else -1
    for gname in sub.target_gene.astype(str).str.upper():
        gene_expected.setdefault(gname, []).append(expected)
union_genes = sorted(gene_expected)
resp = {}
for gname, exps in gene_expected.items():
    ok = any((gname in de_up) if e==1 else (gname in de_dn) for e in exps)
    resp[gname] = int(ok)

# --- HPA v23 predictors ---
imm = pd.read_csv('sources/services/hpa/rna_immune_cell.tsv.zip', sep='\t', compression='zip')
imm['g'] = imm['Gene name'].astype(str).str.upper()
blood_nx = imm.groupby('g')['nTPM'].max()
cons = pd.read_csv('sources/services/hpa/rna_tissue_consensus.tsv.zip', sep='\t', compression='zip')
cons['g'] = cons['Gene name'].astype(str).str.upper()
heart = cons[cons['Tissue']=='heart muscle'].groupby('g')['nTPM'].max()
rows = []
for gname in union_genes:
    if gname not in blood_nx.index and gname not in heart.index: continue
    b = float(blood_nx.get(gname, 0.0)); h = float(heart.get(gname, 0.0))
    rows.append({'gene': gname, 'replication': resp[gname],
                 'blood_abundance': np.log1p(b),
                 'cardiac_spec': np.log2((h+1)/(b+1))})
D = pd.DataFrame(rows)
X = sm.add_constant(D[['blood_abundance','cardiac_spec']])
fit = sm.Logit(D['replication'], X).fit(disp=0)
auc = roc_auc_score(D['replication'], fit.predict(X))
out = {'n_genes': len(D), 'n_replicated': int(D['replication'].sum()),
       'coef': fit.params.to_dict(), 'wald_p': fit.pvalues.to_dict(),
       'descriptive_auc': float(auc),
       'hpa_pin': 'HPA v23 (v23.proteinatlas.org), files dated 2023-06-17',
       'hpa_files_sha256': {
         'rna_immune_cell.tsv.zip': 'ac4a219a92fac1d27c7610f35d0053b0502b1aedefd3b95c4b37a42021028ba2',
         'rna_tissue_consensus.tsv.zip': 'b9ac6bbdf8152524ff845767638f4c92389b2933c5ef02060dee68b873176095'},
       'substitution': 'rna_blood_cell.tsv.zip 404 on v22-v25; used HPA blood-cell dataset file rna_immune_cell.tsv.zip (disclosed)',
       'response_rule': 'any-consistent across targeting candidates (disclosed disambiguation, fixed before results)',
       'reading_locked': 'ADDENDUM_4 D1: blood_abundance coef>0 & p<0.05 supports compartment specificity; null/negative downgrades 9.3b'}
json.dump(out, open('results/d1_compartment_specificity.json','w'), indent=2)
D.to_csv('results/d1_gene_level.csv', index=False)
print(fit.summary2().tables[1])
print('n', len(D), 'replicated', int(D['replication'].sum()), 'AUC', round(auc,3))
