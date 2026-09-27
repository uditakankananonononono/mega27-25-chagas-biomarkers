#!/usr/bin/env python3
"""E2 per ADDENDUM_5 (locked 2026-09-27T09:23): leave-one-cohort-out
biological validation + cross-compartment concordance.
DirectionScore per candidate = fraction of independent tissue legs agreeing
with the predicted direction. Legs (locked): serum GSE299582 (own direction
in CCC-vs-control DE), hiPSC-CM GSE203525 (target-program enrichment),
whole blood GSE244827 (same), heart LV-wall GSE191081 sevCCC-vs-CTRL
(matrix retrieved 2026-09-27, GEO supplementary, sha256 logged).
Enrichment legs use the hypergeometric analytic equivalent of gate-(b)'s
size-preserving without-replacement null (same distribution, disclosed;
feasible at the 10,000-set scale where per-test permutations are not).
Null: 10,000 random 6-sets matched member-wise to CORE6 on miRTarBase
validated-target count (x0.5-2), serum mean log2CPM (x0.5-2) and
publication count (x0.5-2; Europe PMC TITLE_ABS counts - disclosed
substitution for the locked esearch fallback, esearch backend down at
run time). READING locked: concordance unusual only if CORE6 mean
DirectionScore >= 95th percentile of the matched null."""
import gzip, json, time
import numpy as np, pandas as pd
from scipy.stats import ttest_ind, hypergeom

SEED = 20260927; NSET = 10000
rng = np.random.default_rng(SEED)
t0 = time.time()
LEGACY_ALIASES = {'MIR-375-3P': ['MIR-375']}
CORE6 = ['miR-1-3p','miR-122-5p','miR-192-5p','miR-30c-5p','miR-145-5p','miR-194-5p']

# ---- serum panel + DE ----
xw = pd.read_csv('sources/GSE299582_sample_crosswalk.csv')
xw['ch'] = xw['characteristics'].map(lambda s: json.loads(s.replace('""','"')) if s else {})
xw['sev'] = xw['ch'].map(lambda c: c.get('chronic chagas_cardiomyopathy_severity','?'))
xw['prefix'] = xw['title'].str.split(' - ').str[0].str.strip()
with gzip.open('sources/matrices/GSE299582_normalized_counts.csv.gz','rt') as f:
    M = pd.read_csv(f, index_col=0)
od = xw[xw.sev.isin(['mild','moderate','severe']) | (xw.label=='control')].copy()
od['y'] = od.apply(lambda r: 0 if r['label']=='control' else ['mild','moderate','severe'].index(r['sev'])+1, axis=1)
Xs = M[od['prefix']]
Ls = np.log2(Xs.div(Xs.sum(axis=0), axis=1)*1e6 + 1)
names = np.array([m.replace('hsa-','') for m in M.index])
abund = Ls.values.mean(axis=1)
y = od['y'].values
case_idx = np.where(y>0)[0]; ctrl_idx = np.where(y==0)[0]
t_s, p_s = ttest_ind(Ls.values[case_idx].T if False else Ls.values[:,case_idx], Ls.values[:,ctrl_idx], axis=1, equal_var=False)
fc_s = Ls.values[:,case_idx].mean(axis=1) - Ls.values[:,ctrl_idx].mean(axis=1)
p_s = np.nan_to_num(p_s, nan=1.0)
print('serum DE done', flush=True)

# ---- external cohort DE (heart added) ----
def de_stats(genes, X, case_idx, ctrl_idx):
    lib = X.sum(axis=0); L = np.log2(X/lib*1e6+1)
    t,p = ttest_ind(L[:,case_idx], L[:,ctrl_idx], axis=1, equal_var=False)
    fc = L[:,case_idx].mean(axis=1) - L[:,ctrl_idx].mean(axis=1)
    return genes.values if hasattr(genes,'values') else genes, fc, np.nan_to_num(p, nan=1.0)

bmap = pd.read_csv('sources/services/ensembl_biomart/ensg_symbol_map.tsv', sep='\t')
e2s = dict(zip(bmap['Gene stable ID'], bmap['Gene name'].astype(str).str.upper()))

with gzip.open('sources/matrices/GSE244827_CHAVArawcounts.txt.gz','rt',errors='replace') as f:
    df = pd.read_csv(f, sep='\t')
g1 = df.iloc[:,0].astype(str).str.upper().map(lambda x: e2s.get(x, x))
X1 = df.iloc[:,6:].apply(pd.to_numeric, errors='coerce').fillna(0).values
lmap = pd.read_csv('sources/GSE244827_column_label_map.csv'); b2l = dict(zip(lmap.b_code, lmap.label))
labels = [b2l[c] for c in df.columns[6:]]
case = [i for i,l in enumerate(labels) if l=='case']; ctrl = [i for i,l in enumerate(labels) if l=='control']
g1, fc1, p1 = de_stats(g1, X1, case, ctrl)
print('blood DE done', flush=True)

with gzip.open('sources/matrices/GSE203525_Counts.txt.gz','rt',errors='replace') as f:
    df2 = pd.read_csv(f, sep='\t')
g2 = df2.iloc[:,0].astype(str).str.upper()
X2 = df2.iloc[:,1:].apply(pd.to_numeric, errors='coerce').fillna(0).values
cols = df2.columns[1:]
cc = [i for i,c in enumerate(cols) if c.startswith('CC') and '0hpi' in c]
ind = [i for i,c in enumerate(cols) if c.startswith('IND') and '0hpi' in c]
g2, fc2, p2 = de_stats(g2, X2, cc, ind)
print('hipsc DE done', flush=True)

with gzip.open('sources/matrices/GSE191081_Count_table.txt.gz','rt',errors='replace') as f:
    df3 = pd.read_csv(f, sep='\t')
g3 = df3.iloc[:,0].astype(str).map(lambda x: e2s.get(x, x.upper() if isinstance(x,str) else x)).str.upper()
X3 = df3.iloc[:,1:].apply(pd.to_numeric, errors='coerce').fillna(0).values
cols3 = list(df3.columns[1:])
ccc = [i for i,c in enumerate(cols3) if c.startswith('sevCCC')]
hctrl = [i for i,c in enumerate(cols3) if c.startswith('CTRL')]
g3, fc3, p3 = de_stats(g3, X3, ccc, hctrl)
print('heart DE done', flush=True)

cohorts = {'serum_GSE299582': (names, fc_s, p_s, True),
           'hipsc_GSE203525': (g2, fc2, p2, False),
           'blood_GSE244827': (g1, fc1, p1, False),
           'heart_GSE191081': (g3, fc3, p3, False)}
# precomputed per-cohort: gene->idx, and per expected-direction DE index sets
COH = {}
for lbl, (gg, ff, pp, is_serum) in cohorts.items():
    if is_serum: continue
    gene_idx = {g: i for i, g in enumerate(gg)}
    de_up = set(np.where((np.sign(ff) > 0) & (pp <= 0.05))[0])
    de_dn = set(np.where((np.sign(ff) < 0) & (pp <= 0.05))[0])
    COH[lbl] = (gene_idx, de_up, de_dn, len(gg))

# ---- miRTarBase targets per miRNA (panel-wide) ----
tg = pd.read_csv('sources/services/mirtarbase/validated_targets.tsv', sep='\t')
tg['key'] = tg.mirna.str.upper()
def _targets_for(mir):
    key = mir.upper()
    sub = tg[tg.key.str.contains(key, regex=False)]
    if len(sub) == 0 and key in LEGACY_ALIASES:
        sub = tg[tg.key.isin([f'HSA-{a}' for a in LEGACY_ALIASES[key]])]
    return sub.target_gene.astype(str).str.upper().unique().tolist() if len(sub) else []
TG_MAP = {n: _targets_for(n) for n in names}
def targets_for(mir):
    return TG_MAP.get(mir, [])

# ---- per-miRNA legs ----
def leg_scores(mir, d_sign):
    """d_sign: +1 up-in-severe, -1 down-in-severe (from screen). Returns dict leg->bool."""
    out = {}
    j = np.where(names == mir)[0]
    if len(j):
        j = j[0]
        out['serum_GSE299582'] = bool(np.sign(fc_s[j]) == d_sign and p_s[j] <= 0.05)
    tgts = targets_for(mir)
    expected = 1 if d_sign < 0 else -1   # canonical repression
    for lbl, (gene_idx, de_up, de_dn, G) in COH.items():
        tidx = [gene_idx[t] for t in tgts if t in gene_idx]
        if len(tidx) < 5:
            out[lbl] = None; continue
        de_set = de_up if expected == 1 else de_dn
        D = len(de_set); k = len(tidx)
        H = sum(1 for i in tidx if i in de_set)
        pv = hypergeom.sf(H-1, G, D, k) if D > 0 else 1.0
        out[lbl] = bool(pv <= 0.05)
    return out

cand = pd.read_csv('results/h2_gate_a_passing.csv')
cand = cand[cand['class']=='CANDIDATE']
d_map = dict(zip(cand['mirna'].str.replace('hsa-','', regex=False), np.sign(cand['d_sev_mild'])))
name2idx = {n: i for i, n in enumerate(names)}
pub = json.load(open('sources/services/europepmc/mirna_pubcounts.json'))
tc_map = {n: len(targets_for(n)) for n in CORE6}
profiles = {}
for m in CORE6:
    profiles[m] = {'tc': tc_map[m], 'ab': abund[name2idx[m]], 'pc': pub.get(m)}

def matched_pool(prof):
    pool = []
    for i, n in enumerate(names):
        if n in CORE6: continue
        tcn = len(targets_for(n))
        if tcn < 5: continue
        pcn = pub.get(n)
        if pcn is None or prof['pc'] is None: continue
        if not (prof['tc']/2 <= tcn <= prof['tc']*2): continue
        if not (prof['ab']/2 <= abund[i] <= prof['ab']*2): continue
        if not (max(prof['pc'],1)/2 <= pcn <= max(prof['pc'],1)*2): continue
        pool.append(n)
    return pool

pools = {m: matched_pool(profiles[m]) for m in CORE6}
print({m: len(p) for m, p in pools.items()}, flush=True)

# observed CORE6 concordance
obs_legs = {m: leg_scores(m, d_map[m]) for m in CORE6}
def ds(legs):
    vals = [v for v in legs.values() if v is not None]
    return sum(vals)/len(vals) if vals else np.nan
obs_ds = {m: ds(obs_legs[m]) for m in CORE6}
obs_mean = float(np.nanmean(list(obs_ds.values())))

# matched null
null_means = np.full(NSET, np.nan)
for b in range(NSET):
    draws = []
    ok = True
    for m in CORE6:
        if len(pools[m]) == 0: ok = False; break
        draws.append(pools[m][rng.integers(0, len(pools[m]))])
    if not ok: continue
    dss = []
    for n in draws:
        # random miRNA direction: its own severe-mild sign is meaningless for
        # a non-candidate; use its serum CCC-vs-control sign as its 'direction'
        j = name2idx[n]
        sgn = 1 if fc_s[j] >= 0 else -1
        dss.append(ds(leg_scores(n, sgn)))
    null_means[b] = np.nanmean(dss)
    if (b+1) % 500 == 0: print(f'set {b+1}/{NSET} ({time.time()-t0:.0f}s)', flush=True)
valid = ~np.isnan(null_means)
pct = float((null_means[valid] < obs_mean).mean()*100) if valid.any() else np.nan
out = {'core6_direction_scores': obs_ds, 'core6_legs': obs_legs,
       'core6_mean_ds': obs_mean,
       'null_n': int(valid.sum()), 'null_median': float(np.nanmedian(null_means)),
       'null_q95': float(np.nanquantile(null_means, 0.95)),
       'core6_percentile_vs_matched_null': pct,
       'pool_sizes': {m: len(p) for m, p in pools.items()},
       'candidates_ds': {},
       'spec': 'ADDENDUM_5 E2; unusual only if >=95th percentile of matched null',
       'runtime_s': round(time.time()-t0, 1)}
for _, r in cand.iterrows():
    m = r['mirna'].replace('hsa-','')
    if m in CORE6: continue
    out['candidates_ds'][m] = ds(leg_scores(m, d_map[m]))
json.dump(out, open('results/e2_cross_tissue_concordance.json','w'), indent=2)
print(json.dumps(out, indent=2))
