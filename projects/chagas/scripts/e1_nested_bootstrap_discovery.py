#!/usr/bin/env python3
"""E1 per ADDENDUM_5 (locked 2026-09-27T09:23): nested bootstrap discovery.
B=1000 bootstrap resamples of the 146 graded GSE299582 samples; inside EACH
bootstrap, full discovery pipeline on the bootstrap sample only: KW screen
(BH FDR<=0.05 & |d|>=0.5) -> frozen 11-marker exclusion -> gate-(b)
direction-consistent miRTarBase enrichment on BOTH external cohorts (1,000
perms/test, disclosed reduction from 10k) -> AND intersection -> signed
module (bootstrap sign weights + z-params) -> score OUT-OF-BOOTSTRAP.
Readings locked in ADDENDUM_5 before this run. Tie-corrected KW is
vectorized (identical statistic to scripts/d2_oof_module_reconstruction.py,
faster implementation)."""
import gzip, json, time, itertools
import numpy as np, pandas as pd
from scipy.stats import rankdata, chi2, ttest_ind

SEED = 20260927; B = 1000; INNER_PERM = 1000
EXCL = set("""miR-143-3p miR-223-3p miR-486-5p miR-3960 miR-6734-5p miR-1285-5p
miR-10527-5p miR-1228-5p miR-30c-3p miR-146a miR-208a""".split())
LEGACY_ALIASES = {'MIR-375-3P': ['MIR-375']}
t0 = time.time()

# ---- discovery cohort (identical construction to h1prime/d2) ----
xw = pd.read_csv('sources/GSE299582_sample_crosswalk.csv')
xw['ch'] = xw['characteristics'].map(lambda s: json.loads(s.replace('""','"')) if s else {})
xw['sev'] = xw['ch'].map(lambda c: c.get('chronic chagas_cardiomyopathy_severity','?'))
xw['prefix'] = xw['title'].str.split(' - ').str[0].str.strip()
with gzip.open('sources/matrices/GSE299582_normalized_counts.csv.gz','rt') as f:
    M = pd.read_csv(f, index_col=0)
od = xw[xw.sev.isin(['mild','moderate','severe']) | (xw.label=='control')].copy()
od['y'] = od.apply(lambda r: 0 if r['label']=='control' else ['mild','moderate','severe'].index(r['sev'])+1, axis=1)
X = M[od['prefix']]
L = np.log2(X.div(X.sum(axis=0), axis=1)*1e6 + 1).T.values
y = od['y'].values
names = np.array([m.replace('hsa-','') for m in M.index])
excl_mask = np.isin(names, list(EXCL))
print(f'cohort {L.shape}, setup {time.time()-t0:.0f}s', flush=True)

def kw_screen(V, yv):
    """Tie-corrected KW across 4 groups, vectorized. V: samples x features."""
    N = len(yv)
    R = np.apply_along_axis(rankdata, 0, V)
    H = np.zeros(V.shape[1])
    for g in range(4):
        mask = yv == g
        H += R[mask].sum(axis=0)**2 / mask.sum()
    H = 12/(N*(N+1)) * H - 3*(N+1)
    # vectorized tie correction: encode (value,column) pairs, one np.unique
    F = V.shape[1]
    code = (np.round(V * 1e6).astype(np.int64) * F + np.arange(F)).ravel()
    _, inv, cnt = np.unique(code, return_inverse=True, return_counts=True)
    col_ids = np.arange(code.size) % F
    tsum = np.bincount(col_ids, weights=(cnt[inv]**3 - cnt[inv]).astype(float),
                       minlength=F)
    tie = 1 - tsum/(N**3 - N)
    Hc = np.where(tie > 0, H/np.maximum(tie, 1e-12), 0.0)
    return chi2.sf(Hc, 3)

def medians_by_group(V, yv):
    return np.stack([np.median(V[yv==g], axis=0) for g in range(4)])

# ---- external cohort DE (computed ONCE; frozen machinery from h2_gate_b) ----
def load_counts(path, meta_cols):
    with gzip.open(path,'rt',errors='replace') as f:
        df = pd.read_csv(f, sep='\t')
    genes = df.iloc[:,0].astype(str).str.upper()
    Xc = df.iloc[:,meta_cols:].apply(pd.to_numeric, errors='coerce').fillna(0).values
    return genes, Xc

def de_stats(genes, Xc, case_idx, ctrl_idx):
    lib = Xc.sum(axis=0); Lc = np.log2(Xc/lib*1e6+1)
    t,p = ttest_ind(Lc[:,case_idx], Lc[:,ctrl_idx], axis=1, equal_var=False)
    fc = Lc[:,case_idx].mean(axis=1) - Lc[:,ctrl_idx].mean(axis=1)
    return genes, fc, np.nan_to_num(p, nan=1.0)

g1, X1 = load_counts('sources/matrices/GSE244827_CHAVArawcounts.txt.gz', 6)
with gzip.open('sources/matrices/GSE244827_CHAVArawcounts.txt.gz','rt') as f:
    mat_cols = f.readline().rstrip('\n').split('\t')[6:]
lmap = pd.read_csv('sources/GSE244827_column_label_map.csv')
b2l = dict(zip(lmap.b_code, lmap.label))
labels = [b2l[c] for c in mat_cols]
case = [i for i,l in enumerate(labels) if l=='case']; ctrl = [i for i,l in enumerate(labels) if l=='control']
bmap = pd.read_csv('sources/services/ensembl_biomart/ensg_symbol_map.tsv', sep='\t')
e2s = dict(zip(bmap['Gene stable ID'], bmap['Gene name'].astype(str).str.upper()))
g1 = g1.map(lambda x: e2s.get(x, x))
g1, fc1, p1 = de_stats(g1, X1, case, ctrl)
g2m, X2 = load_counts('sources/matrices/GSE203525_Counts.txt.gz', 1)
cols = pd.read_csv('sources/matrices/GSE203525_Counts.txt.gz', sep='\t', nrows=0).columns[1:]
cc = [i for i,c in enumerate(cols) if c.startswith('CC') and '0hpi' in c]
ind = [i for i,c in enumerate(cols) if c.startswith('IND') and '0hpi' in c]
g2, fc2, p2 = de_stats(g2m, X2, cc, ind)
cohorts = {}
for lbl, gg, ff, pp in [('blood_GSE244827', g1, fc1, p1), ('hipsc_GSE203525', g2, fc2, p2)]:
    gene_idx = {g: i for i, g in enumerate(gg)}
    cohorts[lbl] = (gene_idx, np.sign(ff), pp, len(gg))
print(f'external DE ready {time.time()-t0:.0f}s', flush=True)

# ---- miRTarBase targets per miRNA ----
tg = pd.read_csv('sources/services/mirtarbase/validated_targets.tsv', sep='\t')
tg['key'] = tg.mirna.str.upper()
def targets_for(mir):
    key = mir.upper()
    sub = tg[tg.key.str.contains(key, regex=False)]
    if len(sub) == 0 and key in LEGACY_ALIASES:
        sub = tg[tg.key.isin([f'HSA-{a}' for a in LEGACY_ALIASES[key]])]
    return sub.target_gene.tolist() if len(sub) else None

def enrich_p(mir, expected, rng):
    """1,000-perm size-preserving enrichment p per cohort; None if <5 mapped."""
    tgts = targets_for(mir)
    if not tgts: return None
    out = {}
    for lbl, (gene_idx, sgn, pp, G) in cohorts.items():
        tidx = [gene_idx[t.upper()] for t in tgts if t.upper() in gene_idx]
        if len(tidx) < 5: return None
        tidx = np.array(tidx)
        de = (expected * sgn > 0) & (pp <= 0.05)
        obs = de[tidx].mean()
        draws = rng.integers(0, G, size=(INNER_PERM, len(tidx)))
        null = de[draws].mean(axis=1)
        out[lbl] = (1 + (null >= obs).sum())/(INNER_PERM+1)
    return out

def cindex(yy, s):
    c = t = 0
    for i, j in itertools.combinations(range(len(yy)), 2):
        if yy[i] == yy[j]: continue
        t += 1; c += (s[i]-s[j])*(yy[i]-yy[j]) > 0
    return c/t if t else np.nan

# ---- bootstrap loop ----
rng = np.random.default_rng(SEED)
N = len(y)
disc_count = pd.Series(0, index=names, dtype=int)   # full-pipeline discoveries
oob_c = []; member_hist = []; screen_pass_hist = []
for b in range(B):
    idx = rng.integers(0, N, N)
    oob = np.setdiff1d(np.arange(N), np.unique(idx))
    yb = y[idx]
    if len(np.unique(yb)) < 4 or len(oob) < 10:
        oob_c.append(np.nan); member_hist.append(0); screen_pass_hist.append(0); continue
    V = L[idx]
    p = kw_screen(V, yb)
    o = np.argsort(p); bh = p[o]*len(p)/(np.arange(len(p))+1)
    fdr = np.empty_like(bh); fdr[o] = np.minimum.accumulate(bh[::-1])[::-1]
    gm = medians_by_group(V, yb); d = gm[3] - gm[1]
    sel = (fdr <= 0.05) & (np.abs(d) >= 0.5) & ~excl_mask
    screen_pass_hist.append(int(sel.sum()))
    # gate-(b) both cohorts
    tests = []
    for j in np.where(sel)[0]:
        r = enrich_p(names[j], 1 if d[j] < 0 else -1, rng)
        if r: tests.append((j, r))
    if not tests:
        oob_c.append(np.nan); member_hist.append(0); continue
    ps = np.array([pval for _, r in tests for pval in r.values()])
    o = np.argsort(ps); bhv = ps[o]*len(ps)/(np.arange(len(ps))+1)
    fdrs = np.minimum.accumulate(bhv[::-1])[::-1][np.argsort(o)]
    k = 0; members = []
    for j, r in tests:
        pb, ph = fdrs[k], fdrs[k+1]; k += 2
        if pb <= 0.05 and ph <= 0.05:
            members.append(j)
    member_hist.append(len(members))
    for j in members: disc_count[names[j]] += 1
    if not members or len(oob) < 10:
        oob_c.append(np.nan); continue
    members = np.array(members)
    w = np.sign(d[members])
    mu = V[:, members].mean(axis=0); sd = V[:, members].std(axis=0)
    sd = np.where(sd == 0, 1.0, sd)
    s_oob = (w * ((L[oob][:, members] - mu)/sd)).sum(axis=1)
    oob_c.append(cindex(y[oob], s_oob))
    if (b+1) % 100 == 0:
        print(f'bootstrap {b+1}/{B} ({time.time()-t0:.0f}s)', flush=True)

oob_c = np.array(oob_c, dtype=float)
valid = ~np.isnan(oob_c)
core6 = ['miR-1-3p','miR-122-5p','miR-192-5p','miR-30c-5p','miR-145-5p','miR-194-5p']
out = {'B': B, 'inner_perm': INNER_PERM,
       'oob_cindex_median': float(np.nanmedian(oob_c)),
       'oob_cindex_q025': float(np.nanquantile(oob_c, 0.025)),
       'oob_cindex_q975': float(np.nanquantile(oob_c, 0.975)),
       'n_valid_bootstraps': int(valid.sum()),
       'core6_discovery_frequency': {m: int(disc_count.get(m, 0))/B for m in core6},
       'top20_discovery_frequency': disc_count.sort_values(ascending=False).head(20).div(B).to_dict(),
       'member_count_median': float(np.median(member_hist)),
       'screen_pass_median': float(np.median(screen_pass_hist)),
       'spec': 'ADDENDUM_5 E1 locked 2026-09-27T09:23; readings fixed both ways',
       'runtime_s': round(time.time()-t0, 1)}
json.dump(out, open('results/e1_nested_bootstrap_discovery.json','w'), indent=2)
pd.DataFrame({'bootstrap': range(B), 'oob_cindex': oob_c, 'n_members': member_hist,
              'n_screen_pass': screen_pass_hist}).to_csv('results/e1_bootstrap_detail.csv', index=False)
print(json.dumps(out, indent=2))
