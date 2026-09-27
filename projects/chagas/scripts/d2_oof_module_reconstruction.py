#!/usr/bin/env python3
"""D2 per ADDENDUM_4 (locked 2026-09-27T08:58, judge-02 Q2): OOF module
reconstruction inside the locked outer 5-fold splits (seed 20260926, same
StratifiedKFold construction as scripts/h1prime_ordinal.py), plus a
1,000-permutation label null. Per outer fold, on train only: KW screen
(BH FDR<=0.05 AND |d|>=0.5) -> frozen 11-marker exclusion -> intersect
with the committed both-tissue CORE6 filter (external cohorts) -> train
sign weights + train z-params -> score the held-out fold. Claim stands
only if empirical permutation p < 0.05 (reading locked before the run).
Vectorized KW with tie correction (ranks are label-invariant per
permutation; implementation detail, not design)."""
import gzip, json, itertools, time
import numpy as np, pandas as pd
from sklearn.model_selection import StratifiedKFold
from scipy.stats import rankdata, chi2, kruskal

SEED = 20260926; NPERM = 1000
CORE6 = ['miR-1-3p','miR-122-5p','miR-192-5p','miR-30c-5p','miR-145-5p','miR-194-5p']
EXCL = set("""miR-143-3p miR-223-3p miR-486-5p miR-3960 miR-6734-5p miR-1285-5p
miR-10527-5p miR-1228-5p miR-30c-3p miR-146a miR-208a""".split())

xw = pd.read_csv('sources/GSE299582_sample_crosswalk.csv')
xw['ch'] = xw['characteristics'].map(lambda s: json.loads(s.replace('""','"')) if s else {})
xw['sev'] = xw['ch'].map(lambda c: c.get('chronic chagas_cardiomyopathy_severity','?'))
xw['prefix'] = xw['title'].str.split(' - ').str[0].str.strip()
with gzip.open('sources/matrices/GSE299582_normalized_counts.csv.gz','rt') as f:
    M = pd.read_csv(f, index_col=0)
od = xw[xw.sev.isin(['mild','moderate','severe']) | (xw.label=='control')].copy()
od['y'] = od.apply(lambda r: 0 if r['label']=='control' else ['mild','moderate','severe'].index(r['sev'])+1, axis=1)
X = M[od['prefix']]
L = np.log2(X.div(X.sum(axis=0), axis=1)*1e6 + 1).T.values   # samples x mirnas
y = od['y'].values
names = np.array([m.replace('hsa-','') for m in M.index])
core6_mask = np.isin(names, CORE6)
excl_mask = np.isin(names, list(EXCL))

def kw_screen(V, yv):
    """Vectorized KW with tie correction across 4 groups. V: samples x features."""
    N = len(yv)
    R = np.apply_along_axis(rankdata, 0, V)          # ranks per column
    H = np.zeros(V.shape[1]); tie = np.ones(V.shape[1])
    for g in range(4):
        mask = yv == g
        H += R[mask].sum(axis=0)**2 / mask.sum()
    H = 12/(N*(N+1)) * H - 3*(N+1)
    # tie correction per column: 1 - sum(t^3-t)/(N^3-N)
    for j in range(V.shape[1]):
        _, counts = np.unique(V[:,j], return_counts=True)
        tie[j] = 1 - (counts**3 - counts).sum()/(N**3 - N)
    Hc = np.where(tie > 0, H/np.maximum(tie, 1e-12), 0.0)
    p = chi2.sf(Hc, 3)
    return p

def medians_by_group(V, yv):
    return np.stack([np.median(V[yv==g], axis=0) for g in range(4)])

def cindex(yy, s):
    c = t = 0
    for i, j in itertools.combinations(range(len(yy)), 2):
        if yy[i] == yy[j]: continue
        t += 1; c += (s[i]-s[j])*(yy[i]-yy[j]) > 0
    return c/t if t else float('nan')

def oof_scores(yv, folds):
    oof = np.zeros(len(yv)); member_counts = []
    for tr, te in folds:
        p = kw_screen(L[tr], yv[tr])
        o = np.argsort(p); bh = p[o]*len(p)/(np.arange(len(p))+1)
        fdr = np.empty_like(bh); fdr[o] = np.minimum.accumulate(bh[::-1])[::-1]
        gm = medians_by_group(L[tr], yv[tr])
        d = gm[3] - gm[1]
        sel = (fdr <= 0.05) & (np.abs(d) >= 0.5) & ~excl_mask & core6_mask
        member_counts.append(int(sel.sum()))
        if sel.sum() == 0:
            continue
        w = np.sign(d[sel])
        mu = L[tr][:, sel].mean(axis=0); sd = L[tr][:, sel].std(axis=0)
        sd = np.where(sd == 0, 1.0, sd)
        oof[te] = (w * ((L[te][:, sel] - mu)/sd)).sum(axis=1)
    return oof, member_counts

t0 = time.time()
folds_true = list(StratifiedKFold(5, shuffle=True, random_state=SEED).split(L, y))
oof, mcounts = oof_scores(y, folds_true)
c_obs = cindex(y, oof)
kw_obs = kruskal(*[oof[y==g] for g in range(4)]).pvalue
print(f'observed OOF c-index {c_obs:.4f}, KW p {kw_obs:.3g}, members/fold {mcounts}', flush=True)

rng = np.random.default_rng(SEED)
null = np.zeros(NPERM)
for b in range(NPERM):
    yp = rng.permutation(y)
    folds_p = list(StratifiedKFold(5, shuffle=True, random_state=SEED).split(L, yp))
    oofp, _ = oof_scores(yp, folds_p)
    null[b] = cindex(yp, oofp)
    if (b+1) % 100 == 0:
        print(f'perm {b+1}/{NPERM} ({time.time()-t0:.0f}s)', flush=True)
emp_p = (1 + (null >= c_obs).sum())/(NPERM+1)
out = {'oof_cindex': c_obs, 'oof_kw_p': float(kw_obs), 'per_fold_members': mcounts,
       'n_perm': NPERM, 'perm_p': float(emp_p),
       'null_median': float(np.median(null)), 'null_q975': float(np.quantile(null, 0.975)),
       'in_sample_cindex_reference': 0.7557,
       'spec': 'ADDENDUM_4 D2 locked 2026-09-27T08:58; claim stands only if perm_p<0.05',
       'runtime_s': round(time.time()-t0, 1)}
json.dump(out, open('results/d2_oof_module_reconstruction.json','w'), indent=2)
print(json.dumps(out, indent=2))
