#!/usr/bin/env python3
"""F7 per ADDENDUM_6 (compute-only, no causal claims): prioritization
ranking of the 20 frozen-gate candidates. Four sub-scores, all from
committed artifacts: (1) cohort direction-consistency (gate-b FDR passes,
results/h2_gate_b_enrichment.csv); (2) network centrality (committed
CORE6 STRING graph: mean target degree of member's strong-support
targets; computable for CORE6 only - disclosed NA elsewhere);
(3) tissue agreement (E2 DirectionScore, results/e2_*.json);
(4) disease-pathway proximity (fraction of strong-support targets in the
pinned-library Chagas/cardiac term union, E13 sets). Composite = mean of
available sub-scores. Heuristic only; causal wording banned."""
import json, re, csv
import pandas as pd

CORE6 = ['miR-1-3p','miR-122-5p','miR-192-5p','miR-30c-5p','miR-145-5p','miR-194-5p']

# candidate list + serum direction from gate-a
ga = pd.read_csv('results/h2_gate_a_passing.csv')
ga['bare'] = ga.mirna.str.replace('hsa-','',regex=False)
cands = ga.bare.tolist()  # 21 rows: 20 candidates + replication class? filter to gate-b tested
gb = pd.read_csv('results/h2_gate_b_enrichment.csv')
tested = sorted(gb.mirna.unique())
assert len(tested) == 20, len(tested)

# sub-score 1: cohort pass fraction
s1 = {}
for m in tested:
    sub = gb[gb.mirna == m]
    s1[m] = float((sub.fdr <= 0.05).sum()) / len(sub)

# sub-score 2: STRING centrality (CORE6 only)
net = pd.read_csv('results/string_core6_network.tsv', sep='\t')
from collections import Counter
deg = Counter()
for _, r in net.iterrows():
    deg[r.preferredName_A.upper()] += 1
    deg[r.preferredName_B.upper()] += 1
tg = pd.read_csv('sources/services/mirtarbase/validated_targets.tsv', sep='\t')
tg['bare'] = tg.mirna.str.upper().str.replace('HSA-','',regex=False)
strong = tg[tg.support == 'Functional MTI']
s2raw = {}
for m in CORE6:
    genes = set(g.upper() for g in strong[strong.bare == m.upper()].target_gene.unique())
    s2raw[m] = sum(deg[g] for g in genes) / len(genes) if genes else 0.0
mx = max(s2raw.values()) or 1.0
s2 = {m: (s2raw[m]/mx if m in s2raw else None) for m in tested}

# sub-score 3: E2 DirectionScore
e2 = json.load(open('results/e2_cross_tissue_concordance.json'))
ds = dict(e2['candidates_ds']); ds.update(e2['core6_direction_scores'])
s3 = {m: ds.get(m) for m in tested}

# sub-score 4: disease-pathway proximity (E13 pinned-library disease terms)
pat = re.compile(r'chagas|cardiomyopath|heart failure|cardiac|myocard|fibros|hypertroph', re.I)
dunion = set()
for libf in ['KEGG_2021_Human','WikiPathways_2024_Human','GO_Biological_Process_2025']:
    libj = json.load(open(f'results/enrichr/{libf}.json'))
    for row in libj.get(libf, libj):
        if isinstance(row, list) and len(row) > 5 and isinstance(row[5], list):
            term, genes = row[1], row[5]
        elif isinstance(row, dict):
            term, genes = row.get('term',''), row.get('geneSymbols',[])
        else:
            continue
        if pat.search(term):
            dunion.update(x.upper() for x in genes)
s4raw = {}
for m in tested:
    genes = set(g.upper() for g in strong[strong.bare == m.upper()].target_gene.unique())
    s4raw[m] = len(genes & dunion) / len(genes) if genes else 0.0
mx4 = max(s4raw.values()) or 1.0
s4 = {m: s4raw[m]/mx4 for m in tested}

rows = []
for m in tested:
    parts = [s1[m]] + ([s2[m]] if s2[m] is not None else []) + \
            ([s3[m]] if s3[m] is not None else []) + [s4[m]]
    rows.append({
        'mirna': m,
        'cohort_pass_frac': round(s1[m], 3),
        'string_centrality_norm': (round(s2[m], 3) if s2[m] is not None else 'NA'),
        'e2_direction_score': (round(s3[m], 3) if s3[m] is not None else 'NA'),
        'disease_proximity_norm': round(s4[m], 3),
        'n_subscores': len(parts),
        'composite': round(sum(parts)/len(parts), 4),
        'core6': m in CORE6,
    })
R = pd.DataFrame(rows).sort_values('composite', ascending=False).reset_index(drop=True)
R.insert(0, 'rank', R.index + 1)
R.to_csv('results/f7_prioritization_ranking.csv', index=False)
json.dump({'spec': 'ADDENDUM_6 F7; prioritization heuristic only; causal wording banned',
           'disease_union_genes': len(dunion), 'string_nodes': len(deg),
           'top5': R.head(5).to_dict('records')},
          open('results/f7_prioritization_ranking.json','w'), indent=2)
print(R.to_string(index=False))
