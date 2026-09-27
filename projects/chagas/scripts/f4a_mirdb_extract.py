#!/usr/bin/env python3
"""F4 extraction step (committed for blind reproducibility): from the
pinned miRDB v6.0 file, keep the 20 candidates' human predictions at
score >= 80 -> sources/services/mirdb/mirdb_candidates_ge80.json.
Regenerate with: python3 scripts/f4a_mirdb_extract.py (needs the 59MB
pin present)."""
import gzip, json, csv
cands = sorted(set(r['mirna'] for r in csv.DictReader(open('results/h2_gate_b_enrichment.csv'))))
want = set('hsa-' + c for c in cands)
hits = {('hsa-' + c): [] for c in cands}
with gzip.open('sources/services/mirdb/miRDB_v6.0_prediction_result.txt.gz', 'rt', errors='replace') as f:
    for line in f:
        m = line.split('\t', 1)[0]
        if m in want:
            p = line.rstrip('\n').split('\t')
            if float(p[2]) >= 80.0:
                hits['hsa-' + c].append((p[1], float(p[2])))
json.dump(hits, open('sources/services/mirdb/mirdb_candidates_ge80.json', 'w'))
print('extracted', sum(len(v) for v in hits.values()))
