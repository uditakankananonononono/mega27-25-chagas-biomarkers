#!/usr/bin/env python3
"""E2 support retrieval (ADDENDUM_5 E2, locked 2026-09-27T09:23): per-miRNA
publication counts for the 2,632-serum-panel matching null. Locked fallback
was NCBI PubMed esearch; the esearch backend was DOWN at run time
(500/Cannot connect to SOLR, 2026-09-27 ~09:39 IST, retried) so the counts
come from Europe PMC (ledgered service) TITLE_ABS queries over the same
PubMed corpus - DISCLOSED SUBSTITUTION. Cached JSON + sha256."""
import gzip, json, time, hashlib, urllib.request, urllib.parse
import pandas as pd
with gzip.open('sources/matrices/GSE299582_normalized_counts.csv.gz','rt') as f:
    M = pd.read_csv(f, index_col=0)
names = [m.replace('hsa-','') for m in M.index]
out = {}
for i, mir in enumerate(names):
    q = urllib.parse.quote(f'TITLE_ABS:"{mir}"')
    url = f'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={q}&format=json'
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=20) as r:
                out[mir] = json.loads(r.read())['hitCount']
            break
        except Exception as e:
            time.sleep(2*(attempt+1))
    else:
        out[mir] = None
    if (i+1) % 250 == 0:
        print(f'{i+1}/{len(names)}', flush=True)
        time.sleep(1)
    else:
        time.sleep(0.15)
raw = json.dumps(out, sort_keys=True).encode()
open('sources/services/europepmc/mirna_pubcounts.json','wb').write(raw)
print('sha256', hashlib.sha256(raw).hexdigest(), 'n', len(out),
      'nulls', sum(1 for v in out.values() if v is None))
