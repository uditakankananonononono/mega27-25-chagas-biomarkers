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
from concurrent.futures import ThreadPoolExecutor
def fetch(mir):
    q = urllib.parse.quote(f'TITLE_ABS:"{mir}"')
    url = f'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={q}&format=json'
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=25) as r:
                return mir, json.loads(r.read())['hitCount']
        except Exception:
            time.sleep(2*(attempt+1))
    return mir, None
done = 0
with ThreadPoolExecutor(max_workers=8) as ex:
    for mir, cnt in ex.map(fetch, names):
        out[mir] = cnt
        done += 1
        if done % 250 == 0:
            print(f'{done}/{len(names)}', flush=True)
import os
os.makedirs('sources/services/europepmc', exist_ok=True)
raw = json.dumps(out, sort_keys=True).encode()
open('sources/services/europepmc/mirna_pubcounts.json','wb').write(raw)
print('sha256', hashlib.sha256(raw).hexdigest(), 'n', len(out),
      'nulls', sum(1 for v in out.values() if v is None))
