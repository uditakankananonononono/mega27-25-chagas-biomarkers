#!/usr/bin/env python3
"""Re-execution revalidation of the live-service CORE6 analyses (2026-09-27).

Re-runs, against the live services, the queries behind three committed
results files and compares headline numbers to the committed values:
  - mygene.info symbol validation   (results/mygene_core6_symbol_validation.json)
  - Open Targets tractability overlay (results/druggability_overlay.json)
  - STRING network + PPI enrichment (results/string_core6_*.tsv)

Input symbol set: CORE6 miRNAs' strong-support ("Functional MTI") unique
targets from sources/services/mirtarbase/validated_targets.tsv (352 symbols;
asserted equal to the committed mygene validation query set).

Output: results/service_revalidation_2026-09-27.json with per-service
status = match | drift | error, committed vs observed values. Raw responses
under sources/services/revalidation_2026-09-27/. This is a re-USE of
already-ledgered services (mygene=36, Open Targets=8/33, STRING=10/34),
not new tools for the 40-service gate.
"""
import json, sys, time, urllib.request, urllib.parse
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUTDIR = ROOT / 'sources/services/revalidation_2026-09-27'
OUTDIR.mkdir(parents=True, exist_ok=True)
CORE6 = ['hsa-miR-1-3p','hsa-miR-122-5p','hsa-miR-192-5p','hsa-miR-30c-5p','hsa-miR-145-5p','hsa-miR-194-5p']

t = pd.read_csv(ROOT/'sources/services/mirtarbase/validated_targets.tsv', sep='\t')
symbols = sorted(t[(t.mirna.isin(CORE6)) & (t.support=='Functional MTI')].target_gene.unique())
assert len(symbols) == 352, len(symbols)
committed_mg = json.load(open(ROOT/'results/mygene_core6_symbol_validation.json'))
assert sorted(r['query'] for r in committed_mg) == symbols, 'mygene query set mismatch'

report = {'date': '2026-09-27', 'input_symbols': 352, 'services': {}}

def post_json(url, data, timeout=30):
    req = urllib.request.Request(url, data=urllib.parse.urlencode(data).encode(),
        headers={'User-Agent': 'chagas-biomarker-revalidation/1.0'})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode()

def get_json(url, timeout=60):
    req = urllib.request.Request(url, headers={'User-Agent': 'chagas-biomarker-revalidation/1.0'})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode()

# --- A) mygene batch validation --------------------------------------------
try:
    raw = post_json('https://mygene.info/v3/query', {
        'q': ','.join(symbols), 'scopes': 'symbol',
        'fields': 'symbol,entrezgene,taxid', 'species': 'human'})
    (OUTDIR/'mygene_query.json').write_text(raw)
    hits = json.loads(raw)
    matched = [h for h in hits if '_id' in h]
    notfound = sorted(h['query'] for h in hits if '_id' not in h)
    committed_nf = sorted(r['query'] for r in committed_mg if '_id' not in r)
    report['services']['mygene'] = {
        'status': 'match' if (len(matched)==344 and notfound==committed_nf) else 'drift',
        'committed': {'matched': 344, 'notfound': committed_nf},
        'observed': {'matched': len(matched), 'notfound': notfound}}
except Exception as e:
    report['services']['mygene'] = {'status': 'error', 'error': repr(e)}
print('mygene done:', report['services']['mygene']['status'], flush=True)

# --- B) Open Targets tractability overlay ----------------------------------
try:
    m = pd.read_csv(ROOT/'sources/services/ensembl_biomart/ensg_symbol_map.tsv', sep='\t')
    m.columns = ['ensg', 'symbol']
    sym2ensg = dict(zip(m.symbol, m.ensg))
    mapped = {s: sym2ensg[s] for s in symbols if s in sym2ensg}
    q = ('query($id: String!){ target(ensemblId:$id){ approvedSymbol '
         'tractability { modality label value } } }')
    tp, ad, ot_seen = {}, [], 0
    for i, (sym, ensg) in enumerate(sorted(mapped.items())):
        body = json.dumps({'query': q, 'variables': {'id': ensg}}).encode()
        req = urllib.request.Request('https://api.platform.opentargets.org/api/v4/graphql',
            data=body, headers={'Content-Type': 'application/json',
                                'User-Agent': 'chagas-biomarker-revalidation/1.0'})
        with urllib.request.urlopen(req, timeout=30) as r:
            resp = json.loads(r.read().decode())
        tgt = resp.get('data', {}).get('target')
        if tgt:
            ot_seen += 1
            buckets = sorted(x['label'] for x in tgt['tractability'] if x['value'])
            if buckets:
                tp[sym] = buckets
                if 'Approved Drug' in buckets:
                    ad.append(sym)
        if i % 50 == 0:
            print(f'OT {i}/{len(mapped)}', flush=True)
        time.sleep(0.05)
    (OUTDIR/'opentargets_tractability.json').write_text(json.dumps(tp))
    committed_ov = json.load(open(ROOT/'results/druggability_overlay.json'))
    ok = (len(mapped)==344 and len(tp)==332 and len(ad)==49
          and set(tp)==set(committed_ov['tractability_positive'])
          and sorted(ad)==sorted(committed_ov['approved_drug']))
    report['services']['opentargets'] = {
        'status': 'match' if ok else 'drift',
        'committed': {'mapped': 344, 'tractability_positive': 332, 'approved_drug': 49},
        'observed': {'mapped': len(mapped), 'ot_seen': ot_seen,
                     'tractability_positive': len(tp), 'approved_drug': len(ad)}}
except Exception as e:
    report['services']['opentargets'] = {'status': 'error', 'error': repr(e)}
print('opentargets done:', report['services']['opentargets']['status'], flush=True)

# --- C) STRING network + PPI enrichment ------------------------------------
try:
    idlist = '%0d'.join(sorted(mapped))  # same 344-symbol input set
    base = 'https://string-db.org/api/tsv'
    net = get_json(f'{base}/network?identifiers={idlist}&species=9606'
                   f'&required_score=700&caller_identity=chagas_revalidation')
    (OUTDIR/'string_network.tsv').write_text(net)
    enr = get_json(f'{base}/ppi_enrichment?identifiers={idlist}&species=9606'
                   f'&required_score=700&caller_identity=chagas_revalidation')
    (OUTDIR/'string_ppi_enrichment.tsv').write_text(enr)
    n_edges = max(0, len(net.strip().splitlines()) - 1)
    ep = pd.read_csv(__import__('io').StringIO(enr), sep='\t').iloc[0]
    ok = (n_edges==1051 and int(ep['number_of_nodes'])==299
          and int(ep['expected_number_of_edges'])==418 and float(ep['p_value'])<1e-16)
    report['services']['string'] = {
        'status': 'match' if ok else 'drift',
        'committed': {'edges': 1051, 'nodes': 299, 'expected_edges': 418, 'p': '<1e-16'},
        'observed': {'edges': n_edges, 'nodes': int(ep['number_of_nodes']),
                     'expected_edges': int(ep['expected_number_of_edges']),
                     'p': float(ep['p_value'])}}
except Exception as e:
    report['services']['string'] = {'status': 'error', 'error': repr(e)}
print('string done:', report['services']['string']['status'], flush=True)

out = ROOT/'results/service_revalidation_2026-09-27.json'
out.write_text(json.dumps(report, indent=2))
print('WROTE', out, flush=True)
