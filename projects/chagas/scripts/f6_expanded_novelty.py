import json, time, urllib.request, urllib.parse, csv, re

def epmc(path, params, tries=3):
    url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/' + path + '?' + urllib.parse.urlencode(params)
    for a in range(tries):
        try:
            return json.load(urllib.request.urlopen(url, timeout=40))
        except Exception as e:
            print('retry', path, a, e, flush=True); time.sleep(4)
    return {}

cands = sorted(set(r['mirna'] for r in csv.DictReader(open('/home/sandbox/work/repo/projects/chagas/results/h2_gate_b_enrichment.csv'))))
out = {'per_candidate': {}, 'citation_expansion': {}}

# Stage 1: per-candidate Chagas searches - title/abstract and open full text
for c in cands:
    bare = c  # e.g. miR-1-3p
    q_ta = f'("{bare}" OR "hsa-{bare}") AND Chagas'
    r1 = epmc('search', {'query': q_ta, 'format': 'json', 'pageSize': 25})
    hits1 = r1.get('resultList', {}).get('result', [])
    q_ft = f'("{bare}" OR "hsa-{bare}") AND Chagas AND OPEN_ACCESS:Y'
    r2 = epmc('search', {'query': q_ft, 'format': 'json', 'pageSize': 25})
    hits2 = r2.get('resultList', {}).get('result', [])
    out['per_candidate'][c] = {
        'tiab_count': r1.get('hitCount', 0),
        'tiab_pmids': [h.get('pmid', h.get('id','')) for h in hits1],
        'tiab_titles': [h.get('title','') for h in hits1],
        'oa_count': r2.get('hitCount', 0),
        'oa_pmids': [h.get('pmid', h.get('id','')) for h in hits2],
    }
    print('cand', c, out['per_candidate'][c]['tiab_count'], out['per_candidate'][c]['oa_count'], flush=True)
    time.sleep(0.4)

# Stage 2: citation expansion of the 121-seed
seed = json.load(open('/home/sandbox/work/repo/projects/chagas/prereg/pubmed_screen_ids.json'))['esearchresult']['idlist']
pat = re.compile(r'\b(' + '|'.join(re.escape(c) for c in cands) + r')\b', re.I)
citing = {}
for pmid in seed:
    r = epmc(f'MED/{pmid}/citations', {'format': 'json', 'pageSize': 1000})
    for cit in r.get('citationList', {}).get('citation', []):
        cid = cit.get('pmid', cit.get('id', ''))
        if cid and cid not in citing:
            citing[cid] = (cit.get('title', '') + ' ' + cit.get('abstractText', ''))
    time.sleep(0.25)
print('citing corpus:', len(citing), flush=True)
exp_hits = {c: [] for c in cands}
for cid, text in citing.items():
    for m in pat.finditer(text):
        mm = m.group(1)
        for c in cands:
            if c.lower() == mm.lower():
                exp_hits[c].append(cid)
out['citation_expansion'] = {
    'seed_pmids': len(seed), 'citing_papers_screened': len(citing),
    'hits_per_candidate': {c: sorted(set(v)) for c, v in exp_hits.items() if v},
}
json.dump(out, open('/tmp/f6_expanded_novelty.json', 'w'), indent=1)
print('DONE', flush=True)
