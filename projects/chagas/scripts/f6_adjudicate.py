import json, re, time, urllib.request, urllib.parse

d = json.load(open('/tmp/f6_expanded_novelty.json'))
out = {}
def fetch_abs(pmid):
    url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?' + urllib.parse.urlencode(
        {'query': f'EXT_ID:{pmid} AND SRC:MED', 'format': 'json', 'resultType': 'core'})
    try:
        r = json.load(urllib.request.urlopen(url, timeout=30))
        res = r.get('resultList', {}).get('result', [])
        return (res[0].get('title','') + ' ' + res[0].get('abstractText','')) if res else ''
    except Exception as e:
        print('fetchfail', pmid, e, flush=True); return ''

for c, v in d['per_candidate'].items():
    pmids = v['tiab_pmids'][:25]
    exact = []
    # exact token: miR-199b-5p or hsa-miR-199b-5p, boundaries; also arm-less exact "miR-199b" NOT counted (family ambiguity)
    pat = re.compile(r'\b(?:hsa-)?' + re.escape(c) + r'\b', re.I)
    for pmid in pmids:
        if not pmid: continue
        text = fetch_abs(pmid)
        if pat.search(text):
            chagas = bool(re.search(r'chagas|cruzi', text, re.I))
            exact.append({'pmid': pmid, 'chagas_context': chagas})
        time.sleep(0.2)
    out[c] = exact
    print(c, len(exact), [e['pmid'] for e in exact if e['chagas_context']], flush=True)
json.dump(out, open('/tmp/f6_adjudication.json','w'), indent=1)
print('DONE')
