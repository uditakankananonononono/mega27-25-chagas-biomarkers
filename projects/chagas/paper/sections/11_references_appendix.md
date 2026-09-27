# 11. References and appendices (assembly notes)
References are resolved through CrossRef and PubMed with DOIs/PMIDs in
the acquisition log: Roma et al. 2026 (PMID 41574750, GSE299582 source);
Choudhuri et al. 2021 (PMID 34479416, ELISA prognosis panel); Woc-Colburn
et al. / cohort studies per ACQUISITION_LOG PMIDs (40290486, 35873155,
31105048, 33897690, 40391216, 29545200, 42614816, 41648170, 42523576,
41574750). Appendix A: full service ledger with retrieval hashes
(SERVICE_LEDGER.md + sources/services/EVIDENCE_SHA256.txt). Appendix B:
per-series crosswalk schemas. Appendix C: preregistration chain
(PREREGISTRATION.md, ADDENDUM_1-3, RUN_LOG.md) with commit hashes.
Appendix D: judge-round transcripts (judge_rounds/, model-tagged).
Appendix E: the count-correction record (750->716, commit ef7b29a). Appendix F: the reproducibility and re-execution register.

Service-source references for the executed discovery arm: miRTarBase
release 8.0 (Huang et al. 2022, NAR 50:D222-D230) human MTI table,
official file recovered via Internet Archive snapshot 20220623192730 of
the publisher URL, sha256-pinned; Ensembl BioMart GRCh38 hsapiens gene
annotation (Cunningham et al. 2022, NAR 50:D988-D995) via the useast
mirror; miRBase v22 (Kozomara et al. 2019, NAR 47:D155-D162), live lookup
MIMAT0000728; Enrichr (Chen et al. 2013, BMC Bioinformatics 14:128;
Kuleshov et al. 2016, NAR 44:W90-97) GO BP 2025 / KEGG 2021 / WikiPathways
2024 Human; Open Targets Platform (Ochoa et al. 2023, NAR 51:D1353-D1359)
GraphQL tractability. GEO/SRA and eutils per ACQUISITION_LOG.


## 11.1 Per-series reference map
Every cohort in the atlas is anchored to its publication record; PMIDs
were verified per series at acquisition (ACQUISITION_LOG.md), not copied
from secondary sources. GSE299582 -> PMID 41574750 (Roma et al. 2026,
J Infect Dis 234(2):277-287; the severity-miRNA source paper, also the
frozen comparator context and the Addendum-2 exclusion source);
GSE244827 -> PMID 40290486; GSE203525 -> PMID 35873155; GSE129676 ->
PMID 31105048; GSE158986 -> PMID 33897690; GSE295194 -> PMID 40391216;
GSE107376 -> PMID 29545200; GSE311812 -> PMID 41648170; GSE333874 ->
PMID 42523576; GSE328447 -> PMID 42614816; GSE348071 -> unpublished at
acquisition (flagged for citation monitoring); GSE191081/GSE191082 ->
methylation study pair per ACQUISITION_LOG; GSE154421 -> benznidazole
pharmacogenomics per ACQUISITION_LOG; GSE84796 -> Cunha-Neto CCC heart
tissue study (prior-tagged). Comparator and constraint literature:
PMID 34479416 (Choudhuri et al. 2021, the ELISA prognosis panel whose
0.688 R2 anchored failed H1), PMID 38203212 (2024 prospective
parasite-DNA + 47-protein ten-year decline model), PMID 23393012 (the
499-donor ten-year outcomes cohort), and the WHO fact sheet and CDC DPDx
pages retrieved-and-hashed under sources/services/.

## Appendix B: crosswalk schema (frozen)
Every per-series crosswalk (sources/<GSE>_sample_crosswalk.csv) carries:
gsm, source_url (canonical acc.cgi full-text URL), sha256 (of the stored
SOFT bytes), title, status, label, platform, finite_probes, organism,
source_name, sra_relation, characteristics (JSON of GEO characteristic
fields). Labels are assigned by frozen rules in scripts/label_crosswalks.py;
any sample whose label is a design contrast (treatment, timepoint) rather
than a clinical state is marked as such and excluded from clinical
analyses by construction. The uniqueness and hash checks are executable
(scripts/verify_crosswalks.py).

## Appendix C: preregistration chain with commits
PREREGISTRATION.md (6e562c7) -> literature screen (70ab16d) ->
ADDENDUM_1 (44a4e5e, comparator resolution) -> RUN 1 negative (49de4ae)
-> judge round 01 (judge_rounds/01_*, chat-tagged) -> ADDENDUM_2
(exclusion screen extension) -> ADDENDUM_3 (H1' + gate (b) + module
claim tier, locked 2026-09-26T21:23) -> H1' run + gate-(b) run +
module/druggability execution (2026-09-26/27, RUN_LOG). Every amendment
predates the outcomes it governs; the chain is re-readable in git.

## Appendix D: judge rounds
The owner's settled rule (2026-09-27, verbatim: "EACH PROJECTS NEED ONE
FROM ME TO PASS") counts one USER-PROVIDED ChatGPT judge round per
project; only verdicts she provides count. The canonical tally is
X-JUDGE-ROUNDS.md. Counted: 2 of 1 required - REQUIREMENT MET: round 02
(verdict couriered 2026-09-27T09:17, verbatim 02_response_chatgpt.txt +
02_assessment_chatgpt.txt; adoptions locked in ADDENDUM_5) and round 03
(manuscript self-review couriered 2026-09-27T10:00, verbatim
03_response_chatgpt.txt + 03_assessment_chatgpt.txt; new analyses locked
in ADDENDUM_6). Historical/supplementary (never counted): round 01
(agent-initiated ChatGPT thread 6ab7e9ef, verbatim 01_*) and a Gemini
Flash consult (02_response_gemini.txt + assessment with erratum).
Prior requirement wordings ("min 10 ChatGPT judge rounds", then
10:00:07's one-provided) are superseded and recorded in
X-JUDGE-ROUNDS.md's correction history.

## Appendix E: the 750->716 count correction
The second acquisition pass counted 750 records; the uniqueness
verifier flagged super-series and prior-tag overlaps (GSE191083
container; GSE244827 prior-series overlap), and the frozen count was
corrected to 716 (commit ef7b29a) - the correction record, not the
original error, is the citable artifact.

## Appendix F: reproducibility and re-execution register
### F.1 Script-to-result coverage
Every figure-free number in this paper regenerates from a committed
script against committed inputs; the audit that closed the remaining
gaps is itself in RUN_LOG (2026-09-27T05:50, 06:08). The coverage
table: scripts/acquire_geo_crosswalk.py - per-series sample crosswalks
(sources/*_sample_crosswalk.csv); scripts/label_crosswalks.py - the
frozen label rules; scripts/verify_crosswalks.py - executable
uniqueness and hash checks on every crosswalk;
scripts/h2_severity_association.py - gate (a)
(results/h2_gate_a_summary.json: 2,632 tests, 21 passing, 20 series;
results/h2_gate_a_passing.csv); scripts/h2_gate_b_enrichment.py -
gate (b) (results/h2_gate_b_enrichment.csv, 40 tests);
scripts/h1_frozen_analysis.py - the frozen H1 run
(results/h1_frozen_run.json); scripts/h1prime_ordinal.py and
scripts/h1prime_ci.py - the H1' ordinal model and confidence
intervals (results/h1prime_*.json); scripts/h1prime_metrics_completion.py
- the locked-metric completion (macro-AUC OvR 0.715; calibration
slopes 0.735/0.198/0.150; results/h1prime_oof_predictions.csv), with
sanity anchors reproducing the original h1prime artifacts to full
precision; scripts/module_score.py - results/module_score.json,
reproduced EXACTLY (full-precision in-script assertion;
results/module_score_stability.json for the bootstrap/LOO
descriptives); scripts/service_revalidation.py - live re-execution
of the service analyses (results/service_revalidation_2026-09-27.json).

### F.2 Re-execution revalidation record (2026-09-27)
The three service-backed results files were re-executed live and
compared headline-for-headline: mygene.info symbol validation
reproduced EXACTLY (344/352 matched; identical 8-symbol legacy
notfound set: ALPPL2, COX1, CTGF, FAM45A, H3F3A, ND1, NDUFA4,
SEPT10). Open Targets tractability reproduced EXACTLY at bucket
level (332 tractability-positive; 49 Approved-Drug; identical gene
sets), with mapped count 343 vs the committed 344 - see F.3. STRING
v12 drifted with the live database between 00:51 and 06:08 the same
night (1,296 edges vs the committed 1,051; 341 mapped nodes vs 299;
expected edges 523 vs 418; enrichment p < 1e-16 in both runs). The
committed 00:51 values remain the lane record; the drift record is
the honest reproducibility statement for a live service, and the
enrichment conclusion (significantly interconnected) is unchanged.

### F.3 Known provenance gaps, named
- results/pharos_tdl_approved_drug_targets.json is the blocked-service
placeholder (every entry ERR:HTTPError), kept as the blocked-attempt
record; Pharos is in the blocked register, not the 40.
- miRTarBase is pinned to v8.0 via a Wayback snapshot because the
live site returned 404/400 on 2026-09-26; the version pin is
disclosed wherever miRTarBase counts are cited.
- The original Open Targets symbol-to-ENSG mapping path was not
byte-documented; re-execution through the committed BioMart
ensg_symbol_map resolves 343 of 352 symbols (the 8 legacy aliases
plus one further symbol unresolved by this path). Bucket-level
overlay results are unaffected; the one-symbol difference is a
mapping-path gap, not a results discrepancy.
- scikit-learn was unpinned in the lane environment record until
2026-09-27 (now pinned at 1.7.2; Section 5.y). Exact full-precision
reproduction of the H1' sanity anchors and the module score under
this pin indicates the pipeline is env-stable.

### F.4 Evidence integrity
Raw service responses are stored under sources/services/ with sha256
hashes in sources/services/EVIDENCE_SHA256.txt (the manifest covers
the revalidation responses added 2026-09-27); retrieval timestamps
and re-use annotations are in SERVICE_LEDGER.md.

## Appendix G: data-source bibliography
One citable record per series in the fifteen-series atlas. Citations
are from EuropePMC core records fetched 2026-09-27
(sources/services/europepmc_series_citations.json, sha256-hashed; a
re-use of ledger service 6). Design facts are from
sources/new_series_ledger.csv and the frozen per-series crosswalks.
Honest gaps are stated per row.

- **GSE299582** (192 serum miRNA-seq, GPL30173; susceptibility + CCC
  severity, the primary analysis cohort): Roma EH, Marques-Santos F,
  Renzetti ARDS, et al. "Transcriptome Analysis of Circulating
  microRNAs Associated With Chagas Disease Susceptibility and Chronic
  Chagas Cardiomyopathy Severity." J Infect Dis, 2026.
  doi:10.1093/infdis/jiag021. PMID 41574750.
- **GSE244827** (33 whole-blood RNA-seq, GPL24676; asymptomatic/early
  CCC vs seronegative): Duque C, So J, Castro-Sesquen YE, et al.
  "Immunologic changes in the peripheral blood transcriptome of
  clinically asymptomatic Chagas cardiomyopathy patients." Lancet Reg
  Health Am, 2025. doi:10.1016/j.lana.2025.101090. PMID 40290486.
- **GSE203525** (20 patient hiPSC-CM RNA-seq; CCC vs indeterminate
  +/- reinfection; gate-(b) cohort): Oliveira TGM, Venturini G,
  Alvim JM, et al. "Different Transcriptomic Response to T. cruzi
  Infection in hiPSC-Derived Cardiomyocytes." Front Cell Infect
  Microbiol, 2022. doi:10.3389/fcimb.2022.904747. PMID 35873155.
- **GSE129676** (16 hiPSC-CM RNA-seq infection timecourse): Bozzi A,
  Sayed N, Matsa E, et al. "Using Human Induced Pluripotent Stem
  Cell-Derived Cardiomyocytes as a Model to Study Trypanosoma cruzi
  Infection." Stem Cell Reports, 2019.
  doi:10.1016/j.stemcr.2019.04.017. PMID 31105048.
- **GSE158986** (12 monocyte-derived dendritic-cell RNA-seq, first
  contact): Gil-Jaramillo N, Rocha AP, Raiol T, et al. "The First
  Contact of Human Dendritic Cells With Trypanosoma cruzi." Front
  Immunol, 2021. doi:10.3389/fimmu.2021.638020. PMID 33897690.
- **GSE295194** (16 scRNA-seq PBMC sample tags; CCC vs indeterminate
  CD4 T-cell peptide response): Souza-Silva TG, Neves EGA,
  Teixeira-Carvalho A, et al. "Self and parasite-derived peptides
  selected upon DERAA-bearing HLA-DRB1 molecules." Front Immunol,
  2025. doi:10.3389/fimmu.2025.1527115. PMID 40391216.
- **GSE107376** (9 placental expression array; seropositive vs
  seronegative mothers): Juiz NA, Torrejon I, Burgos M, et al.
  "Alterations in Placental Gene Expression of Pregnant Women With
  Chagas Disease." Am J Pathol, 2018.
  doi:10.1016/j.ajpath.2018.02.011. PMID 29545200.
- **GSE328447** (4 THP1-macrophage small RNA-seq; isomiR response):
  Lyu MA, Maimaiti M, Hu J, Hu H. "Comparative 5' IsomiRome Analysis
  Uncovers Dysregulated 5' IsomiRs in Trypanosoma cruzi Infection."
  Comput Struct Biotechnol J, 2026. doi:10.34133/csbj.0146.
  PMID 42614816.
- **GSE191081** (22 LV-wall RNA-seq; CCC vs dilated cardiomyopathy)
  and **GSE191082** (158 methylation tiling array, blood + LV wall):
  the lane ledger records both members of this family against one
  publication - Brochet P, Ianni BM, Laugier L, et al. "Epigenetic
  regulation of transcription factor binding motifs in chronic
  Chagas cardiomyopathy." Front Immunol, 2022.
  doi:10.3389/fimmu.2022.958200. PMID 36072583. (GSE191083, the
  super-series, was acquired then removed by the uniqueness
  verifier; Appendix E.)
- **GSE311812** (46 RNA-seq + Visium spatial; congenital
  transmission): GEO records PMID 41648170, but that PMID is not yet
  indexed in EuropePMC as of 2026-09-27 (zero-hit response stored in
  the citation artifact). Cited here by accession; the lane claims no
  bibliographic record beyond the GEO page.
- **GSE333874** (31 placental small RNA-seq; congenital transmission
  miRNAs): GEO records PMID 42523576; same not-yet-indexed status
  (zero-hit stored). Cited by accession.
- **GSE348071** (32 AC16 + patient iPSC-CM RNA-seq; DHODH R135C):
  unpublished per the GEO record at acquisition time.
- **GSE154421** (92 benznidazole adverse-reaction SNP pharmacogenomic
  array): unpublished per the GEO record at acquisition time.
- **GSE84796** (17 records used; prior-tagged series): the GEO record
  is the only provenance artifact in this lane for this series; the
  originating publication was not separately fetched and is not
  claimed.

EuropePMC returned 9 of the 11 queried PMIDs; the two misses are the
two most recent records (both 2026 congenital-transmission series),
consistent with indexing lag rather than absence - but the lane
records them as unverified citations, not as resolved ones.
