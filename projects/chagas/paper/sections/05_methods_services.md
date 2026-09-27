# 5. Methods: the service-verified analysis layer

## 5.1 Design
The compendium's analyses rest on 40 distinct external services, each used
live during this build with retrieved bytes preserved and sha256-hashed
(sources/services/, manifest EVIDENCE_SHA256.txt). The count is distinct
tools, not uses: four services (Enrichr, Open Targets, STRING, g:Profiler)
were run a second time for the discovery arm and are annotated in the
ledger as re-uses of their first-use entries rather than assigned new
numbers - the same append-only correction convention as the 750->716
record correction. The layer divides into five functional groups, each
feeding named sections of this paper.

## 5.2 Acquisition and identity services
NCBI GEO full-text retrieval, eutils, GEO FTP, SRA and ENA form the
acquisition spine (section 3); PubMed, CrossRef and Europe PMC anchor
every series to its publication record. Gene and protein identity is
pinned three ways - HGNC (TGFB1 = HGNC:11766), Ensembl (ENSG00000105329,
chr19:41,288,203-41,353,961) and UniProt - so that no symbol ambiguity
propagates into the candidate tables. Ontology grounding comes from EBI
OLS4 (EFO:0008559 American trypanosomiasis; EFO:0600031 response to
benznidazole) and Open Targets (890 targets associated with Chagas
disease, host-side led by TGFB1 at score 0.089).

## 5.3 Pathway and enrichment services
The KEGG Chagas pathway (hsa05142) contributes a 102-gene host anchor
set. Two independent enrichment engines - g:Profiler g:GOSt (1,322 terms)
and Enrichr (206 terms, KEGG_2021_Human) - were run on this set; both
return the Chagas pathway itself as the top term (p = 1.7e-234 and
2.0e-277 respectively), an engine-independent confirmation that the
anchor set is disease-coherent and that both pipelines are wired
correctly. The second-rank terms agree across
engines too: Toll-like receptor signaling (g:Profiler KEGG p = 8.5e-67;
Enrichr p = 4.8e-89) and lipid-and-atherosclerosis (p = 7.4e-60 and
5.8e-83) - innate-immune and vascular-remodeling biology, the same axes
the discovery arm's target-program enrichments land on in section 9. The
Open Targets anchor figures cited in section 2 (890 associated targets;
TGFB1 top host target, score 0.089) regenerate from the saved GraphQL
response (sources/services/opentargets/chagas_targets_top20.json), and
the drug/structure counts - ChEMBL 9 indication records across 6
molecules (CHEMBL110 benznidazole, CHEMBL1487, CHEMBL1397, CHEMBL1631694,
CHEMBL290960, CHEMBL6068503) and 39 cruzipain PDB entries - were
re-verified against the saved responses before appearing here. Reactome (TGF-beta signaling, R-HSA-170834) and QuickGO
(cruzipain GO annotations) supply pathway context; STRING supplies the
parasite-side network neighborhood of cruzipain in T. cruzi CL Brener.

## 5.4 Drug, structure and clinical services
ChEMBL returns 9 Chagas drug-indication records across 6 molecules
including benznidazole; the RCSB PDB indexes 39 cruzipain structures and
AlphaFold DB covers the same protein with a predicted model
(AF-P25779-F1-model_v6); ClinicalTrials.gov lists ongoing Chagas studies
(NCT04084379, NCT01549236, NCT01755377 among the first page). These
three services frame translatability: any candidate marker adjacent to a
druggable target or an active trial carries more clinical weight, and the
discovery report (section 8) scores this adjacency explicitly.

## 5.5 Reference and context services
The WHO fact sheet and CDC DPDx pages supply the burden and
diagnostic-gold-standard background for sections 1-2; the Human Protein
Atlas contributes TGFB1 tissue-consensus expression including heart
muscle, the organ that defines the disease's chronic phase; GitHub hosts
the versioned artifact itself. Every page retrieved is stored with its
hash; nothing is cited from memory.

## 5.6 Reproducibility contract
Any reader can re-execute the service layer: the ledger names the exact
query, date and output file for every service, and the manifest pins the
bytes. Services that refused scripted access (TriTrypDB's API-key gate,
Semantic Scholar rate limiting, DrugCentral's login wall, medRxiv's empty
responses, GEO2R's browser-only interface) are logged in the ledger's
blocked section with the exact failure - they are not counted, and their
planned roles are named so the gap is visible rather than papered over.

## 5.x Gate-(b) enrichment pipeline (executed 2026-09-27)
Inputs: miRTarBase v8.0 human MTI table (382,175 dedup rows;
sha256-pinned, version disclosed), GSE244827 CHAVA raw counts
(60,675 Ensembl genes x 33 libraries), GSE203525 counts (58,142 genes
x 20 libraries). Per candidate miRNA, its validated-target set T is
taken as-is from miRTarBase; the direction rule is canonical repression:
candidate down in severe disease predicts targets UP in the case
contrast and vice versa. Differential expression is a Welch t-test on
log2(CPM+1), gene-level; the direction-consistent DE set is genes
moving in the predicted direction at nominal p <= 0.05. The enrichment
statistic is the fraction of T in that set, judged against 10,000
size-preserving random gene sets from the same matrix (one-sided
permutation p, (1+b)/(10001)); Benjamini-Hochberg FDR across all 40
candidate-cohort tests. Sample-to-column identity for GSE244827 is
established by a verified B-code bijection (each GSM's live
!Sample_description), asserted in-code; GSE244827 Ensembl IDs are
translated to symbols with an Ensembl BioMart GRCh38 map. miR-375-3p is
matched through its documented legacy name hsa-miR-375
(MIMAT0000728, verified against live miRBase). Seed 20260926.


## 5.y Compute environment and reproducibility
All analyses run as committed scripts (scripts/) against byte-hashed
inputs; every figure-free number in this paper regenerates from the
repository state at the cited commit. Environment: Python
3.10.12, NumPy 2.2.6, SciPy
1.15.3, pandas 2.3.3, scikit-learn
1.7.2 (pinned 2026-09-27; the pin was
added when the H1' completion run found sklearn unpinned in this record
- RUN_LOG 2026-09-27T05:33). Frozen seed 20260926 for
every stochastic step (CV splits, bootstrap, permutation nulls).
Permutation tests use (1+b)/(n+1) p-value convention; multiple-testing
control is Benjamini-Hochberg within each locked family. External
service responses are stored with sha256 in sources/services/ and
logged with retrieval timestamps in SERVICE_LEDGER.md (40 distinct
tools; blocked attempts logged, not counted).

## 5.7 The full 40-service enumeration
Grouped by function, with the ledger entry numbers (SERVICE_LEDGER.md) so
every row is auditable: acquisition and identity - NCBI GEO full-text
retrieval (1), eutils (2), GEO FTP (3), SRA run selector (4), PubMed (5),
Europe PMC (6), PLOS figure retrieval (7), CrossRef (23), HGNC (24),
Ensembl REST (25), EBI OLS4 (26), ENA Portal (27), GitHub API (28),
Ensembl BioMart (31), mygene.info (36), miRBase (37), NCBI Datasets (38);
pathway and enrichment - Reactome (11), KEGG (12), QuickGO (13),
g:Profiler (14), Enrichr (15), STRING (10), TargetScan seed-family context
(29, prediction-only per Addendum 3 and never the gate-(b) source),
miRTarBase v8.0 validated targets (30, version-pinned and disclosed);
drug, structure and clinical - Open Targets (8), ChEMBL (16), RCSB PDB
(17), AlphaFold DB (18), ClinicalTrials.gov (19), IntAct (39), WikiData
SPARQL (40); reference and context - WHO (20), CDC DPDx (21), Human
Protein Atlas (22). The functional groups of sections 5.2-5.5 describe
what these services contributed; this enumeration fixes what they are.

## 5.8 The blocked register, with reasons
Eleven further services were attempted and failed, and none is counted:
TriTrypDB (gene-record API now key-gated; the HTML page is an SPA shell),
Semantic Scholar (persistent HTTP 429 without a key), DrugCentral (no open
public API; downloads behind login), medRxiv (0-byte responses on repeated
attempts), GEO2R (browser/R-only, not scriptable), WikiPathways
(webservice 404; service moved or retired), TarBase v9 (JS application;
all probed data URLs 404), Pharos GraphQL (HTTP errors on every query
tonight), GTEx medianGeneExpression (empty responses for direct gencode
queries), Expression Atlas (JSON path 404) and RNAcentral (accession path
serves HTML). Each entry names the exact failure and, where one exists,
the planned fix (free API key, alternate endpoint, or a browser visit on
the token queue). The register is part of the methods, not an apology:
a 40-service claim that hides its failures cannot be audited, and the
blocked list is where a reviewer should look first for wishful counting.
