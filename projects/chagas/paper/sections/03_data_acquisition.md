# 3. Data acquisition and the provenance model

## 3.1 Design principle
Every record in this compendium is individually retrievable, byte-verified,
and re-checkable by an independent reader with nothing more than the public
URL recorded against it. We rejected the two common shortcuts in
secondary-analysis papers - citing series-level accessions without
sample-level evidence, and trusting prior download caches - because both
fail silently when a source record changes or was never what the analyst
assumed. The unit of provenance here is the individual GEO sample (GSM),
not the series (GSE).

## 3.2 Canonical retrieval
For each GSM we fetched the NCBI GEO full-text record from its canonical
URL (acc.cgi?acc=<GSM>&targ=self&form=text&view=full) and stored the exact
response bytes under sources/soft/. The sha256 of those bytes is recorded
in the per-series crosswalk (sources/<GSE>_sample_crosswalk.csv), one row
per sample with columns gsm, source_url, sha256, title, status, label,
platform, organism, and SRA relation where present. Series-level records
were fetched and hashed the same way. Nothing enters the compendium from
memory, screenshots, or third-party mirrors.

## 3.3 Hermetic verification
A committed verifier (scripts/verify_crosswalks.py) re-reads the committed
evidence and checks, without any network access, that: (i) every crosswalk
sha256 matches the bytes of its stored SOFT file; (ii) every GSM identifier
is unique across the compendium - the check that catches super-series
double-counting (a lesson applied from the leishmaniasis lane's GSE191083
duplicate and, this lane, the GSE244827 prior-series overlap in section
3.6); (iii) no acquired series collides with the program's prior-tagged
exclusion set; and (iv) every sample carries a non-empty label or an
explicitly documented honest-label limitation (section 3.5). Live
re-verification additionally re-fetches a random sample of records
(5-10%, new seed per run) and byte-compares against the recorded sha256.

## 3.4 Discovery sweep and relevance screening
Candidate series were found with live NCBI eutils queries (db=gds, Chagas
keyword, Homo sapiens, entry type GSE; 25 hits on 2026-09-26). Keyword
matches are not evidence of relevance: every candidate was screened against
its series title and design text, and non-Chagas matches were rejected with
the reason logged in ACQUISITION_LOG.md (e.g. GSE78975 anxiety methylome;
GSE27353/GSE27054 thymocyte hormone studies). Relevance screening, not
query recall, is what makes the compendium trustworthy.

## 3.5 Honest labels
Where GEO encodes per-sample case/control or severity attributes, labels
were taken verbatim. Where it does not - for example pooled-donor designs
or stimulation timecourses - we labelled the actual treatment contrast and
recorded the limitation rather than forcing clinical labels the source
does not support. Forced labels are the quiet failure mode of
secondary-analysis cohorts; we prefer an honest "treated vs control" to a
borrowed "case vs control".

## 3.6 De-duplication and the count-correction record
Program record counts are computed against the prior-tagged manifest on
the main branch (results/dataset_manifest.csv plus results/series/), never
from this lane's status text. Two double-counts were caught and corrected
by live re-verification on 2026-09-26: GSE244827 (33 GSM + 1 GSE) was
already prior-tagged under the P31 work and had been re-acquired in the
expansion; the lane total was corrected 750 -> 716 with exact-match proof
(commit ef7b29a). The standing rule that produced the catch - no count is
reported from a status file without a live recomputation - is now applied
to every claim in this paper. The corrected, live-verified total is
716 unique records: 15 GSE series, 700 GSM samples, and 1 other record
(Open Targets disease entry EFO_0008559), against a program floor of 120.

## 3.7 Record inventory
The 700 crosswalk GSM rows span 16 crosswalk files covering 15 distinct series
(the sixteenth file is the prior provenance file for GSE84796,
GSE84796_used_sample_crosswalk.csv, not a separate series); 650 rows are new acquisitions and 50 are
prior-tagged rows retained as provenance. Modalities: bulk RNA-seq (blood,
hiPSC-cardiomyocyte, dendritic cell, placenta), small/miRNA-seq (serum,
placenta, macrophage), single-cell RNA-seq (PBMC), spatial transcriptomics
(Visium, placenta), expression microarray (placenta, cardiac), DNA
methylation (blood), and SNP pharmacogenomics. Organism: Homo sapiens
throughout; parasite-side context is supplied by service retrievals
(section 5), not by mixing organisms into the human cohort atlas.

## 3.8 Verification statistics and the spot-check record
The first acquisition pass preserved 411 per-sample SOFT files plus 11
series-level files, every one sha256-hashed at fetch time. The hermetic
verifier re-hashes all of them on every run and checks SOFT
well-formedness, GSM uniqueness across series, disjointness from the six
previously tagged series, non-empty labels, and record accounting; all
checks pass at the current tip. Hash-matching alone only proves we kept
what we downloaded, not that the download was ever faithful, so the
protocol adds live spot re-verification: after the first pass, 21 GSMs
drawn at random (seed 25, about 5% of 411) were re-fetched independently
and all 21 refetched byte strings reproduced their recorded sha256 with
zero mismatches; a separate single-record reproduction at acquisition time
(GSM9683415) also matched its recorded sha256 exactly. One structural
caveat is logged rather than smoothed over: series-matrix files exist on
the GEO FTP mirror for 10 of 11 first-pass series, while GSE158986
(dual-organism RNA-seq) ships supplementary count files only, so its
downstream processing path differs from the series-matrix path and is
documented as such.

## 3.9 Label-assignment machinery
Labels are not read off accession titles. scripts/label_crosswalks.py
applies explicit per-series rules over the frozen characteristics and
titles, and the verifier fails if any row is left unmapped; the current
corpus has zero unmapped rows. The vocabulary is deliberately split
between clinical cohorts (case/control, with severity and disease form
retained verbatim in characteristics) and mechanism contrasts
(infected/treated/variant/reference - section 4.6b). A worked example of
why this machinery exists: the GSE311812 series design text states five
transmitter blood samples, but the deposited sample titles record six; we
treated the titles as the record of truth (25 case / 21 control) and
logged the discrepancy in ACQUISITION_LOG.md rather than silently
choosing either number.

## 3.10 The GSE244827 matrix-column bijection
The GSE244827 expression matrix column headers (B052..H754, 33 libraries)
are laboratory codes that do not describe themselves, and the acquisition
crosswalk rows carry no B-code field, so column-to-label assignment could
not be taken from either artifact alone. We resolved it against live GEO:
each GSM's SOFT record carries its CHAVA B-code in !Sample_description
(e.g. GSM7830424 -> B052). All 33 records were re-fetched on 2026-09-27
00:35 IST; the resulting map (sources/GSE244827_column_label_map.csv) is
bijective with the 33 matrix columns and confirms that acquisition-time
GSM order matches matrix column order exactly. Labels then come from the
frozen crosswalk (10 case / 23 control), not from the live page. One
honest caveat is recorded in prereg/RUN_LOG.md: the re-fetched SOFT page
bytes do not hash-match the acquisition-time crosswalk hashes because the
page carries dynamic content, but the B-code and label metadata lines are
stable and mutually consistent with the matrix header. The analysis
script (scripts/h2_gate_b_enrichment.py) asserts the bijection and the
per-group sizes before computing anything and fails loudly on mismatch,
so a future matrix-or-map drift breaks the run instead of silently
mislabelling samples. The same pattern - assert group structure on the
live header rather than assume it - is applied to GSE203525 (6 CCC + 6
indeterminate at 0hpi verified against the header).

## 3.11 Rejected-series register
Relevance screening rejected more candidates than it accepted, and every
rejection is logged with its reason in ACQUISITION_LOG.md so the corpus
boundary is auditable. Concrete rejections from the 2026-09-26 sweep:
GSE78975 (anxiety-disorder methylome; keyword match only, no Chagas
content), GSE27353 and GSE27054 (thymocyte hormone studies; no Chagas
content), and GSE7047 (a 2007 infected-cell-line array whose platform
GPL1053 is already implicated in the program's skip_labels exclusions).
GSE191083 was acquired and then removed by the hermetic verifier itself:
its uniqueness check showed the series is the super-series container of
GSE191081 and GSE191082, so all 180 of its GSMs duplicate records already
held through the component series. That removal is exactly the failure
the uniqueness gate exists to catch, and it is why compendium counts are
always recomputed from the frozen crosswalks rather than quoted from
status text (section 3.6).

## 3.12 The acquisition timeline, in passes
The compendium was built in two passes plus corrections, each logged
with its own verification record. Pass one (2026-09-26) closed the
67-record shortfall to the 120-record floor: 11 new series, 411 new GSM
records, chosen by the live eutils sweep (25 hits, relevance-screened)
with the six previously tagged series excluded; the hermetic verifier
passed and 21/21 spot re-fetches matched their recorded hashes. Pass two
(same day, on direction to acquire the reserved leads) added 3 series
and 271 GSMs - the methylation pair GSE191081/GSE191082 and the
benznidazole pharmacogenomics series GSE154421 - and immediately
produced the lane's two teaching errors: GSE191083 was acquired and
removed within the pass when the uniqueness check exposed it as the
super-series container of the methylation pair, and the expansion total
was corrected 750 -> 716 when live recomputation showed GSE244827 had
been re-acquired despite already sitting in the prior manifest (the
exclusion list had not been updated before the pass). Both corrections
are recorded append-only (ACQUISITION_LOG.md, commit ef7b29a) and both
are why every count in this paper is recomputed from the frozen
crosswalks at claim time: the errors were caught precisely because the
verifier, not the narrative, keeps the books.
