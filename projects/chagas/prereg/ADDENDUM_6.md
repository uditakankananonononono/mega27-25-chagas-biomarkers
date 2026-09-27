# ADDENDUM_6 - locked 2026-09-27T10:02 IST BEFORE any F-run.
Origin: counted judge round 03 (ChatGPT, user-couriered manuscript review,
judge_rounds/03_response_chatgpt.txt + 03_assessment_chatgpt.txt). Judge text
advisory; these are OUR locks. Round-02 overlap items (E1-E13) are NOT
duplicated here; this addendum locks only the NEW round-03 analyses.

## F1. Provenance-Aware Omics Discovery Framework benchmark (new method claim)
Formalize the acquisition/provenance pipeline as a measurable method.
Benchmark table: traditional GEO workflow vs this workflow on the 15 series:
metadata errors detected, sample mismatches found, missing annotations
recovered, reproducibility failures prevented - each row tied to a committed
artifact (RUN_LOG corrections, crosswalk fixes, hash mismatches). No invented
counts: every cell cites a committed artifact or is marked not-tracked.
READING: if fewer than 3 measurable categories have committed evidence, the
framework claim is scoped down to "documented practice", not a benchmark.

## F2. Extended null framework (extends D2)
On the H1' pipeline, same folds/seed: (a) label shuffle - DONE in D2 style
for the module; add for the H1' model itself; (b) shuffled miRNA identities
(feature-name permutation); (c) random 100-feature sets, 1,000 draws each.
Report observed OOF c-index vs each null distribution with empirical p.
READING both ways: observed must exceed the 95th percentile of each null to
keep the "not random" sentence in the paper; otherwise that sentence is cut.

## F3. Permutation importance on the locked H1' model
Permutation importance of the top-100 features within each outer fold
(1,000 feature permutations, OOF c-index drop); reported alongside E1
selection frequencies. Descriptive; no claim change either way.

## F4. Cross-database target confirmation (extends D3)
miRDB and miRWalk arms for the 20 candidates, same gate-(b) machinery,
non-gate. Retrieval pinned + sha256; if a database is not scriptably
retrievable on free tiers, that arm is documented infeasible, not faked.
Validated (miRTarBase) evidence always reported separately from predicted.

## F5. Standard-pipeline comparison
Same severity screen re-done as: (a) plain differential expression
(no FDR+effect double gate), (b) PCA clustering, (c) RF feature selection,
(d) limma-style workflow. Compare candidate sets against the frozen 21:
overlap counts only. READING: framework superiority claimed only where the
frozen pipeline's candidates are not recovered by the plain alternatives;
overlap is reported verbatim either way.

## F6. Expanded novelty screen
Europe PMC full-text where open, citation-network expansion from the 121-
abstract seed corpus, automated title/abstract similarity against the 20
candidates. Boundaries: same exclusion rule as the frozen screen; any new
hit downgrades that candidate from "novel" to "reported", disclosed.
No candidate is added by this screen alone.

## F7. Causal-prioritization ranking (no causal claims)
Rank the 20 candidates by: direction consistency across cohorts, network
centrality (committed STRING graph), tissue agreement (E2), disease-pathway
proximity (E13 sets). Output is a ranked table with the four sub-scores.
Causal wording banned; this is a prioritization heuristic, labeled as such.

## F8. Blind machine reproduction
Clean-clone rerun of the pipeline from the repo alone (fresh clone, fresh
venv, pinned deps): compare outputs, figures and sha256 hashes against
committed artifacts. Human-independent arm infeasible solo - documented.
READING: reproducibility claim graded by fraction of artifacts that
reproduce byte-identical; mismatches listed, not fixed silently.

## F9. Central-claim restructure (paper edit, D4 batch)
Central claim becomes: "a provenance-first computational framework improves
reliability and reproducibility of public omics biomarker discovery";
compendium, severity model, module and druggability become its validation.
Title/abstract/intro/discussion restructured to match; no new science claim.

Execution order (frozen): F2 -> F3 -> F5 -> F7 (compute-only) -> F1 (ledger
mining) -> F4 (retrieval-dependent) -> F6 (retrieval-dependent) -> F8 ->
F9 (paper restructure, last). Every run in RUN_LOG; honest negatives
reported; no gate claims until verifiably met.
