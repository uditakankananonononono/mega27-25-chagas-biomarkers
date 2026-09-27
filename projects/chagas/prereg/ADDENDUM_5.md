# ADDENDUM_5 - locked 2026-09-27T09:23 IST BEFORE any E-run.

Origin: counted judge round 02 (ChatGPT, user-couriered; judge_rounds/02_response_chatgpt.txt,
02_assessment_chatgpt.txt) + her directive "ADD ALL THIS TO IT AND WHEN COMPLETE I WILL
CHECK AGAIN" (WhatsApp 09:17:56, wamid...M0VCMEMwMzdDRUQwQjE1RkRBOEJFNwA=). The judge text
is advisory; adoptions below are OUR locks, readings fixed both ways before running.
Overlap map: D1 (ADDENDUM_4) still runs as locked (HPA predictors of blood replication);
D2 done; D3 still runs (TargetScan sensitivity); D4 edits absorb the framing adoptions.

## E1. Nested bootstrap discovery (CENTERPIECE; extends D2 to full pipeline)
B=1000 bootstrap resamples of the 146 graded samples (rng seed 20260927). Inside EACH
bootstrap, on the bootstrap sample only: (a) KW screen, BH FDR<=0.05 AND |d|>=0.5;
(b) frozen 11-marker literature exclusion (verbatim list, scripts/h2_severity_association.py);
(c) gate-(b) direction-consistent enrichment per surviving candidate on BOTH external
cohorts (GSE203525, GSE244827) using the committed miRTarBase v8.0 target sets and the
frozen direction rule - COMPUTE-FAIR ARM LOCKED: 1,000 permutations per test inside
bootstraps (not 10k; disclosed reduction, same size-preserving scheme); (d) AND
intersection across cohorts; (e) signed module from bootstrap survivors (train=bootstrap
members, sign weights + z-params from bootstrap sample); score OUT-OF-BOOTSTRAP samples.
Record: per-miRNA full-pipeline discovery frequency; OOB c-index per bootstrap; member
counts. Empty-module bootstraps: score NaN, counted and reported.
READING (locked both ways): discovery claim STRENGTHENS if all 6 CORE6 members recur at
>=50% AND the OOB c-index 95% percentile interval excludes 0.5; any CORE6 member <50%
recurrence is reported as instability against that member; median OOB below D2's 0.704
is reported as the honest generalization figure. No threshold is edited after the run.

## E2. Leave-one-cohort-out biological validation + cross-compartment concordance
Per the 20 candidates: DirectionScore = (# independent tissues agreeing with predicted
direction)/(# tissues tested) over serum (GSE299582), hiPSC-CM (GSE203525), whole blood
(GSE244827), heart-tissue series among the 15 (as available, disclosed per series).
Null: 10,000 random miRNA sets matched on (i) miRTarBase validated-target count,
(ii) mean log2CPM in GSE299582, (iii) publication count from the committed novelty
corpus metadata (121-abstract frozen corpus; if per-miRNA counts are not recoverable
from committed artifacts, PubMed esearch counts via the ledgered NCBI service, re-use
logged, sha256 where files). READING: concordance "unusual" only if observed
DirectionScore >= 95th percentile of the matched null; otherwise compartment language
stays hypothesis-level and 9.3b wording is corrected (not silently).

## E3. miRTarBase-bias negative controls
Control A: random miRNAs matched on validated-target count, expression, publication
count (same matching as E2). Control B: random validated targets from the same miRNA
families. Control C: TargetScan/miRDB predicted-target arm (= ADDENDUM_4 D3 machinery,
non-gate). READING: the cardiac 18/20 enrichment claim holds only if it exceeds all
three controls at the locked FDR; failure -> gate-(b) claim downgraded to
"database-bias-consistent", disclosed in paper 9.x.

## E4. Blood-signature pathway correlation (non-gate, descriptive)
Disease-pathway enrichment (Chagas/heart-failure sets, Enrichr re-use ledgered) of
blood DE genes among candidate targets; FDR<=0.05 descriptive only.

## E5. Influence diagnostics (severe n=30)
Cook's distance per sample on the locked H1' OOF logistic fits; leave-5%-out (7-sample)
stability; leave-severe-out sensitivity (refit on 116, report control-vs-mild/moderate
c-index + KW without severe). READING: signal attributed to a true gradient only if
leave-severe-out KW stays <0.05 AND no single sample's Cook's D exceeds 4/n rule
without disclosure; outlier-driven findings downgrade the gradient claim.

## E6. Post-hoc recalibration (disclosed post-hoc)
Isotonic and Platt scaling on OOF predictions of the locked H1' model; report
calibration slopes post-recalibration; paper keeps "NOT decision-grade" language.

## E7. Decision-curve analysis
Net benefit of the OOF ordinal score across threshold range vs treat-all/none;
descriptive, no clinical-use claim.

## E8. Endpoint hierarchy (prereg note + paper edit)
PRIMARY claim (one): an ordinal serum-miRNA severity model beats internal comparators
(age/sex, best single miRNA) on locked OOF folds. SECONDARY: 6-miRNA both-tissue core
(module, with D2/E1 stability caveats); TERTIARY/exploratory: everything else
(enrichment, druggability, concordance). Paper restructured to match.

## E9. Controlled ML baselines (same folds, seed 20260926)
Elastic net, random forest, gradient boosting (XGBoost if importable in this sandbox,
else sklearn HistGradientBoosting - substitution disclosed), RBF-SVM, ordinal logistic
- identical outer 5 folds and inner 3 folds, OOF ordinal c-index. READING: H1'
advantage stated only over models it actually beats on OOF; any baseline >= H1' is
reported as parity/loss, not hidden.

## E10. Batch robustness
If a batch/plate/lane covariate exists in GSE299582 metadata: batch-stratified KW +
ComBat sensitivity. If absent: documented infeasible (no fake batch arm).

## E11. Confounder adjustment
Age/sex-augmented ordinal logistic on the locked folds (BMI/comorbidity NOT in the
series metadata -> documented); compare OOF c-index vs unadjusted.

## E12. qPCR-style reduced panel simulation
Top 3-6 markers by E1 discovery frequency; OOF eval on the locked folds; report
c-index vs full module (assay-transfer evidence, simulation-level only).

## E13. Disease-specific enrichment (non-gate)
Enrichment against Chagas/heart-failure gene sets (Enrichr re-use, ledgered) for
module targets; replaces generic-only pathway language where supported.

## Deferred-infeasible (documented, not hidden)
- Longitudinal trajectory: all 15 series cross-sectional; no public follow-up cohort.
- Mendelian randomization: no miRNA GWAS instruments available; causal language banned.

## Execution order (frozen): E1 -> E2 -> E3 (with D3 arm) -> E5 -> E9 -> E11 -> E6 ->
E7 -> E12 -> E4 -> E13 -> E10 (metadata check) -> E8/D4 paper edits. D1 (ADDENDUM_4)
runs interleaved as compute allows; HPA retrieval first. Every run logged in RUN_LOG
with its script; results committed; honest negatives reported; no gate claims until met.
