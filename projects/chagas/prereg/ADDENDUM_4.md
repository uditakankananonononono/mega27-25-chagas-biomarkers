# PREREGISTRATION ADDENDUM 4 - locked 2026-09-27T08:58 IST, BEFORE any
D1/D2/D3 run. Origin: judge round 02 (SUPPLEMENTARY Gemini consult,
judge_rounds/02_*, thread app/2613ca9536a087c2 - not counted toward the
ChatGPT-only 10-round gate; X-JUDGE-ROUNDS.md). Adopted by us via this
lock; the judge's text itself is advisory.

## D1. Compartment-specificity decisive test (frozen)
Judge-02 Q1 ruled the cardiac-18/20 vs blood-6/20 "compartment
specificity" reading a post-hoc patch until tested. Frozen test:
- Target sets: the same per-candidate miRTarBase validated-target sets
  used by the committed gate-(b) run (no re-curation).
- Response (gene level): direction-consistent DE in the GSE244827
  early-CCC-vs-seronegative contrast, binary per gene, same DE
  machinery as scripts/h2_gate_b_enrichment.py (Welch t on log2(CPM+1),
  predicted direction per the frozen direction rule).
- Predictors from HPA (ledgered service; new file retrievals logged as
  re-uses, sha256-hashed): blood_abundance = log1p(max NX across blood
  cell types, HPA blood-cell dataset); cardiac_spec =
  log2((heart_muscle_NX+1)/(max_blood_cell_NX+1)). If the HPA blood-cell
  file is not scriptably retrievable, the fallback is HPA consensus
  tissue NX with blood-related tissues where present; the fallback is
  disclosed in the run record.
- Model: logistic regression, gene level over the union of the 20
  candidates' mapped targets: replication ~ blood_abundance +
  cardiac_spec. Report coefficients, Wald p, descriptive AUC.
- READING, locked both ways BEFORE the run: blood_abundance coefficient
  positive with p<0.05 supports compartment specificity; a null or
  negative result DOWNGRADES the paper's 9.3b tissue reading (edit
  recorded as a correction, not silently removed). Candidate-level
  (6-pass vs 14-fail vs median blood abundance) is descriptive only
  (n=20), never inferential.

## D2. OOF module reconstruction + permutation null (frozen)
Judge-02 Q2: the module's in-sample circularity stands until the score
is reconstructed out-of-fold.
- Folds: the same outer 5-fold stratified splits of the 146 graded
  samples, seed 20260926, as the locked H1' pipeline.
- Per outer fold, on the 4/5 training split ONLY: KW per miRNA across
  the four grades -> BH FDR<=0.05 AND |d(severe-mild)|>=0.5 -> the
  frozen 11-marker literature exclusion -> intersect survivors with the
  COMMITTED gate-(b) both-tissue filter (external cohorts; no serum-
  fold leakage) -> sign weights from training medians only -> z-score
  parameters (mean, sd per member) estimated on training samples only.
- Score the held-out fold; pool all out-of-fold scores; report pooled
  OOF ordinal c-index and KW p.
- Null: 1,000 permutations of the grade labels (permuted before the
  whole per-fold pipeline; expression ranks are label-invariant, so the
  implementation may use vectorized rank sums - implementation detail,
  not design). Empirical p = (1 + #{null >= observed}) / 1001.
- READING, locked both ways BEFORE the run: the claim "the module
  tracks severity beyond its in-sample construction" stands ONLY IF
  permutation p < 0.05. The OOF c-index magnitude is reported as-is,
  with the gap to the in-sample 0.756 stated. A failure is documented
  in the paper (sections 9.6/9b.1), not hidden.

## D3. TargetScan sensitivity (frozen, NON-GATE)
Judge-02 Q3: prove the gate-(b) enrichment is not a miRTarBase-v8.0
annotation artifact.
- Retrieval: TargetScan 8.0 human predicted-target table (version
  pinned, sha256-hashed, ledgered as a re-use of service 29).
- Procedure: the committed gate-(b) pipeline, unchanged, with TargetScan
  predicted targets (family mapping via miR_Family_Info) substituted for
  miRTarBase validated sets, in both cohorts, same permutation null and
  BH family.
- READING: this is a SENSITIVITY, never a gate: gate-(b) verdicts stand
  on miRTarBase validated targets per ADDENDUM_3 C2 and judge-01's
  demotion of bare predicted-target databases. Overlap, divergences and
  CORE6 membership under TargetScan enrichment are reported as-is.

## D4. Presentation adoptions (no gate, no run)
- D4a: present the CORE6 as the output of the two-stage sequential
  filtration (2,632 -> 20 serum candidates -> 6 both-tissue core)
  wherever membership is described (9.1/9.3/9.7).
- D4b: publish the discordance audit for the 12 cardiac-passing
  candidates excluded from the core (gate-(a) directions x gate-(b)
  blood results; descriptive).
- D4c: abstract/intro framing audit toward judge-02's ranking:
  framing (c) primary, (b) hedged until D1+D2 complete, (a) never.

## Execution order (frozen): D2 (pure compute) -> D1 (HPA retrieval) ->
D3 (TargetScan retrieval) -> D4 edits. Every run logged in RUN_LOG with
its script; results files committed; honest negatives reported.
