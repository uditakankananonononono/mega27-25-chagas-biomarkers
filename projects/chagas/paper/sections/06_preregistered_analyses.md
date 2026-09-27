# 6. Preregistered analyses and the first gated run

## 6.1 What was locked, and when
Two hypotheses were frozen in prereg/PREREGISTRATION.md (commit 6e562c7)
before any outcome value was computed: H1, a benchmark-beat claim on the
192-sample serum miRNA severity cohort GSE299582, and H2, a cross-modal
convergence discovery gate requiring an effect, an orthogonal replication,
and absence from a frozen literature screen. The literature screen itself
- 121 PubMed abstracts, the full relevance-sorted set, with a curated
named-marker exclusion list - was committed next (commit 70ab16d). When
the frozen comparator rule resolved to an empty marker subset (the
published panel is protein-based; the cohort is miRNA), the resolution
was locked as Addendum 1 (commit 44a4e5e), again before outcomes: the
primary contrast became symptomatic CCC versus asymptomatic seropositive
- the same contrast the published prognosis arm reports - with the
published best Cox & Snell R2 of 0.688 as the benchmark value.

## 6.1b The frozen hypotheses, verbatim in structure
H1's lock named its endpoint (CCC severity ordinal grade,
control<mild<moderate<severe), its primary metric (macro-averaged ROC-AUC,
one-vs-rest, with quadratic weighted kappa secondary), its comparator
anchor (PMID 34479416's panel as the cross-sectional benchmark; the
longitudinal parasite-DNA + 47-protein model, PMID 38203212, frozen as
secondary context, explicitly not the cross-sectional benchmark), its
pipeline (70/30 stratified split, seed 20260926, train-only univariate
FDR<=0.05 selection, elastic-net ordinal logistic with 5-fold inner CV,
single frozen-test fit), and its beat criterion (primary metric exceeding
the comparator by >=0.05 absolute with the bootstrap CI of the difference
excluding zero; ties are negatives). H2's lock named its three gates: (a)
differential severity association in GSE299582 (FDR<=0.05 and |effect|>=0.5
log2 fold-equivalent), (b) direction-consistent replication in the
orthogonal blood or cardiomyocyte cohorts, (c) absence from the frozen
literature screen - plus the survival clause that a feature becomes a
"new discovery" only when all three gates hold AND it survives the
judge rounds' novelty critiques; anything less is logged as replication
or negative. The frozen dataset identifiers named GSE299582 as primary
(the only acquired cohort with graded severity at n>100), GSE244827 as
orthogonal secondary, the remaining series as context, and the
prior-tagged exclusion set that the count-correction later enforced. The
integrity rules completed the lock: no outcome value read into the repo
before the preregistration and screen were committed; splits, seeds and
code committed before first run; every run logged with commit hashes.

## 6.2 Run 1: an honest negative on the primary arm
The frozen pipeline (70/30 stratified split, seed 20260926, train-only
univariate feature selection, elastic-net logistic regression) was
executed exactly as locked. On the primary binary arm (n=150: 104 graded
CCC, 46 indeterminate), the frozen model scored R2 = -55.7 on the frozen
test set against the 0.688 benchmark - worse than the null model, a clear
NEGATIVE under the locked tiers (commit 49de4ae, results JSON and run log
in prereg/RUN_LOG.md). The frozen pipeline as specified does not beat the
published protein panel on this contrast.

## 6.3 The signal that did appear
The secondary ordinal arm told a different story: four-class severity
(control, mild, moderate, severe) separated with macro-AUC 0.779 under
the same frozen split discipline. Severity-graded information is present
in the serum miRNA signal; the frozen binary pipeline's failure is a
pipeline-and-contrast failure, not evidence of an empty cohort. This
asymmetry - negative on the locked primary, positive descriptive signal
on the secondary - is exactly what preregistration is for: it keeps the
negative honest and prevents the positive from being silently promoted
into a claim it was never registered as.

## 6.4 What happens next (and what does not)
Per the standing pivot rule, the negative triggers a documented
redirection: a ChatGPT consult on pipeline redesign, then Amendment 2 -
locked before any rerun - with candidate fixes named in the run log
(log1p transform, class weighting, revised feature-count rule, and
re-examination of whether R2 is a fair cross-study metric at n=150
against a published n=60). Run 1 is not edited, re-run silently, or
removed; it stands as the first entry of the run log. Whatever Amendment
2 produces will be measured against the same locked benchmark value and
reported under the same tier rules.

## 6.4b Addendum 2: the source-paper blind spot, locked
Two minutes after Run 1 was logged, a second amendment closed a hole the
original screen could not see. The frozen 121-abstract screen never
covered the GSE299582 SOURCE paper itself (Roma et al. 2026, J Infect Dis
234(2):277-287, PMID 41574750): its abstract names nine miRNAs by
contrast direction (miR-143-3p, miR-223-3p up in ChD; miR-486-5p,
miR-3960 down in ChD; miR-6734-5p, miR-1285-5p, miR-10527-5p, miR-1228-5p
up severe vs mild; miR-30c-3p down severe vs mild) and reports 40
differentially expressed miRNAs in total, but the full list is
paywalled. Addendum 2 (B1) added all abstract-named miRNAs to the
replication class and froze the residual-risk rule: candidates matching
the unseen 40-DEM list cannot be excluded, so every H2 candidate is
additionally checked against the abstract contrast directions, and the
contamination risk stands as a logged limitation rather than being waved
away. The same addendum (B2) operationalized gate (a) exactly as later
executed - per-miRNA Kruskal-Wallis across the four ordinal groups
(n=146), BH FDR<=0.05, and |median log2(CPM+1) severe - mild| >= 0.5 -
while explicitly leaving gate-(b) routing PENDING the redirection
consult, with the anticipation that target mapping, if chosen, would use
TargetScan or miRDB with logged service evidence. Addendum 3's C2 later
overrode that anticipation with the stricter miRTarBase-only rule; the
sequence - anticipate, then tighten by amendment - is preserved because
each step is its own timestamped commit.

## 6.5 Addendum 3: the judge-01 redesign, locked before any rerun
Judge round 01 (judge_rounds/01_*, advisory text; adoption is ours)
produced Addendum 3, locked 2026-09-26T21:23 IST before any new outcome
run. Three clauses. C1: H1 closed as a documented negative; H1'
registered - ordinal severity prediction on the 146 graded samples with
immediate-threshold logistic regression, nested 5x3 CV, seed 20260926,
top-100 inner-fold feature selection, ordinal c-index primary, and
molecule-matched benchmarks (clinical age/sex model; best single miRNA
chosen on inner folds). The beat criterion requires exceeding BOTH with
bootstrap CIs of the difference excluding zero. C2: gate (b) routed to
experimentally validated targets only (miRTarBase; prediction-only
databases explicitly rejected), frozen direction rule, two orthogonal
cohorts, 10,000 size-preserving permutations, BH FDR <= 0.05; same-
molecule placenta replication and monotonic-trend "internal replication"
rejected as replication. C3: the discovery deliverable frozen as a
multi-omic severity module score resting only on gate-(b)-passing
programs plus a druggability overlay, with the claim wording frozen as
"a conserved regulatory module tracks transition toward cardiac disease"
- never "miR-X predicts CCC".

## 6.6 Execution record against the amendments
H1' ran as locked (section 8): CLEAR BEAT on both locked comparators.
Gate (b) ran 2026-09-27 exactly as C2 specifies - miRTarBase validated
targets, frozen direction rule, both named cohorts, the locked
permutation null - with two execution-time disclosures made in the run
log rather than silently absorbed: the retrieval pinned miRTarBase v8.0
(the live site's download paths had restructured; the official file was
recovered byte-hashed from an archive snapshot of the publisher URL),
and the GSE244827 sample-to-column identity, assumed positional at
staging, was re-derived as a verified B-code bijection from live GEO
records and asserted in code. Both disclosures strengthen the run; the
design did not move. The module score and druggability overlay followed
per C3 (section 9): the both-tissue core satisfies every C3 tier; the
in-sample nature of the score is stated wherever the number appears.

## 6.7 Judge round 01: the record behind Addendum 3
Round 01 ran 2026-09-26T21:21 IST via the text-paste route (ChatGPT Free
thread 6ab7e9ef; prompt and full response preserved verbatim in
judge_rounds/01_prompt.txt and 01_response.txt and re-read in full before
the consult ended). The persona then was a simulated competition-style
critique; the framing has since been standardized to an expert research
reviewer (competition framing dropped), and the advisory status is
unchanged: judge content is input, never authority - every adopted change
was locked by us in Addendum 3 before any run used it. The verdict:
prediction weak, biology promising, novelty moderate, with the standing
warning that "the biggest mistake would be trying to resurrect the binary
classifier." Its three routings, and what we did with them. Q1 (H1
redesign): abandon the benchmark-beat as primary - "structurally
unsalvageable", the -55.7 R2 being catastrophic instability at n=45 test,
not a tuning problem - and re-register ordinal severity prediction with
nested CV and pre-specified metrics against clinical and single-marker
baselines; post-hoc metric rescue, ordinal-vs-DEM benchmarking and
protein-subset benchmarking were each rejected with stated reasons.
Adopted as C1. Q2 (H2 routing): first choice miRNA-to-validated-target
direction-predicted enrichment (miRTarBase, CLIP-supported preferred,
explicitly not bare prediction-database fishing), second gene-level
convergence reframing, third placenta miRNA (biologically poor match),
and a standing demotion of monotonic-trend evidence to dose-response
support, never replication. Adopted as C2. Q3 (novelty): the multi-omic
severity module score with a pharmacogenomic druggability overlay, with
the claim sentence that now anchors section 9. Adopted as C3. The
round's throughline - let the failed classifier stay failed and rebuild
the claim on what the data actually support - is the reason this paper
has a section 8 at all.

## 6.8 Endpoint hierarchy (E8, ADDENDUM_5, locked 2026-09-27)
To bound researcher degrees of freedom in interpretation, the claims of
this report are ranked once, here, and every other section follows the
ranking:

PRIMARY claim (one only): an ordinal serum-miRNA severity model beats
the internal comparators - a clinical age/sex baseline and the best
single miRNA - on locked out-of-fold predictions (section 8: c-index
0.787; +0.167 over age/sex, bootstrap CI [+0.087, +0.242]; +0.074 over
the best single miRNA, CI [+0.009, +0.136]).

SECONDARY claim: the six-miRNA both-tissue core (the module), carried
with its stability caveats - D2 decorrelated-selection OOF c-index
0.704, E1 full-pipeline bootstrap median 0.615 with 3 of 6 members
stably re-discovered (9.5/9.6). Module membership is the output of the
two-stage sequential filtration (9.1) and nothing stronger.

TERTIARY / exploratory: everything else - target-program enrichment
(downgraded to database-bias-consistent, 9.5), the corrected F10
survivor test, druggability overlay, STRING coherence, concordance
scores, pathway signatures. None of these supports the primary or
secondary claim; they are reported for the record with their controls.
