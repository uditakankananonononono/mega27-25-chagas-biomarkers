# 6b. Numbered formulas (program gate: 10+)
(1) log2 CPM: x_ij' = log2( x_ij / sum_k(x_kj) * 10^6 + 1 )
(2) Kruskal-Wallis H: H = 12/(N(N+1)) * sum_g( R_g^2/n_g ) - 3(N+1)
(3) Benjamini-Hochberg FDR: q_(i) = p_(i) * m / i (step-up, monotone)
(4) Effect size: d = median(x'_severe) - median(x'_mild); gate |d| >= 0.5
(5) Ordinal concordance index: c = #{(i,j): sign(s_i-s_j)=sign(y_i-y_j)} /
     #{(i,j): y_i != y_j}
(6) Immediate-threshold ordinal model: P(y >= k | x) = sigma(a_k + b'x),
     score s = 1 + sum_k P(y >= k | x)
(7) Cox & Snell R2: R2 = 1 - exp( -2/n * (LL_1 - LL_0) )
(8) Elastic-net objective: min_b [ -LL(b)/n + lambda( alpha||b||_1 +
     (1-alpha)||b||_2^2/2 ) ], alpha = 0.5
(9) Direction-predicted enrichment p: p = (1 + #{b: f_b >= f_obs}) / 10001,
     f over 10,000 size-preserving random gene sets
(10) Severity module score: S_j = sum_m( w_m * z(x'_mj) ), w_m = sign(d_m),
     z = per-gene standardization over samples
(11) Bootstrap CI of delta c-index: delta = c_model - c_bench over B=1000
     resamples; beat requires the 95% percentile interval to exclude 0.

## 6b.1 Where each formula is used
Formulas (1)-(4) define the gate-(a) candidate screen on GSE299582:
per-miRNA Kruskal-Wallis across the four ordinal grades (2) on
log2(CPM+1) values (1), Benjamini-Hochberg correction across the tested
set (3), and the mild-to-severe median effect-size gate (4); the 20
passing candidates are listed in results/h2_gate_a_passing.csv. Formulas
(5)-(8) define the benchmark arm of section 8: the ordinal concordance
index (5) is the primary H1' metric, the immediate-threshold score (6) is
the locked model family, Cox & Snell R2 (7) is the failed H1 benchmark
metric of Run 1, and the elastic-net objective (8) belongs to the frozen
Run-1 pipeline. Formula (9) is the gate-(b) enrichment test of section 7:
the direction-predicted p-value over 10,000 size-preserving random gene
sets, run separately in GSE244827 (blood) and GSE203525 (hiPSC-CM) with
BH correction. Formulas (10)-(11) define the discovery deliverable and
its benchmark discipline: the signed severity module score (10) computed
per sample from the both-tissue core, and the bootstrap CI of the c-index
difference (11) whose exclusion of zero operationalizes the CLEAR BEAT
criterion of Addendum 3.

## 6b.2 Worked numeric examples (one per formula, from the artifacts)
(1) log2 CPM: sample OM1 (severe) carries miR-223-3p at raw value
2448.93 in the normalized-counts matrix against a library size of
234,682: log2(2448.93/234682 x 10^6 + 1) = 13.3493. (The matrix ships
pre-normalized non-integer counts; the transform is applied on top,
exactly as the locked script does.)
(2) Kruskal-Wallis: for miR-223-3p across control/mild/moderate/severe
(n = 42/37/37/30), H = 50.294, p = 6.92e-11 - the screen's strongest
passing association.
(3) BH: the screen's rank-1 p (miR-629-5p, p = 2.508e-15) gets
q = p x 2632/1 = 6.60e-12, matching the artifact exactly; it then
fails the effect gate (d = 0.03), which is why it is not a candidate.
(4) Effect size: miR-1-3p medians 7.494 (severe) - 6.926 (mild)
= 0.568 >= 0.5 - a pass; miR-629-5p's 0.029 is a fail at the same
gate.
(5) Concordance: the H1' model's out-of-fold c = 0.7872
(pooled; mean-of-folds 0.7872 to the same 4 decimals -
results/h1prime_locked_metrics_completion.json).
(6) Immediate-threshold model: the locked metric completion reports
per-class OvR AUCs 0.969 (control), 0.655 (mild), 0.576 (moderate),
0.661 (severe), macro 0.715, and calibration slopes 0.735 / 0.198 /
0.150 for the three thresholds - the worked example of why ranking
and decision-grade probability are different claims (section 8.5).
(7) Cox & Snell R2: the failed Run-1 binary arm recorded test
R2(C&S) = -55.72 out-of-sample (RUN_LOG) - the negative number that
retired H1's original metric arm and forced the Addendum-3 redesign.
(8) Elastic-net: the frozen Run-1 fit under this objective produced
the frozen ordinal macro-AUC 0.779 on the full signal (section 8) -
the number the H1' redesign had to be honest against.
(9) Permutation p: miR-1-3p's blood-cohort gate-(b) test observed
zero random sets at or above the observed fraction, so
p = (1+0)/10001 = 9.999e-5 - the minimum attainable at 10,000
permutations, reported as 0.0001, never as p = 0.
(10) Module score: the signed CORE6 score S reaches KW
p = 1.70e-11 and c-index 0.7557 across the four grades
(results/module_score.json); leave-one-out values 0.728-0.767 are in
results/module_score_stability.json.
(11) Bootstrap delta: against the clinical baseline the recorded
result is +0.167 c-index with 95% CI [+0.087, +0.242] (section 8) -
the interval excludes zero, which is what CLEAR BEAT means
operationally in Addendum 3.
