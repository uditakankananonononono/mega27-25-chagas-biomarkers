# 1. Abstract
Chagas disease kills through a slow, silent progression from asymptomatic
infection to chronic cardiomyopathy (CCC), and the field still lacks
validated markers that track its severity spectrum. This project builds a
provenance-first compendium of 716 individually byte-verified public
records across 15 human Chagas series - serum miRNA, blood and tissue
RNA-seq, single-cell, spatial, methylation and pharmacogenomics - and runs
two preregistered analyses on the only severity-graded serum miRNA cohort
(GSE299582, n=192). First, a locked ordinal severity model
(immediate-threshold logistic, nested cross-validation) tracks the
control-to-severe gradient at out-of-fold concordance 0.787, beating both
a clinical age/sex baseline (+0.167, bootstrap CI [+0.087, +0.242]) and the
best single miRNA (+0.074, CI [+0.009, +0.136]). Second, a frozen
novelty-screened association analysis isolates 20 severity-linked miRNAs
absent from the screened Chagas biomarker literature, with a
muscle-lineage thread (miR-1, miR-206, miR-145, miR-199b) consistent with
progressive cardiomyocyte injury. The registered orthogonal replication
has now run: direction-predicted enrichment of miRTarBase-validated
target programs passes for 18/20 candidates in patient-derived
cardiomyocytes but only 6/20 in peripheral blood - an honest
compartment split, and one that matched-popularity controls show is
consistent with database bias rather than candidate-specific validation
(interpretation downgraded accordingly, 9.5). Six miRNAs (miR-1-3p, miR-122-5p, miR-192-5p,
miR-30c-5p, miR-145-5p, miR-194-5p) replicate in both tissues - the
output of a two-stage sequential filtration (2,632 serum miRNAs tested,
20 frozen-gate candidates, 6 both-tissue core members) - and form
a signed severity module that tracks the gradient monotonically
(in-sample c-index 0.756; full-pipeline bootstrap median 0.615, 3 of 6 members stably re-discovered - 9.6); their 352 strong-support validated targets
are significantly interconnected (STRING, p < 1e-16) and 49 are
approved-drug targets (Open Targets). Every count, hash,
run and negative - including one failed classifier kept as a documented
audit - is committed and re-verifiable end to end.
