# 9b. Discussion

## 9b.1 What was actually found, restated without the machinery
In the endpoint hierarchy of 6.8, the primary claim is the severity
model: an ordinal serum-miRNA model separates the control-to-severe
gradient of chronic Chagas cardiomyopathy at out-of-fold c-index
0.787, beating age+sex (+0.167, 95% CI excludes zero) and the best
single miRNA (0.714), on a cohort nobody in this program collected.
The secondary claim is the six-miRNA both-tissue core (miR-1-3p,
miR-122-5p, miR-192-5p, miR-30c-5p, miR-145-5p, miR-194-5p): the
output of a two-stage sequential filtration (2,632 tested, 20
frozen-gate candidates, 6 both-tissue members), tracked against the
severity spectrum with its stability caveats printed (9.5/9.6:
bootstrap median 0.615, 3 of 6 members stably re-discovered). Every
member passed a frozen association gate, a frozen novelty
screen, and a direction-predicted validated-target enrichment test
in two orthogonal cohorts. The module is the discovery carrier, not
any single molecule: no leave-one-out subset drops below 0.72, one
member is mildly dilutive and one redundant, and both facts are
printed in the paper (9.6, 9.7). Everything else - enrichment,
druggability, concordance - is tertiary and exploratory.

## 9b.2 The biological reading, held at the right strength
The candidate set carries a muscle-lineage thread (miR-1-3p,
miR-206, miR-145-5p, miR-199b-5p among the 21 passers) that is
consistent with progressive cardiomyocyte injury leaking into
serum, and the module's heaviest load-bearer (miR-122-5p, LOO
0.728) is liver-enriched in the literature - consistent with the
systemic metabolic involvement chronic Chagas is known to have,
and equally consistent with confounding by hepatic state, which
this cohort cannot exclude. The gate-(b) asymmetry (section 9.3b)
says the candidates' validated target programs are visibly engaged
in cardiomyocytes (18/20) and only sometimes detectable in blood
(6/20): exactly what a tissue-of-origin marker should look like,
read descriptively. The four-of-six non-monotone median shapes
(7.8) warn against any simple "more damage, more marker" story;
the ordinal model absorbs that shape rather than denying it. None
of these readings is a mechanism claim. The lane's own
falsifiability note (2.3b) applies to this paragraph first.

## 9b.3 What the replication cross-check is worth
The single REPLICATION row (miR-223-3p) is the strongest passing
association in the screen, and it is a molecule the GSE299582
depositors themselves named among their headline upregulated
markers (4.9). Two independent readings of the same cohort - the
source study's and this lane's frozen screen - agree on its most
visible signal. That agreement does not validate the novel
candidates, but it does say the screen is reading the same biology
the depositors read, which is the strongest internal sanity check
available without a new cohort.

## 9b.4 What stands between this and a usable marker
The gap is not statistical any more; it is clinical and practical,
and the paper says so in its limitations. In order: (i) external
validation in an independent severity-graded serum cohort - none
exists in the atlas today, and the rejected-series register (3.11)
records why the candidates for one failed; (ii) longitudinal
sampling - every dataset here is cross-sectional, so "progression"
is inferred from staged severity, not observed, and the language of
the claims keeps that distinction; (iii) calibration - the model
ranks well but its upper-threshold probabilities are badly shrunk
(8.5), so any decision-grade use needs a recalibration step the
lane has not done; (iv) assay transfer - miRNA-seq medians are not
a qPCR panel, and the zero-inflation named in 7.8 will look
different on a targeted assay; (v) the dilutive member - a
preregistered-panel decision about miR-192-5p belongs to a future
amendment, not to post-hoc editing (9.7). Each step is a named
piece of future work, not a caveat waved at the reader.

## 9b.5 What this project demonstrates even where the biology is young
Independent of whether the module survives external validation, the
lane demonstrates a working discipline for public-data biomarker
claims: 716 byte-verified records with per-series crosswalks, a
preregistration chain where every amendment predates the outcomes
it governs, a novelty screen that removes the strongest passing
result rather than let it anchor the paper, an honest-negative
first run that stayed in the record and forced the redesign, a
revalidation script that caught a live-service drift within hours,
and precision corrections (750 to 716; the 7.9 "strongest"
misstatement) that were published as corrections rather than edited
silently. The discovery claim stands on that discipline; the
discipline does not stand on the discovery claim.
