# 10. Limitations and the honest-negative register
This project's claims are bounded by design. The severity model is
cross-validated on one 146-sample cohort, not prospectively validated;
no second severity-graded serum miRNA cohort exists to hold out. The
benchmarks it beats are internal comparators (clinical covariates,
single markers), fairly chosen but not external champions. The frozen
binary benchmark-beat failed on Run 1 (test R2 -55.7 vs 0.688; commit
49de4ae) and is kept as a documented audit, not deleted - the successful
ordinal model answers a different, registered question. The literature
exclusion screen covers 121 PubMed abstracts plus abstract-named markers
from the source paper; the source paper's full 40-DEM list is paywalled,
so a residual chance remains that a "candidate" appears in that list
(logged in Addendum 2). The methylation and pharmacogenomic cohorts are
context, not endpoints. Acute-phase sampling is absent from every public
human series we could verify. Gate (b) tests experimentally validated
targets only; candidates without such support drop out rather than
borrowing weaker evidence. ChatGPT/DeepSeek/Gemini judge rounds are
simulated critiques used for redesign, not peer review; every adopted
change is locked by us in a timestamped amendment.

## 10.1 Register additions from the gate-(b) execution (2026-09-27)
Gate (b) replication is compartment-asymmetric: strong in the
patient-derived cardiomyocyte cohort (18/20 candidates) but partial in
peripheral blood (6/20; three candidates show zero predicted-direction
signal). The discovery module therefore rests on a 6-miRNA both-tissue
core, not the full candidate set. The module score is computed in-sample
on the same cohort that selected its members and weights - it describes,
it does not validate. The 18-member superset fails as a signed module
(c-index 0.509) and is reported as a negative. miRTarBase was recovered
as v8.0 (2022), not the current release: validated-target sets may miss
post-2022 evidence. GSE244827 labels ride on a B-code bijection verified
against live GEO SOFT records, but re-fetched SOFT pages hash differently
than acquisition-time records (dynamic page content); the metadata lines
used are stable and mutually consistent. Target-program druggability
(Open Targets tractability) does not imply the miRNAs themselves are
druggable. Enrichr pathway context is post-hoc and descriptive, outside
the preregistered gates.

## 10.2 Compendium-level limitations
The 716 records are records, not independent studies: they span 15
independent series, and the GSM samples are nested within their GSE
parents. Any reader who treats the record count as a cohort count will
overstate the evidence base by two orders of magnitude; the per-disease
gate counts records, and we report both numbers together for exactly this
reason. Provenance is single-source: every record comes from NCBI GEO, so
a GEO-side metadata error that is internally consistent would survive our
checks - byte-hashing proves fidelity to the source, not truth of the
source. De-duplication trusts GEO's super-series structure: the
GSE191083 catch worked because the container was declared, and an
undeclared overlap of the same samples across unrelated series would be
harder to detect at series level. The 750->716 correction (section 3.6)
is itself evidence that expansion-time counts are fragile; the standing
rule - recompute from frozen crosswalks, never quote status text - exists
because this failure already happened once in this lane. Live
re-verification samples 5-10% of records per run, not the whole corpus;
an unsampled record could in principle have changed upstream, though the
hermetic hash check still bounds what we use to what we recorded.

## 10.3 Statistical limitations
Every performance number in sections 7-9 is computed inside one cohort
family. Nested 5x3 cross-validation with feature selection inside the
outer fold controls the obvious optimism, but the out-of-fold predictions
remain pooled across resamples of the same 146 donors; the bootstrap
confidence intervals on their differences measure resampling stability,
not transportability. The severe grade is 30 samples, so any
severe-specific subgroup claim is underpowered by construction. The
gate-(b) enrichment null is size-preserving over random gene sets; it
does not preserve the correlation structure of co-expressed genes, so its
p-values are calibrated for the question "would an arbitrary set of this
size do as well" and not for "would a correlated set do as well". The
module score's members and weights come from the same cohort that scores
them (in-sample, descriptive); the KW p=1.7e-11 and c-index 0.756 for
CORE6 quantify separation inside GSE299582 only. Kruskal-Wallis across
four grades is a heterogeneity test; monotonicity is assessed separately
and reported as dose-response support, not as a replication event.

## 10.4 Service and evidence-currency limitations
Three services carry version or access caveats that bound their evidence.
miRTarBase is pinned to v8.0 (2022) via a Wayback snapshot of the
publisher URL because the live site returns 404/400; validated-target
sets may miss post-2022 evidence, and the pin is disclosed in the
judge-round record. Enrichr results are tied to a specific userListId and
library versions at run time; reruns can drift as libraries update. Open
Targets tractability and target-disease association scores are
live-service outputs captured as JSON snapshots; the service's underlying
data releases move, so the druggability overlay is a dated snapshot, not
a standing fact. WHO/CDC burden figures are retrieved-and-hashed web
documents whose numbers update on the source's schedule.

## 10.5 What would close each gap
In order of leverage: (i) an external severity-graded serum miRNA cohort
would convert the section-8 internal beat into a transport claim - none
is public at sweep time; (ii) a longitudinal seropositive cohort with
adjudicated cardiac outcomes would unlock the progression-marker question
that section 2.3b shows is currently out of reach; (iii) a miRTarBase
version refresh, once the live service is reachable, would tighten gate
(b) without changing its rule; (iv) same-molecule placenta or blood
miRNA cohorts would let the discovery arm test replication in the strict
sense Addendum 3 currently routes around. Each item is a named, checkable
next step, not a rhetorical hedge: if any becomes available, the lane's
standing rules require it to be acquired, hashed, and registered before
outcomes are examined.

## 10.6 Live-service numbers drift; the lane's answer is a revalidation record, not a frozen screenshot
Three results in this paper depend on live external services, and one
of them demonstrably moved during a single night: the STRING v12
network, re-executed six hours after the committed run, returned
1,296 edges where the committed artifact holds 1,051 (341 mapped
nodes vs 299; expected edges 523 vs 418; enrichment p < 1e-16 in
both). The mygene validation and the Open Targets tractability
buckets reproduced exactly under re-execution (344/352 with the same
8 legacy notfounds; 332 tractability-positive and 49 Approved-Drug
with identical gene sets), but "reproduced today" is not "will
reproduce next quarter" for any live service. The lane's mitigations
are structural rather than rhetorical: raw responses are stored with
sha256 (sources/services/, EVIDENCE_SHA256.txt), every committed
number names its service version or retrieval timestamp where the
service exposes one, the re-execution script is committed
(scripts/service_revalidation.py) so any reviewer can re-run the
comparison, and drift is logged in RUN_LOG rather than silently
absorbed. Where a live number is load-bearing for a claim (the
STRING interconnectivity sentence), the claim is worded so that both
the committed and the re-executed values support it.
