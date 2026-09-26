# 2. Introduction: the burden and the biomarker gap
## 2.1 The disease
Chagas disease (American trypanosomiasis, Trypanosoma cruzi) affects
roughly 6-7 million people, mostly in Latin America (WHO fact sheet,
retrieved 2026-09-26, sources/services/who/). After an often-unnoticed
acute phase, infection persists for life; 20-30% of those infected
eventually develop chronic Chagas cardiomyopathy, the form that kills.
The tragedy is mechanical: patients feel well while fibrosis and
conduction damage accumulate, and by the time symptoms declare, the
myocardium is already remodeled. The CDC's diagnostic guidance
(sources/services/cdc/) confirms the tools that exist - serology, PCR,
imaging - answer "infected?" and "damaged?" but not "progressing?".
## 2.2 The gap this project attacks
A marker that tracks progression would change triage: who needs annual
echocardiography, who needs treatment escalation, who is safe to watch.
The published record offers fragments - protein panels (galectin-3,
BNP/NT-proBNP, MMPs; the frozen screen in prereg/), a prognostic ELISA
panel (PMID 34479416), and the first severity-graded serum miRNA survey
(Roma et al. 2026, PMID 41574750, whose data this project re-analyzes
under registration) - but no validated multivariate severity model, and
no cross-modal mechanistic bridge from circulating miRNAs to cardiac
biology. Open Targets associates 890 targets with the disease
(sources/services/opentargets/), yet translation into a progression
marker has not happened.
## 2.2b What "already named" means here
A discovery claim is only as honest as its exclusion list. Before any
candidate was examined, we froze a literature screen (PubMed
"Chagas"[Title/Abstract] AND "biomarker"[Title/Abstract], relevance-sorted,
121 records screened in full on 2026-09-26; abstracts and IDs preserved in
prereg/) and hand-curated the markers it names: protein and clinical
markers (galectin-3, MMP-2, MMP-9, the TIMP family, BNP, NT-proBNP,
troponin I/T, CK-MB, ApoA1, IL-6, IL-10, IL-17A, TGF-beta, CCL2, CCL17,
MPO, sST2, CRP), miRNAs already named in Chagas contexts (miR-146a,
miR-208a), parasite-side markers (kDNA PCR load, cruzipain-derived and
F29 antigens) and cellular/serological profiles (CD4/CD8 responses, IgG1).
The source paper of our primary cohort adds nine further miRNAs
(miR-143-3p, miR-223-3p, miR-486-5p, miR-3960, miR-6734-5p, miR-1285-5p,
miR-10527-5p, miR-1228-5p, miR-30c-3p) to the replication class; its full
40-marker list is paywalled, and that residual risk is logged rather than
ignored (Addendum 2, B1). Any candidate matching this register - by exact
symbol or unambiguous synonym - is classified as replication, never as
new discovery. The screen bounds "already named"; it cannot prove global
novelty, and we say so in its own header.

## 2.3 The approach
Three commitments distinguish this work. Provenance first: every one of
the 716 records is individually retrievable and byte-hashed, because
secondary analyses fail silently when their inputs are assumed rather
than verified (section 3). Registration before outcomes: hypotheses,
splits, seeds, metrics, comparator rules and exclusion screens were
committed before any outcome was computed, and every amendment is
timestamped with its reason (section 6). Honest negatives: one frozen
classifier failed hard and stands in the record; the successful ordinal
model was designed and locked after that failure was documented, not
instead of documenting it (sections 6-8). The discovery arm then asks the
question that matters: do the severity-linked miRNAs point, through
experimentally validated targets, at cardiac remodeling programs visible
in independent cohorts (section 9).


## 2.3b What would refute the contribution
The strongest version of this project is a progression marker with
measured added value over existing parasite-load, immune-protein and
clinical predictors in baseline seropositive asymptomatic patients. That
version is currently refuted as a claim, and the refutation is part of
the design: it would require training only on baseline asymptomatic
participants whose future cardiac outcome is adjudicated, an untouched
longitudinal cohort, and head-to-head comparison against parasite-load
and inflammation-protein prognostic models on the same people and
endpoint. No public cohort of that shape exists in this compendium; the
33 cross-sectional GSE244827 blood libraries cannot substitute for a
time-to-event outcome. The published record constrains ambition further:
a 2024 prospective study already found baseline parasite DNA and 47
immune proteins associated with ten-year decline (21 progressors, 31
controls; PMID 38203212), and a 499-person seropositive donor cohort
established ten-year cardiac outcomes without an omics crosswalk
(PMID 23393012). What remains - and what this paper actually claims - is
narrower: a provenance-complete public compendium, a locked ordinal
severity model that beats its molecule-matched internal comparators
(section 8), and a cross-modal mechanistic bridge tested under
registration (section 9). NOVELTY_PLAN.md states the standing verdict
plainly: no new named biological discovery, no validated new method, no
demonstrated same-task advantage - an unproven research direction with
its falsification criteria written down.

## 2.4 Host biology: the pathway the parasite hijacks
The KEGG Chagas disease pathway (hsa05142, retrieved with hash into
sources/services/kegg/) frames the cardiac mechanism this project reads
in molecular data: T. cruzi invades cardiomyocytes, activates Ca2+
signaling through cruzipain, oligopeptidase B and trans-sialidase,
escapes the parasitophorous vacuole via TcTOX, replicates in the cytosol,
and drives cardiomyocyte hypertrophy while disturbing T-cell responses.
The pathway's human gene set (102 genes in the retrieved flat file) is
enriched for calcium-handling, immune and remodeling terms - the same
axes on which the severity-linked miRNA programs of section 9 land
(PI3K-Akt, focal adhesion, VEGFA-VEGFR2 among their validated targets).
This grounding matters for claim discipline: a serum miRNA module that
tracks severity is biologically plausible exactly to the extent its
target programs intersect these host processes, and that intersection -
not the association alone - is what the discovery arm tests.

## 2.5 Scope within the wider program
This paper is the Chagas project of a ten-disease program (MEGA-PROGRAM-27)
whose shared core manuscript covers the cross-disease machinery and the
negative-results register. The program's gates are per-disease: 120
records, 40 distinct external services, 50 substantive pages, a
benchmark-beat and a discovery result, each met within the disease
project, not pooled across projects. Where this paper quotes a count -
716 records, 40 services - it is the Chagas lane's own verified number,
computed from this repository's frozen artifacts (section 3.6); shared-core
counts do not transfer, and no gate is claimed met until it is verifiably
met here. The compendium, the preregistered analyses and this manuscript
live together in one public repository so that every number in these
pages can be traced to the artifact that produced it.
