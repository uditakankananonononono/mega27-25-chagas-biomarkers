# 4. Cohort atlas: fifteen annotated series

Each entry: accession, n samples, assay, biological question, label status,
and the role the cohort plays in the preregistered analyses (prereg/
PREREGISTRATION.md). PMIDs were verified per series in ACQUISITION_LOG.md.

## 4.1 Severity-graded blood cohorts (analysis-grade)
- GSE299582 (n=192, serum miRNA-seq). The largest acquired cohort and the
  primary severity dataset: chronic Chagas cardiomyopathy (CCC) graded
  mild/moderate/severe plus controls. Anchors preregistered hypothesis H1
  (benchmark-beat severity classification). PMID 41574750. Composition
  from the frozen per-sample characteristics (all 192 rows carry age,
  gender, disease-form and severity fields): 104 CCC samples graded mild
  (37), moderate (37) and severe (30); 46 indeterminate-form samples; 42
  non-chagasic controls. Sex split 103 female / 89 male; recorded ages
  span 30-83 years. The near-balance of mild and moderate grades and the
  smaller severe grade shape every severity-graded analysis downstream:
  ordinal tests have reasonable support across the full grade ladder, but
  pairwise severe-vs-control contrasts run on 30 cases and are reported
  with that power limitation attached.
- GSE244827 (n=33, whole-blood RNA-seq; prior-tagged under P31, retained
  as provenance). Asymptomatic/early-CCC vs seronegative; orthogonal
  replication cohort for H2 cross-modal convergence. PMID 40290486.
  Crosswalk composition: cardiac-stage field CARD 16 / non-CARD 17 with
  sex 19F/14M recorded per donor; the label column (10 case / 23 control)
  follows the frozen crosswalk, and the matrix-column identity is the
  verified B-code bijection of section 3.10.

## 4.2 Cardiac tissue and cardiomyocyte models
- GSE84796 (n=17, expression array; prior-tagged). Cunha-Neto CCC heart
  tissue study; the comparator anchor for the leish-lane audit pattern and
  the historical heart-vs-blood sign analysis whose fragility on
  leave-one-person-out is documented in NOVELTY_PLAN.md.
- GSE203525 (n=20, RNA-seq, patient hiPSC-derived cardiomyocytes). CCC vs
  indeterminate patient lines, with and without T. cruzi reinfection;
  graded by donor clinical status. PMID 35873155. Crosswalk composition:
  12 CCC-donor and 8 indeterminate-donor cardiomyocyte samples, split
  10 non-infected / 10 infected with T. cruzi Y strain. This is the
  gate-(b) cardiac cohort; the analysis contrast (CCC vs indeterminate
  at 0hpi-equivalent, 6 vs 6) is asserted against the live header
  before every run (section 3.10).
- GSE129676 (n=16, RNA-seq, hiPSC-CM). Chagas-patient vs control
  cardiomyocyte infection timecourse. Timecourse design labels are
  honest treatment contrasts, not clinical labels. PMID 31105048.
  Crosswalk composition: 8 Chagas-disease and 8 control donor lines,
  four libraries each at 0h, 24h, 48h and 72h post-infection - the
  only time-resolved cardiomyocyte design in the atlas.
- GSE348071 (n=32, RNA-seq, AC16 cells + patient iPSC-CM). DHODH R135C
  mitochondrial vulnerability in CCC; stimulation contrasts
  (e.g. IFN-gamma) labelled as such. Unpublished at acquisition time;
  flagged for citation monitoring. Crosswalk composition: 16 AC16 and
  16 patient iPSC-CM libraries across three genotype arms - 8 reference
  C/C, 16 heterozygous p.Arg135Cys (C/T), 8 CRISPR-corrected C/C - and
  three 48-hour treatments (16 untreated, 8 IFN-gamma, 8 IFN-gamma +
  TNF-alpha). The CRISPR-corrected arm is the design's strength: it
  isolates the variant's effect from donor background.

## 4.3 Congenital transmission
- GSE311812 (n=46, RNA-seq + Visium spatial). Maternal blood and placenta
  with transmitter/non-transmitter contrast; the only spatial dataset in
  the compendium. PMID 41648170. Crosswalk composition: 24 peripheral
  blood, 20 central-placenta (fetal side) and 2 central chorionic-villi
  records. The transmitter count discrepancy (6 in titles vs 5 in the
  design text) is resolved in section 3.9; labels follow the titles.
- GSE333874 (n=31, small RNA-seq, placenta). Congenital-transmission
  miRNAs; tissue-restricted complement to the serum miRNA severity cohort.
  PMID 42523576.
- GSE107376 (n=9, expression array, middle-section placenta). Seropositive
  vs seronegative mothers; smallest cohort, retained for transmission-theme
  completeness with its power limitation stated. PMID 29545200. Crosswalk
  composition: 6 seropositive mothers (3 PCR-positive, 3 PCR-negative) and
  3 seronegative mothers - serostatus and PCR status are recorded
  separately, so the tiny n does not force the two into one axis.

## 4.4 Innate immune response models
- GSE158986 (n=12, RNA-seq, monocyte-derived dendritic cells). Human
  first-contact response to T. cruzi; treatment-contrast labels.
  PMID 33897690. Dual-platform deposit (6 libraries on GPL16791, 6 on
  GPL29219, section 4.6) with donor and replicate fields preserved per
  library; the only series in the atlas without a GEO series-matrix
  file (supplementary counts only, section 3.8).
- GSE328447 (n=4, small RNA-seq, THP1 macrophages). isomiR response in an
  infection model; retained as an exploratory isomiR lead with explicit
  small-n caution. PMID 42614816. The recorded treatment contrast is
  verbatim: miR-1246+1 mimic transfection vs control mimic transfection
  (2+2) - a mimic-intervention design, labelled treated/control and used
  as nothing more than an exploratory lead.
- GSE295194 (n=16, scRNA-seq PBMC with sample tags). CCC vs indeterminate
  CD4 T-cell peptide response; single-cell modality. PMID 40391216.
  Crosswalk composition: 16 wild-type PBMC libraries with batch, cell-type
  and genotype fields; all 16 carry the case label with the
  donor-status/peptide-response contrast inside the sample records
  (section 4.6b).

## 4.5 Methylation and pharmacogenomics (context modalities)
- GSE191081 (n=22, heart LV-wall RNA-seq) and GSE191082 (n=158,
  methylation tiling array: 144 blood, 14 heart). The CCC methylation
  study pair; 180 GSMs acquired with the super-series GSE191083
  deliberately dropped after the uniqueness check flagged it as a
  container that would double-count its children. GSE191081's group
  field is three-way - 8 chronic chagasic cardiomyopathy, 8 dilated
  cardiomyopathy, 6 non-chagasic control - and the label column
  collapses it to case/control (CCC vs the rest), so any disease-
  specificity question must read the group field, not the label.
  GSE191082 labels 104 case / 54 control across its two tissues.
  GSE154421 (n=92, SNP pharmacogenomics) covers benznidazole-response
  genotypes - a treatment-response modality orthogonal to every
  expression cohort; its real contrast, adverse reaction to
  benznidazole (63 yes / 29 no), lives in the characteristics, and all
  92 donors are Chagas patients. These three series are context for
  discussion, not expression endpoints.

## 4.6 Atlas-level properties
Disease spectrum: indeterminate/asymptomatic, graded CCC (mild to severe),
congenital transmission, and in-vitro infection models - the full natural
history except acute-phase sampling, which no public human series offered
at sweep time (gap logged). Tissue breadth: heart, blood, serum, placenta,
PBMC, DC, macrophage, cardiomyocyte lines. Technology breadth: five assay families across nine
platforms, enumerated from the frozen crosswalks rather than from series
text: GPL14550 (GSE84796 array), GPL16791 (GSE107376, GSE328447 and the
human arm of GSE158986), GPL17301 (GSE129676), GPL18573 (GSE203525),
GPL21145 (GSE191082 methylation), GPL24676 (GSE191081, GSE244827,
GSE295194, GSE311812, GSE333874), GPL28868 (GSE154421 SNP array),
GPL29219 (the second arm of GSE158986's dual-organism design, 6 samples
against 6 on GPL16791), and GPL30173 (GSE299582, GSE348071). The atlas is deliberately heterogeneous:
the preregistered convergence test (H2) uses that heterogeneity as the
replication filter rather than treating it as noise to be normalized away.

## 4.6b Label vocabulary by design family
Labels are assigned per series by scripts/label_crosswalks.py under
explicit per-series rules, with zero unmapped rows at verification time.
Clinical-cohort series use case/control (with the underlying severity or
form retained in characteristics); in-vitro mechanism series use the
honest treatment vocabulary instead: infected/control (GSE158986),
treated/control (GSE328447), variant/reference (GSE348071). Three series
carry their primary contrast inside the sample record rather than in the
label column: GSE154421 (all 92 rows are benznidazole-treated Chagas
patients; the adverse-reaction yes/no field is the real contrast),
GSE203525 (20 CCC/indeterminate patient lines; the line-status and
reinfection contrasts live in characteristics), and GSE333874 (all 31
rows placental; the transmitter/non-transmitter contrast lives in
characteristics). A reader who groups naively by the label column will
misread these three series; the crosswalk characteristics are the
authoritative record.

## 4.7 Prior-tagged series (provenance only, not re-counted)
GSE84796 and GSE244827 remain in the crosswalk set so that a reader can
verify the prior 50 samples byte-for-byte; per the count-correction record
(commit ef7b29a) they contribute to the 716-record total exactly once,
through the prior manifest.

## 4.8 Machine-checked per-series table (generated from the frozen crosswalks)
| Series | n GSM | Platform | Label counts |
|---|---|---|---|
| GSE107376 | 9 | GPL16791 | case 6, control 3 |
| GSE129676 | 16 | GPL17301 | case 8, control 8 |
| GSE154421 | 92 | GPL28868 | case 92 (genotype contrast within) |
| GSE158986 | 12 | GPL16791 + GPL29219 | infected 6, control 6 |
| GSE191081 | 22 | GPL24676 | case 8, control 14 |
| GSE191082 | 158 | GPL21145 | case 104, control 54 |
| GSE203525 | 20 | GPL18573 | case 20 (line-level contrasts within) |
| GSE244827 | 33 | GPL24676 | case 10, control 23 |
| GSE295194 | 16 | GPL24676 | case 16 (donor-status tags within) |
| GSE299582 | 192 | GPL30173 | case 150, control 42 |
| GSE311812 | 46 | GPL24676 | case 25, control 21 |
| GSE328447 | 4 | GPL16791 | treated 2, control 2 |
| GSE333874 | 31 | GPL24676 | case 31 (transmitter contrast within) |
| GSE348071 | 32 | GPL30173 | variant 16, reference 16 |
| GSE84796 | 17 | GPL14550 | case 10, control 7 |
Every row regenerates from sources/<GSE>_sample_crosswalk.csv
(sha256-hashed per-GSM GEO SOFT records); the table is a view, not a
source of truth.

## 4.9 What the source studies themselves concluded
The depositors' own findings, paraphrased tightly from the summary
field of each series record (sources/new_series_ledger.csv, verbatim
GEO text; GSE84796 from the prior-tag record). This is the
related-work thread the lane builds on: where a source study already
named a molecule, the lane's novelty screen classifies it as
REPLICATION, not discovery.

- **GSE299582** (Roma et al., J Infect Dis 2026): 40 differentially
  expressed miRNAs between Chagas disease patients and healthy
  controls (miR-199b-5p, miR-153-3p, miR-143-3p, miR-223-3p
  upregulated; miR-150-3p, miR-4508, miR-486-5p, miR-3960
  downregulated), with severity-trending candidates named
  (miR-6734-5p, miR-1285-5p, miR-10527-5p, miR-31-5p, miR-5187-5p,
  miR-6515-5p higher in severe; miR-30c-2-3p lower). Direct
  cross-check: the lane's single REPLICATION row is miR-223-3p, one
  of the source study's headline upregulated markers - the frozen
  novelty screen agrees with the depositors' own list.
- **GSE244827** (Duque et al., Lancet Reg Health Am 2025): early CCC
  is associated with peripheral downregulation of immune-response
  genes (reduced antigen presentation and T-cell activation),
  distinct from early cardiomyopathy in Chagas-negative patients -
  the study that established blood-based early-CCC signal exists.
- **GSE311812** (Bolivian congenital-transmission cohort): integrative
  blood + placenta transcriptomics comparing infected
  transmitter vs non-transmitter mothers; the transmitter-specific
  maternal blood signature is the depositors' core finding (title
  audit in section 3.x found 6 transmitter blood samples vs 5 stated
  - recorded honestly in ACQUISITION_LOG).
- **GSE333874**: nine candidate placental miRNAs associated with
  congenital transmission status (DESeq2, BH FDR) - a placenta-side
  miRNA lead set orthogonal to this lane's blood-severity focus.
- **GSE348071**: a rare heterozygous DHODH variant (p.Arg135Cys)
  sensitizes cardiomyocytes to IFN-gamma-driven mitochondrial
  dysfunction - a gene-environment mechanism for why only ~30% of
  infected individuals develop CCC.
- **GSE203525** (Oliveira et al., Front Cell Infect Microbiol 2022):
  patient-derived hiPSC cardiomyocytes from CCC vs indeterminate
  donors show different transcriptional responses to T. cruzi
  reinfection - the deposit behind gate (b)'s cardiomyocyte cohort.
- **GSE129676** (Bozzi et al., Stem Cell Reports 2019): hiPSC-CMs as
  a T. cruzi infection model; human cardiomyocyte response
  timecourse.
- **GSE158986** (Gil-Jaramillo et al., Front Immunol 2021): 468
  differentially expressed genes in human dendritic cells at 12h
  first contact with infective T. cruzi forms.
- **GSE295194** (Souza-Silva et al., Front Immunol 2025):
  DERAA-motif HLA-DRB1 alleles (*0103, *0402, *1301, *1302) linked
  to severe cardiomyopathy; CD4+ T-cell peptide response profiled by
  scRNA-seq.
- **GSE107376** (Juiz et al., Am J Pathol 2018): placental gene
  expression differences between seropositive and seronegative
  mothers (pooled RNA-seq, 9 pools of 2).
- **GSE328447** (Lyu et al., Comput Struct Biotechnol J 2026):
  5'-isomiR dysregulation in the T. cruzi-macrophage context; the
  deposited design is a miR-1246+1 mimic transfection contrast
  (section 4.4) - exploratory lead only.
- **GSE154421**: SNP scan for benznidazole adverse-reaction
  pharmacogenomics; unpublished at acquisition time.
- **GSE191081 / GSE191082** (Brochet et al., Front Immunol 2022):
  bulk RNA-seq + EPIC methylation in left-ventricular tissue;
  epigenetic regulation of transcription-factor binding motifs in
  CCC vs dilated cardiomyopathy.
- **GSE84796** (Cunha-Neto CCC heart array, prior-tagged): carried
  for provenance continuity with the parent program; the lane holds
  the crosswalk, not a fresh analysis, and claims nothing new from
  it.
