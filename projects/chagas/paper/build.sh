#!/bin/bash
# Assemble sections into a working TNR-style A4 PDF for the 50+ text-body
# page measurement. Minimal template (template_min.tex) - this TeX install
# lacks xcolor/booktabs/hyperref deps; toprule etc. are aliased to \hline.
# Working measurement build, not the final typeset.
set -e
cd "$(dirname "$0")"
ORDER="01_abstract 02_introduction 03_data_acquisition 04_cohort_atlas 05_methods_services 06_preregistered_analyses 06b_formulas 07_results_severity 08_benchmark 09_discovery 10_limitations 11_references_appendix"
for s in $ORDER; do cat "sections/$s.md"; printf '\n\\clearpage\n\n'; done > /tmp/chagas_paper_all.md
pandoc /tmp/chagas_paper_all.md -o manuscript_working.pdf --pdf-engine=pdflatex \
  --template=template_min.tex --toc \
  -M title="A provenance-first public-data compendium and preregistered comparator analysis for chronic Chagas disease biomarkers" \
  -M author="MEGA-PROGRAM-27 chagas project - working manuscript" \
  -M date="$(date +%Y-%m-%d) - not yet complete"
pdfinfo manuscript_working.pdf | grep -E 'Pages|Page size'
