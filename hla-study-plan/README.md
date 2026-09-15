# HLA pangenome study plan — v1.0.0

Read [the PDF](plan.pdf) or [the LaTeX source](plan.tex).

## Build

Requires a TeX installation with `latexmk`, `pdflatex`, and the packages declared
in `plan.tex`.

```sh
cd hla-study-plan
make
```

`make clean` removes auxiliary files and preserves the PDF.

## Document versioning

`VERSION` is the single source for the version printed in the document.
Semantic versioning applies to the study specification:

- **Major:** changes to the principal question, endpoint, or evaluation design
  that change the interpretation of planned comparisons.
- **Minor:** additional compatible experiments, analyses, or substantive detail.
- **Patch:** corrections or editorial changes without a design change.

Update `VERSION` and `CHANGELOG.md`, rebuild the PDF, and commit them together.
Document releases use namespaced Git tags, beginning with `hla-plan-v1.0.0`.

## Review basis and decisions

Reviewed the tracked workflow at repository commit `9b298d8` and the local,
untracked `hla-viz/` analysis snapshot on 15 September 2026. The latter is active
work and is external to this document release. Relevant local files were:

- `hla-viz/hla_calls.tsv`: 610 unique haplotype identifiers, grouped by cohort.
- `hla-viz/figures/CAPTIONS.md`: extraction/annotation provenance, assembly QC
  observations, and existing bundle analyses. Its biological conclusions are
  prior analysis claims requiring their underlying validation in the study.
- `hla-viz/gene_bundles/HLA-A.ctg.summary.tsv`: example bundle summaries.
- `workflows/variant-calling/main-vg.cwl` and `vg-giraffe.cwl`: KMC-based
  sampling and Giraffe v1.54.0, with projected BAM output.

Local tables contain APR 106, HPRC r2 464, Japanese JaSaPaGe 20,
Saudi JaSaPaGe 18, and reference 2, totalling 610 haplotypes. Robert confirmed
that adding Korea and CRC brings the intended panel above 700. These additions
require reconciliation with the local accession manifest. HPRC includes Asian
and non-Asian donors; quantify the final ancestry composition per donor. The
study emphasises Asian diversity within this mixed-ancestry panel.

The principal design refinements are:

1. Separate whole-donor subsets from the cited paper's blockwise haplotype
   sampling, and preserve graph-native alignments for inference.
2. Exclude test donors and relatives; verify experimental truth provenance and
   resolution. DDBJ files were not inspected for this document.
3. Preserve mates and unmapped reads, assess hg38 recruitment loss, and use a
   fixed genomic background to assess off-target alignment.
4. Separate haplotype reconstruction, sequence annotation, and each typer's
   native read pipeline; account for internal realignment and supported inputs.
5. Test real long reads on matched individuals, with explicit coverage and
   technology strata, and distinguish gene-level from extended-MHC phasing.
6. Use label-blind bundle clustering before HLA-label comparison; develop
   functional interpretation through binding profiles and independent evidence.

Citations and the complete experiment specification are in the LaTeX document.
This release contains a study plan; the proposed experiments have not been run.
