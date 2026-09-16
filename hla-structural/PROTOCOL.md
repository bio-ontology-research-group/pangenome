# Structural MHC typing: prospective technical pilot

Authorised 2026-09-16. Starting panel and completed HLA analysis remain unchanged.

Primary endpoints: diploid component count and unordered pair of structural
signatures, evaluated separately at RCCX and DRA–DRB1. The DRB endpoint includes
annotated expressed DRB paralogs; it is not a complete pseudogene or breakpoint
classification. RCCX homolog names are not functional CYP21/TNX assignments.

Extract graph walks at assembly-derived coordinates, verify exact sequence
agreement and exclude gaps. Retain sequence hashes, path steps, donor/family,
source QC, annotation limits and available matched CRAMs. Labels are assembly
comparators, not experimental truth. Missing annotations remain provisional.

Prototype: infer pairs of observed local paths from quality-filtered canonical
31-mer read counts, retaining within-path multiplicity. Compare against binary
path features and independent copy-number marker estimates. A flat sequence
catalogue with identical markers must produce identical predictions: do not claim
that graph storage alone improves inference. The supplied graph was built using
all donors; reference paths and marker eligibility must exclude query families.
Graph topology remains a fixed-panel representation, not an independently rebuilt
alignment. Marker definitions depend only on retained training sequences at each
fold. No phenotype data will be used.

Use deterministic selection before viewing read results. Include donors with
both gap-free haplotypes and local recruited CRAMs; sample across diploid RCCX
counts, plus the known DRB5 absence control. Describe enrichment and report all
selected donors, including failed/no-call cases. Original WGS recruitment and
whole-genome paralog specificity are outside this pilot's guarantees.

Represented structures and structures absent from training are separate strata.
No forced exact arrangement when the read features are non-identifying. Score
margins are not calibrated probabilities. The phased pair is unordered: no claim
of phase relative to distant HLA loci or independent chromosome assignment.

Deliverables: frozen catalogue, graph/sequence concordance audit, reproducible
read measurements, family-excluded predictions and baselines, ambiguity and
failure analysis, report and recommendation for further assay validation.

## Locked extension before prediction inspection

The initial 18 selected donors remain the pilot set. All remaining 88 eligible
donors are an additional evaluation set, chosen before computing any typing
predictions. The same marker filters, normalization, optimization and independent
baselines apply without outcome-driven tuning. Per-set results must be reported.
`source/pre_evaluation_freeze.json` records the implementation snapshot.

## Pilot-trained secondary model, frozen before validation predictions

The 18-donor pilot showed systematic RCCX overdosage. Retain the unchanged
primary benchmark. A secondary model estimates one multiplicative dosage factor
as the median pilot conserved-C4 estimate divided by assembly C4 copy number,
then constrains path pairs to the rounded corrected dosage. Its evaluation is
restricted to the remaining 88 donors. Exclude any query-family pilot record
from calibration. `source/dosage_calibration.json` freezes the factor, records,
code hash and timing. This is supervised technical calibration, not novel SV
discovery, and its results must be reported separately from the original method.
