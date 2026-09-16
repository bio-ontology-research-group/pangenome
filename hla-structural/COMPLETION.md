# Structural typing experiment completed

The authorised technical objective was to construct and benchmark structural
typing at RCCX and DRB, identify recoverable copy-number/arrangement features,
compare simpler approaches, and connect structure to classical HLA labels.

| Requirement | Evidence | Result / boundary |
|---|---|---|
| Sequence-backed catalogue | `source/catalogue.tsv`, `source/graph_paths.jsonl.gz` | 1,182 intervals, exact graph/assembly concordance, 591 haplotypes at two loci |
| Additional HLA dimension | `results/extended_HLA_panel.tsv`, `results/same_hla_structural_contrasts.tsv` | Classical calls plus DRB/RCCX descriptors; two full-HLA-matched RCCX contrasts |
| Read-based technical evaluation | `results/predictions.tsv`, `results/split_summary.tsv` | 106 donors, query-family exclusion, fixed original graph, assembly comparators |
| Prospective model calibration | `source/dosage_calibration.json`, `results/calibrated_predictions.tsv` | One factor trained on 18 donors; frozen before 88-donor validation predictions |
| Simpler baselines | `results/depth_baseline.tsv`, `results/calibrated_depth_baseline.tsv`, `results/summed_dosage_baseline.tsv` | Calibrated ordinary depth ties 88/88 RCCX copy-number accuracy; graph-specific benefit not established |
| Arrangement limits | `FAILURE_ANALYSIS.md`, `results/arrangement_errors.tsv` | 73/88 full signatures, 74/88 C4-only arrangements; no validated unknown-structure rejection |
| Verification/provenance | `results/checks.json`, `source/input_manifest.json`, `results/artifact_manifest.json` | Tests pass; 699 source input hashes, scripts, archived measurements and logs |

DRB gene-content labels match 106/106 assemblies. RCCX copy number is recoverable
after calibration over the validation range of three to five total copies;
ordinary reference depth performs equally well. The tested path model remains
imperfect for C4 form assignment and phase. Its flat-sequence equivalent produces
the same predictions, so it does not establish a benefit from graph storage or
alignment itself. The report retains the initially unsuccessful uncalibrated
method and explicitly marks later fairness diagnostics.

This completes the defined prototype and benchmark, not clinical validation or
complete SV calling. Functional CYP21/TNX identity, exact breakpoints, WGS
recruitment, whole-genome specificity, robust novel-structure detection and disease
association are not established. Dedicated A/B and long/short junction markers,
independent long-read validation and a new untouched evaluation set are the next
assay-development stage. The 88 examined validation donors cannot serve as a new
blind test for further tuning.
