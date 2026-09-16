# Structural MHC typing

Follow-up to the completed HLA analysis. Read [REPORT.md](REPORT.md) for results
and [PROTOCOL.md](PROTOCOL.md) for the prospective experiment and limits.

The custom prototype asks whether observed local haplotype paths plus marker
multiplicity recover assembly-labelled copy number and structural signatures from
MHC-recruited short reads. It is not PanGenie or a replacement for a validated SV
caller. Graph walks and flat sequence catalogues are equivalent inputs to this
particular inference method.

## Reproduce from archived measurements

Dependencies: Python, NumPy, Matplotlib; g++ with C++17 for read-counting tests.

```sh
g++ -O3 -std=c++17 hla-structural/count_markers.cpp -o hla-structural/count_markers
python3 hla-structural/test_contracts.py
OPENBLAS_NUM_THREADS=2 python3 hla-structural/benchmark.py
OPENBLAS_NUM_THREADS=2 python3 hla-structural/calibrated_benchmark.py
python3 hla-structural/failure_analysis.py
python3 hla-structural/report.py
python3 hla-structural/plot_results.py
python3 hla-structural/check_results.py
```

`prepare_catalogue.py` regenerates interval requests from the previous frozen QC
panel and annotation files. `extract_paths.py`, run on DDBJ by `extract.sbatch`,
extracts actual GFA walks and verifies their sequence against the source assemblies.
`prepare_markers.py` regenerates the deterministic 31-mer vocabulary, sequence
profiles and donor lists. Marker vocabulary is a measurement union: each evaluation
fold separately excludes query families and requires independent training support.

`compile.sbatch` builds the exact C++ counter on a compute node. `background.sbatch`
counts markers across whole source MHC assemblies to screen off-locus hits.
`count.sbatch` name-collates recruited CRAM reads, filters records and Q20 marker
windows, deduplicates mate support, counts markers, and records a reference-depth
baseline. Its default donor list is the 18 pilot donors; set
`DONOR_LIST=validation_donors.txt` for the remaining 88. Slurm arrays use zero-based
indices. Decode reference, source directories and job IDs are in
`source/remote_jobs.json`. SSH uses alias `ddbj`, then the existing DDBJ key for
`a001`; no credentials are copied into the analysis.

Primary outputs:

- `source/catalogue.tsv`: 1,182 interval definitions, sequence hashes and QC.
- `source/graph_paths.jsonl.gz`: exact oriented node traversals, with clipped end
  offsets and repeated nodes retained.
- `results/structural_types.tsv`: structural signature catalogue.
- `results/extended_HLA_panel.tsv`: classical HLA plus structural descriptors for
  each eligible assembled haplotype, with explicit phase and evidence limits.
- `results/same_hla_structural_contrasts.tsv`: extra structure within matching HLA.
- `results/predictions.tsv`: all attempted donor/locus/method predictions.
- `results/dosage_baseline.tsv`, `results/depth_baseline.tsv`: simpler comparisons.
- `results/split_summary.tsv`, `results/uncertainty.tsv`: performance and uncertainty.
- `results/calibrated_predictions.tsv`: the pilot-trained secondary model, scored
  only on the other 88 donors. `source/dosage_calibration.json` preserves the
  factor, pilot records and code hash frozen before validation predictions.
- `FAILURE_ANALYSIS.md`: remaining errors and the next assay-development target.
- `results/benchmark.png` and `.pdf`: standalone validation-set figure.

The observed assembly paths are the comparator, not independent experimental
structural truth. Full graph construction included the queries; genotyping
references and marker selection exclude their families. Reads are pre-recruited
using linear-reference mappings. Full WGS recruitment, whole-genome specificity,
breakpoint confirmation, functional paralog identity and phenotype association
remain separate validation tasks.

The original pre-prediction model snapshot is retained as
`source/benchmark_frozen_before_predictions.py.txt`; its hash matches the recorded
freeze. The current `benchmark.py` adds only donor-list and output-directory CLI
adapters to that model. The scalar-calibrated reference-depth and joint-form
dosage comparisons are post-pilot fairness diagnostics, not predeclared primary
algorithms. Read results and all selected donors are retained, including errors.
