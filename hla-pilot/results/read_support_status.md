# DRB5 targeted read-support result

**Completed:** DDBJ Slurm job **20604529**, partition `asianhla-c32`, node `asianhla-vm`, 4 CPUs, 12 GB requested. The job log ends in `COMPLETE`. Initial job 20604506 failed because a symlinked CRAM reference was unavailable on the compute node; the successful job used the existing full local reference under `hla/ref1kg/local/`. No cluster jobs from this pilot remain running.

## Test

Map the existing MHC-recruited 1000G reads for HG03139 (candidate) and HG03195 (positive control) against the same competitive reference: both assembled HG03139 MHC haplotypes plus the HG03195 haplotype-2 DRB5 gene and 2 kb flanks. Both paired and unpaired primary reads were retained and aligned. Aligner: minimap2 2.28-r1209, `-ax sr`; samtools 1.21.

The control segment is 16,882 bp long. The annotated DRB5 gene occupies positions 2001–14882 inclusive (12,882 bp) in this reverse-complemented segment. Depth thresholds are mapping quality ≥20 and base quality ≥20.

| Sample | Mean gene depth | Median gene depth | Gene bases at ≥5 reads | Gene bases with zero depth |
|---|---:|---:|---:|---:|
| HG03139 | 0.00 | 0 | 0% | 100% |
| HG03195 | 20.10 | 20 | 99.86% | 0% |

One high-mapping-quality HG03139 alignment reaches the **flanking** part of the control segment; none covers the annotated DRB5 gene at these thresholds. The positive control has 2,282 reads in the complete control segment according to `samtools coverage`. Both candidate MHC haplotypes otherwise receive substantial read alignment.

This supports **DRB5 absence in HG03139** together with the assembly comparison. It does not establish an exact deletion breakpoint, validate chromosome-wide phase, prove a dosage/expression effect, or demonstrate a new structural class. DRB1*15:03 without DRB5 is already described in the literature.

## Limits

- Input CRAMs were recruited through GRCh38/MHC alternate/HLA contigs, so recruitment bias remains possible; absence of alignments alone is not definitive evidence of absence.
- The extra DRB5 reference uses one control allele. An exceptionally divergent alternative could map poorly. A broader allele panel or donor-matched long reads would strengthen the inference.
- This competitive reference covers the MHC rather than the whole genome, and does not constitute a general genome-wide copy-number caller.
- The control is not a fully matched whole-genome background. It tests that the assay recovers a known DRB5 sequence from the same kind of read data.

## Artifacts

Local summaries, depth files, logs and figures: `read_support/` beside this file. Remote BAMs, indexes, FASTQs and competitive reference: `/home/leechuck/hla/codex-pilot/drb5/` on DDBJ. The job script is `hla-pilot/drb5_read_support.sbatch` in the local repository and is also saved directly in the remote `drb5/` directory.

The candidate and control assembly block alignments are saved locally as `DR_block_alignment.paf` and `HG03139_DRB5_alignment.paf` in the parent results folder. They do not by themselves define precise breakpoints.
