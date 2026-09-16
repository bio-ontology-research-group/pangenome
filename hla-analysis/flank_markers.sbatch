#!/bin/bash
#SBATCH --job-name=hla_flank_support
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=00:30:00
#SBATCH --output=/home/leechuck/hla/codex-analysis/flank-validation/job-%A_%a.log
set -euo pipefail
export PATH=/home/leechuck/hla/mm/envs/hla/bin:$PATH
work=/home/leechuck/hla/codex-analysis/flank-validation
decode_ref=/home/leechuck/hla/ref1kg/local/GRCh38_full_analysis_set_plus_decoy_hla.fa
cd "$work"
donor=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" read_marker_donors.txt)
test -n "$donor"
mkdir -p "$donor"
cram=/home/asianhla/data/upload/1000G_MHC/cram/${donor}.mhc.cram
samtools collate -@1 -u -O --reference "$decode_ref" "$cram" "$donor/collate" \
 | samtools fastq -@1 -F 0xF00 -n -1 "$donor/r1.fq.gz" -2 "$donor/r2.fq.gz" \
     -0 "$donor/other.fq.gz" -s "$donor/single.fq.gz" -
python3 count_read_markers.py read_markers.tsv "$donor" "$donor/support.tsv" \
 "$donor/r1.fq.gz" "$donor/r2.fq.gz" "$donor/other.fq.gz" "$donor/single.fq.gz"
date -u
echo COMPLETE
