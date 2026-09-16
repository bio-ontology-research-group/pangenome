#!/bin/bash
#SBATCH --job-name=mhc_struct_counts
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/home/leechuck/hla/codex-structural/count-%A_%a.log
set -euo pipefail
export PATH=/home/leechuck/hla/mm/envs/hla/bin:$PATH
cd /home/leechuck/hla/codex-structural
mkdir -p read_counts
donor=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "${DONOR_LIST:-donors.txt}")
test -n "$donor"
ref=/home/leechuck/hla/ref1kg/local/GRCh38_full_analysis_set_plus_decoy_hla.fa
samtools collate -@1 -u -O --reference "$ref" "/home/asianhla/data/upload/1000G_MHC/cram/${donor}.mhc.cram" "read_counts/${donor}.tmp" \
 | samtools fastq -n -F 0xF04 - \
 | ./count_markers fastq markers.txt - "read_counts/${donor}.bin"
samtools view -T "$ref" -u -F 0xF04 "/home/asianhla/data/upload/1000G_MHC/cram/${donor}.mhc.cram" chr6:31960000-32110000 \
 | samtools depth -a -q20 -Q0 -b rccx_depth_targets.bed - \
 | gzip > "read_counts/${donor}.depth.tsv.gz"
echo COMPLETE "$donor"
