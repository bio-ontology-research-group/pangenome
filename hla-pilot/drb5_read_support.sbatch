#!/bin/bash
#SBATCH --job-name=codex_drb5_support
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --output=/home/leechuck/hla/codex-pilot/drb5/read-support-%j.log
set -euo pipefail
export PATH=/home/leechuck/hla/mm/envs/hla/bin:$PATH
export OPENBLAS_NUM_THREADS=1
work=/home/leechuck/hla/codex-pilot/drb5
decode_ref=/home/leechuck/hla/ref1kg/local/GRCh38_full_analysis_set_plus_decoy_hla.fa
cd "$work"
date -u
hostname
minimap2 --version
samtools --version | head -1
# Competitive diploid MHC reference plus a DRB5-containing control segment.
# This is targeted read support, not a genome-wide CNV assay.
cat /home/leechuck/hla/mhc_all/HG03139_1.mhc.fa \
    /home/leechuck/hla/mhc_all/HG03139_2.mhc.fa DRB5_control.fa > competitive.fa
samtools faidx competitive.fa
for sample in HG03139 HG03195; do
    mkdir -p "$sample"
    cram=/home/asianhla/data/upload/1000G_MHC/cram/${sample}.mhc.cram
    samtools collate -@1 -u -O --reference "$decode_ref" "$cram" "$sample/collate" \
      | samtools fastq -@1 -F 0x900 -n -1 "$sample/r1.fq.gz" -2 "$sample/r2.fq.gz" \
          -0 "$sample/other.fq.gz" -s "$sample/single.fq.gz" -
    minimap2 -t2 -ax sr competitive.fa "$sample/r1.fq.gz" "$sample/r2.fq.gz" \
      | samtools sort -@1 -m 1G -o "$sample/paired.bam" -
    # Keep and align singleton/other primary reads as well as pairs.
    cat "$sample/single.fq.gz" "$sample/other.fq.gz" > "$sample/unpaired.fq.gz"
    minimap2 -t2 -ax sr competitive.fa "$sample/unpaired.fq.gz" \
      | samtools sort -@1 -m 1G -o "$sample/unpaired.bam" -
    samtools merge -@1 "$sample/competitive.bam" "$sample/paired.bam" "$sample/unpaired.bam"
    samtools index "$sample/competitive.bam"
    samtools flagstat "$sample/competitive.bam" > "$sample/flagstat.txt"
    samtools coverage -q 20 -Q 20 "$sample/competitive.bam" > "$sample/coverage.tsv"
    samtools depth -aa -q 20 -Q 20 -r 'HG03195#2#HLA-DRB5' "$sample/competitive.bam" > "$sample/DRB5.depth.tsv"
    samtools view -c -F 0x904 -q 20 "$sample/competitive.bam" 'HG03195#2#HLA-DRB5' > "$sample/DRB5.primary_mapq20_count.txt"
done
date -u
echo COMPLETE
