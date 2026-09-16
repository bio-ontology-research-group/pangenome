#!/bin/bash
#SBATCH --job-name=rccx_read_depth
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=/home/leechuck/hla/codex-analysis/rccx/depth-%j.log
set -euo pipefail
export PATH=/home/leechuck/hla/mm/envs/hla/bin:$PATH
cd /home/leechuck/hla/codex-analysis/rccx
ref=/home/leechuck/hla/ref1kg/local/GRCh38_full_analysis_set_plus_decoy_hla.fa
for donor in NA21144 HG03195; do
    samtools view -T "$ref" -u -F 0xF04 \
        /home/asianhla/data/upload/1000G_MHC/cram/${donor}.mhc.cram chr6:31960000-32110000 \
        | samtools depth -a -q20 -Q0 -b rccx_depth_targets.bed - > "${donor}.depth.tsv"
done
echo COMPLETE
