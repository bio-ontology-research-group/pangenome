#!/bin/bash
#SBATCH --job-name=hla_rccx_screen
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/home/leechuck/hla/codex-analysis/rccx/map-%j.log
set -euo pipefail
export PATH=/home/leechuck/hla/mm/envs/hla/bin:$PATH
cd /home/leechuck/hla/codex-analysis/rccx
mkdir -p paf
minimap2 --version
while read -r entry; do
    minimap2 -t2 -x asm20 -c --secondary=yes -N50 -p0.1 \
        "/home/leechuck/hla/mhc_all/${entry}.mhc.fa" rccx_prototypes.fa \
        > "paf/${entry}.paf" 2> "paf/${entry}.log"
    echo "$entry"
done < rccx_panel_ids.txt
echo COMPLETE
