#!/bin/bash
#SBATCH --job-name=mhc_struct_background
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/home/leechuck/hla/codex-structural/background-%j.log
set -euo pipefail
cd /home/leechuck/hla/codex-structural
mkdir -p background_counts
while read -r entry; do
 ./count_markers fasta markers.txt "/home/leechuck/hla/mhc_all/${entry}.mhc.fa" "background_counts/${entry}.bin"
done < assembly_entries.txt
echo COMPLETE
