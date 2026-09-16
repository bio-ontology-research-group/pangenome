#!/bin/bash
#SBATCH --job-name=mhc_struct_paths
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --output=/home/leechuck/hla/codex-structural/extract-%j.log
set -euo pipefail
cd /home/leechuck/hla/codex-structural
/home/leechuck/hla/mm/envs/hla/bin/python3 extract_paths.py
