#!/bin/bash
#SBATCH --job-name=mhc_struct_provenance
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=/home/leechuck/hla/codex-structural/provenance-%j.log
set -euo pipefail
cd /home/leechuck/hla/codex-structural
/home/leechuck/hla/mm/envs/hla/bin/python3 freeze_inputs.py
