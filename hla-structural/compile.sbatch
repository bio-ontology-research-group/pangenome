#!/bin/bash
#SBATCH --job-name=mhc_struct_compile
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=/home/leechuck/hla/codex-structural/compile-%j.log
set -euo pipefail
cd /home/leechuck/hla/codex-structural
g++ -O3 -std=c++17 count_markers.cpp -o count_markers
g++ --version
