#!/bin/bash
#SBATCH --job-name=ctxm14_sim_assembly
#SBATCH --output=logs/ctxm14_sim_assembly.%j.out
#SBATCH --error=logs/ctxm14_sim_assembly.%j.err
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=1000G

set -euo pipefail

echo "Started: $(date --iso-8601=seconds)"
echo "Host: $(hostname)"
echo "Working directory: $(pwd)"

echo "Loading Conda shell support"
source /hpc/opt/conda/etc/profile.d/conda.sh

environment_path="/hpc/users/marvoj/.conda/envs/assembly_simulator"
export PATH="${environment_path}/bin:${PATH}"

echo "Using environment: ${environment_path}"
echo "Snakemake: $(command -v snakemake)"
echo "Nextflow: $(command -v nextflow)"

mkdir -p logs benchmarks

snakemake \
  --snakefile Snakefile \
  --cores 64 \
  --use-conda \
  --rerun-incomplete \
  --printshellcmds \
  --resources mem_mb=1000000
