#!/bin/bash
#SBATCH --job-name=ctxm14_sim_assembly
#SBATCH --output=logs/ctxm14_sim_assembly.%j.out
#SBATCH --error=logs/ctxm14_sim_assembly.%j.err
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=80
#SBATCH --mem=1500G

set -euo pipefail

source /hpc/opt/conda/etc/profile.d/conda.sh
conda activate /hpc/projects/2025/multibiomine/well_trodden_paths/conda_envs/magician_base

snakemake \
  --snakefile Snakefile \
  --cores "${SLURM_CPUS_PER_TASK}" \
  --use-conda \
  --rerun-incomplete \
  --printshellcmds \
  --resources mem_mb=1000000
