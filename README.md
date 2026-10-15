# Assembly simulator

This workflow performs two analysis steps:

1. CAMISIM 2.0.1 simulates paired-end reads for all samples using the custom
   NovaSeq 6000 ART profiles.
2. metaSPAdes assembles every sample and retains its complete output directory.

CAMISIM itself uses Nextflow. Snakemake invokes one CAMISIM run and then runs
the independent metaSPAdes assemblies.

## Input

`input/abundance.tsv` is a wide, tab-separated abundance table:

```text
genomes  seq_type  sample1  sample2
/absolute/genome1.fna  chromosome  1  0
/absolute/plasmid1.fna plasmid     1  0
/absolute/genome2.fna  chromosome  0  1
```

Values are non-negative relative genome-copy abundances. The preparation rule
converts this table into CAMISIM's genome-location, metadata, and per-sample
`distribution_<index>.txt` files. `work/camisim_input/sample_map.tsv` records
the relationship between CAMISIM's numerical indices and the original sample
names.

## Install CAMISIM

Clone the current release into the location configured in `config/config.yaml`:

```bash
mkdir -p tools
git clone --branch 2.0.1 --depth 1 \
  https://github.com/CAMI-challenge/CAMISIM.git tools/CAMISIM
```

CAMISIM 2.0.1 uses Nextflow and creates its own process environments. The
workflow pins Nextflow 23.10.1 because CAMISIM's configuration syntax is not
accepted by Nextflow 26.

## Create the workflow environment

The included `requirements.yaml` can be used to create a standalone outer
environment. On the current HPC, `run.slurm` instead activates the existing
`snakemake-8.30.0` environment, which is already known to work on compute
nodes. Snakemake creates the CAMISIM and metaSPAdes rule environments from
`envs/`.

To create the standalone environment elsewhere:

```bash
conda env create \
  --prefix /hpc/projects/2025/multibiomine/well_trodden_paths/conda_envs/assembly_simulator \
  --file requirements.yaml
```

This environment supplies Snakemake, Nextflow, Mamba, and Git.

## Run

From `assembly_simulator/`:

```bash
conda activate snakemake-8.30.0
mkdir -p logs benchmarks
snakemake --cores 80 --use-conda --resources mem_mb=1000000
```

On Slurm:

```bash
mkdir -p logs benchmarks
sbatch run.slurm
```

CAMISIM may use the full 80-core allocation during read simulation. Each
metaSPAdes assembly uses 20 threads and 120 GB, so an 80-core run executes at
most four assemblies at once.

## Outputs

- `results/camisim/`: the complete CAMISIM output, including reads, mappings,
  source genomes, distributions, and seeds retained by CAMISIM.
- `results/metaspades/<sample>/`: the complete metaSPAdes directory for each
  sample.
- `logs/` and `benchmarks/`: workflow and assembly logs and resource metrics.

The default simulation uses 0.96 Gbp per sample, 150 bp paired reads, a
350 +/- 35 bp fragment size, and random seed 14. These parameters can be
changed in `config/config.yaml`.
