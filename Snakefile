import csv
import pathlib


configfile: "config/config.yaml"


workflow_dir = pathlib.Path(workflow.current_basedir)


def resolve_path(path):
    path = pathlib.Path(path)
    if path.is_absolute():
        return str(path)
    return str((workflow_dir / path).resolve())


abundance_table = resolve_path(config["abundance_table"])
camisim_dir = resolve_path(config["camisim_dir"])
error_profile = resolve_path(config["error_profile"])


with open(abundance_table, newline="") as abundance_file:
    header = next(csv.reader(abundance_file, delimiter="\t"))

if header[:2] != ["genomes", "seq_type"]:
    raise ValueError("The abundance table must begin with: genomes, seq_type")

samples = header[2:]
sample_index = {sample: index for index, sample in enumerate(samples)}
distribution_files = expand(
    "work/camisim_input/distribution_{index}.txt",
    index=range(len(samples)),
)
camisim_reads_1 = expand(
    "results/camisim/sample_{index}/reads/fastq/sample_{index}_01.fq.gz",
    index=range(len(samples)),
)
camisim_reads_2 = expand(
    "results/camisim/sample_{index}/reads/fastq/sample_{index}_02.fq.gz",
    index=range(len(samples)),
)


rule all:
    input:
        expand("results/metaspades/{sample}", sample=samples),
        "results/camisim/.complete",


rule prepare_camisim_input:
    input:
        abundance_table,
    output:
        genome_locations="work/camisim_input/genome_locations.tsv",
        metadata="work/camisim_input/metadata.tsv",
        sample_map="work/camisim_input/sample_map.tsv",
        distributions=distribution_files,
        config="work/camisim_input/camisim.config",
    params:
        camisim_dir=camisim_dir,
        error_profile=error_profile,
        output_dir=lambda wildcards: str((workflow_dir / "results" / "camisim").resolve()),
        sample_size_gb=config["sample_size_gb"],
        read_length=config["read_length"],
        fragment_mean=config["fragment_mean"],
        fragment_sd=config["fragment_sd"],
        seed=config["seed"],
    script:
        "scripts/prepare_camisim_input.py"


rule simulate_reads:
    input:
        camisim_config=rules.prepare_camisim_input.output.config,
        genome_locations=rules.prepare_camisim_input.output.genome_locations,
        metadata=rules.prepare_camisim_input.output.metadata,
        distributions=rules.prepare_camisim_input.output.distributions,
    output:
        reads_1=camisim_reads_1,
        reads_2=camisim_reads_2,
        done=touch("results/camisim/.complete"),
    log:
        "logs/camisim.log",
    threads: config["camisim_cores"]
    resources:
        mem_mb=config["camisim_mem_mb"],
    conda:
        "envs/camisim.yaml"
    params:
        camisim_dir=camisim_dir,
        work_dir=lambda wildcards: str((workflow_dir / "work" / "nextflow").resolve()),
    shell:
        """
        nextflow run "{params.camisim_dir}/main.nf" \
          --config "{input.camisim_config}" \
          -work-dir "{params.work_dir}" \
          -resume > "{log}" 2>&1
        """


rule assemble_metaspades:
    input:
        reads_1=lambda wildcards: (
            f"results/camisim/sample_{sample_index[wildcards.sample]}/reads/fastq/"
            f"sample_{sample_index[wildcards.sample]}_01.fq.gz"
        ),
        reads_2=lambda wildcards: (
            f"results/camisim/sample_{sample_index[wildcards.sample]}/reads/fastq/"
            f"sample_{sample_index[wildcards.sample]}_02.fq.gz"
        ),
    output:
        assembly=directory("results/metaspades/{sample}"),
    log:
        out="logs/metaspades/{sample}.out",
        err="logs/metaspades/{sample}.err",
    benchmark:
        "benchmarks/metaspades/{sample}.tsv"
    threads: config["metaspades_threads"]
    resources:
        mem_mb=config["metaspades_mem_mb"],
    conda:
        "envs/metaspades.yaml"
    params:
        memory_gb=lambda wildcards, resources: max(1, resources.mem_mb // 1000),
        kmers=",".join(str(kmer) for kmer in config["metaspades_kmers"]),
    shell:
        """
        metaspades.py \
          -1 "{input.reads_1}" \
          -2 "{input.reads_2}" \
          -o "{output.assembly}" \
          -t {threads} \
          --memory {params.memory_gb} \
          -k {params.kmers} \
          > "{log.out}" 2> "{log.err}"
        """
