import csv
from pathlib import Path


def read_abundance_table(path):
    with open(path, newline="") as input_file:
        reader = csv.DictReader(input_file, delimiter="\t")
        if reader.fieldnames is None or reader.fieldnames[:2] != ["genomes", "seq_type"]:
            raise ValueError("The abundance table must begin with: genomes, seq_type")

        samples = reader.fieldnames[2:]
        rows = list(reader)

    if not samples or not rows:
        raise ValueError("The abundance table needs genomes and at least one sample")

    return rows, samples


def validate_rows(rows, samples):
    paths = set()

    for row_number, row in enumerate(rows, start=2):
        genome_path = row["genomes"]
        if not genome_path or genome_path in paths:
            raise ValueError(f"Missing or duplicate genome path on row {row_number}")
        paths.add(genome_path)

        for sample in samples:
            try:
                abundance = float(row[sample])
            except ValueError as error:
                raise ValueError(
                    f"Invalid abundance for {sample} on row {row_number}"
                ) from error
            if abundance < 0:
                raise ValueError(f"Negative abundance for {sample} on row {row_number}")

    for sample in samples:
        if sum(float(row[sample]) for row in rows) == 0:
            raise ValueError(f"Sample {sample} has no genomes with positive abundance")


def validate_inputs(rows):
    missing_genomes = [row["genomes"] for row in rows if not Path(row["genomes"]).is_file()]
    if missing_genomes:
        preview = "\n".join(missing_genomes[:5])
        raise FileNotFoundError(f"Genome files were not found:\n{preview}")

    profile_base = Path(snakemake.params.error_profile)
    # CAMISIM appends 1.txt and 2.txt to this profile prefix.
    missing_profiles = [
        f"{profile_base}{read_number}.txt"
        for read_number in (1, 2)
        if not Path(f"{profile_base}{read_number}.txt").is_file()
    ]
    if missing_profiles:
        raise FileNotFoundError(
            "ART error-profile files were not found:\n" + "\n".join(missing_profiles)
        )

    camisim_main = Path(snakemake.params.camisim_dir) / "main.nf"
    if not camisim_main.is_file():
        raise FileNotFoundError(f"CAMISIM entry point was not found: {camisim_main}")


def write_camisim_files(rows, samples):
    genome_ids = [f"genome_{index:04d}" for index in range(1, len(rows) + 1)]

    with open(snakemake.output.genome_locations, "w") as output_file:
        for genome_id, row in zip(genome_ids, rows):
            output_file.write(f"{genome_id}\t{row['genomes']}\n")

    with open(snakemake.output.metadata, "w") as output_file:
        output_file.write("genome_ID\tOTU\tNCBI_ID\tnovelty_category\n")
        for genome_id in genome_ids:
            output_file.write(f"{genome_id}\t{genome_id}\t2\tknown_strain\n")

    with open(snakemake.output.sample_map, "w") as output_file:
        output_file.write("sample_index\tsample\n")
        for sample_index, sample in enumerate(samples):
            output_file.write(f"{sample_index}\t{sample}\n")

    for sample_index, sample in enumerate(samples):
        distribution_path = Path(snakemake.output.distributions[sample_index])
        with open(distribution_path, "w") as output_file:
            for genome_id, row in zip(genome_ids, rows):
                output_file.write(f"{genome_id}\t{float(row[sample]):g}\n")


def write_config(genome_count, sample_count):
    config_text = f'''includeConfig "{snakemake.params.camisim_dir}/pipelines/metagenomic/config/metagenomic.config"

params.outdir = "{snakemake.params.output_dir}"
params.size = {snakemake.params.sample_size_gb}
params.type = "art"
params.number_of_samples = {sample_count}
params.gsa = false
params.pooled_gsa = []
params.anonymization = false
params.seed = {snakemake.params.seed}
params.genome_locations_file = "{Path(snakemake.output.genome_locations).resolve()}"
params.metadata_file = "{Path(snakemake.output.metadata).resolve()}"
params.distribution_files = "{Path(snakemake.output.distributions[0]).resolve().parent}/distribution_*.txt"
params.genomes_total = {genome_count}
params.genomes_real = {genome_count}
params.profile_read_length = {snakemake.params.read_length}
params.fragment_size_mean = {snakemake.params.fragment_mean}
params.fragment_size_sd = {snakemake.params.fragment_sd}
params.base_profile_name = "{snakemake.params.error_profile}"
'''

    Path(snakemake.output.config).write_text(config_text)


rows, samples = read_abundance_table(snakemake.input[0])
validate_rows(rows, samples)
validate_inputs(rows)
write_camisim_files(rows, samples)
write_config(len(rows), len(samples))
