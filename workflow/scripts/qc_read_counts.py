import subprocess
import concurrent.futures
import pandas as pd


def count_fasta_reads(file):
    command = f"grep -c '^>' {file}"
    process = subprocess.Popen(
        command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    output, error = process.communicate()
    if process.returncode == 0:
        try:
            count = int(output.strip())
            return count
        except ValueError:
            print("Error: Unable to parse the output of grep -c")
    else:
        print(f"Error: {error.decode().strip()}")
    return 0


def count_fastq_reads(file):
    command = f"zcat {file} | wc -l"
    process = subprocess.Popen(
        command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    output, error = process.communicate()
    if process.returncode == 0:
        try:
            count = int(output.strip())
            return int(count / 4)
        except ValueError:
            print("Error: Unable to parse the output of wc -l")
    else:
        print(f"Error: {error.decode().strip()}")
    return 0


def parse_snakemake(snakemake):
    infiles, annotations, lambdas = list(), list(), list()

    def inner(files, ann, f):
        for file in files:
            infiles.append(file)
            annotations.append(ann)
            lambdas.append(f)

    n_samples = len(snakemake.input.raw_reads)
    assert n_samples % 2 == 0, "Expecting paired samples"
    middle = n_samples // 2

    inner(snakemake.input.raw_reads[0:middle], "Forward raw reads", count_fastq_reads)
    inner(
        snakemake.input.raw_reads[middle:n_samples],
        "Reverse raw reads",
        count_fastq_reads,
    )
    inner(
        snakemake.input.trimmed_reads[0:middle],
        "Forward trimmed reads",
        count_fastq_reads,
    )
    inner(
        snakemake.input.trimmed_reads[middle:n_samples],
        "Reverse trimmed reads",
        count_fastq_reads,
    )
    inner(
        snakemake.input.SSU_sortmerna[0:middle], "Forward SSU reads", count_fastq_reads
    )
    inner(
        snakemake.input.SSU_sortmerna[middle:n_samples],
        "Reverse SSU reads",
        count_fastq_reads,
    )
    inner(
        snakemake.input.not_LSU_sortmerna[0:middle],
        "Forward not rRNA reads",
        count_fastq_reads,
    )
    inner(
        snakemake.input.not_LSU_sortmerna[middle:n_samples],
        "Reverse not rRNA reads",
        count_fastq_reads,
    )

    inner([snakemake.input.rRNA], "rRNA assembly", count_fasta_reads)
    inner([snakemake.input.mRNA], "mRNA assembly", count_fasta_reads)
    inner([snakemake.input.filtered_mRNA], "Filtered mRNA assembly", count_fasta_reads)
    return infiles, annotations, lambdas


infiles, annotations, lambdas = parse_snakemake(snakemake)
outfile = str(snakemake.output)
with concurrent.futures.ThreadPoolExecutor(max_workers=snakemake.threads) as executor:
    # Create a list to store the future objects
    futures = []
    # Iterate over the forward and reverse reads
    for file, f in zip(infiles, lambdas):
        future = executor.submit(f, file)
        futures.append(future)

    counts = [future.result() for future in concurrent.futures.as_completed(futures)]
    df = pd.DataFrame({"Annotation": annotations, "Count": counts, "File": infiles})
    df.to_csv(outfile, index=False)
