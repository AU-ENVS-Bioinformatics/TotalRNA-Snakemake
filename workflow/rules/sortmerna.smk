idx_SSU_directory = config.get("SORTMERNA_SSU_DATABASE_INDEX")
idx_LSU_directory = config.get("SORTMERNA_LSU_DATABASE_INDEX")


rule sortmerna_ssu:
    input:
        fasta=["results/trimmed/{sample}_R1.fq.gz", "results/trimmed/{sample}_R2.fq.gz"],
        database=ancient(config.get("SORTMERNA_SSU_REF_DATABASE")),
        database_index=idx_SSU_directory,
    shadow:
        "minimal"
    output:
        aligned=[
            protected("results/rrna/{sample}_fwd.fq.gz"),
            protected("results/rrna/{sample}_rev.fq.gz"),
        ],
        not_aligned=[
            protected("results/sortmerna/not_SSU/{sample}_fwd.fq.gz"),
            protected("results/sortmerna/not_SSU/{sample}_rev.fq.gz"),
        ],
        stats="results/rrna/{sample}.log",
    params:
        extra=" ".join(config.get("sortmerna", "")),
    log:
        "logs/sortmerna/{sample}_SSU.log",
    conda:
        "../envs/sortmerna.yaml"
    threads: config["threads"]["sortmerna"]
    script:
        "../scripts/sortmerna.py"


rule sortmerna_lsu:
    input:
        fasta=[
            "results/sortmerna/not_SSU/{sample}_fwd.fq.gz",
            "results/sortmerna/not_SSU/{sample}_rev.fq.gz",
        ],
        database=ancient(config.get("SORTMERNA_LSU_REF_DATABASE")),
        database_index=idx_LSU_directory,
    shadow:
        "minimal"
    output:
        aligned=[
            protected("results/sortmerna/LSU/{sample}_fwd.fq.gz"),
            protected("results/sortmerna/LSU/{sample}_rev.fq.gz"),
        ],
        not_aligned=[
            protected("results/sortmerna/not_LSU/{sample}_fwd.fq.gz"),
            protected("results/sortmerna/not_LSU/{sample}_rev.fq.gz"),
        ],
        stats="results/sortmerna/LSU/{sample}.log",
    params:
        extra=" ".join(config.get("sortmerna", "")),
    log:
        "logs/sortmerna/{sample}_LSU.log",
    conda:
        "../envs/sortmerna.yaml"
    threads: config["threads"]["sortmerna"]
    script:
        "../scripts/sortmerna.py"
