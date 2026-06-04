def get_sortmerna_input(wildcards):
    """Get input based on rRNA type"""
    if wildcards.rrna_type == "SSU":
        return {
            "fasta": [
                "results/trim_galore/{sample}_fwd.fq.gz",
                "results/trim_galore/{sample}_rev.fq.gz",
            ]
        }
    elif wildcards.rrna_type == "LSU":
        return {
            "fasta": [
                "results/sortmerna/not_SSU/{sample}_fwd.fq.gz",
                "results/sortmerna/not_SSU/{sample}_rev.fq.gz",
            ]
        }
    else:
        raise ValueError(f"Unknown rRNA type: {wildcards.rrna_type}")


rule sortmerna:
    input:
        unpack(get_sortmerna_input),
        database=lambda wc: config.get(
            f"SORTMERNA_{wc.rrna_type}_REF_DATABASE",
            rules.databases_sortmerna_idx.output.fasta.format(RNA_TYPE=wc.rrna_type),
        ),
        database_index=lambda wc: config.get(
            f"SORTMERNA_{wc.rrna_type}_DATABASE_INDEX",
            rules.databases_sortmerna_idx.output.idx_dir.format(RNA_TYPE=wc.rrna_type),
        ),
    output:
        aligned=[
            "results/sortmerna/{rrna_type}/{sample}_fwd.fq.gz",
            "results/sortmerna/{rrna_type}/{sample}_rev.fq.gz",
        ],
        not_aligned=[
            "results/sortmerna/not_{rrna_type}/{sample}_fwd.fq.gz",
            "results/sortmerna/not_{rrna_type}/{sample}_rev.fq.gz",
        ],
        stats="results/sortmerna/{rrna_type}/{sample}.log",
    shadow:
        "minimal"
    params:
        extra=" ".join(config.get("sortmerna", [])),
        aligned_prefix="results/sortmerna/{rrna_type}/{sample}",
        not_aligned_prefix="results/sortmerna/not_{rrna_type}/{sample}",
    log:
        "logs/sortmerna/{rrna_type}/{sample}.log",
    benchmark:
        "results/benchmarks/sortmerna_{rrna_type}_{sample}.txt",
    conda:
        "../envs/sortmerna.yaml"
    threads: config["threads"]["sortmerna"]
    shell:
        """
        sortmerna -ref {input.database} \
        --idx-dir {input.database_index} \
        --workdir . \
        --threads {threads} \
        {params.extra} \
        --aligned {params.aligned_prefix} \
        --other {params.not_aligned_prefix} \
        --reads {input.fasta[0]} --reads {input.fasta[1]} \
        > {log} 2>&1
        """