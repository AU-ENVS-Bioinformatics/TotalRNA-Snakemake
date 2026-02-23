def get_sortmerna_input(wildcards):
    """Get input based on rRNA type"""
    if wildcards.rrna_type == "SSU":
        return {
            "fasta": [
                "results/trim_galore/{sample}_R1.fq.gz",
                "results/trim_galore/{sample}_R2.fq.gz",
            ]
        }
    else:  # LSU
        return {
            "fasta": [
                "results/sortmerna/not_SSU/{sample}_fwd.fq.gz",
                "results/sortmerna/not_SSU/{sample}_rev.fq.gz",
            ]
        }


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
    shadow:
        "minimal"
    output:
        aligned=[
            ("results/sortmerna/{rrna_type}/{sample}_fwd.fq.gz"),
            ("results/sortmerna/{rrna_type}/{sample}_rev.fq.gz"),
        ],
        not_aligned=[
            ("results/sortmerna/not_{rrna_type}/{sample}_fwd.fq.gz"),
            ("results/sortmerna/not_{rrna_type}/{sample}_rev.fq.gz"),
        ],
        stats="results/sortmerna/{rrna_type}/{sample}.log",
    params:
        extra=" ".join(config.get("sortmerna", [])),
    log:
        "logs/sortmerna/{rrna_type}/{sample}.log",
    conda:
        "../envs/sortmerna.yaml"
    threads: config["threads"]["sortmerna"]
    script:
        "../scripts/sortmerna.py"
