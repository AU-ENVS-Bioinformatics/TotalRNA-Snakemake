def get_sortmerna_input(wildcards):
    """Get input based on rRNA type"""
    if wildcards.rrna_type == "SSU":
        return {
            "fasta": [
                "results/trim_galore/{sample}_R1.fq.gz",
                "results/trim_galore/{sample}_R2.fq.gz",
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
    shadow:
        "minimal"
    output:
        aligned=[
            "results/sortmerna/{rrna_type}/{sample}_fwd.fq.gz",
            "results/sortmerna/{rrna_type}/{sample}_rev.fq.gz",
        ],
        not_aligned=[
            "results/sortmerna/not_{rrna_type}/{sample}_fwd.fq.gz",
            "results/sortmerna/not_{rrna_type}/{sample}_rev.fq.gz",
        ],
    params:
        extra=" ".join(config.get("sortmerna", [])),
    log:
        "logs/sortmerna/{rrna_type}/{sample}.log",
    conda:
        "../envs/sortmerna.yaml"
    threads: config["threads"]["sortmerna"]
    shell:
        """
        sortmerna -ref {input.database} \
        --idx-dir {input.database_index} \
        --workdir . \
        --threads {threads} \
        {params.extra} --log \
        --aligned aligned \
        --other not_aligned \
        --reads {input.fasta[0]} --reads {input.fasta[1]} \
        > {log} 2>&1
        
        mkdir -p $(dirname {output.aligned[0]})
        mkdir -p $(dirname {output.not_aligned[0]})

        cp aligned_fwd.fq.gz {output.aligned[0]} > {log} 2>&1
        cp aligned_rev.fq.gz {output.aligned[1]} > {log} 2>&1
        cp not_aligned_fwd.fq.gz {output.not_aligned[0]} > {log} 2>&1
        cp not_aligned_rev.fq.gz {output.not_aligned[1]} > {log} 2>&1
        """
