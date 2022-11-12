TRIMMED_READS_FILEPATH = config.get("TRIMMED_READS_FILEPATH", "trimmed/")
SORTMERNA_FILEPATH = config.get("SORTMERNA_FILEPATH", "sortmerna/")
RRNA_FILEPATH = config.get("RRNA_FILEPATH", "rrna/")
mRNA_FILEPATH = config.get("RRNA_FILEPATH", "mRNA/")
AVAILABLE_THREADS = int(workflow.cores * 0.5)


rule sortmerna_ssu:
    input:
        R1=f"results/trimmed/{{sample}}_R1_val_1.fq.gz",
        R2=f"results/trimmed/{{sample}}_R2_val_2.fq.gz",
        database_ref_ssu=config.get("SORTMERNA_SSU_REF_DATABASE"),
    shadow:
        "minimal"
    output:
        protected(f"results/rrna/{{sample}}_fwd.fq.gz"),
        protected(f"results/rrna/{{sample}}_rev.fq.gz"),
        protected(f"results/sortmerna/not_SSU/{{sample}}_fwd.fq.gz"),
        protected(f"results/sortmerna/not_SSU/{{sample}}_rev.fq.gz"),
    params:
        extra=" ".join(config.get("sortmerna", "")),
        aligned=lambda wildcards, output: output[0][:-10],
        other=lambda wildcards, output: output[2][:-10],
        outdir=lambda wildcards, output: output[2][:-10].replace("not_SSU/", ""),
    log:
        "logs/sortmerna/{sample}.log",
    conda:
        "../envs/sortmerna.yaml"
    threads: AVAILABLE_THREADS
    shell:
        "sortmerna -ref {input.database_ref_ssu} "
        "--threads {threads} "
        "--workdir {params.outdir} "
        "{params.extra} "
        "--aligned {params.aligned} "
        "--other {params.other} "
        "--reads {input.R1} --reads {input.R2} "
        ">> {log} 2>&1 "

rule sortmerna_LSU:
    input:
        R1=f"results/rrna/{{sample}}_fwd.fq.gz",
        R2=f"results/rrna/{{sample}}_rev.fq.gz",
        database_ref=config.get("SORTMERNA_LSU_REF_DATABASE"),
    shadow:
        "minimal"
    output:
        protected(f"results/mRNA/{{sample}}_fwd.fq.gz"),
        protected(f"results/mRNA/{{sample}}_rev.fq.gz"),
        protected(f"results/sortmerna/not_LSU/{{sample}}_fwd.fq.gz"),
        protected(f"results/sortmerna/not_LSU/{{sample}}_rev.fq.gz"),
    params:
        extra=" ".join(config.get("sortmerna", "")),
        aligned=lambda wildcards, output: output[0][:-10],
        other=lambda wildcards, output: output[2][:-10],
        outdir=lambda wildcards, output: output[2][:-10].replace("not_LSU/", ""),
    log:
        "logs/sortmerna_LSU/{sample}.log",
    conda:
        "../envs/sortmerna.yaml"
    threads: AVAILABLE_THREADS
    shell:
        "sortmerna -ref {input.database_ref} "
        "--threads {threads} "
        "--workdir {params.outdir} "
        "{params.extra} "
        "--aligned {params.aligned} "
        "--other {params.other} "
        "--reads {input.R1} --reads {input.R2} "
        ">> {log} 2>&1 "