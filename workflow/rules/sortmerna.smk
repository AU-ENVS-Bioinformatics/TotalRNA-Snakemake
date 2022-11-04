TRIMMED_READS_FILEPATH = config.get("TRIMMED_READS_FILEPATH", "trimmed/")
SORTMERNA_FILEPATH = config.get("SORTMERNA_FILEPATH", "sortmerna/")
RRNA_FILEPATH = config.get("RRNA_FILEPATH", "rrna/")
mRNA_FILEPATH = config.get("RRNA_FILEPATH", "mRNA/")
AVAILABLE_THREADS = int(workflow.cores * 0.75)


rule sortmerna:
    input:
        R1=f"{DEFAULT_DEST_FILEPATH}{TRIMMED_READS_FILEPATH}{{sample}}_R1_val_1.fq.gz",
        R2=f"{DEFAULT_DEST_FILEPATH}{TRIMMED_READS_FILEPATH}{{sample}}_R2_val_2.fq.gz",
        database_ref=config.get("SORTMERNA_REF_DATABASE"),
    shadow:
        "minimal"
    output:
        protected(f"{DEFAULT_DEST_FILEPATH}{RRNA_FILEPATH}{{sample}}_fwd.fq.gz"),
        protected(f"{DEFAULT_DEST_FILEPATH}{RRNA_FILEPATH}{{sample}}_rev.fq.gz"),
        protected(
            f"{DEFAULT_DEST_FILEPATH}{SORTMERNA_FILEPATH}not_SSU/{{sample}}_fwd.fq.gz"
        ),
        protected(
            f"{DEFAULT_DEST_FILEPATH}{SORTMERNA_FILEPATH}not_SSU/{{sample}}_rev.fq.gz"
        ),
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
        "sortmerna -ref {input.database_ref} "
        "--threads {threads} "
        "--workdir {params.outdir} "
        "{params.extra} "
        "--aligned {params.aligned} "
        "--other {params.other} "
        "--reads {input.R1} --reads {input.R2} "
        ">> {log} 2>&1 "

rule sortmerna_LSU:
    input:
        R1=f"{DEFAULT_DEST_FILEPATH}{SORTMERNA_FILEPATH}not_SSU/{{sample}}_fwd.fq.gz",
        R2=f"{DEFAULT_DEST_FILEPATH}{SORTMERNA_FILEPATH}not_SSU/{{sample}}_rev.fq.gz",
        database_ref=config.get("SORTMERNA_LSU_REF_DATABASE"),
    shadow:
        "minimal"
    output:
        protected(f"{DEFAULT_DEST_FILEPATH}mRNA/{{sample}}_fwd.fq.gz"),
        protected(f"{DEFAULT_DEST_FILEPATH}mRNA/{{sample}}_rev.fq.gz"),
        protected(
            f"{DEFAULT_DEST_FILEPATH}{SORTMERNA_FILEPATH}LSU/{{sample}}_fwd.fq.gz"
        ),
        protected(
            f"{DEFAULT_DEST_FILEPATH}{SORTMERNA_FILEPATH}LSU/{{sample}}_rev.fq.gz"
        ),
    params:
        extra=" ".join(config.get("sortmerna", "")),
        aligned=lambda wildcards, output: output[0][:-10],
        other=lambda wildcards, output: output[2][:-10],
        outdir=lambda wildcards, output: output[2][:-10].replace("not_SSU/", ""),
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