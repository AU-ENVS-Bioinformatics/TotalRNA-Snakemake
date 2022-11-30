AVAILABLE_THREADS = int(config.get("SORTMERNA-THREADS", 50))
idx_SSU_directory = config.get("SORTMERNA_SSU_DATABASE_INDEX", "")
idx_LSU_directory = config.get("SORTMERNA_LSU_DATABASE_INDEX", "")

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
        idx_flag = f"--idx-dir {idx_SSU_directory}" if idx_SSU_directory else ""
    log:
        "logs/sortmerna/{sample}.log",
    benchmark:
        "benchmarks/sortmerna/{sample}.log"
    conda:
        "../envs/sortmerna.yaml"
    threads: AVAILABLE_THREADS
    shell:
        "sortmerna -ref {input.database_ref_ssu} "
        "--threads {threads} "
        "--workdir {params.outdir} "
        "{params.extra} "
        "{params.idx_flag} "
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
        protected(f"results/sortmerna/LSU/{{sample}}_fwd.fq.gz"),
        protected(f"results/sortmerna/LSU/{{sample}}_rev.fq.gz"),
        protected(f"results/sortmerna/not_LSU/{{sample}}_fwd.fq.gz"),
        protected(f"results/sortmerna/not_LSU/{{sample}}_rev.fq.gz"),
    params:
        extra=" ".join(config.get("sortmerna", "")),
        aligned=lambda wildcards, output: output[0][:-10],
        other=lambda wildcards, output: output[2][:-10],
        outdir=lambda wildcards, output: output[2][:-10].replace("not_LSU/", ""),
        idx_flag = f"--idx-dir {idx_LSU_directory}" if idx_LSU_directory else ""
    log:
        "logs/sortmerna_LSU/{sample}.log",
    benchmark:
        "benchmarks/sortmerna_LSU/{sample}.log"
    conda:
        "../envs/sortmerna.yaml"
    threads: AVAILABLE_THREADS
    shell:
        "sortmerna -ref {input.database_ref} "
        "--threads {threads} "
        "--workdir {params.outdir} "
        "{params.extra} "
        "{params.idx_flag} "
        "--aligned {params.aligned} "
        "--other {params.other} "
        "--reads {input.R1} --reads {input.R2} "
        ">> {log} 2>&1 "
