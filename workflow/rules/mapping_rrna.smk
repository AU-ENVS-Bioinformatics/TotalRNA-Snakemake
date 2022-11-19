RRNA_FILEPATH = config.get("RRNA_FILEPATH", "rrna/")
METARIB_FILEPATH = config.get("METARIB_FILEPATH", "MetaRib/")
OTU_FILEPATH = config.get("OTU_FILEPATH", "mapped_reads_to_contigs.tsv")
AVAILABLE_THREADS = int(workflow.cores * 0.5)


rule prepare_mapping_rrna:
    input:
        R1=expand(
            f"results/MetaRib/data/{{sample}}.1.fq",
            sample=unique_samples,
        ),
        R2=expand(
            f"results/MetaRib/data/{{sample}}.2.fq",
            sample=unique_samples,
        ),
    output:
        R1=expand(
            f"results/MetaRib/SSU_fastq/{{sample}}_R1.fastq",
            sample=unique_samples,
        ),
        R2=expand(
            f"results/MetaRib/SSU_fastq/{{sample}}_R2.fastq",
            sample=unique_samples,
        ),
        readsdir=directory(f"results/MetaRib/SSU_fastq/"),
    log:
        "logs/mapping_rrna/make_symlinks.log",
    conda:
        "../envs/base_python.yaml"
    script:
        "../scripts/make_symlinks.py"


rule map_reads_to_contigs:
    input:
        R1=expand(
            f"results/MetaRib/SSU_fastq/{{sample}}_R1.fastq",
            sample=unique_samples,
        ),
        R2=expand(
            f"results/MetaRib/SSU_fastq/{{sample}}_R2.fastq",
            sample=unique_samples,
        ),
        readsdir=f"results/MetaRib/SSU_fastq/",
        filtered=f"results/MetaRib/MetaRib/Abundance/all.dedup.filtered.fasta",
    output:
        f"results/MetaRib/SSU_fastq/mapped_reads_to_contigs.tsv",
    params:
        script="workflow/scripts/CoMW/scripts/map_reads_to_contigs.py",
    threads: AVAILABLE_THREADS
    conda:
        "../envs/CoMW.yaml"
    log:
        "logs/mapping_rrna/map_reads_to_contigs.log",
    shell:
        "python {params.script} "
        "-f {input.filtered} "
        "-i {input.readsdir} "
        "-o {output} "
        "-t {threads} "
        ">> {log} 2>&1"
