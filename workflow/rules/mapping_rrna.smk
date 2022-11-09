RRNA_FILEPATH = config.get("RRNA_FILEPATH", "rrna/")
METARIB_FILEPATH = config.get("METARIB_FILEPATH", "MetaRib/")
OTU_FILEPATH = config.get("OTU_FILEPATH", "mapped_reads_to_contigs.tsv")
AVAILABLE_THREADS = int(workflow.cores * 0.75)


rule prepare_mapping_rrna:
    input: 
        R1=expand(
            f"{DEFAULT_DEST_FILEPATH}{METARIB_FILEPATH}data/{{sample}}.1.fq",
            sample=unique_samples,
        ),
        R2=expand(
            f"{DEFAULT_DEST_FILEPATH}{METARIB_FILEPATH}data/{{sample}}.2.fq",
            sample=unique_samples,
        )
    output: 
        R1=expand(
            f"{DEFAULT_DEST_FILEPATH}{METARIB_FILEPATH}SSU_fastq/{{sample}}_R1.fastq",
            sample=unique_samples,
        ),
        R2=expand(
            f"{DEFAULT_DEST_FILEPATH}{METARIB_FILEPATH}SSU_fastq/{{sample}}_R2.fastq",
            sample=unique_samples,
        ),
        readsdir = directory(f"{DEFAULT_DEST_FILEPATH}{METARIB_FILEPATH}SSU_fastq/")
    log: "logs/mapping_rrna/make_symlinks.log"
    conda:
        "../envs/base_python.yaml"
    script:
        "../scripts/make_symlinks.py"

rule map_reads_to_contigs:
    input: 
        R1=expand(
            f"{DEFAULT_DEST_FILEPATH}{METARIB_FILEPATH}SSU_fastq/{{sample}}_R1.fastq",
            sample=unique_samples,
        ),
        R2=expand(
            f"{DEFAULT_DEST_FILEPATH}{METARIB_FILEPATH}SSU_fastq/{{sample}}_R2.fastq",
            sample=unique_samples,
        ),
        readsdir = f"{DEFAULT_DEST_FILEPATH}{METARIB_FILEPATH}SSU_fastq/",
        filtered = f"{DEFAULT_DEST_FILEPATH}{METARIB_FILEPATH}MetaRib/Abundance/all.dedup.filtered.fasta"
    output: 
        f"{DEFAULT_DEST_FILEPATH}{METARIB_FILEPATH}SSU_fastq/{OTU_FILEPATH}"
    params:
        script = "workflow/scripts/CoMW/scripts/map_reads_to_contigs.py",
    threads: AVAILABLE_THREADS,
    conda:
        "../envs/CoMW.yaml"
    log: "logs/mapping_rrna/map_reads_to_contigs.log"
    shell: 
        "python {params.script} "
        "-f {input.filtered} "
        "-i {input.readsdir} "
        "-o {output} "
        "-t {threads} "
        ">> {log} 2>&1"