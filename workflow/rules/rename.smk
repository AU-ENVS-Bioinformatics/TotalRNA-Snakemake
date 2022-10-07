rule rename:
    output:
        expand(
            f"{DEFAULT_DEST_FILEPATH}{RENAMED_READS_FILEPATH}{{sample}}_{{read}}.fastq.gz",
            zip,
            sample=sample,
            read=read,
        ),
    log:
        "logs/renaming_files/{wildcards.sample}_{wildcards.sample}.log",
    conda:
        "../envs/base_python.yaml"
    script:
        "../scripts/rename_files.py"