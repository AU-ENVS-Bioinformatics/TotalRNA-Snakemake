RRNA_FILEPATH = config.get("RRNA_FILEPATH", "rrna/")

rule prepare_metarib:
    input:
        R1= expand(f"{DEFAULT_DEST_FILEPATH}{RRNA_FILEPATH}{{sample}}_fwd.fq.gz", sample = sample),
        R2= expand(f"{DEFAULT_DEST_FILEPATH}{RRNA_FILEPATH}{{sample}}_rev.fq.gz", sample = sample),
    output:
        R1 = f"{DEFAULT_DEST_FILEPATH}{RRNA_FILEPATH}all.1.fq.gz",   
        R2 = f"{DEFAULT_DEST_FILEPATH}{RRNA_FILEPATH}all.2.fq.gz",
        sample_list = report(
            f"{DEFAULT_DEST_FILEPATH}{RRNA_FILEPATH}samples_list.txt",
            caption = "report/prepare_metarib.rst"
            ),
    log:
        "logs/sortmerna/concatenate_all.log",
    conda:
        "../envs/base_python.yaml"
    threads: AVAILABLE_THREADS
    params:
        samples_names = "\n".join(sorted(set(sample))),
    shell:
        "cat {input.R1} > {output.R1} && "
        "echo 'Forward files were successfully concatenated' >> {log} && "
        "cat {input.R2} > {output.R2} && "
        "echo 'Reverse files were successfully concatenated' >> {log} && "
        "echo '{params.samples_names}' > {output.sample_list} && "
        "echo 'Copy samples names into samples_list.txt' >> {log}"
