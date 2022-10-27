DEFAULT_DEST_FILEPATH = config.get("DEFAULT_DEST_FILEPATH", "results/")
METARIB_FILEPATH = config.get("METARIB_FILEPATH", "MetaRib/")
QUAST_FILEPATH = config.get("QUAST_FILEPATH", "quast/")
AVAILABLE_THREADS = int(workflow.cores * 0.75)


rule quast:
    input:
        fasta=f"{DEFAULT_DEST_FILEPATH}{METARIB_FILEPATH}MetaRib/Abundance/all.dedup.filtered.fasta",
        R1 = f"{DEFAULT_DEST_FILEPATH}{METARIB_FILEPATH}data/all.1.fq",   
        R2 = f"{DEFAULT_DEST_FILEPATH}{METARIB_FILEPATH}data/all.2.fq",
    output:
        outdir=directory(f"{DEFAULT_DEST_FILEPATH}{QUAST_FILEPATH}"),
        report_txt=report(
            f"{DEFAULT_DEST_FILEPATH}{QUAST_FILEPATH}report.txt",
            category="Genome assembly",
        ),
        report_tsv=report(
            f"{DEFAULT_DEST_FILEPATH}{QUAST_FILEPATH}report.tsv",
            category="Genome assembly",
        ),
        report_html=report(
            f"{DEFAULT_DEST_FILEPATH}{QUAST_FILEPATH}report.html",
            caption="report/quast_html.rst",
            category="Genome assembly",
        ),
    log:
        "logs/quast.log",
    conda:
        "../envs/quast.yaml"
    params:
        extra=" ".join(config.get("quast", "")),
    threads: AVAILABLE_THREADS
    shell:
        "quast {params.extra} "
        "--threads {threads} "
        "-1 {input.R1} -2 {input.R2} "
        "-o {output.outdir} "
        "{input.fasta} "
        ">> {log} 2>&1 "