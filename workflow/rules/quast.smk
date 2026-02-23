rule quast:
    input:
        fasta="results/metarib/all.dedup.filtered.fasta",
        R1="results/metarib/data/all.1.fq",
        R2="results/metarib/data/all.2.fq",
    output:
        outdir=directory("qc/quast/metarib"),
        report_txt="qc/quast/metarib/report.txt",
        report_tsv="qc/quast/metarib/report.tsv",
        report_html="qc/quast/metarib/report.html",
    log:
        "logs/quast-metarib.log",
    conda:
        "../envs/quast.yaml"
    params:
        extra=" ".join(config.get("quast", "")),
    threads: config["threads"]["quast"]
    shell:
        "quast {params.extra} "
        "--threads {threads} "
        "-1 {input.R1} -2 {input.R2} "
        "-o {output.outdir} "
        "{input.fasta} "
        ">> {log} 2>&1 "
