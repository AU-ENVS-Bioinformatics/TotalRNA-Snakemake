rule quast:
    input:
        fasta="results/MetaRib/all.dedup.filtered.fasta",
        R1="results/MetaRib/data/all.1.fq",
        R2="results/MetaRib/data/all.2.fq",
    output:
        outdir=directory("qc/quast/MetaRib"),
        report_txt="qc/quast/MetaRib/report.txt",
        report_tsv="qc/quast/MetaRib/report.tsv",
        report_html="qc/quast/MetaRib/report.html",
    log:
        "logs/quast-MetaRib.log",
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
