AVAILABLE_THREADS = int(config.get("QUAST-THREADS", 50))


rule quast:
    input:
        fasta=f"results/MetaRib/MetaRib/Abundance/all.dedup.filtered.fasta",
        R1=f"results/MetaRib/data/all.1.fq",
        R2=f"results/MetaRib/data/all.2.fq",
    output:
        outdir=directory(f"results/quast/"),
        report_txt=report(
            f"results/quast/report.txt",
            category="Genome assembly",
        ),
        report_tsv=report(
            f"results/quast/report.tsv",
            category="Genome assembly",
        ),
        report_html=report(
            f"results/quast/report.html",
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
