rule quast:
    input:
        fasta="results/MetaRib/all.dedup.filtered.fasta",
        pe1="results/MetaRib/data/all.1.fq",
        pe2="results/MetaRib/data/all.2.fq",
    output:
        multiext(
            "qc/quast/all.dedup.filtered.fasta/report.",
            "html",
            "tex",
            "txt",
            "pdf",
            "tsv",
        ),
        multiext(
            "qc/quast/all.dedup.filtered.fasta/transposed_report.", "tex", "txt", "tsv"
        ),
        multiext(
            "qc/quast/all.dedup.filtered.fasta/basic_stats/",
            "cumulative_plot.pdf",
            "GC_content_plot.pdf",
            "gc.icarus.txt",
            "genome_GC_content_plot.pdf",
            "NGx_plot.pdf",
            "Nx_plot.pdf",
        ),
        multiext(
            "qc/quast/all.dedup.filtered.fasta/contigs_reports/",
            "all_alignments_genome.tsv",
            "contigs_report_genome.mis_contigs.info",
            "contigs_report_genome.stderr",
            "contigs_report_genome.stdout",
        ),
        "qc/quast/all.dedup.filtered.fasta/contigs_reports/minimap_output/genome.coords_tmp",
        "qc/quast/all.dedup.filtered.fasta/icarus.html",
        "qc/quast/all.dedup.filtered.fasta/icarus_viewers/contig_size_viewer.html",
        "qc/quast/all.dedup.filtered.fasta/quast.log",
    log:
        "logs/rrna/quast.log",
    threads: config["threads"]["quast"]
    params:
        extra=" ".join(config.get("quast", "")),
    wrapper:
        "v2.6.0/bio/quast"
