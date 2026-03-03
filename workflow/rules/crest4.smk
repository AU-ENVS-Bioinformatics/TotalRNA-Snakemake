rule CREST4:
    input:
        fasta="results/metarib/final_contigs.fasta",
        otu="results/rRNA/mapped_reads_to_contigs.tsv",
    output:
        "results/crest4/assignments.txt",
        "results/crest4/otus_by_rank.tsv",
        "results/crest4/otus_cumulative.tsv",
        "results/crest4/search.hits",
    params:
        CREST4_DIR=config.get("CREST4_DIR", "~/.crest4/"),
        extra=" ".join(config.get("crest4", "")),
    conda:
        "../envs/crest4.yaml"
    log:
        "logs/crest4/crest4.log",
    threads: config["threads"]["crest4"]
    shell:
        "export CREST4_DIR={params.CREST4_DIR} ; "
        "crest4 -f {input.fasta} "
        "-t {threads} "
        "-u {input.otu} "
        "-o results/crest4 "
        "{params.extra} "
        "> {log} 2>&1 || true"


rule process_crest4:
    input:
        infile_assignments="results/crest4/assignments.txt",
        infile_counts="results/rRNA/mapped_reads_to_contigs.tsv",
    output:
        outfile_otu="results/crest4/mapped_reads_to_contigs.tsv",
        outfile_gg="results/crest4/mapped_reads_to_contigs_gg.tsv",
        outphyseq="results/crest4/physeq.Rds",
    conda:
        "../envs/phyloseq.yaml"
    log:
        "logs/crest4/process.log",
    script:
        "../scripts/process_crest4_phyloseq.R"


rule edit_taxonomy:
    input:
        "results/crest4/mapped_reads_to_contigs_gg.tsv",
    output:
        "results/crest4/mapped_reads_to_contigs.tsv.edited",
    log:
        "logs/crest4/edit_taxonomy.log",
    conda:
        "../envs/base_python.yaml"
    priority: 10
    shell:
        "cp {input} {output} && "
        "echo 'Just a copy of {input} to keep consistency with previous versions' > {log}"
