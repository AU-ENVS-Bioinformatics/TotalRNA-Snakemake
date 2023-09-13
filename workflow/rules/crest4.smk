AVAILABLE_THREADS = int(config.get("BLAST-THREADS", 50))


rule CREST4:
    input:
        fasta="results/MetaRib/all.dedup.filtered.fasta",
        otu="results/rrna/mapped_reads_to_contigs.tsv",
    output:
        "results/crest4_results/assignments.txt",
        "results/crest4_results/otus_by_rank.tsv",
        "results/crest4_results/otus_cumulative.tsv",
        "results/crest4_results/search.hits",
    params:
        CREST4_DIR=config.get("CREST4_DIR", "~/.crest4/"),
    conda:
        "../envs/crest4.yaml"
    log:
        "logs/crest4/crest4.log",
    threads: AVAILABLE_THREADS
    shell:
        "export CREST4_DIR={params.CREST4_DIR} ; "
        "crest4 -f {input.fasta} "
        "-t {threads} "
        "-u {input.otu} "
        "-o results/crest4_results "
        "> {log} 2>&1 || true"
        # Crest4 sends an error signal, so we need to make it always work
