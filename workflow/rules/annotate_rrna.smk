AVAILABLE_THREADS = int(config.get("BLAST-THREADS", 50))


rule CREST4:
    input:
        fasta="results/MetaRib/MetaRib/Abundance/all.dedup.filtered.fasta",
        otu="results/MetaRib/mapped_reads_to_contigs.tsv",
    output:
        "results/CREST_Results/assignments.txt",
    params:
        CREST4_DIR=config.get("CREST4_DIR", "~/.crest4/"),
    conda:
        "../envs/CREST.yaml"
    log:
        "logs/mapping_rrna/crest.log",
    threads: AVAILABLE_THREADS
    benchmark:
        "benchmarks/mapping_rrna/crest.log"
    shell:
        "export CREST4_DIR={params.CREST4_DIR} ; "
        "crest4 -f {input.fasta} "
        "-t {threads} "
        "-u {input.otu} "
        "-o results/CREST_Results "
        "> {log} 2>&1 || true"


rule add_taxa_mapped:
    input:
        taxa="results/CREST_Results/assignments.txt",
        otu="results/MetaRib/mapped_reads_to_contigs.tsv",
    output:
        "results/CREST_Results/mapped_reads_to_contigs.tsv",
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/mapping_rrna/add_taxa.log",
    script:
        "../scripts/add_taxa_mapped_contigs.py"
