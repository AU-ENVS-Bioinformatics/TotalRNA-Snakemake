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
        "-o results/crest4_results "
        "{params.extra} "
        "> {log} 2>&1 || true"
        # Crest4 sends an error signal


rule add_taxa_to_mapped_contigs:
    input:
        taxa="results/crest4_results/assignments.txt",
        otu="results/rrna/mapped_reads_to_contigs.tsv",
    output:
        "results/crest4_results/mapped_reads_to_contigs.tsv",
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/crest4/pandas.log",
    script:
        "../scripts/add_taxa_mapped_contigs.py"


rule edit_taxonomy:
    input:
        "results/crest4_results/mapped_reads_to_contigs.tsv",
    output:
        "results/crest4_results/mapped_reads_to_contigs.tsv.edited",
    log:
        "logs/crest4/edit_taxonomy.log",
    conda:
        "../envs/base_python.yaml"
    shell:
        "python3 workflow/external_scripts/crest4_taxonomy_edit.py "
        "{input} {output} >> {log} 2>&1"
