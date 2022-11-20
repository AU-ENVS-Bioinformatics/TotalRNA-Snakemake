rule taxonomy_edit:
    input:
        f"results/CREST_Results/mapped_reads_to_contigs.tsv",
    output:
        f"results/CREST_Results/mapped_reads_to_contigs.tsv.edited",
    log:
        "logs/mapping_rrna/Taxonomy_edits.log",
    benchmark:
        "benchmarks/mapping_rrna/Taxonomy_edits.log"
    conda:
        "../envs/base_python.yaml"
    shell:
        "python3 workflow/scripts/TaxonomyEdits.py "
        "{input} {output} >> {log} 2>&1"
