rule taxonomy_edit:
    input: 
        f"{DEFAULT_DEST_FILEPATH}CREST_Results/mapped_reads_to_contigs.tsv"
    output: 
        f"{DEFAULT_DEST_FILEPATH}CREST_Results/mapped_reads_to_contigs.tsv.edited"
    log: 
        "logs/mapping_rrna/Taxonomy_edits.log"
    conda:
        "../envs/base_python.yaml"
    shell: 
        "python3 workflow/scripts/TaxonomyEdits.py "
        "{input} {output} >> {log} 2>&1"