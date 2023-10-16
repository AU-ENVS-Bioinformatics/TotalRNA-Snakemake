rule annotree:
    input:
        trinity="results/mRNA/abundance_filtered/Trinity_contigs_ncrna_filtered.fasta",
        count_table="results/mRNA/abundance_filtered/mapped_reads_to_contigs.tsv",
        database=config["ANNOTREE"]["DATABASE"],
        mapping=config["ANNOTREE"]["MAPPING"],
        brite=config["BRITE"],
    output:
        diamond="results/mRNA/diamond/annotree.tsv",
        brite_ann="results/mRNA/diamond/annotree_brite_ann.tsv",
        taxonomy="results/mRNA/diamond/annotree_taxonomy_ann.tsv",
        metabolism="results/mRNA/diamond/annotree_metabolism_ann.tsv",
    log:
        notebook="notebooks/annotree.ipynb",
    threads: config["threads"]["diamond"]
    conda:
        "../envs/notebook_diamond.yaml"
    notebook:
        "../notebooks/AnnoTree.py.ipynb"
