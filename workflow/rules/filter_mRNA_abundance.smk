rule find_contigs_to_keep:
    input:
        "results/mRNA/mapped_reads_to_contigs.tsv",
    output:
        "results/mRNA/abundance_filtered/contigs_to_keep.txt",
    params:
        minimum=config.get("minimum_abundance", 1),
    conda:
        "../envs/vegan.yaml"
    log:
        "logs/mRNA_abundance_filter/vegan_table.log",
    script:
        "../scripts/vegan_filter_abundance.R"


rule filter_contigs:
    input:
        fasta="results/mRNA/Trinity_contigs_ncrna_filtered.fasta",
        contigs="results/mRNA/abundance_filtered/contigs_to_keep.txt",
    output:
        "results/mRNA/abundance_filtered/Trinity_contigs_ncrna_filtered.fasta",
    conda:
        "../envs/seqtk.yaml"
    log:
        "logs/mRNA_abundance_filter/seqtk_subset.log",
    shell:
        "seqtk subseq {input.fasta} {input.contigs} > {output} 2> {log}"
