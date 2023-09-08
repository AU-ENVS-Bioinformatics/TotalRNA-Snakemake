rule find_contigs_to_keep:
    input:
        "results/mRNA/mapped_reads_to_contigs.tsv",
    output:
        "results/mRNA/abundance_filtered/mapped_reads_to_contigs.tsv",
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
        table="results/mRNA/abundance_filtered/mapped_reads_to_contigs.tsv",
    output:
        "results/mRNA/abundance_filtered/Trinity_contigs_ncrna_filtered.fasta",
    conda:
        "../envs/seqtk.yaml"
    log:
        "logs/mRNA_abundance_filter/seqtk_subset.log",
    shell:
        "awk '{{print $1}}' {input.table} | seqtk subseq {input.fasta} - > {output} 2> {log}"
