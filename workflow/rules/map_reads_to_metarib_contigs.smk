rule metarib_bwa_index:
    input:
        "results/MetaRib/all.dedup.filtered.fasta",
    output:
        idx=multiext(
            "results/MetaRib/all.dedup.filtered.fasta",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa",
        ),
    log:
        "logs/bwa/index/metarib.log",
    params:
        algorithm="bwtsw",
    wrapper:
        "v2.6.0/bio/bwa/index"


rule bwa_contig_rrna:
    input:
        contig=multiext(
            "results/MetaRib/all.dedup.filtered.fasta",
            "",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa",
        ),
        read="results/rrna/{sample_dir}.fq.gz",
    output:
        temp("results/rrna/bwa/{sample_dir}.sam"),
    log:
        "logs/bwa/{sample_dir}_rrna.log",
    threads: 50
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        bwa mem -t {threads} {input.contig[0]} {input.read} > {output} 2> {log}
        """


rule sort_bwa_contig_rRNA:
    input:
        "results/rrna/bwa/{sample_dir}.sam",
    output:
        "results/rrna/bwa/{sample_dir}_sorted.bam",
    log:
        "logs/samtools/{sample_dir}_sort_rrna.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -b -S {input} | samtools sort -o {output} -
        """


rule samtools_index_rrna:
    input:
        "results/rrna/bwa/{sample_dir}_sorted.bam",
    output:
        "results/rrna/bwa/{sample_dir}_sorted.bam.bai",
    log:
        "logs/samtools/index/{sample_dir}_rrna.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v2.6.0/bio/samtools/index"


rule samtools_idxstats_rrna:
    input:
        bam="results/rrna/bwa/{sample_dir}_sorted.bam",
        idx="results/rrna/bwa/{sample_dir}_sorted.bam.bai",
    output:
        "results/rrna/bwa/{sample_dir}_sorted.bam.idxstats",
    log:
        "logs/samtools/idxstats/{sample_dir}_rrna.log",
    params:
        extra="",  # optional params string
    wrapper:
        "v2.6.0/bio/samtools/idxstats"


rule processing_mapped_contigs_rrna:
    input:
        expand(
            "results/rrna/bwa/{sample}_{dir}_sorted.bam.idxstats",
            sample=unique_samples,
            dir=["fwd", "rev"],
        ),
    output:
        "results/rrna/mapped_reads_to_contigs.tsv",
    log:
        "logs/samtools/mapped_reads_table_rrna.log",
    params:
        extension_fwd="_fwd",
        extension_rev="_rev",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/process_mapped_reads_contigs.py"
