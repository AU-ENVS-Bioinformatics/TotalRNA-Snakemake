rule trinity_bwa_index:
    input:
        "results/mRNA/trinity.Trinity.fasta",
    output:
        idx=multiext(
            "results/mRNA/trinity.Trinity.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"
        ),
    log:
        "logs/bwa_index/trinity.log",
    params:
        algorithm="bwtsw",
    wrapper:
        "v2.6.0/bio/bwa/index"


rule bwa_contig:
    input:
        contig=multiext(
            "results/mRNA/trinity.Trinity.fasta",
            "",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa",
        ),
        read="results/sortmerna/not_LSU/{sample_dir}.fq.gz",
    output:
        temp("results/mRNA/bwa/{sample_dir}.sam"),
    log:
        "logs/mRNA_bwa/{sample_dir}.log",
    threads: 50
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        bwa mem -t {threads} {input.contig[0]} {input.read} > {output} 2> {log}
        """


rule sort_bwa_contig:
    input:
        "results/mRNA/bwa/{sample_dir}.sam",
    output:
        "results/mRNA/bwa/{sample_dir}_sorted.bam",
    log:
        "logs/mRNA_samtools/{sample_dir}_sort.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -b -S {input} | samtools sort -o {output} -
        """


rule samtools_index:
    input:
        "results/mRNA/bwa/{sample_dir}_sorted.bam",
    output:
        "results/mRNA/bwa/{sample_dir}_sorted.bam.bai",
    log:
        "logs/samtools/index/{sample_dir}.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v2.6.0/bio/samtools/index"


rule samtools_idxstats:
    input:
        bam="results/mRNA/bwa/{sample_dir}_sorted.bam",
        idx="results/mRNA/bwa/{sample_dir}_sorted.bam.bai",
    output:
        "results/mRNA/bwa/{sample_dir}_sorted.bam.idxstats",
    log:
        "logs/samtools/idxstats/{sample_dir}.log",
    params:
        extra="",  # optional params string
    wrapper:
        "v2.6.0/bio/samtools/idxstats"


rule map_reads_to_trinity_contigs:
    input:
        left=expand(
            "results/sortmerna/not_LSU/{sample}_fwd.fq.gz", sample=unique_samples
        ),
        right=expand(
            "results/sortmerna/not_LSU/{sample}_rev.fq.gz", sample=unique_samples
        ),
        contigs="results/mRNA/trinity.Trinity.fasta",
    output:
        bam="results/mRNA/trinity.Trinity.bam",
        bai="results/mRNA/trinity.Trinity.bam.bai",
    threads: 50
    log:
        "logs/trinity/trinity_bam.log",
    script:
        "scripts/map_reads_to_trinity_contigs.py"
