rule prepare_assemble_reads:
    input:
        R1=expand(
            f"results/sortmerna/not_LSU/{{sample}}_fwd.fq.gz",
            sample=unique_samples,
        ),
        R2=expand(
            f"results/sortmerna/not_LSU/{{sample}}_rev.fq.gz",
            sample=unique_samples,
        ),
    output:
        R1=expand(
            f"results/mRNA/renamed/{{sample}}_R1.fastq",
            sample=unique_samples,
        ),
        R2=expand(
            f"results/mRNA/renamed/{{sample}}_R2.fastq",
            sample=unique_samples,
        ),
        readsdir=directory(f"results/mRNA/renamed/"),
    log:
        "logs/assemble_mRNA/pigz.log",
    benchmark:
        "benchmarks/assemble_mRNA/pigz.log"
    threads: int(config.get("PIGZ-THREADS", 50))
    conda:
        "../envs/base_python.yaml"
    script:
        "../scripts/pigz_reads.py"


rule assemble_reads:
    input:
        R1=expand(f"results/mRNA/renamed/{{sample}}_R1.fastq", sample=unique_samples),
        R2=expand(
            f"results/mRNA/renamed/{{sample}}_R2.fastq",
            sample=unique_samples,
        ),
        readsdir=f"results/mRNA/renamed/",
    output:
        outdir=directory(f"results/mRNA/trinity/"),
        outfile="results/mRNA/trinity/Trinity.fasta",
    params:
        script=config.get("CoMW_REPOSITORY", "workflow/scripts/CoMW/")
        + "scripts/assemble_reads.py",
        extra=" ".join(config.get("assemble_reads", "")),
    threads: int(config.get("assemble_reads-THREADS", 50))
    conda:
        "../envs/CoMW.yaml"
    log:
        "logs/assemble_mRNA/assemble_reads.log",
    benchmark:
        "benchmarks/assemble_mRNA/assemble_reads.log"
    shell:
        "python {params.script} "
        "-i {input.readsdir} "
        "-o {output.outdir} "
        "-c {threads} "
        "{params.extra} "
        ">> {log} 2>&1"


rule filter_non_coding_rna:
    input:
        fasta="results/mRNA/trinity/Trinity.fasta",
    output:
        f"results/mRNA/trinity/contigs_ncrna_filtered.fasta",
    params:
        script=config.get("CoMW_REPOSITORY", "workflow/scripts/CoMW/")
        + "scripts/filter_ncRNA_edited.py",
        extra=" ".join(config.get("filter_ncRNA", "")),
    threads: int(config.get("filter_non_coding_rna-THREADS", 50))
    conda:
        "../envs/biopython.yaml"
    log:
        "logs/assemble_mRNA/filter_non_coding_rna.log",
    benchmark:
        "benchmarks/assemble_mRNA/filter_non_coding_rna.log"
    shell:
        "python {params.script} "
        "-f {input.fasta} "
        "-o {output} "
        "-t {threads} "
        "{params.extra} "
        ">> {log} 2>&1"


rule map_reads_to_contigs_mRNA:
    input:
        fasta=f"results/mRNA/trinity/contigs_ncrna_filtered.fasta",
        indir=f"results/mRNA/renamed/",
    output:
        f"results/mRNA/mapped_reads_to_contigs.tsv",
    params:
        script=config.get("CoMW_REPOSITORY", "workflow/scripts/CoMW/")
        + "scripts/map_reads_to_contigs.py",
        extra=" ".join(config.get("map_reads_to_contigs_mRNA", "")),
    threads: int(config.get("map_reads_to_contigs-THREADS", 50))
    conda:
        "../envs/CoMW.yaml"
    log:
        "logs/assemble_mRNA/map_reads_to_contigs.log",
    benchmark:
        "benchmarks/assemble_mRNA/map_reads_to_contigs.log"
    shell:
        "python {params.script} "
        "-f {input.fasta} "
        "-i {input.indir} "
        "-o {output} "
        "-t {threads} "
        "{params.extra} "
        ">> {log} 2>&1"


rule filter_table_by_abundance:
    input:
        fasta=f"results/mRNA/trinity/contigs_ncrna_filtered.fasta",
        indir=f"results/mRNA/mapped_reads_to_contigs.tsv",
    output:
        f"results/mRNA/mapped_reads_to_contigs_AbundanceFiltered.tsv",
    params:
        script=config.get("CoMW_REPOSITORY", "workflow/scripts/CoMW/")
        + "scripts/filter_table_by_abundance.py",
        extra=" ".join(config.get("filter_table_by_abundance", "")),
    conda:
        "../envs/vegan.yaml"
    log:
        "logs/assemble_mRNA/filter_table_by_abundance.log",
    benchmark:
        "benchmarks/assemble_mRNA/filter_table_by_abundance.log"
    shell:
        "cd results/mRNA && "
        "python {params.script} "
        "-i ./mapped_reads_to_contigs.tsv "
        "-f ./trinity/contigs_ncrna_filtered.fasta "
        "{params.extra} "
        "-o ./mapped_reads_to_contigs "
        ">> ../../{log} 2>&1"


rule align_contigs_to_database:
    input:
        fasta=f"results/mRNA/trinity/contigs_ncrna_filtered.fasta",
    output:
        f"results/mRNA/ML_SWORD_result.tsv",
    params:
        script=config.get("CoMW_REPOSITORY", "workflow/scripts/CoMW/")
        + "scripts/align_contigs_to_database.py",
        extra=" ".join(config.get("align_contigs_to_database", "")),
    threads: int(config.get("align_contigs_to_database-THREADS", 50))
    conda:
        "../envs/pyfasta.yaml"
    log:
        "logs/assemble_mRNA/align_contigs_to_database.log",
    benchmark:
        "benchmarks/assemble_mRNA/align_contigs_to_database.log"
    shell:
        "python {params.script} "
        "-f {input.fasta} "
        "-o {output} "
        "-t {threads} "
        "{params.extra} "
        ">> {log} 2>&1"
