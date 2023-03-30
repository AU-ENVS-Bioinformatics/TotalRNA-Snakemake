AVAILABLE_THREADS = int(config.get("BLAST-THREADS", 50))


rule blast_silvamod:
    input:
        query=f"results/MetaRib/MetaRib/Abundance/all.dedup.filtered.fasta",
        blastdb=multiext(
            config.get("BLAST_DATABASE", ""),
            ".nhr",
            ".nin",
            ".nsq",
        ),
    output:
        f"results/ml_rRNA/ml_rRNA_silvamod.xml",
    log:
        "logs/ml_rRNA_silvamod.blast.log",
    benchmark:
        "benchmarks/ml_rRNA_silvamod.blast.log"
    threads: AVAILABLE_THREADS
    params:
        # Usable options and specifiers for the different output formats are listed here:
        # https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/blast/blastn.html.
        format="5",
        extra=" ".join(config.get("blastn", "")),
    wrapper:
        "v1.18.3/bio/blast/blastn"


rule LCAClassifier:
    input:
        silva=f"results/ml_rRNA/ml_rRNA_silvamod.xml",
        otu=f"results/MetaRib/mapped_reads_to_contigs.tsv",
    output:
        outdir=directory(f"results/CREST_Results"),
        otu=f"results/CREST_Results/mapped_reads_to_contigs.tsv",
    params:
        script=config.get("CREST_LCAClassifier_BINARY", ""),
        tmp=f"results/tmp_crest",
    conda:
        "../envs/base_python.yaml"
    log:
        "logs/mapping_rrna/crest.log",
    benchmark:
        "benchmarks/mapping_rrna/crest.log"
    shell:
        "{params.script} "
        "-t {input.otu} "
        "{input.silva} "
        "-o {params.tmp} >> {log} 2>&1 && "
        "mv {params.tmp}/* {output.outdir} && "
        "rm -rf {params.tmp}"
