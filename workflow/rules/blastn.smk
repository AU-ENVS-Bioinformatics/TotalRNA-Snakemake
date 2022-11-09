AVAILABLE_THREADS = int(workflow.cores * 0.75)
rule blast_silvamod128:
    input:
        query = f"{DEFAULT_DEST_FILEPATH}{METARIB_FILEPATH}MetaRib/Abundance/all.dedup.filtered.fasta",
        blastdb=multiext(config.get("BLAST_DATABASE", ""),
            ".nhr",
            ".nin",
            ".nsq",
        )
    output:
        f"{DEFAULT_DEST_FILEPATH}ml_rRNA/ml_rRNA_silvamod.xml"
    log:
        "logs/ml_rRNA_silvamod.blast.log"
    threads:
        AVAILABLE_THREADS
    params:
        # Usable options and specifiers for the different output formats are listed here:
        # https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/blast/blastn.html.
        format="5",
        extra=" ".join(config.get("blastn", ""))
    wrapper:
        "v1.18.3/bio/blast/blastn"

rule crest4:
    input: 
        silva = f"{DEFAULT_DEST_FILEPATH}ml_rRNA/ml_rRNA_silvamod.xml",
        otu = f"{DEFAULT_DEST_FILEPATH}{METARIB_FILEPATH}SSU_fastq/{OTU_FILEPATH}"
    output:
        outdir = directory(f"{DEFAULT_DEST_FILEPATH}CREST_Results"),
        otu = f"{DEFAULT_DEST_FILEPATH}CREST_Results/mapped_reads_to_contigs.tsv"
    params:
        script = config.get("CREST_LCAClassifier_BINARY", ""),
        tmp = f"{DEFAULT_DEST_FILEPATH}tmp_crest"
    conda:
        "../envs/base_python.yaml"
    log: "logs/mapping_rrna/crest.log"
    shell: 
        "{params.script} "
        "-t {input.otu} "
        "{input.silva} "
        "-o {params.tmp} >> {log} 2>&1 && "
        "mv {params.tmp}/* {output.outdir} && " 
        "rm -rf {params.tmp}"