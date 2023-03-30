DIAMOND_DATABASES_DIR = config.get("DIAMOND_DATABASES_DIR")
rule diamond_blastx:
    input:
        fname_fastq = ancient(choose_filepath(config)),
        fname_db = ancient(f"{DIAMOND_DATABASES_DIR}{{database}}.dmnd")
    output:
        fname = "results/mRNA/diamond/{database}.tsv.gz"
    log:
        "logs/diamond/{database}.log"
    benchmark:
        "benchmarks/diamond/{database}.log"
    params:
        extra=" ".join(config.get("diamond", "--compress 1")),
    threads: 50
    wrapper:    
        "v1.21.4/bio/diamond/blastx"
