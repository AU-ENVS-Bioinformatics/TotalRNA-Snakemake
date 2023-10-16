rule trinity:
    input:
        left=expand(
            "results/sortmerna/not_LSU/{sample}_fwd.fq.gz", sample=unique_samples
        ),
        right=expand(
            "results/sortmerna/not_LSU/{sample}_rev.fq.gz", sample=unique_samples
        ),
    output:
        dir=temp(directory("results/mRNA/trinity")),
        fas="results/mRNA/trinity.Trinity.fasta",
        map="results/mRNA/trinity.Trinity.fasta.gene_trans_map",
    log:
        "logs/trinity/trinity.log",
    params:
        extra=" ".join(config.get("trinity", "")),
    threads: config["threads"]["trinity"]
    resources:
        mem_gb=int(config["memory"]["trinity"]),
    wrapper:
        "v2.6.0/bio/trinity"
