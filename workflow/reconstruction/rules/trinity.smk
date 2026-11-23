rule trinity_coassembly:
    conda:
        "../envs/trinity.yaml"
    wildcard_constraints:
        target=COASSEMBLY_TARGET_PATTERNS["trinity"]
    message:
        "[Trinity] Co-assembling {wildcards.target} reads across all samples"
    input:
        r1=f"{ASSEMBLY_DIR}/coassembly/{{target}}/concatenated/{{target}}_1.fastq.gz",
        r2=f"{ASSEMBLY_DIR}/coassembly/{{target}}/concatenated/{{target}}_2.fastq.gz"
    output:
        transcripts=f"{ASSEMBLY_DIR}/coassembly/{{target}}/trinity/Trinity.fasta",
        link_transcript=f"{ASSEMBLY_DIR}/coassembly/{{target}}/{{target}}_transcripts.fasta"
    log:
        stdout=f"{ASSEMBLY_DIR}/coassembly/{{target}}/logs/{{target}}_trinity.log"
    benchmark:
        f"{ASSEMBLY_DIR}/coassembly/{{target}}/benchmarks/{{target}}_trinity_benchmark.txt"
    params:
        outdir=f"{ASSEMBLY_DIR}/coassembly/{{target}}/trinity",
        max_memory_gb=config["reconstruction"]["trinity"]["max_memory_gb"]
    threads:
        config["reconstruction"]["trinity"].get("threads", 12)
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {params.outdir}) $(dirname {log.stdout})
        Trinity \
            --seqType fq \
            --left {input.r1} \
            --right {input.r2} \
            --CPU {threads} \
            --max_memory {params.max_memory_gb}G \
            --output {params.outdir} \
            > {log.stdout} 2>&1
        if [[ ! -s "{output.transcripts}" ]]; then
            echo "ERROR: Trinity.fasta is missing or empty." >> {log.stdout}
            exit 1
        fi

        ln -sfn $(realpath {output.transcripts}) {output.link_transcript}
        """