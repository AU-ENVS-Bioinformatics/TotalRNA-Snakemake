################################################################################
# 1. SAMPLE SPECIFIC rRNA assembly - primarily useful for phylogenetic analysis
################################################################################
rule rnaspades_sample_specific:
    conda:
        "../envs/rnaspades.yaml"
    wildcard_constraints:
        target=SAMPLE_TARGET_PATTERN
    message:
        "[rnaSPAdes] Reconstructing {wildcards.target} for sample {wildcards.sample}"
    input:
        r1=lambda wc: rnaspades_sample_input(
            wc,
            read=1,
            results_dir=RESULTS_DIR,
        ),
        r2=lambda wc: rnaspades_sample_input(
            wc,
            read=2,
            results_dir=RESULTS_DIR,
        ),
    output:
        transcripts=f"{ASSEMBLY_DIR}/{{sample}}/{{target}}/rnaspades/transcripts.fasta"
    log:
        stdout=f"{ASSEMBLY_DIR}/{{sample}}/{{target}}/logs/{{sample}}_{{target}}_rnaspades.log",
        time=f"{ASSEMBLY_DIR}/{{sample}}/{{target}}/benchmarks/{{sample}}_{{target}}_rnaspades_time.txt"
    benchmark:
        f"{ASSEMBLY_DIR}/{{sample}}/{{target}}/benchmarks/{{sample}}_{{target}}_rnaspades_benchmark.txt"
    params:
        outdir=f"{ASSEMBLY_DIR}/{{sample}}/{{target}}/rnaspades",
        logdir=f"{ASSEMBLY_DIR}/{{sample}}/{{target}}/logs",
        memory=config["reconstruction"]["rnaspades"].get("memory","")
    threads:
        config["reconstruction"]["rnaspades"].get("threads",12)
    shell:
        r"""
        set -euo pipefail

        mkdir -p {params.outdir}

        /usr/bin/time \
            -v \
            -o {log.time} \
            rnaspades.py \
                -1 {input.r1} \
                -2 {input.r2} \
                -t {threads} \
                -o {params.outdir} \
                --memory {params.memory} \
                > {log.stdout} 2>&1
        """