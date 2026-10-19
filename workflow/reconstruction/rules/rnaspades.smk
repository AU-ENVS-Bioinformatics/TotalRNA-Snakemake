################################################################################
# 1. SAMPLE SPECIFIC rRNA assembly
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
        options=config["reconstruction"]["rnaspades"].get("options","")
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
                {params.options} \
                > {log.stdout} 2>&1
        """

################################################################################
# 2. CROSS-SAMPLE NON-rRNA CO-ASSEMBLY
################################################################################

rule rnaspades_cross_sample:
    conda:
        "../envs/rnaspades.yaml"
    wildcard_constraints:
        target=CROSS_SAMPLE_TARGET_PATTERN
    message:
        "[rnaSPAdes] Co-assembling non-rRNA reads across all samples"
    input:
        r1=f"{ASSEMBLY_DIR}/coassembly/{{target}}/concatenated/nonrRNA_1.fastq.gz",
        r2=f"{ASSEMBLY_DIR}/coassembly/{{target}}/concatenated/nonrRNA_2.fastq.gz"
    output:
        transcripts=f"{ASSEMBLY_DIR}/coassembly/{{target}}/rnaspades/transcripts.fasta"
    log:
        stdout=f"{ASSEMBLY_DIR}/coassembly/{{target}}/logs/nonrRNA_rnaspades.log",
        time=f"{ASSEMBLY_DIR}/coassembly/{{target}}/benchmarks/nonrRNA_rnaspades_time.txt",
    benchmark:
        f"{ASSEMBLY_DIR}/coassembly/{{target}}/benchmarks/nonrRNA_rnaspades_benchmark.txt"
    params:
        outdir=f"{ASSEMBLY_DIR}/coassembly/{{target}}/rnaspades",
        options=config["reconstruction"]["rnaspades"].get("options", "")
    threads:
        config["reconstruction"]["rnaspades"].get("threads", 12)
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
                {params.options} \
                > {log.stdout} 2>&1
        """