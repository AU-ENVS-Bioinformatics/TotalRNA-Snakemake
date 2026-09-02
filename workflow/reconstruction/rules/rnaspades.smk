################################################################################
# RNA SPADES RECONSTRUCTION
################################################################################

RECONSTRUCTION_METHOD = config["reconstruction"]["method"] #e.g. rnaspades
SECONDARY_RNA_METHOD = config["RNA"].get("secondary_method", "none")

SUPPORTED_SECONDARY_METHODS = {
    "none",
    "ITS",
}

################################################################################
# SAMPLE-SPECIFIC RECONSTRUCTION TARGETS
################################################################################

# SSU is always reconstructed.
SAMPLE_RECONSTRUCTION_TARGETS = ["SSU"]

if SECONDARY_RNA_METHOD not in SUPPORTED_SECONDARY_METHODS:
    raise ValueError(
        f"Unsupported RNA.secondary_method: {SECONDARY_RNA_METHOD}. Supported values are: " + ", ".join(sorted(SUPPORTED_SECONDARY_METHODS))
    )

# ITS candidates are only reconstructed when the secondary ITS branch is active.
if SECONDARY_RNA_METHOD == "ITS":
    SAMPLE_RECONSTRUCTION_TARGETS.append("ITS_candidates") #ensure it fits with the folder name

# SSU is always reconstructed, but the ITS can be added
SAMPLE_TARGET_PATTERN = "|".join(SAMPLE_RECONSTRUCTION_TARGETS)

# Non-rRNA is always reconstructed.
CROSS_SAMPLE_TARGET_PATTERN = "non_rRNA"

################################################################################
# INPUT ROUTING FOR SAMPLE-SPECIFIC ASSEMBLIES
################################################################################

def rnaspades_sample_input(
    wildcards,
    read,
):
    """
    Return the correct FASTQ input for an SSU or ITS_candidates assembly.
    """

    sample = wildcards.sample
    target = wildcards.target

    if target == "SSU":
        return f"{RNA_CLASSIFIED_DIR}/{sample}/SSU/{sample}_SSU_{read}.fastq.gz"

    if target == "ITS_candidates":
        return f"{RNA_INTERMEDIATE_DIR}/{sample}/ITS_candidates/{sample}_ITS_candidates_{read}.fastq.gz"

    raise ValueError(
        f"Unsupported sample reconstruction target: {target}"
    )


################################################################################
# 1. SAMPLE-SPECIFIC SSU AND OPTIONAL ITS ASSEMBLY
################################################################################

if RECONSTRUCTION_METHOD in ["spades", "rnaspades"]:
    rule rnaspades_sample_specific:
        conda:
            "../envs/rnaspades.yaml"
        wildcard_constraints:
            target=SAMPLE_TARGET_PATTERN
        message:
            "[rnaSPAdes] Reconstructing {wildcards.target} for sample {wildcards.sample}"
        input:
            r1=lambda wc: rnaspades_sample_input(wc, read=1),
            r2=lambda wc: rnaspades_sample_input(wc, read=2),
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