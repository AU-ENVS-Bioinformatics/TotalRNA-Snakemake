################################################################################
# RECONSTRUCTION CONCATENATION SETTINGS
################################################################################

CONCATENATION_TARGETS = ["SSU"]

if SECONDARY_RNA_METHOD == "ITS":
    CONCATENATION_TARGETS.append("ITS")

CONCATENATION_TARGET_PATTERN = "|".join(CONCATENATION_TARGETS)


################################################################################
# INPUT ROUTING FOR CROSS-SAMPLE CONCATENATION
################################################################################

def prefixed_reconstructed_fastas(wildcards):
    """
    Return all prefixed per-sample FASTAs for SSU or ITS.
    """

    if wildcards.target == "SSU":
        return expand(
            f"{ASSEMBLY_DIR}/{{sample}}/SSU/rnaspades/transcripts_prefixed.fasta",
            sample=SAMPLES,
        )

    if wildcards.target == "ITS":
        return expand(
            f"{ASSEMBLY_DIR}/{{sample}}/ITS_candidates/ITSx/{{sample}}_ITSx.full_prefixed.fasta",
            sample=SAMPLES,
        )

    raise ValueError(
        f"Unsupported concatenation target: {wildcards.target}"
    )


################################################################################
# 1. PREFIX RECONSTRUCTED SEQUENCE IDS
################################################################################

rule prefix_reconstructed_sequence_ids:
    message:
        "[Assembly postprocessing] Prefixing sequence IDs in {input.fasta}"

    input:
        fasta=f"{ASSEMBLY_DIR}/{{sequence_path}}.fasta"

    output:
        fasta=f"{ASSEMBLY_DIR}/{{sequence_path}}_prefixed.fasta"

    wildcard_constraints:
        sequence_path=(
            r"[^/]+/(?:"
            r"SSU/rnaspades/transcripts|"
            r"ITS_candidates/ITSx/[^/]+_ITSx\.full"
            r")"
        )

    params:
        sample=lambda wc: wc.sequence_path.split("/")[0]

    shell:
        r"""
        set -euo pipefail

        awk -v sample="{params.sample}" '
            /^>/ {{
                sub(/^>/, ">" sample "__")
            }}
            {{ print }}
        ' "{input.fasta}" > "{output.fasta}"
        """


################################################################################
# 2. COMBINE PREFIXED RECONSTRUCTIONS ACROSS SAMPLES
################################################################################

rule combine_reconstructed_sequences:
    wildcard_constraints:
        target=CONCATENATION_TARGET_PATTERN
    message:
        "[Assembly postprocessing] Combining reconstructed {wildcards.target} sequences across all samples"
    input:
        fasta=prefixed_reconstructed_fastas
    output:
        fasta=f"{ASSEMBLY_DIR}/concatenated/{{target}}/cross_sample_{{target}}.fasta"
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.fasta})

        cat {input.fasta:q} > "{output.fasta}"
        """