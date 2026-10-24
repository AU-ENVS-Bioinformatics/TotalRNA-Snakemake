################################################################################
# PREFIX RECONSTRUCTED SEQUENCE IDS
################################################################################

rule prefix_reconstructed_sequence_ids:
    message:
        "[Assembly postprocessing] Prefixing sequence IDs in {input.fasta}"
    input:
        fasta=f"{ASSEMBLY_DIR}/cross_assembly/{{sequence_path}}.fasta"
    output:
        fasta=f"{ASSEMBLY_DIR}/cross_assembly/{{sequence_path}}_prefixed.fasta"
    wildcard_constraints:
        sequence_path=r"[^/]+/(?:SSU/rnaspades/transcripts|ITS_candidates/ITSx/[^/]+_ITSx\.full)"
    params:
        sample=lambda wc: wc.sequence_path.split("/")[0]
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.fasta})

        awk -v sample={params.sample} '
            /^>/ {{
                sub(/^>/, ">" sample "__")
            }}
            {{ print }}
        ' {input.fasta} > {output.fasta}
        """


################################################################################
# COMBINE PREFIXED RECONSTRUCTIONS ACROSS SAMPLES
################################################################################

rule combine_reconstructed_sequences:
    wildcard_constraints:
        target=CONCATENATION_TARGET_PATTERN
    message:
        "[Assembly postprocessing] Combining reconstructed {wildcards.target} sequences across all samples"
    input:
        fasta=lambda wc: prefixed_reconstructed_fastas(wc, results_dir=RESULTS_DIR, samples=SAMPLES)
    output:
        fasta=f"{ASSEMBLY_DIR}/cross_assembly/concatenated/{{target}}/cross_sample_{{target}}.fasta"
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.fasta})

        cat {input.fasta} > {output.fasta}
        """