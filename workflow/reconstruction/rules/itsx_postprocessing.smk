################################################################################
# FILTER ITSx REGIONS BY LENGTH
################################################################################

rule filter_itsx_regions:
    conda:
        "../envs/seqkit.yaml"
    wildcard_constraints:
        region=r"full|ITS1|5_8S|ITS2"
    message:
        "[ITSx postprocessing] Filtering {wildcards.region} sequences for {wildcards.sample}"
    input:
        fasta=lambda wc: (
            f"{ASSEMBLY_DIR}/cross_assembly/{wc.sample}/ITS_candidates/ITSx/{wc.sample}_ITSx.full_and_partial.fasta"
            if wc.region == "full"
            else f"{ASSEMBLY_DIR}/cross_assembly/{wc.sample}/ITS_candidates/ITSx/{wc.sample}_ITSx.{wc.region}.full_and_partial.fasta"
        )
    output:
        fasta=f"{ASSEMBLY_DIR}/cross_assembly/{{sample}}/ITS_candidates/ITSx/filtered/{{sample}}_ITSx.{{region}}.full_and_partial.fasta"
    log:
        stdout=f"{ASSEMBLY_DIR}/cross_assembly/{{sample}}/ITS_candidates/ITSx/filtered/logs/{{sample}}_{{region}}.log"
    benchmark:
        f"{ASSEMBLY_DIR}/cross_assembly/{{sample}}/ITS_candidates/ITSx/filtered/benchmarks/{{sample}}_{{region}}.txt"
    params:
        minlength=lambda wc: config["reconstruction"]["itsx"]["minlength"][wc.region]
    threads:
        1
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.fasta}) $(dirname {log.stdout})

        if [[ -s {input.fasta} ]]; then
            seqkit seq \
                -m {params.minlength} \
                {input.fasta} \
                > {output.fasta} 2> {log.stdout}
        else
            : > {output.fasta}
            echo "[ITSx postprocessing] Empty input FASTA: {input.fasta}" > {log.stdout}
        fi
        """
