############################################
# FILTER CONTIGS FOR VIRAL ANALYSIS
############################################

rule viral_filter_contigs:
    message:
        "[Seqkit seq] Retain reconstructed non-rRNA contigs that meet the minimum sequence-length requirement for viral candidate analysis."
    conda:
        "../envs/seqkit.yaml"
    input:
        nonrRNA_transcripts=f"{ASSEMBLY_DIR}/coassembly/non_rRNA/rnaspades/transcripts.fasta"
    output:
        viral_filter_contigs=f"{VIRAL_ANALYSIS_DIR}/transcript/transcripts.filtered.fasta",
    params:
        options=config["viral_sequence_analysis"]["seqkit"].get("options",""),
    log:
        stdout=f"{VIRAL_ANALYSIS_DIR}/transcript/logs/filter_contigs.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.viral_filter_contigs})

        seqkit seq \
            {params.options} \
            {input.nonrRNA_transcripts} \
            > {output.viral_filter_contigs} \
            2> {log.stdout}
            """

rule viral_contigs_stats:
    message:
        "[Seqkit seq] Retain reconstructed non-rRNA contigs that meet the minimum sequence-length requirement for viral candidate analysis."
    conda:
        "../envs/seqkit.yaml"
    input:
        viral_filter_contigs=rules.viral_filter_contigs.output.viral_filter_contigs
    output:
        stats=f"{VIRAL_ANALYSIS_DIR}/transcript/transcripts.filtered.stats.tsv"
    log:
        stdout=f"{VIRAL_ANALYSIS_DIR}/transcript/logs/stats_contigs.log"
    shell:
        r"""
        mkdir -p {params.output_dir}

        seqkit stats \
            --tabular \
            {input.viral_filter_contigs} \
            > {output.stats}
        """