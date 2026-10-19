rule geNomad:
    message:
        "[geNomad] Identifying RNA-SPAdes transcripts with viral evidence."
    conda:
        "../envs/geNomad.yaml"
    input:
        viral_filter_contigs=f"{VIRAL_ANALYSIS_DIR}/transcript/transcripts.filtered.fasta"
    output:
        plasmid_summary=f"{VIRAL_ANALYSIS_DIR}/identification/genomad/transcripts.filtered_summary/transcripts.filtered_plasmid_summary.tsv",
        virus_summary=f"{VIRAL_ANALYSIS_DIR}/identification/genomad/transcripts.filtered_summary/transcripts.filtered_virus_summary.tsv",
        viral_fasta=f"{VIRAL_ANALYSIS_DIR}/identification/genomad/transcripts.filtered_summary/transcripts.filtered_virus.fna",
        viral_genes=f"{VIRAL_ANALYSIS_DIR}/identification/genomad/transcripts.filtered_summary/transcripts.filtered_virus_genes.tsv",
    params:
        output_dir=f"{VIRAL_ANALYSIS_DIR}/identification/genomad",
        db=config["viral_sequence_analysis"]["geNomad"]["database"],
        options=config["viral_sequence_analysis"]["geNomad"].get("options", "")
    threads:
        config["viral_sequence_analysis"]["geNomad"].get("threads", 1)
    log:
        stdout=f"{VIRAL_ANALYSIS_DIR}/identification/genomad/logs/genomad.log"
    benchmark:
        f"{VIRAL_ANALYSIS_DIR}/identification/genomad/benchmarks/genomad.txt"
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.plasmid_summary})

        genomad end-to-end \
            {input.viral_filter_contigs} \
            {params.output_dir} \
            {params.db} \
            --threads {threads} \
            {params.options} \
            > {log.stdout} 2>&1
        """