rule viral_checkv:
    message:
        "[CheckV] Evaluating geNomad viral candidates for completeness, contamination and host-derived sequence."
    conda:
        "../envs/checkv.yaml"
    input:
        viral_fasta=f"{VIRAL_ANALYSIS_DIR}/identification/genomad/transcripts_summary/transcripts_virus.fna"
    output:
        quality=f"{VIRAL_ANALYSIS_DIR}/quality/checkv/quality_summary.tsv",
        completeness=f"{VIRAL_ANALYSIS_DIR}/quality/checkv/completeness.tsv",
        contamination=f"{VIRAL_ANALYSIS_DIR}/quality/checkv/contamination.tsv"
    params:
        output_dir=f"{VIRAL_ANALYSIS_DIR}/quality/checkv",
        db=config["viral_sequence_analysis"]["checkv"]["db"]
    threads:
        config["viral_sequence_analysis"]["checkv"].get("threads", 8)
    log:
        stdout=f"{VIRAL_ANALYSIS_DIR}/quality/checkv/logs/checkv.log"
    benchmark:
        f"{VIRAL_ANALYSIS_DIR}/quality/checkv/benchmarks/checkv.txt"
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.quality})

        if ! grep -q '^>' {input.viral_fasta}; then
            printf \
                "contig_id\tcontig_length\tprovirus\tproviral_length\tgene_count\tviral_genes\thost_genes\tcheckv_quality\tmiuvig_quality\tcompleteness\tcompleteness_method\tcontamination\tkmer_freq\twarnings\n" \
                > {output.quality}

            : > {output.completeness}
            : > {output.contamination}

            echo "No geNomad viral candidates were available for CheckV." \
                > {log.stdout}

            exit 0
        fi

        checkv end_to_end \
            {input.viral_fasta} \
            {params.output_dir} \
            -d {params.db} \
            -t {threads} \
            > {log.stdout} 2>&1
        """