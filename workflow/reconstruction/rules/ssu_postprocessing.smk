################################################################################
# SSU: RETAIN WHOLE CONTIGS WITH >=90% SSU COVERAGE
################################################################################

rule coassembly_filter_ssu:
    message:
        "[SSU filtering] Retaining whole contigs with sufficient SSU coverage"
    input:
        fasta=f"{ASSEMBLY_DIR}/coassembly/SSU/postprocessing/SSU_dereplicated.fasta",
        tblout=f"{ASSEMBLY_DIR}/coassembly/SSU/postprocessing/SSU_rfam_hits.tbl",
    output:
        fasta=f"{ASSEMBLY_DIR}/coassembly/SSU/postprocessing/SSU_rfam_filtered.fasta",
        removed=f"{ASSEMBLY_DIR}/coassembly/SSU/postprocessing/SSU_rfam_removed.fasta",
        report=f"{ASSEMBLY_DIR}/coassembly/SSU/postprocessing/SSU_rfam_filter_report.tsv"
    params:
        script=f"{workflow.basedir}/reconstruction/scripts/filter_ssu.py",
        coverage=config["reconstruction"]["coassembly_postprocessing"]["SSU"].get("min_ssu_contig_coverage", 0.90)
    log:
        stdout=f"{ASSEMBLY_DIR}/coassembly/SSU/logs/SSU_filter_rfam.log"
    threads:
        1
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.fasta}) $(dirname {log.stdout})

        python3 {params.script} \
            --fasta {input.fasta} \
            --tblout {input.tblout} \
            --output {output.fasta} \
            --removed_out {output.removed} \
            --report {output.report} \
            --min-ssu-coverage {params.coverage} \
            > {log.stdout} 2>&1
        """

################################################################################
# SSU: ORIENTATION AND REFERENCE-BASED CHIMERA SCREENING
################################################################################

rule coassembly_ssu_chimera:
    conda:
        "../envs/vsearch.yaml"
    message:
        "[VSEARCH] Orienting SSU contigs and screening for apparent chimeras"
    input:
        fasta=f"{ASSEMBLY_DIR}/coassembly/SSU/postprocessing/SSU_rfam_filtered.fasta"
    output:
        oriented=f"{ASSEMBLY_DIR}/coassembly/SSU/postprocessing/SSU_oriented.fasta",
        unoriented=f"{ASSEMBLY_DIR}/coassembly/SSU/postprocessing/SSU_unoriented.fasta",
        chimeras=f"{ASSEMBLY_DIR}/coassembly/SSU/postprocessing/SSU_chimeras.fasta",
        borderline=f"{ASSEMBLY_DIR}/coassembly/SSU/postprocessing/SSU_chimera_borderline.fasta",
        report=f"{ASSEMBLY_DIR}/coassembly/SSU/postprocessing/SSU_chimera_report.tsv",
        nonchimeras=f"{ASSEMBLY_DIR}/coassembly/SSU/final/transcript.fasta"
    params:
        database=config["RNA_classification"]["databases"]["sortmeRNA_ssu"],
        options=config["reconstruction"]["coassembly_postprocessing"]["SSU"]["chimera"].get("options", ""),
        minlength=config["reconstruction"]["coassembly_postprocessing"]["SSU"].get("minlength") or 1,
        maxlength=config["reconstruction"]["coassembly_postprocessing"]["SSU"].get("maxlength") or 2147481646
    log:
        stdout=f"{ASSEMBLY_DIR}/coassembly/SSU/logs/SSU_chimera.log"
    threads:
        config["reconstruction"]["coassembly_postprocessing"]["SSU"]["chimera"].get("threads", 12)
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.nonchimeras}) $(dirname {log.stdout})

        # Initialize outputs, including when the input is empty.
        for file in {output}; do
            : > "$file"
        done

        : > {log.stdout}

        if [ -s {input.fasta} ]; then
            vsearch \
                --orient {input.fasta} \
                --db {params.database} \
                --fastaout {output.oriented} \
                --notmatched {output.unoriented} \
                --threads {threads} \
                --minseqlength {params.minlength} \
                --maxseqlength {params.maxlength} \
                {params.options} \
                >> {log.stdout} 2>&1
        fi

        if [ -s {output.oriented} ]; then
            vsearch \
                --uchime_ref {output.oriented} \
                --db {params.database} \
                --nonchimeras {output.nonchimeras} \
                --chimeras {output.chimeras} \
                --borderline {output.borderline} \
                --uchimeout {output.report} \
                --threads {threads} \
                --minseqlength {params.minlength} \
                --maxseqlength {params.maxlength} \
                {params.options} \
                >> {log.stdout} 2>&1
        fi
        """