################################################################################
# SSU: RETAIN WHOLE CONTIGS WITH SUFFICIENT SSU COVERAGE
################################################################################

rule assembly_filter_ssu:
    wildcard_constraints:
        assembly="coassembly|cross_assembly"
    message:
        "[SSU filtering] Filtering {wildcards.assembly}/SSU by Rfam coverage"
    input:
        fasta=f"{ASSEMBLY_DIR}/{{assembly}}/SSU/postprocessing/SSU_dereplicated.fasta",
        tblout=f"{ASSEMBLY_DIR}/{{assembly}}/SSU/postprocessing/SSU_rfam_hits.tbl"
    output:
        fasta=f"{ASSEMBLY_DIR}/{{assembly}}/SSU/postprocessing/SSU_rfam_filtered.fasta",
        removed=f"{ASSEMBLY_DIR}/{{assembly}}/SSU/postprocessing/SSU_rfam_removed.fasta",
        report=f"{ASSEMBLY_DIR}/{{assembly}}/SSU/postprocessing/SSU_rfam_filter_report.tsv"
    params:
        script=f"{workflow.basedir}/reconstruction/scripts/filter_ssu.py",
        coverage=config["reconstruction"]["coassembly_postprocessing"]["SSU"].get("min_ssu_contig_coverage", 0.80)
    log:
        stdout=f"{ASSEMBLY_DIR}/{{assembly}}/SSU/logs/SSU_filter_rfam.log"
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

################################################################################
# SSU: SPLIT RETAINED CONTIGS INTO DOMAINS
################################################################################

rule assembly_split_ssu_domains:
    conda:
        "../envs/seqkit.yaml"
    wildcard_constraints:
        assembly="cross_assembly",
        domain="archaea|bacteria|eukaryotes"
    message:
        "[SSU grouping] Extracting {wildcards.domain} from {wildcards.assembly}/SSU"
    input:
        fasta=f"{ASSEMBLY_DIR}/{{assembly}}/SSU/postprocessing/SSU_rfam_filtered.fasta",
        report=f"{ASSEMBLY_DIR}/{{assembly}}/SSU/postprocessing/SSU_rfam_filter_report.tsv"
    output:
        fasta=f"{ASSEMBLY_DIR}/{{assembly}}/SSU/domains/{{domain}}.fasta",
        ids=f"{ASSEMBLY_DIR}/{{assembly}}/SSU/domains/{{domain}}.ids"
    log:
        stdout=f"{ASSEMBLY_DIR}/{{assembly}}/SSU/logs/{{domain}}_split.log"
    threads:
        1
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.fasta}) $(dirname {log.stdout})

        awk -F '\t' -v domain={wildcards.domain} '
            BEGIN {{
                groups["SSU_rRNA_archaea"] = "archaea"
                groups["SSU_rRNA_bacteria"] = "bacteria"
                groups["SSU_rRNA_eukarya"] = "eukaryotes"
                groups["SSU_rRNA_microsporidia"] = "eukaryotes"
                groups["SSU_trypano_mito"] = "eukaryotes"
            }}
            $12 == "retained" && groups[$3] == domain {{
                print $1
            }}
        ' {input.report} > {output.ids}

        if [ -s {output.ids} ]; then
            seqkit grep \
                --threads {threads} \
                --pattern-file {output.ids} \
                {input.fasta} \
                > {output.fasta} \
                2> {log.stdout}
        else
            : > {output.fasta}
            printf 'No retained contigs assigned to {wildcards.domain}.\n' > {log.stdout}
        fi
        """

################################################################################
# SAMPLE-DERIVED SSU: DOMAIN-SPECIFIC CHIMERA SCREENING
################################################################################

rule cross_assembly_ssu_domain_chimera:
    conda:
        "../envs/vsearch.yaml"
    wildcard_constraints:
        domain="archaea|bacteria|eukaryotes"
    message:
        "[VSEARCH] Orienting and screening cross-assembly SSU: {wildcards.domain}"
    input:
        fasta=f"{ASSEMBLY_DIR}/cross_assembly/SSU/domains/{{domain}}.fasta"
    output:
        oriented=f"{ASSEMBLY_DIR}/cross_assembly/SSU/domains/{{domain}}/postprocessing/SSU_oriented.fasta",
        unoriented=f"{ASSEMBLY_DIR}/cross_assembly/SSU/domains/{{domain}}/postprocessing/SSU_unoriented.fasta",
        chimeras=f"{ASSEMBLY_DIR}/cross_assembly/SSU/domains/{{domain}}/postprocessing/SSU_chimeras.fasta",
        borderline=f"{ASSEMBLY_DIR}/cross_assembly/SSU/domains/{{domain}}/postprocessing/SSU_chimera_borderline.fasta",
        report=f"{ASSEMBLY_DIR}/cross_assembly/SSU/domains/{{domain}}/postprocessing/SSU_chimera_report.tsv",
        nonchimeras=f"{ASSEMBLY_DIR}/cross_assembly/SSU/domains/{{domain}}/final/transcript.fasta"
    params:
        database=lambda wc: config["reconstruction"]["ssu_postprocessing"]["chimera"]["references"][wc.domain],
        options=config["reconstruction"]["ssu_postprocessing"]["chimera"].get("options", ""),
        minlength=config["reconstruction"]["coassembly_postprocessing"]["SSU"].get("minlength") or 1,
        maxlength=config["reconstruction"]["coassembly_postprocessing"]["SSU"].get("maxlength") or 2147481646
    log:
        stdout=f"{ASSEMBLY_DIR}/cross_assembly/SSU/logs/{{domain}}_chimera.log"
    threads:
        config["reconstruction"]["ssu_postprocessing"]["chimera"].get("threads", 12)
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.oriented}) $(dirname {output.nonchimeras}) $(dirname {log.stdout})

        # Initialize all outputs, including for empty domain FASTAs.
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
        else
            printf 'No contigs assigned to {wildcards.domain}.\n' >> {log.stdout}
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