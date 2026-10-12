################################################################################
# 1. EXTRACT AND LENGTH-FILTER PREDICTED SSU REGIONS
################################################################################

rule filter_ssu_length:
    conda:
        "../envs/seqkit.yaml"
    wildcard_constraints:
        kingdom=BARRNAP_KINGDOM_PATTERN
    message:
        "[SSU postprocessing] Length-filtering SSU sequences: {wildcards.kingdom}"
    input:
        fasta=f"{ASSEMBLY_DIR}/SSU_candidates/barrnap/{{kingdom}}/cross_sample_SSU_{{kingdom}}.barrnap.fasta"
    output:
        fasta=f"{ASSEMBLY_DIR}/SSU_candidates/postprocessing/{{kingdom}}/SSU_{{kingdom}}_filtered.fasta"
    log:
        stdout=f"{ASSEMBLY_DIR}/SSU_candidates/postprocessing/{{kingdom}}/logs/length_filter.log"
    params:
        min_length=config["reconstruction"]["ssu_postprocessing"].get("min_length", 1000)
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.fasta})

        seqkit seq \
            --min-len {params.min_length} \
            {input.fasta} \
            > {output.fasta} \
            2> {log.stdout}
        """

################################################################################
# 2. EXACT DEREPLICATION WITHIN EACH KINGDOM
################################################################################

rule dereplicate_ssu:
    conda:
        "../envs/vsearch.yaml"
    wildcard_constraints:
        kingdom=BARRNAP_KINGDOM_PATTERN
    message:
        "[SSU postprocessing] Exact dereplication: {wildcards.kingdom}"
    input:
        fasta=f"{ASSEMBLY_DIR}/SSU_candidates/postprocessing/{{kingdom}}/SSU_{{kingdom}}_filtered.fasta"
    output:
        fasta=f"{ASSEMBLY_DIR}/SSU_candidates/postprocessing/{{kingdom}}/SSU_{{kingdom}}_dereplicated.fasta",
        uc=f"{ASSEMBLY_DIR}/SSU_candidates/postprocessing/{{kingdom}}/SSU_{{kingdom}}_dereplicated.uc"
    log:
        stdout=f"{ASSEMBLY_DIR}/SSU_candidates/postprocessing/{{kingdom}}/logs/dereplication.log"
    threads:
        config["reconstruction"]["ssu_postprocessing"]["dereplication"].get("threads", 1)
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.fasta})

        if [ ! -s {input.fasta} ]; then
            : > {output.fasta}
            : > {output.uc}
            printf 'No SSU regions passed extraction and length filtering.\n' > {log.stdout}
        else
            vsearch \
                --derep_fulllength {input.fasta} \
                --output {output.fasta} \
                --uc {output.uc} \
                --strand plus \
                --minseqlength 1 \
                --sizeout \
                --maxseqlength 1000000 \
                --threads {threads} \
                --fasta_width 0 \
                > {log.stdout} 2>&1
        fi
        """


################################################################################
# 3. REFERENCE-BASED CHIMERA SCREENING WITHIN EACH KINGDOM
################################################################################

rule screen_ssu_chimeras:
    conda:
        "../envs/vsearch.yaml"
    wildcard_constraints:
        kingdom=BARRNAP_KINGDOM_PATTERN
    message:
        "[SSU postprocessing] Reference-based chimera screening: {wildcards.kingdom}"
    input:
        fasta=f"{ASSEMBLY_DIR}/SSU_candidates/postprocessing/{{kingdom}}/SSU_{{kingdom}}_dereplicated.fasta",
        reference=lambda wc: config["reconstruction"]["ssu_postprocessing"]["chimera"]["references"][wc.kingdom]
    output:
        nonchimeras=f"{ASSEMBLY_DIR}/SSU_candidates/postprocessing/{{kingdom}}/SSU_{{kingdom}}_nonchimeras.fasta",
        chimeras=f"{ASSEMBLY_DIR}/SSU_candidates/postprocessing/{{kingdom}}/SSU_{{kingdom}}_chimeras.fasta",
        borderline=f"{ASSEMBLY_DIR}/SSU_candidates/postprocessing/{{kingdom}}/SSU_{{kingdom}}_borderline.fasta",
        report=f"{ASSEMBLY_DIR}/SSU_candidates/postprocessing/{{kingdom}}/SSU_{{kingdom}}_uchime.tsv"
    log:
        stdout=f"{ASSEMBLY_DIR}/SSU_candidates/postprocessing/{{kingdom}}/logs/chimera_screening.log"
    params:
        options=config["reconstruction"]["ssu_postprocessing"]["chimera"].get("options", "")
    threads:
        config["reconstruction"]["ssu_postprocessing"]["chimera"].get("threads", 12)
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.nonchimeras})

        if [ ! -s {input.fasta} ]; then
            : > {output.nonchimeras}
            : > {output.chimeras}
            : > {output.borderline}
            : > {output.report}
            printf 'No dereplicated SSU sequences to screen.\n' > {log.stdout}
        else
            vsearch \
                --uchime_ref {input.fasta} \
                --db {input.reference} \
                --nonchimeras {output.nonchimeras} \
                --chimeras {output.chimeras} \
                --borderline {output.borderline} \
                --uchimeout {output.report} \
                --minseqlength 1 \
                --maxseqlength 1000000 \
                --threads {threads} \
                --fasta_width 0 \
                {params.options} \
                > {log.stdout} 2>&1
        fi
        """

################################################################################
# 4. CONCATENATE NON-CHIMERA AND BORDERLINE TO CREATE A CONTIG CATALOGUE FOR 
#    DOWNSTREAM PROCESSING (PHYLOGENY, NOVEL IDENTIFICATION ETC.)
################################################################################

rule concatenate_ssu:
    conda:
        "../envs/seqkit.yaml"
    wildcard_constraints:
        kingdom=BARRNAP_KINGDOM_PATTERN
    message:
        "[SSU postprocessing] Concatenating non-chimeras and borderline sequences: {wildcards.kingdom}"
    input:
        nonchimeras=f"{ASSEMBLY_DIR}/SSU_candidates/postprocessing/{{kingdom}}/SSU_{{kingdom}}_nonchimeras.fasta",
        borderline=f"{ASSEMBLY_DIR}/SSU_candidates/postprocessing/{{kingdom}}/SSU_{{kingdom}}_borderline.fasta"
    output:
        fasta=f"{ASSEMBLY_DIR}/SSU_candidates/postprocessing/{{kingdom}}/SSU_{{kingdom}}_catalogue.fasta"
    log:
        stdout=f"{ASSEMBLY_DIR}/SSU_candidates/postprocessing/{{kingdom}}/logs/concatenation.log"
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.fasta})

        seqkit concat \
            {input.nonchimeras} \
            {input.borderline} \
            > {output.fasta} \
            2> {log.stdout}
        """