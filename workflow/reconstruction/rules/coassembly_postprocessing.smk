################################################################################
# ASSEMBLY LENGTH FILTERING
################################################################################

rule assembly_length_filter:
    conda:
        "../envs/seqkit.yaml"
    wildcard_constraints:
        assembly="coassembly|cross_assembly",
        target=ASSEMBLY_POSTPROCESSING_TARGET_PATTERN
    message:
        "[SeqKit] Length filtering {wildcards.assembly}/{wildcards.target}"
    input:
        transcripts=lambda wc: assembly_postprocessing_input(wc, results_dir=RESULTS_DIR, config=config)
    output:
        filtered=f"{ASSEMBLY_DIR}/{{assembly}}/{{target}}/postprocessing/{{target}}_length_filtered.fasta"
    params:
        minlength=lambda wc: config["reconstruction"]["coassembly_postprocessing"][wc.target]["minlength"],
        maxlength=lambda wc: str(config["reconstruction"]["coassembly_postprocessing"][wc.target]["maxlength"])
    log:
        stdout=f"{ASSEMBLY_DIR}/{{assembly}}/{{target}}/logs/{{target}}_length_filter.log"
    threads:
        1
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.filtered}) $(dirname {log.stdout})

        length_options=(--min-len {params.minlength})

        if [[ {params.maxlength} != None ]]; then
            length_options+=(--max-len {params.maxlength})
        fi

        seqkit seq \
            --threads {threads} \
            "${{length_options[@]}}" \
            {input.transcripts} \
            > {output.filtered} \
            2> {log.stdout}
        """


################################################################################
# ASSEMBLY EXACT DEREPLICATION
################################################################################

rule assembly_dereplication:
    conda:
        "../envs/vsearch.yaml"
    wildcard_constraints:
        assembly="coassembly|cross_assembly",
        target=ASSEMBLY_POSTPROCESSING_TARGET_PATTERN
    message:
        "[VSEARCH] Exact dereplication of {wildcards.assembly}/{wildcards.target}"
    input:
        fasta=f"{ASSEMBLY_DIR}/{{assembly}}/{{target}}/postprocessing/{{target}}_length_filtered.fasta"
    output:
        fasta=f"{ASSEMBLY_DIR}/{{assembly}}/{{target}}/postprocessing/{{target}}_dereplicated.fasta",
        uc=f"{ASSEMBLY_DIR}/{{assembly}}/{{target}}/postprocessing/{{target}}_dereplicated.uc"
    params:
        minlength=lambda wc: config["reconstruction"]["coassembly_postprocessing"][wc.target]["minlength"],
        maxlength=lambda wc: config["reconstruction"]["coassembly_postprocessing"][wc.target].get("maxlength") or 2147481646,
        options=config["reconstruction"]["coassembly_postprocessing"].get("dereplicate_options", "")
    log:
        stdout=f"{ASSEMBLY_DIR}/{{assembly}}/{{target}}/logs/{{target}}_dereplication.log"
    threads:
        1
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.fasta}) $(dirname {log.stdout})

        if [ ! -s {input.fasta} ]; then
            : > {output.fasta}
            : > {output.uc}
            printf 'No contigs passed length filtering.\n' > {log.stdout}
        else
            vsearch \
                --derep_fulllength {input.fasta} \
                --output {output.fasta} \
                --uc {output.uc} \
                --minseqlength {params.minlength} \
                --maxseqlength {params.maxlength} \
                {params.options} \
                > {log.stdout} 2>&1
        fi
        """


################################################################################
# ASSEMBLY RFAM SEARCH
################################################################################

rule assembly_cmsearch:
    conda:
        "../envs/infernal.yaml"
    wildcard_constraints:
        assembly="coassembly|cross_assembly",
        target=ASSEMBLY_POSTPROCESSING_TARGET_PATTERN
    message:
        "[Infernal] Searching {wildcards.assembly}/{wildcards.target} against Rfam"
    input:
        fasta=f"{ASSEMBLY_DIR}/{{assembly}}/{{target}}/postprocessing/{{target}}_dereplicated.fasta"
    output:
        tblout=f"{ASSEMBLY_DIR}/{{assembly}}/{{target}}/postprocessing/{{target}}_rfam_hits.tbl"
    params:
        database=lambda wc: config["reconstruction"]["coassembly_postprocessing"][wc.target]["rfam_database"],
        options=config["reconstruction"]["coassembly_postprocessing"]["cmsearch"].get("options", "--cut_ga")
    log:
        stdout=f"{ASSEMBLY_DIR}/{{assembly}}/{{target}}/logs/{{target}}_cmsearch.log"
    threads:
        config["reconstruction"]["coassembly_postprocessing"]["cmsearch"].get("threads", 16)
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.tblout}) $(dirname {log.stdout})

        if [ ! -s {input.fasta} ]; then
            printf '# No sequences available for Rfam screening.\n' > {output.tblout}
            printf 'Skipping cmsearch: input FASTA is empty.\n' > {log.stdout}
        else
            cmsearch \
                --cpu {threads} \
                --noali \
                {params.options} \
                --tblout {output.tblout} \
                {params.database} \
                {input.fasta} \
                > {log.stdout} 2>&1
        fi
        """