################################################################################
# CO-ASSEMBLY LENGTH FILTERING
################################################################################

rule coassembly_length_filter:
    conda:
        "../envs/seqkit.yaml"
    wildcard_constraints:
        target="(?:" + "|".join(re.escape(target) for target in config["reconstruction"]["coassembly_targets"]) + ")"
    message:
        "[SeqKit] Filtering {wildcards.target} co-assembly by sequence length"
    input:
        transcripts=f"{ASSEMBLY_DIR}/coassembly/{{target}}/{{target}}_transcripts.fasta"
    output:
        filtered=f"{ASSEMBLY_DIR}/coassembly/{{target}}/postprocessing/{{target}}_length_filtered.fasta"
    params:
        minlength=lambda wc: config["reconstruction"]["coassembly_postprocessing"][wc.target]["minlength"],
        maxlength=lambda wc: str(config["reconstruction"]["coassembly_postprocessing"][wc.target]["maxlength"])
    log:
        stdout=f"{ASSEMBLY_DIR}/coassembly/{{target}}/logs/{{target}}_length_filter.log"
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
# CO-ASSEMBLY EXACT DEREPLICATION
################################################################################

rule coassembly_dereplication:
    conda:
        "../envs/vsearch.yaml"
    wildcard_constraints:
        target="(?:" + "|".join(re.escape(target) for target in config["reconstruction"]["coassembly_targets"]) + ")"
    message:
        "[VSEARCH] Exact dereplication of {wildcards.target} co-assembly"
    input:
        fasta=f"{ASSEMBLY_DIR}/coassembly/{{target}}/postprocessing/{{target}}_length_filtered.fasta"
    output:
        fasta=f"{ASSEMBLY_DIR}/coassembly/{{target}}/postprocessing/{{target}}_dereplicated.fasta",
        uc=f"{ASSEMBLY_DIR}/coassembly/{{target}}/postprocessing/{{target}}_dereplicated.uc"
    params:
        minlength=lambda wc: config["reconstruction"]["coassembly_postprocessing"][wc.target]["minlength"],
        maxlength=lambda wc: config["reconstruction"]["coassembly_postprocessing"][wc.target].get("maxlength") or 2147481646,
        options=lambda wc: config["reconstruction"]["coassembly_postprocessing"].get("dereplicate_options", "")
    log:
        stdout=f"{ASSEMBLY_DIR}/coassembly/{{target}}/logs/{{target}}_dereplication.log"
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
# CO-ASSEMBLY RFAM SEARCH
################################################################################

rule coassembly_cmsearch:
    conda:
        "../envs/infernal.yaml"
    wildcard_constraints:
        target="(?:" + "|".join(re.escape(target) for target in config["reconstruction"]["coassembly_targets"]) + ")"
    message:
        "[Infernal] Searching {wildcards.target} co-assembly against Rfam models"
    input:
        fasta=f"{ASSEMBLY_DIR}/coassembly/{{target}}/postprocessing/{{target}}_dereplicated.fasta",
        database=lambda wc: config["reconstruction"]["coassembly_postprocessing"][wc.target]["rfam_database"]
    output:
        tblout=f"{ASSEMBLY_DIR}/coassembly/{{target}}/postprocessing/{{target}}_rfam_hits.tbl"
    params:
        options=config["reconstruction"]["coassembly_postprocessing"]["cmsearch"].get("options", "--cut_ga")
    log:
        stdout=f"{ASSEMBLY_DIR}/coassembly/{{target}}/logs/{{target}}_cmsearch.log"
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
                {input.database} \
                {input.fasta} \
                > {log.stdout} 2>&1
        fi
        """