HUMANN_MERGED_DIR = f"{FUNCTION_DIR_READS}/humann_merged"
HUMANN_JOIN_DIR = f"{HUMANN_MERGED_DIR}/humann_link_input"


################################################################################
# 1. LINK SAMPLE TABLES
################################################################################

rule humann_link_data:
    wildcard_constraints:
        metric="genefamilies|pathabundance|pathcoverage"
    input:
        table=f"{FUNCTION_DIR_READS}/{{sample}}/humann/{{sample}}_{{metric}}.tsv"
    output:
        table=f"{HUMANN_JOIN_DIR}/{{sample}}_{{metric}}.tsv"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname "{output.table}")"

        ln -sf "$(realpath "{input.table}")" "{output.table}"
        """


################################################################################
# 2. JOIN ALL SAMPLES FOR EACH METRIC
################################################################################

rule humann_join_tables:
    conda:
        "../envs/humann.yaml"
    wildcard_constraints:
        metric="genefamilies|pathabundance|pathcoverage"
    input:
        tables=expand(
            f"{HUMANN_JOIN_DIR}/{{sample}}_{{metric}}.tsv",
            sample=SAMPLES,
            allow_missing=True,
        )
    output:
        table=f"{HUMANN_MERGED_DIR}/merged_{{metric}}.tsv"
    params:
        input_dir=f"{HUMANN_JOIN_DIR}"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname "{output.table}")"

        humann_join_tables \
            --input "{params.input_dir}" \
            --file_name "{wildcards.metric}" \
            --output "{output.table}"
        """


################################################################################
# 3. SPLIT STRATIFIED AND UNSTRATIFIED TABLES
################################################################################

rule humann_split_stratified:
    conda:
        "../envs/humann.yaml"
    wildcard_constraints:
        metric="genefamilies|pathabundance|pathcoverage"
    input:
        table=f"{HUMANN_MERGED_DIR}/merged_{{metric}}.tsv"
    output:
        stratified=f"{HUMANN_MERGED_DIR}/merged_{{metric}}_stratified.tsv",
        unstratified=f"{HUMANN_MERGED_DIR}/merged_{{metric}}_unstratified.tsv"
    params:
        output_dir=f"{HUMANN_MERGED_DIR}"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "{params.output_dir}"

        humann_split_stratified_table \
            --input "{input.table}" \
            --output "{params.output_dir}"
        """


################################################################################
# 4. NORMALIZE ABUNDANCE TABLES
################################################################################

rule humann_renorm:
    conda:
        "../envs/humann.yaml"
    wildcard_constraints:
        metric="genefamilies|pathabundance"
    input:
        table=f"{HUMANN_MERGED_DIR}/merged_{{metric}}_unstratified.tsv"
    output:
        table=f"{HUMANN_MERGED_DIR}/merged_{{metric}}_relab.tsv"
    shell:
        r"""
        set -euo pipefail

        humann_renorm_table \
            --input "{input.table}" \
            --output "{output.table}" \
            --units relab
        """


################################################################################
# 5. REGROUP UNIREF90 GENE FAMILIES INTO REACTIONS
################################################################################

rule humann_regroup_rxn:
    conda:
        "../envs/humann.yaml"
    input:
        table=f"{HUMANN_MERGED_DIR}/merged_genefamilies_relab.tsv"
    output:
        table=f"{HUMANN_MERGED_DIR}/merged_genefamilies_rxn.tsv"
    shell:
        r"""
        set -euo pipefail

        humann_regroup_table \
            --input "{input.table}" \
            --groups uniref90_rxn \
            --output "{output.table}"
        """