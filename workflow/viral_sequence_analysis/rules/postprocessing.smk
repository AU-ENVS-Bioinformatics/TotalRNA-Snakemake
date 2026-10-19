################################################################################
# 01: COMBINE GENOMAD AND CHECKV
################################################################################

rule viral_candidate_table:
    message:
        "[Viral integration] Combining geNomad and CheckV results."
    input:
        genomad=rules.geNomad.output.virus_summary,
        checkv=rules.viral_checkv.output.quality
    output:
        table=f"{VIRAL_ANALYSIS_DIR}/intermediate/01_candidate_table/viral_candidates.tsv",
        ids=f"{VIRAL_ANALYSIS_DIR}/intermediate/01_candidate_table/viral_candidate_ids.txt"
    params:
        script=f"{VIRAL_SCRIPT_DIR}/01_candidate_table.py"
    log:
        stdout=f"{VIRAL_ANALYSIS_DIR}/intermediate/01_candidate_table/logs/viral_candidate_table.log"
    shell:
        r"""
        set -euo pipefail
        python {params.script} \
            --genomad {input.genomad} \
            --checkv {input.checkv} \
            --output {output.table} \
            --ids {output.ids} \
            --log {log.stdout}
        """

################################################################################
# INDEPENDENT PER-SAMPLE SALMON ABUNDANCE RESULTS
################################################################################

rule viral_salmon_subset:
    message:
        "[Viral abundance] Selecting candidate-parent transcripts for {wildcards.sample}."
    input:
        ids=rules.viral_candidate_table.output.ids,
        quant=f"{FUNCTION_DIR}/salmon/transcripts/{{sample}}/quant.sf"
    output:
        quant=f"{VIRAL_ANALYSIS_DIR}/abundance/{{sample}}/viral_quant.sf"
    params:
        script=f"{VIRAL_SCRIPT_DIR}/salmon_subset.py"
    log:
        stdout=f"{VIRAL_ANALYSIS_DIR}/abundance/logs/{{sample}}.salmon_subset.log"
    shell:
        r"""
        set -euo pipefail
        python {params.script} \
            --ids {input.ids} \
            --quant {input.quant} \
            --output {output.quant} \
            --log {log.stdout}
        """

################################################################################
# 02: TRANSDECODER PROTEIN-TO-TRANSCRIPT MAPPING
################################################################################

rule transdecoder_protein_to_transcript:
    message:
        "[TransDecoder] Creating protein-to-transcript ID mapping."
    input:
        gff=f"{FUNCTION_DIR}/ORF_prediction/final_orf/ORF_genes.gff3"
    output:
        mapping=f"{VIRAL_ANALYSIS_DIR}/intermediate/02_protein_map/protein_to_transcript.tsv"
    params:
        script=f"{VIRAL_SCRIPT_DIR}/02_protein_map.py"
    log:
        stdout=f"{VIRAL_ANALYSIS_DIR}/intermediate/02_protein_map/logs/protein_to_transcript.log"
    shell:
        r"""
        set -euo pipefail
        python {params.script} \
            --gff {input.gff} \
            --output {output.mapping} \
            --log {log.stdout}
        """

################################################################################
# 03: VALIDATE TRANSDECODER AND EGGNOG IDENTIFIERS
################################################################################

rule check_protein_id_mapping:
    message:
        "[Viral integration] Validating TransDecoder and eggNOG identifier overlap."
    input:
        mapping=rules.transdecoder_protein_to_transcript.output.mapping,
        eggnog=f"{FUNCTION_DIR}/eggnog/eggnog.emapper.annotations"
    output:
        overlap=f"{VIRAL_ANALYSIS_DIR}/intermediate/03_mapping_validation/overlapping_protein_ids.txt",
        preview=f"{VIRAL_ANALYSIS_DIR}/intermediate/03_mapping_validation/overlapping_protein_ids.preview.txt",
        summary=f"{VIRAL_ANALYSIS_DIR}/intermediate/03_mapping_validation/protein_id_mapping_summary.tsv"
    params:
        script=f"{VIRAL_SCRIPT_DIR}/03_mapping_validation.py"
    log:
        stdout=f"{VIRAL_ANALYSIS_DIR}/intermediate/03_mapping_validation/logs/check_protein_id_mapping.log"
    shell:
        r"""
        set -euo pipefail
        python {params.script} \
            --mapping {input.mapping} \
            --eggnog {input.eggnog} \
            --overlap {output.overlap} \
            --preview {output.preview} \
            --summary {output.summary} \
            --log {log.stdout}
        """

################################################################################
# 04: SUBSET EXISTING EGGNOG ANNOTATIONS
################################################################################

rule viral_eggnog_subset:
    message:
        "[Viral annotation] Selecting eggNOG annotations on viral candidate parents."
    input:
        ids=rules.viral_candidate_table.output.ids,
        mapping=rules.transdecoder_protein_to_transcript.output.mapping,
        validation=rules.check_protein_id_mapping.output.summary,
        eggnog=f"{FUNCTION_DIR}/eggnog/eggnog.emapper.annotations"
    output:
        annotations=f"{VIRAL_ANALYSIS_DIR}/final_result/viral_eggnog_annotations.tsv"
    params:
        script=f"{VIRAL_SCRIPT_DIR}/04_eggnog_subset.py"
    log:
        stdout=f"{VIRAL_ANALYSIS_DIR}/final_result/logs/viral_eggnog_subset.log"
    shell:
        r"""
        set -euo pipefail
        python {params.script} \
            --ids {input.ids} \
            --mapping {input.mapping} \
            --eggnog {input.eggnog} \
            --output {output.annotations} \
            --log {log.stdout}
        """

################################################################################
# FINAL CANDIDATE SUMMARY
################################################################################

rule viral_annotation_results:
    message:
        "[Viral results] Creating the concise candidate summary."
    input:
        candidates=rules.viral_candidate_table.output.table,
        genes=rules.geNomad.output.viral_genes
    output:
        candidates=f"{VIRAL_ANALYSIS_DIR}/final_result/viral_candidates.tsv"
    params:
        script=f"{VIRAL_SCRIPT_DIR}/viral_results.py"
    log:
        stdout=f"{VIRAL_ANALYSIS_DIR}/final_result/logs/viral_candidates.log"
    shell:
        r"""
        set -euo pipefail
        python {params.script} annotations \
            --candidates {input.candidates} \
            --genomad-genes {input.genes} \
            --candidate-output {output.candidates} \
            --log {log.stdout}
        """
