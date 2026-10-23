################################################################################
# NON-rRNA: REMOVE CONTIGS WITH >90% ncRNA COVERAGE
################################################################################

rule coassembly_filter_ncrna:
    message:
        "[ncRNA filtering] Filtering non-rRNA contigs by ncRNA coverage"
    input:
        fasta=f"{ASSEMBLY_DIR}/coassembly/non_rRNA/postprocessing/non_rRNA_dereplicated.fasta",
        tblout=f"{ASSEMBLY_DIR}/coassembly/non_rRNA/postprocessing/non_rRNA_rfam_hits.tbl"
    output:
        fasta=f"{ASSEMBLY_DIR}/coassembly/non_rRNA/final/transcript.fasta",
        removed=f"{ASSEMBLY_DIR}/coassembly/non_rRNA/postprocessing/non_rRNA_removed_ncrna.fasta",
        removed_report=f"{ASSEMBLY_DIR}/coassembly/non_rRNA/postprocessing/non_rRNA_removed_ncrna_report.tsv"
    params:
        script=f"{workflow.basedir}/reconstruction/scripts/filter_ncrna.py",
        coverage=config["reconstruction"]["coassembly_postprocessing"]["non_rRNA"].get("ncrna_min_contig_coverage", 0.90)
    log:
        stdout=f"{ASSEMBLY_DIR}/coassembly/non_rRNA/logs/non_rRNA_filter_ncrna.log"
    threads:
        1
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.fasta}) $(dirname {output.removed}) $(dirname {log.stdout})

        python3 {params.script} \
            --fasta {input.fasta} \
            --tblout {input.tblout} \
            --output {output.fasta} \
            --removed_out {output.removed} \
            --removed_report {output.removed_report} \
            --max-ncrna-coverage {params.coverage} \
            > {log.stdout} 2>&1
        """