SCRIPTS_DIR = Path(workflow.basedir) / "taxonomy/scripts"

if its_taxonomy_enabled(config):

    ITS_SUMMARY_RANK = config["taxonomy"]["ITS"].get("summary", {}).get(
        "taxonomic_rank", "genus"
    )

    if ITS_SUMMARY_RANK not in {"kingdom", "phylum", "class", "order", "family", "genus", "species"}:
        raise ValueError(
            f"Unsupported ITS summary rank: {ITS_SUMMARY_RANK}"
        )
    
    ################################################################################
    # POST-PROCESS ITS BLASTN RESULTS
    ################################################################################
    
    rule postprocess_its_blast:
        conda:
            "../envs/taxonomy_postprocessing.yaml"
        message:
            "[ITS taxonomy] Filtering {wildcards.region} UNITE hits for {wildcards.sample}"
        input:
            blast=f"{TAXONOMY_DIR}/{{sample}}/ITS/{{region}}/blastn/{{sample}}.{{region}}.UNITE.blastn.tsv",
            query=f"{ASSEMBLY_DIR}/{{sample}}/ITS_candidates/ITSx/{{sample}}_ITSx.{{region}}.fasta"
        output:
            assignments=f"{TAXONOMY_DIR}/{{sample}}/ITS/{{region}}/taxonomy/{{sample}}.{{region}}.UNITE.assignments.tsv"
        log:
            stdout=f"{TAXONOMY_DIR}/{{sample}}/ITS/{{region}}/taxonomy/logs/{{sample}}.{{region}}.UNITE.postprocessing.log"
        benchmark:
            f"{TAXONOMY_DIR}/{{sample}}/ITS/{{region}}/taxonomy/benchmarks/{{sample}}.{{region}}.UNITE.postprocessing.txt"
        params:
            script=workflow.source_path("../scripts/postprocess_unite_blast.py"),
            min_identity=config["taxonomy"]["ITS"]["postprocessing"].get("min_identity", 80),
            min_query_coverage=config["taxonomy"]["ITS"]["postprocessing"].get("min_query_coverage", 80),
            top_bitscore_fraction=config["taxonomy"]["ITS"]["postprocessing"].get("top_bitscore_fraction", 0.99),
            outdir=f"{TAXONOMY_DIR}/{{sample}}/ITS/{{region}}/taxonomy",
            logdir=f"{TAXONOMY_DIR}/{{sample}}/ITS/{{region}}/taxonomy/logs",
            benchmarkdir=f"{TAXONOMY_DIR}/{{sample}}/ITS/{{region}}/taxonomy/benchmarks"
        shell:
            r"""
            set -euo pipefail

            mkdir -p $(dirname {output.assignments})

            python {params.script} \
                --blast {input.blast} \
                --query-fasta {input.query} \
                --output {output.assignments} \
                --sample {wildcards.sample} \
                --region {wildcards.region} \
                --min-identity {params.min_identity} \
                --min-query-coverage {params.min_query_coverage} \
                --top-bitscore-fraction {params.top_bitscore_fraction} \
                > {log.stdout} 2>&1
            """

    ################################################################################
    # COMBINE ITS REGION ASSIGNMENTS PER SAMPLE
    ################################################################################

    rule combine_its_assignments:
        conda:
            "../envs/taxonomy_postprocessing.yaml"
        message:
            "[ITS taxonomy] Combining ITS region assignments for {wildcards.sample}"
        input:
            assembly=f"{ASSEMBLY_DIR}/{{sample}}/ITS_candidates/rnaspades/transcripts.fasta",
            assignments=lambda wc: [
                f"{TAXONOMY_DIR}/{wc.sample}/ITS/{region}/taxonomy/{wc.sample}.{region}.UNITE.assignments.tsv"
                for region in its_blast_regions(config)
            ]
        output:
            assigned=f"{TAXONOMY_DIR}/{{sample}}/ITS/combined/{{sample}}.ITS.assignments.tsv",
            unassigned=f"{TAXONOMY_DIR}/{{sample}}/ITS/combined/{{sample}}.ITS.unassignments.tsv",
            all_assignments=f"{TAXONOMY_DIR}/{{sample}}/ITS/combined/{{sample}}.ITS.allasignments.tsv"
        log:
            stdout=f"{TAXONOMY_DIR}/{{sample}}/ITS/combined/logs/{{sample}}.ITS.combine.log"
        benchmark:
            f"{TAXONOMY_DIR}/{{sample}}/ITS/combined/benchmarks/{{sample}}.ITS.combine.txt"
        shell:
            r"""
            set -euo pipefail

            mkdir -p $(dirname {output.assigned})

            python {SCRIPTS_DIR}/combine_its_assignments.py \
                --assembly-fasta {input.assembly} \
                --assignments {input.assignments} \
                --assigned-output {output.assigned} \
                --unassigned-output {output.unassigned} \
                --all-output {output.all_assignments} \
                --sample {wildcards.sample} \
                > {log.stdout} 2>&1
            """

    ################################################################################
    # COMBINE ITS ASSIGNMENTS ACROSS SAMPLES
    ################################################################################

    rule create_its_summary_matrices:
        conda:
            "../envs/taxonomy_postprocessing.yaml"
        message:
            "[ITS taxonomy] Creating cross-sample ITS count matrices"
        input:
            tables=expand(
                f"{TAXONOMY_DIR}/{{sample}}/ITS/combined/{{sample}}.ITS.allasignments.tsv",
                sample=SAMPLES,
            )
        output:
            assigned=f"{TAXONOMY_DIR}/summary/ITS/matrices/ITS.assignments.{ITS_SUMMARY_RANK}.counts.tsv",
            unassigned=f"{TAXONOMY_DIR}/summary/ITS/matrices/ITS.unassignments.reasons.counts.tsv",
            statuses=f"{TAXONOMY_DIR}/summary/ITS/matrices/ITS.allassignments.status.counts.tsv"
        log:
            stdout=f"{TAXONOMY_DIR}/summary/ITS/matrices/logs/ITS.summary_matrices.log"
        benchmark:
            f"{TAXONOMY_DIR}/summary/ITS/matrices/benchmarks/ITS.summary_matrices.txt"
        params:
            script=workflow.source_path("../scripts/create_its_matrices.py"),
            samples=SAMPLES,
            rank=ITS_SUMMARY_RANK,
            outdir=f"{TAXONOMY_DIR}/summary/ITS/matrices",
            logdir=f"{TAXONOMY_DIR}/summary/ITS/matrices/logs",
            benchmarkdir=f"{TAXONOMY_DIR}/summary/ITS/matrices/benchmarks"
        shell:
            r"""
            set -euo pipefail

            mkdir -p $(dirname {output.assigned})

            python {params.script} \
                --tables {input.tables} \
                --samples {params.samples} \
                --rank {params.rank} \
                --assigned-output {output.assigned} \
                --unassigned-output {output.unassigned} \
                --status-output {output.statuses} \
                > {log.stdout} 2>&1
            """
