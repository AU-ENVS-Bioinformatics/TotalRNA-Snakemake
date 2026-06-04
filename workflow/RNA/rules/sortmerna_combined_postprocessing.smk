from pathlib import Path

SCRIPTS_DIR = Path(workflow.basedir)/"RNA/scripts"

if config["RNA"]["method"] == "sortmerna" and config["RNA"]["refinement"] == "combined":
    rule filter_sortmerna_aligned:
        conda:
            "../envs/rna_python_tools.yaml"
        message:
            "[Sortmerna_postprocessing_filter] for {wildcards.sample} filtering reads of percentage of the read length and at least certain number of nucleotide matches to the reference, and extract the ENA accession number for each read id"
        input:
            sam=rules.sortmerna_combined.output.sam
        output:
            tsv=f"{RESULTS_DIR}/RNA/{{sample}}/sortmerna/aligned/{{sample}}_ssu_lsu.aligned.filtered.tsv",
            sam=f"{RESULTS_DIR}/RNA/{{sample}}/sortmerna/aligned/{{sample}}_ssu_lsu.aligned.filtered.sam"
        log:
            stdout = f"{RESULTS_DIR}/RNA/{{sample}}/logs/sortmerna_postprocessing_filter.log"
        benchmark:
            f"{RESULTS_DIR}/RNA/{{sample}}/benchmarks/sortmerna_postprocessing_filter.txt"
        params:
            fraction=0.8,
            matches=10,
            prefix="SSU_"
        shell:
            """
            python {SCRIPTS_DIR}/sortmerna_read_id_filter.py \
                --sam_in {input.sam} \
                --fraction {params.fraction} \
                --matches {params.matches} \
                --prefix {params.prefix} \
                --read_out {output.tsv} \
                --sam_out {output.sam}
            """
            
    rule extract_sortmerna_accessions:
        conda:
            "../envs/rna_python_tools.yaml"
        message:
            "[Sortmerna_postprocessing_accession_taxonomy] for {wildcards.sample} filtering reads of percentage of the read length and at least certain number of nucleotide matches to the reference, and extract the ENA accession number for each read id"
        input:
            sam=rules.filter_sortmerna_aligned.output.sam,
            taxonomy_table=config["databases"][f"taxonomy_table"]
        output:
            tsv=f"{RESULTS_DIR}/RNA/{{sample}}/sortmerna/aligned/{{sample}}_ssu_lsu.aligned.filtered.accessions_taxonomy.tsv"
        log:
            stdout = f"{RESULTS_DIR}/RNA/{{sample}}/logs/sortmerna_postprocessing_accession_taxonomy.log"
        benchmark:
            f"{RESULTS_DIR}/RNA/{{sample}}/benchmarks/sortmerna_postprocessing_accession_taxonomy.txt"
        params:
            prefix="SSU_",
            email="au811884@uni.au.dk"
        shell:
            """
            python {SCRIPTS_DIR}/extract_read_ena_accesion.py \
                --sam {input.sam} \
                --out {output.tsv} \
                --prefix {params.prefix} \
                --email {params.email} \
                --taxonomy {input.taxonomy_table}
            """

