rule humann:
    conda:
        "/data_2/rasmus/conda-envs/humann_env"
    message:
        "[humann] profiling the presence/absence and abundance of microbial pathways in a community for {wildcards.sample}"
    input:
        nonrRNA_concatenate=rules.concatenate.output.nonrRNA_concatenate,
        metaphlan_profile=rules.metaphlan.output.metaphlan_profile
    output:
        genefamilies=f"{RESULTS_DIR}/nonrRNA/{{sample}}/humann/{{sample}}_genefamilies.tsv",
        abundance=f"{RESULTS_DIR}/nonrRNA/{{sample}}/humann/{{sample}}_pathabundance.tsv",
        coverage=f"{RESULTS_DIR}/nonrRNA/{{sample}}/humann/{{sample}}_pathcoverage.tsv"
    log:
        stdout=f"{RESULTS_DIR}/nonrRNA/{{sample}}/logs/humann.log"
    benchmark:
        f"{RESULTS_DIR}/nonrRNA/{{sample}}/benchmarks/humann.txt"
    params:
        result_dir=f"{RESULTS_DIR}/nonrRNA/{{sample}}/humann/",
        nucleotide_db=config["databases"]["humann_nucleotidedb"],
        protein_db=config["databases"]["humann_proteindb"],
        basename=f"{{sample}}",
        options=config["nonrRNA"]["humann"]["options"]
    threads:
        config["nonrRNA"]["humann"].get("threads", 8)
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.genefamilies})

        humann -i {input.nonrRNA_concatenate} \
            -o {params.result_dir} \
            --threads {threads} \
            --output-basename {params.basename} \
            --taxonomic-profile {input.metaphlan_profile} \
            --nucleotide-database {params.nucleotide_db} \
            --protein-database {params.protein_db} \
            --prescreen-threshold 0 \
            --translated-query-coverage-threshold 50 \
            {params.options} \
            > {log.stdout} 2>&1
        """
#humann -i ANN11_concantenated.fastq.gz -o ../humann/ --threads 24 --taxonomic-profile test_profile.txt --nucleotide-database /data_2/Databases/humann/chocophlan/ --protein-database /data_2/Databases/humann/uniref/ --diamond-options=--fast