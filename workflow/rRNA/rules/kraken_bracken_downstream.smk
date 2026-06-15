ENVS_DIR = Path(workflow.basedir)/"rRNA/envs"
SCRIPTS_DIR = Path(workflow.basedir)/"rRNA/scripts"

rule kraken_biom:
    conda:
        f"{ENVS_DIR}/kraken_biom.yaml"
    message:
        "[Kraken-Biom] converts kraken report to BIOM format for for taxonomic classification"
    input:
        expand(
            f"{RESULTS_DIR}/rRNA/{{sample}}/classification/{{sample}}.k2report",sample=SAMPLES
        )
    output:
        biom=f"{RESULTS_DIR}/rRNA/taxonomy/bracken_species.biom"
    log:
        f"{RESULTS_DIR}/rRNA/taxonomy/logs/kraken_biom.log"
    benchmark:
        f"{RESULTS_DIR}/rRNA/taxonomy/benchmarks/kraken_biom.txt"
    params:
        options=config["rRNA"]["kraken_biom"]["options"]
    shell:
        r"""
        mkdir -p $(dirname {output.biom})

        kraken-biom {input} \
            -o {output.biom} \
            {params.options} \
            > {log} 2>&1
        """

#kraken-biom *.k2report -o bracken_genus.biom --min S --max D
rule biom_table_stats_qualitative:
    conda:
        f"{ENVS_DIR}/kraken_biom.yaml"
    message:
        "[biom summarize-table] Qualitative summary statistics for taxonomic classification"
    input:
        biom = rules.kraken_biom.output.biom
    output:
        biom_qual=f"{RESULTS_DIR}/rRNA/taxonomy/bracken_species_qualitative.txt"
    log:
        f"{RESULTS_DIR}/rRNA/taxonomy/logs/biom_summary_qualitative.log"
    benchmark:
        f"{RESULTS_DIR}/rRNA/taxonomy/benchmarks/biom_summary_qualitative.txt"
    shell:
        r"""
        mkdir -p $(dirname {output.biom_qual})

        biom summarize-table -i {input.biom} \
            -o {output.biom_qual} \
            --qualitative \
            > {log} 2>&1
        """

rule biom_table_stats_observations:
    conda:
        f"{ENVS_DIR}/kraken_biom.yaml"
    message:
        "[biom summarize-table] Observational summary statistics for taxonomic classification"
    input:
        biom = rules.kraken_biom.output.biom
    output:
        biom_obs=f"{RESULTS_DIR}/rRNA/taxonomy/bracken_species_observations.txt"
    log:
        f"{RESULTS_DIR}/rRNA/taxonomy/logs/biom_summary_observations.log"
    benchmark:
        f"{RESULTS_DIR}/rRNA/taxonomy/benchmarks/biom_summary_observations.txt"
    shell:
        r"""
        mkdir -p $(dirname {output.biom_obs})

        biom summarize-table -i {input.biom} \
            -o {output.biom_obs} \
            --observations \
            > {log} 2>&1
        """

rule kraken_biom_convert:
    conda:
        f"{ENVS_DIR}/kraken_biom.yaml"
    message:
        "[biom convert] converts kraken biom object to TSV format for taxonomic classification"
    input:
        biom = rules.kraken_biom.output.biom
    output:
        biom_tsv=f"{RESULTS_DIR}/rRNA/taxonomy/bracken_species.tsv"
    log:
        f"{RESULTS_DIR}/rRNA/taxonomy/logs/biom_convert.log"
    benchmark:
        f"{RESULTS_DIR}/rRNA/taxonomy/benchmarks/biom_convert.txt"
    shell:
        r"""
        mkdir -p $(dirname {output.biom_tsv})

        biom convert -i {input.biom} -o {output.biom_tsv} --to-tsv --header-key taxonomy > {log} 2>&1
        """

rule biom_to_phyloseq:
    conda:
        f"{ENVS_DIR}/r_phyloseq.yaml"
    input:
        biom = rules.kraken_biom.output.biom
    output:
        rds = f"{RESULTS_DIR}/rRNA/taxonomy/phyloseq.rds"
    log:
        f"{RESULTS_DIR}/rRNA/taxonomy/logs/biom_to_phyloseq.log"
    benchmark:
        f"{RESULTS_DIR}/rRNA/taxonomy/benchmarks/biom_to_phyloseq.txt"
    params:
        prefix=f"{RESULTS_DIR}/rRNA/taxonomy/"
    shell:
        r"""
        Rscript {SCRIPTS_DIR}/biom_to_phyloseq.R {input.biom} {output.rds} {params.prefix} > {log} 2>&1
        """

#mamba create -n kraken_biom_test     -c conda-forge     -c bioconda     python=3.10     kraken-biom     biom-format