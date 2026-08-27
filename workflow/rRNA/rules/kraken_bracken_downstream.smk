ENVS_DIR = Path(workflow.basedir)/"rRNA/envs"
SCRIPTS_DIR = Path(workflow.basedir)/"rRNA/scripts"

rule kraken_biom:
    conda:
        f"{ENVS_DIR}/kraken_biom.yaml"
    message:
        "[Kraken-Biom] combines kraken reports across samples against {wildcards.database}"
    input:
        reports=taxonomy_reports_for_database
    output:
        biom=f"{RESULTS_DIR}/rRNA/taxonomy/{{database}}/taxonomy.biom"
    log:
        f"{RESULTS_DIR}/rRNA/taxonomy/{{database}}/logs/kraken_biom.log"
    benchmark:
        f"{RESULTS_DIR}/rRNA/taxonomy/{{database}}/benchmarks/kraken_biom.txt"
    params:
        options=config["rRNA"]["kraken_biom"]["options"]
    shell:
        r"""
        mkdir -p $(dirname {output.biom})

        kraken-biom {input.reports} \
            -o {output.biom} \
            {params.options} \
            > {log} 2>&1
        """

#kraken-biom *.k2report -o bracken_genus.biom --min S --max D
rule biom_table_stats_qualitative:
    conda:
        f"{ENVS_DIR}/kraken_biom.yaml"
    message:
        "[biom summarize-table] qualitative summary statistics across samples against {wildcards.database}"
    input:
        biom=f"{RESULTS_DIR}/rRNA/taxonomy/{{database}}/taxonomy.biom"
    output:
        biom_qual=f"{RESULTS_DIR}/rRNA/taxonomy/{{database}}/taxonomy_qualitative.txt"
    log:
        f"{RESULTS_DIR}/rRNA/taxonomy/{{database}}/logs/biom_summary_qualitative.log"
    benchmark:
        f"{RESULTS_DIR}/rRNA/taxonomy/{{database}}/benchmarks/biom_summary_qualitative.txt"
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
        "[biom summarize-table] Observational summary statistics across samples against {wildcards.database}"
    input:
        biom=f"{RESULTS_DIR}/rRNA/taxonomy/{{database}}/taxonomy.biom"
    output:
        biom_obs=f"{RESULTS_DIR}/rRNA/taxonomy/{{database}}/taxonomy_observations.txt"
    log:
        f"{RESULTS_DIR}/rRNA/taxonomy/{{database}}/logs/biom_summary_observations.log"
    benchmark:
        f"{RESULTS_DIR}/rRNA/taxonomy/{{database}}/benchmarks/biom_summary_observations.txt"
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
        "[biom convert] converts {wildcards.database} database biom object to TSV format for taxonomic classification"
    input:
        biom=f"{RESULTS_DIR}/rRNA/taxonomy/{{database}}/taxonomy.biom"
    output:
        biom_tsv=f"{RESULTS_DIR}/rRNA/taxonomy/{{database}}/taxonomy.tsv"
    log:
        f"{RESULTS_DIR}/rRNA/taxonomy/{{database}}/logs/biom_convert.log"
    benchmark:
        f"{RESULTS_DIR}/rRNA/taxonomy/{{database}}/benchmarks/biom_convert.txt"
    shell:
        r"""
        mkdir -p $(dirname {output.biom_tsv})

        biom convert -i {input.biom} -o {output.biom_tsv} --to-tsv --header-key taxonomy > {log} 2>&1
        """

rule biom_to_phyloseq_raw:
    conda:
        f"{ENVS_DIR}/r_phyloseq.yaml"
    message:
        "[phyloseq] convert BIOM into an unfiltered phyloseq object from {wildcards.database}"
    input:
        biom=f"{RESULTS_DIR}/rRNA/taxonomy/{{database}}/taxonomy.biom"
    output:
        rds=f"{RESULTS_DIR}/rRNA/taxonomy/{{database}}/phyloseq_raw.rds"
    log:
        f"{RESULTS_DIR}/rRNA/taxonomy/{{database}}/logs/biom_to_phyloseq_raw.log"
    benchmark:
        f"{RESULTS_DIR}/rRNA/taxonomy/{{database}}/benchmarks/biom_to_phyloseq_raw.txt"
    shell:
        r"""
        mkdir -p $(dirname {output.rds})

        Rscript -e '
        suppressPackageStartupMessages({{
            library(biomformat)
            library(phyloseq)
        }})

        physeq <- import_biom("{input.biom}")

        saveRDS(
            physeq,
            file="{output.rds}",
            compress="xz"
        )
        ' > {log} 2>&1
        """

rule biom_to_phyloseq_filtered:
    conda:
        f"{ENVS_DIR}/r_phyloseq.yaml"
    message:
        "[phyloseq] convert BIOM into an filtered phyloseq object from {wildcards.database}, by removing multi-cellular organism "
    input:
        biom=f"{RESULTS_DIR}/rRNA/taxonomy/{{database}}/taxonomy.biom"
    output:
        rds=f"{RESULTS_DIR}/rRNA/taxonomy/{{database}}/phyloseq_filtered.rds"
    log:
        f"{RESULTS_DIR}/rRNA/taxonomy/{{database}}/logs/biom_to_phyloseq.log"
    benchmark:
        f"{RESULTS_DIR}/rRNA/taxonomy/{{database}}/benchmarks/biom_to_phyloseq.txt"
    params:
        prefix=f"{RESULTS_DIR}/rRNA/taxonomy/{{database}}/"
    shell:
        r"""
        Rscript {SCRIPTS_DIR}/biom_to_phyloseq.R {input.biom} {output.rds} {params.prefix} > {log} 2>&1
        """

#mamba create -n kraken_biom_test -c conda-forge -c bioconda python=3.10 kraken-biom biom-format