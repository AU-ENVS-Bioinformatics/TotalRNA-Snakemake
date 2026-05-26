rule bracken:
    conda:
        "../envs/bracken.yaml"
    message:
        "[Bracken] resetimates abundance of kraken results for {wildcards.sample} for taxonomic classification"
    input:
        kraken_report=f"{RESULTS_DIR}/rRNA/{{sample}}/classification/{{sample}}.k2report",
    output:
        bracken_output=f"{RESULTS_DIR}/rRNA/{{sample}}/classification/{{sample}}_{{level}}.bracken.tsv",
        bracken_kreport_output=f"{RESULTS_DIR}/rRNA/{{sample}}/classification/{{sample}}_{{level}}.bracken.k2report"
    log:
        stdout=f"{RESULTS_DIR}/rRNA/{{sample}}/logs/bracken.log"
    benchmark:
        f"{RESULTS_DIR}/rRNA/{{sample}}/benchmarks/bracken.txt"
    params:
        db=config["databases"]["kraken_rRNA_db"],
        options=config["rRNA"]["bracken"]["options"]
    wildcard_constraints:
        level="D|P|C|O|F|G|S",
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bracken_output})

        bracken \
          -d {params.db} \
          -i {input.kraken_report} \
          -o {output.bracken_output} \
          -w {output.bracken_kreport_output} \
          > {log.stdout} 2>&1
        """

#kraken-biom *.k2report -o bracken_genus.biom --min S --max D