################################################################################
# PREDICT rRNA REGIONS IN CONCATENATED SSU ASSEMBLIES
################################################################################

rule barrnap_ssu:
    conda:
        "../envs/barrnap.yaml"
    wildcard_constraints:
        kingdom=BARRNAP_KINGDOM_PATTERN
    message:
        "[Barrnap] Predicting rRNA regions using the {wildcards.kingdom} database"
    input:
        fasta=f"{ASSEMBLY_DIR}/concatenated/SSU/cross_sample_SSU.fasta"
    output:
        gff=f"{ASSEMBLY_DIR}/SSU_candidates/barrnap/{{kingdom}}/cross_sample_SSU_{{kingdom}}.barrnap.gff",
        fasta=f"{ASSEMBLY_DIR}/SSU_candidates/barrnap/{{kingdom}}/cross_sample_SSU_{{kingdom}}.barrnap.fasta"
    log:
        stdout=f"{ASSEMBLY_DIR}/SSU_candidates/barrnap/logs/cross_sample_SSU_{{kingdom}}.barrnap.log",
        time=f"{ASSEMBLY_DIR}/SSU_candidates/barrnap/benchmarks/cross_sample_SSU_{{kingdom}}.barrnap_time.txt"
    benchmark:
        f"{ASSEMBLY_DIR}/SSU_candidates/barrnap/benchmarks/cross_sample_SSU_{{kingdom}}.barrnap_benchmark.tsv"
    params:
        options=config["reconstruction"]["barrnap"].get("options", "")
    threads:
        config["reconstruction"]["barrnap"].get("threads", 12)
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.gff})

        /usr/bin/time \
            -v \
            -o {log.time} \
            barrnap \
                --kingdom {wildcards.kingdom} \
                --threads {threads} \
                {params.options} \
                {input.fasta} \
                --outseq {output.fasta} \
                > {output.gff} \
                2> {log.stdout}
        """

################################################################################
# CHECK ALL KINGDOMS AFTER ALL BARRNAP RUNS HAVE COMPLETED
################################################################################

checkpoint check_barrnap_ssu:
    message:
        "[Barrnap] Checking which kingdoms have predicted sequences"
    input:
        fasta=expand(f"{ASSEMBLY_DIR}/SSU_candidates/barrnap/{{kingdom}}/cross_sample_SSU_{{kingdom}}.barrnap.fasta", kingdom=BARRNAP_KINGDOMS)
    output:
        status=f"{ASSEMBLY_DIR}/SSU_candidates/barrnap/barrnap_kingdom_status.tsv"
    run:
        from pathlib import Path

        Path(output.status).parent.mkdir(parents=True, exist_ok=True)

        with open(output.status, "w") as report:
            report.write("kingdom\tstatus\n")

            for kingdom, fasta_path in zip(BARRNAP_KINGDOMS, input.fasta):
                status = "present" if Path(fasta_path).stat().st_size > 0 else "empty"
                report.write(f"{kingdom}\t{status}\n")