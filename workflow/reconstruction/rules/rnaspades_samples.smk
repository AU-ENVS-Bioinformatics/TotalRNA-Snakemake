################################################################################
# SAMPLE-SPECIFIC ASSEMBLY
################################################################################

rule rnaspades_sample_specific:
    conda:
        "../envs/rnaspades.yaml"
    wildcard_constraints:
        target=SAMPLE_TARGET_PATTERN
    message:
        "[rnaSPAdes] Reconstructing {wildcards.target} for sample {wildcards.sample}"
    input:
        r1=lambda wc: rnaspades_sample_input(wc, read=1, results_dir=RESULTS_DIR),
        r2=lambda wc: rnaspades_sample_input(wc, read=2, results_dir=RESULTS_DIR)
    output:
        transcripts=f"{ASSEMBLY_DIR}/cross_assembly/{{sample}}/{{target}}/rnaspades/transcripts.fasta"
    log:
        stdout=f"{ASSEMBLY_DIR}/cross_assembly/{{sample}}/{{target}}/logs/{{sample}}_{{target}}_rnaspades.log",
        time=f"{ASSEMBLY_DIR}/cross_assembly/{{sample}}/{{target}}/benchmarks/{{sample}}_{{target}}_rnaspades_time.txt"
    benchmark:
        f"{ASSEMBLY_DIR}/cross_assembly/{{sample}}/{{target}}/benchmarks/{{sample}}_{{target}}_rnaspades_benchmark.txt"
    params:
        outdir=f"{ASSEMBLY_DIR}/cross_assembly/{{sample}}/{{target}}/rnaspades",
        memory=config["reconstruction"]["rnaspades_sample"]["memory"]
    threads:
        config["reconstruction"]["rnaspades_sample"].get("threads", 10)
    shell:
        r"""
        set -euo pipefail

        mkdir -p {params.outdir} $(dirname {log.stdout}) $(dirname {log.time})

        /usr/bin/time \
            -v \
            -o {log.time} \
            rnaspades.py \
                -1 {input.r1} \
                -2 {input.r2} \
                -t {threads} \
                -o {params.outdir} \
                --memory {params.memory} \
                > {log.stdout} 2>&1
        
        if [ ! -s {output.transcripts} ]; then
            fallback=""
            while IFS= read -r candidate; do
                if [ -s "$candidate" ]; then
                    fallback="$candidate"
                    break
                fi
            done < <(printf '%s\n' {params.outdir}/K[0-9]*/final_contigs.fasta | sort -Vr)
            if [ -z "$fallback" ]; then
                echo "[ERROR] No nonempty transcripts or fallback contigs found." >> {log.stdout}
                exit 1
            fi
            echo "[WARNING] No transcripts produced; linking unresolved contigs from $fallback." >> {log.stdout}
            ln -sf "$(realpath "$fallback")" {output.transcripts}
        fi
        """