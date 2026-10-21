################################################################################
# 2. CROSS-SAMPLE NON-rRNA CO-ASSEMBLY
################################################################################

rule rnaspades_cross_sample:
    conda:
        "../envs/rnaspades.yaml"
    wildcard_constraints:
        target=CROSS_SAMPLE_TARGET_PATTERN
    message:
        "[rnaSPAdes] Co-assembling non-rRNA reads across all samples"
    input:
        r1=f"{ASSEMBLY_DIR}/coassembly/{{target}}/concatenated/nonrRNA_1.fastq.gz",
        r2=f"{ASSEMBLY_DIR}/coassembly/{{target}}/concatenated/nonrRNA_2.fastq.gz"
    output:
        transcripts=f"{ASSEMBLY_DIR}/coassembly/{{target}}/rnaspades/transcripts.fasta"
    log:
        stdout=f"{ASSEMBLY_DIR}/coassembly/{{target}}/logs/nonrRNA_rnaspades.log",
        time=f"{ASSEMBLY_DIR}/coassembly/{{target}}/benchmarks/nonrRNA_rnaspades_time.txt"
    benchmark:
        f"{ASSEMBLY_DIR}/coassembly/{{target}}/benchmarks/nonrRNA_rnaspades_benchmark.txt"
    params:
        outdir=f"{ASSEMBLY_DIR}/coassembly/{{target}}/rnaspades",
        memory=config["reconstruction"]["rnaspades"]["memory"],
        max_memory_gb=250
    threads:
        config["reconstruction"]["rnaspades"].get("threads", 12)
    shell:
        r"""
        set -euo pipefail

        mkdir -p {params.outdir:q}

        read -r memory_flag memory_gb extra <<< {params.memory:q}
        max_memory_gb={params.max_memory_gb:q}

        attempt_log=$(mktemp)
        trap 'rm -f "$attempt_log"' EXIT

        attempt=0
        restart=false

        # Resume an existing SPAdes run when its saved parameters exist.
        if [[ -f {params.outdir:q}/params.txt ]]; then
            restart=true
        fi

        while true; do
            attempt=$((attempt + 1))

            # Cap both the configured memory and subsequent increases.
            if (( memory_gb > max_memory_gb )); then
                memory_gb=$max_memory_gb
            fi

            echo "[rnaSPAdes] Attempt $attempt: memory=$memory_gb GB; restart=$restart." \
                | tee -a {log.stdout:q}

            if [[ "$restart" == true ]]; then
                echo "[rnaSPAdes] Running rnaspades.py --restart-from last --memory $memory_gb." \
                    | tee -a {log.stdout:q}

                if /usr/bin/time -v -a -o {log.time:q} \
                    rnaspades.py \
                        -t {threads} \
                        -o {params.outdir:q} \
                        --restart-from last \
                        --memory "$memory_gb" \
                        > "$attempt_log" 2>&1
                then
                    rc=0
                else
                    rc=$?
                fi
            else
                if /usr/bin/time -v -a -o {log.time:q} \
                    rnaspades.py \
                        -1 {input.r1:q} \
                        -2 {input.r2:q} \
                        -t {threads} \
                        -o {params.outdir:q} \
                        --memory "$memory_gb" \
                        > "$attempt_log" 2>&1
                then
                    rc=0
                else
                    rc=$?
                fi
            fi

            # Append this attempt to the combined log.
            cat "$attempt_log" >> {log.stdout:q}

            if (( rc == 0 )); then
                echo "[rnaSPAdes] Assembly completed successfully with $memory_gb GB." \
                    | tee -a {log.stdout:q}
                break
            fi

            # Check only the latest attempt for retryable errors caused by memory allocation.
            if ! grep -Eiq -- \
                'SIGABRT|SIGABORT|Aborted|OS return value:[[:space:]]*-6([^0-9]|$)|error code:[[:space:]]*-6([^0-9]|$)|error code:[[:space:]]*12([^0-9]|$)|Cannot allocate memory|std::bad_alloc|out of memory' \
                "$attempt_log"
            then
                exit "$rc"
            fi

            if (( memory_gb >= max_memory_gb )); then
                echo "[rnaSPAdes] Failed at the configured maximum of $max_memory_gb GB; stopping." \
                    | tee -a {log.stdout:q}
                exit "$rc"
            fi

            # Increase memory by 1.5 times, rounded upward.
            next_memory_gb=$(( (memory_gb * 3 + 1) / 2 ))

            if (( next_memory_gb > max_memory_gb )); then
                next_memory_gb=$max_memory_gb
            fi
            
            echo "[rnaSPAdes] Retryable failure on attempt $attempt; retrying with --restart-from last --memory $next_memory_gb (previously $memory_gb GB)." \
                | tee -a {log.stdout:q}

            memory_gb=$next_memory_gb
            restart=true
        done
        """