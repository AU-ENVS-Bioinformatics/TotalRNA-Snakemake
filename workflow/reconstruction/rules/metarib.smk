################################################################################
# 1. DECOMPRESS SAMPLE READS
################################################################################

rule metarib_sample_reads:
    wildcard_constraints:
        target=COASSEMBLY_TARGET_PATTERNS["metarib"],
        sample="(?:" + "|".join(re.escape(sample) for sample in SAMPLES) + ")",
        read="1|2"
    input:
        reads=f"{RNA_CLASSIFIED_DIR}/{{sample}}/{{target}}/{{sample}}_{{target}}_{{read}}.fastq.gz"
    output:
        reads=f"{ASSEMBLY_DIR}/coassembly/{{target}}/metarib/data/{{sample}}.{{read}}.fq"
    threads:
        1
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname "{output.reads}")"

        gzip -dc "{input.reads}" > "{output.reads}"
        """


################################################################################
# 2. CONCATENATE READS AND WRITE SAMPLE LIST
################################################################################

rule metarib_prepare_reads:
    wildcard_constraints:
        target=COASSEMBLY_TARGET_PATTERNS["metarib"]
    message:
        "[MetaRib] Preparing pooled {wildcards.target} reads"
    input:
        r1=lambda wc: expand(
            f"{ASSEMBLY_DIR}/coassembly/{wc.target}/metarib/data/{{sample}}.1.fq",
            sample=SAMPLES,
        ),
        r2=lambda wc: expand(
            f"{ASSEMBLY_DIR}/coassembly/{wc.target}/metarib/data/{{sample}}.2.fq",
            sample=SAMPLES,
        )
    output:
        r1=f"{ASSEMBLY_DIR}/coassembly/{{target}}/metarib/data/all.1.fq",
        r2=f"{ASSEMBLY_DIR}/coassembly/{{target}}/metarib/data/all.2.fq",
        samples=f"{ASSEMBLY_DIR}/coassembly/{{target}}/metarib/data/samples.list.txt"
    threads:
        1
    run:
        import shutil

        Path(output.r1).parent.mkdir(parents=True, exist_ok=True)

        for sources, destination in [
            (input.r1, output.r1),
            (input.r2, output.r2),
        ]:
            with open(destination, "wb") as merged:
                for source in sources:
                    with open(source, "rb") as reads:
                        shutil.copyfileobj(reads, merged)

        with open(output.samples, "w") as handle:
            for sample in SAMPLES:
                handle.write(f"{sample}\n")

################################################################################
# 3. WRITE METARIB CONFIGURATION
################################################################################

rule metarib_prepare_config:
    conda:
        "../envs/metarib.yaml"
    wildcard_constraints:
        target=COASSEMBLY_TARGET_PATTERNS["metarib"]
    message:
        "[MetaRib] Creating configuration for {wildcards.target}"
    input:
        reference=config["reconstruction"]["metarib"]["EM_REF"]
    output:
        cfg=f"{ASSEMBLY_DIR}/coassembly/{{target}}/metarib/MetaRib.cfg"
    params:
        project=f"{ASSEMBLY_DIR}/coassembly/{{target}}/metarib",
        sampling_num=config["reconstruction"]["metarib"].get("SAMPLING_NUM", 1000000),
        em_para=config["reconstruction"]["metarib"]["EM_PARA"],
        em_ref=config["reconstruction"]["metarib"]["EM_REF"],
        em_bt=config["reconstruction"]["metarib"]["EM_BT"],
        map_para=config["reconstruction"]["metarib"]["MAP_PARA"],
        cls_para=config["reconstruction"]["metarib"]["CLS_PARA"],
        min_cov=config["reconstruction"]["metarib"].get("MIN_COV", 2),
        min_per=config["reconstruction"]["metarib"].get("MIN_PER", 80)
    threads:
        config["reconstruction"]["metarib"].get("threads", 12)
    log:
        stdout=f"{ASSEMBLY_DIR}/coassembly/{{target}}/logs/metarib/config_file.log"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "{params.project}" "$(dirname "{log.stdout}")"

        project=$(realpath "{params.project}")
        em_ref=$(realpath "{params.em_ref}")
        em_bt=$(realpath -m "{params.em_bt}")
        em_path=$(command -v emirge_amplicon.py)
        bb_path="$CONDA_PREFIX/bin"

        (
            printf '%s\n' '[BASE]'
            printf 'PROJECT_DIR : %s\n' "$project"
            printf 'DATA_DIR : %s/data\n' "$project"
            printf 'SAMPLING_NUM : %s\n' "{params.sampling_num}"
            printf 'THREAD : %s\n' "{threads}"

            printf '%s\n' '[EMIRGE]'
            printf 'EM_PATH : %s\n' "$em_path"
            printf 'EM_PARA : %s\n' "{params.em_para}"
            printf 'EM_REF : %s\n' "$em_ref"
            printf 'EM_BT : %s\n' "$em_bt"

            printf '%s\n' '[BBTOOL]'
            printf 'BB_PATH : %s\n' "$bb_path"
            printf 'MAP_PARA : %s\n' "{params.map_para}"
            printf 'CLS_PARA : %s\n' "{params.cls_para}"

            printf '%s\n' '[FILTER]'
            printf 'MIN_COV : %s\n' "{params.min_cov}"
            printf 'MIN_PER : %s\n' "{params.min_per}"
        ) > "{output.cfg}"

        cp "{output.cfg}" "{log.stdout}"
        """

################################################################################
# 4. RUN METARIB AND PUBLISH THE ASSEMBLY
################################################################################

rule metarib_coassembly:
    conda:
        "../envs/metarib.yaml"
    wildcard_constraints:
        target=COASSEMBLY_TARGET_PATTERNS["metarib"]
    message:
        "[MetaRib] Co-assembling {wildcards.target} reads across all samples"
    input:
        r1=f"{ASSEMBLY_DIR}/coassembly/{{target}}/metarib/data/all.1.fq",
        r2=f"{ASSEMBLY_DIR}/coassembly/{{target}}/metarib/data/all.2.fq",
        samples=f"{ASSEMBLY_DIR}/coassembly/{{target}}/metarib/data/samples.list.txt",
        cfg=f"{ASSEMBLY_DIR}/coassembly/{{target}}/metarib/MetaRib.cfg",
        r1_samples=lambda wc: expand(
            f"{ASSEMBLY_DIR}/coassembly/{wc.target}/metarib/data/{{sample}}.1.fq",
            sample=SAMPLES,
        ),
        r2_samples=lambda wc: expand(
            f"{ASSEMBLY_DIR}/coassembly/{wc.target}/metarib/data/{{sample}}.2.fq",
            sample=SAMPLES,
        ),
        script=f"{workflow.basedir}/reconstruction/scripts/run_MetaRib.py"
    output:
        native_fasta=f"{ASSEMBLY_DIR}/coassembly/{{target}}/metarib/MetaRib/Abundance/all.dedup.filtered.fasta",
        link_transcript=f"{ASSEMBLY_DIR}/coassembly/{{target}}/{{target}}_transcripts.fasta"
    threads:
        config["reconstruction"]["metarib"].get("threads", 12)
    log:
        stdout=f"{ASSEMBLY_DIR}/coassembly/{{target}}/logs/metarib/metarib.log",
        time=f"{ASSEMBLY_DIR}/coassembly/{{target}}/benchmarks/metarib_time.txt"
    benchmark:
        f"{ASSEMBLY_DIR}/coassembly/{{target}}/benchmarks/metarib_benchmark.txt"
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.link_transcript})

        /usr/bin/time -v -o {log.time} \
            python2 {input.script} \
                -cfg {input.cfg} \
                > {log.stdout} 2>&1

        ln -sfn $(realpath {output.native_fasta}) {output.link_transcript}
        """