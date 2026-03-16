rule decompress_rrna:
    input:
        fwd="results/sortmerna/SSU/{sample}_fwd.fq.gz",
        rev="results/sortmerna/SSU/{sample}_rev.fq.gz",
    output:
        r1="results/metarib/data/{sample}.1.fq",
        r2="results/metarib/data/{sample}.2.fq",
    log:
        "logs/metarib/decompress/{sample}.log",
    conda:
        "../envs/pigz.yaml"
    threads: config["threads"]["pigz"]
    shell:
        "pigz -dkf -p{threads} < {input.fwd} > {output.r1} && "
        "echo 'Forward file was successfully decompressed' >> {log} && "
        "pigz -dkf -p{threads} < {input.rev} > {output.r2} && "
        "echo 'Reverse file was successfully decompressed' >> {log} "


emirge_cfg = config["metarib"]["EMIRGE"]
bbtool_cfg = config["metarib"]["BBTOOL"]
iteration_cfg = config["metarib"]["ITERATION"]


WORK_DIR = "results/metarib/work"
LOG_DIR = "logs/metarib"


def iter_dir(iter):
    return f"{WORK_DIR}/iter.{int(iter)}"


def prev_iter_dir(wildcards):
    return iter_dir(max(int(wildcards.iter) - 1, 0))


ITER_DIR = f"{WORK_DIR}/iter.{{iter}}"
ITER0_DIR = f"{WORK_DIR}/iter.0"


ruleorder: metarib_step0 > metarib_step


rule metarib_step0:
    input:
        r1=expand("results/metarib/data/{sample}.1.fq", sample=unique_samples),
        r2=expand("results/metarib/data/{sample}.2.fq", sample=unique_samples),
    output:
        r1=expand(f"{ITER0_DIR}/unmapped_data/{{sample}}.1.fq", sample=unique_samples),
        r2=expand(f"{ITER0_DIR}/unmapped_data/{{sample}}.2.fq", sample=unique_samples),
        contigs=f"{ITER0_DIR}/contigs_derep_next.fasta",
        report=f"{ITER0_DIR}/iteration_report.txt",
        done=f"{ITER0_DIR}/.done",
        output_dir=directory(ITER0_DIR),
        data_dir=directory(f"{ITER0_DIR}/unmapped_data"),
    log:
        f"{LOG_DIR}/iter0_init.log",
    run:
        import os

        os.makedirs(output.output_dir, exist_ok=True)
        os.makedirs(output.data_dir, exist_ok=True)

        for i in range(len(unique_samples)):
            os.symlink(os.path.abspath(input.r1[i]), output.r1[i])
            os.symlink(os.path.abspath(input.r2[i]), output.r2[i])

        with open(output.contigs, "w") as f:
            pass
        with open(output.report, "w") as f:
            pass
        with open(output.done, "w") as f:
            f.write("False")


rule metarib_step:
    input:
        r1=lambda wc: expand(
            f"{prev_iter_dir(wc)}/unmapped_data/{{sample}}.1.fq",
            sample=unique_samples,
        ),
        r2=lambda wc: expand(
            f"{prev_iter_dir(wc)}/unmapped_data/{{sample}}.2.fq",
            sample=unique_samples,
        ),
        contigs=lambda wc: f"{prev_iter_dir(wc)}/contigs_derep_next.fasta",
        prev_report=lambda wc: f"{prev_iter_dir(wc)}/iteration_report.txt",
    output:
        contigs_next=f"{ITER_DIR}/contigs_derep_next.fasta",
        report=f"{ITER_DIR}/iteration_report.txt",
    log:
        f"{LOG_DIR}/iter.{{iter}}_step.log",
    threads: config["threads"]["metarib"]
    conda:
        "../envs/metarib.yaml"
    params:
        nreads=emirge_cfg["SAMPLING_NUM"],
        ref_db=emirge_cfg["EM_REF"],
        bt_idx=emirge_cfg["EM_BT"],
        EM_PARA=emirge_cfg["EM_PARA"],
        MAP_PARA=bbtool_cfg["MAP_PARA"],
        CLS_PARA=bbtool_cfg["CLS_PARA"],
        data_dir=f"{ITER_DIR}",
    shell:
        """
        dirname=$(realpath "$(dirname "{input.r1[0]}")") && \
        source workflow/scripts/subsampling.sh --dirname "$dirname" --num_reads {params.nreads} >> {log} 2>&1

        mv "$dirname"/subsample_R1.fastq {params.data_dir}/subsample_R1.fastq
        mv "$dirname"/subsample_R2.fastq {params.data_dir}/subsample_R2.fastq

        workflow/scripts/metarib_step.sh \
            --r1 {input.r1} \
            --r2 {input.r2} \
            --contigs {input.contigs} \
            --output-dir {params.data_dir} \
            --num-reads {params.nreads} \
            --em-para "$EM_PARA" \
            --map-para '{params.MAP_PARA}' \
            --cls-para '{params.CLS_PARA}' \
            --ref-db {params.ref_db} \
            --bt-idx {params.bt_idx}  >> {log} 2>&1 
        """


checkpoint iteration_check:
    input:
        cur_report=f"{ITER_DIR}/iteration_report.txt",
        prev_report=lambda wc: f"{prev_iter_dir(wc)}/iteration_report.txt",
    output:
        f"{WORK_DIR}/iter.{{iter}}/.done",
    log:
        f"{LOG_DIR}/iter.{{iter}}_check.log",
    params:
        min_reads_threshold=iteration_cfg["MIN_READS_THRESHOLD"],
        convergence_threshold=iteration_cfg["CONVERGENCE_THRESHOLD"],
    threads: 1
    script:
        "../scripts/metarib_convergence.py"


max_iter = iteration_cfg["MAX_ITER"]


def run_iteration(wildcards):
    import warnings

    """
    This function is called after checkpoints run to determine final outputs.
    It walks through iterations until finding one that stopped.
    """
    iter = 1
    while True:
        ck = checkpoints.iteration_check.get(iter=iter)
        done_file = ck.output[0]
        with open(done_file, "r") as f:
            stop = f.read().strip()
        if stop == "True":
            return f"{WORK_DIR}/iter.{iter}/contigs_derep_next.fasta"
        iter += 1
        if iter > max_iter:
            warnings.warn(
                f"Exceeded maximum iterations ({max_iter}) without convergence. Forcing to stop. Check logs for details."
            )
            return f"{WORK_DIR}/iter.{iter}/contigs_derep_next.fasta"


rule metarib:
    input:
        f"{ITER0_DIR}/.done",
        run_iteration,
    output:
        "results/metarib/final_contigs.fasta",
    log:
        "logs/metarib/final_assembly.log",
    threads: 1
    shell:
        "touch {output} && "
        "cat {input[1]} > {output} && "
        "echo 'Final contigs copied to {output}' > {log} 2>&1 && "
        "echo 'All MetaRib reconstructions completed' >> {log} &&"
        "echo 'Final contigs for all samples are available in results/metarib/final_contigs/' >> {log}"
        "echo 'Cleaning up intermediate files...' >> {log} && "
        "rm -rf {WORK_DIR}"
