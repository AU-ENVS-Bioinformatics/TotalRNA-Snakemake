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


rule iteration_0:
    input:
        fwd="results/metarib/data/all.1.fq",
        rev="results/metarib/data/all.2.fq",
    output:
        folder=directory(ITER0_DIR),
        r1=f"{ITER0_DIR}/unmapped_R1_next.fastq",
        r2=f"{ITER0_DIR}/unmapped_R2_next.fastq",
        contigs=f"{ITER0_DIR}/contigs_derep_next.fasta",
        report=f"{ITER0_DIR}/iteration_report.txt",
        done=f"{ITER0_DIR}/.done",
    log:
        f"{LOG_DIR}/iter0_init.log",
    shell:
        """
        mkdir -p {output.folder} && 
        ln -sf $(pwd)/{input.fwd} {output.r1} && 
        ln -sf $(pwd)/{input.rev} {output.r2} &&
        touch {output.contigs} &&
        echo "0" >> {output.report} &&
        r1_reads=$(( $(wc -l < {output.r1}) / 4 )) 
        echo "$r1_reads" >> {output.report} && 
        echo "False" > {output.done} >> {log} 2>&1
        """


rule metarib_step:
    input:
        r1=lambda wc: f"{prev_iter_dir(wc)}/unmapped_R1_next.fastq",
        r2=lambda wc: f"{prev_iter_dir(wc)}/unmapped_R2_next.fastq",
        contigs=lambda wc: f"{prev_iter_dir(wc)}/contigs_derep_next.fasta",
        prev_report=lambda wc: f"{prev_iter_dir(wc)}/iteration_report.txt",
    output:
        contigs_next=f"{ITER_DIR}/contigs_derep_next.fasta",
        r1_next=f"{ITER_DIR}/unmapped_R1_next.fastq",
        r2_next=f"{ITER_DIR}/unmapped_R2_next.fastq",
        report=f"{ITER_DIR}/iteration_report.txt",
        dups=f"{ITER_DIR}/contigs.duplicates.fasta",
    log:
        f"{LOG_DIR}/iter.{{iter}}_step.log",
    # shadow:
    #     "full"
    threads: 8
    conda:
        "../envs/emirge.yaml"
    params:
        n=emirge_cfg["SAMPLING_NUM"],
        ref_db=emirge_cfg["EM_REF"],
        bt_idx=emirge_cfg["EM_BT"],
        EM_PARA=emirge_cfg["EM_PARA"],
        MAP_PARA=bbtool_cfg["MAP_PARA"],
        CLS_PARA=bbtool_cfg["CLS_PARA"],
        output_dir=ITER_DIR,
    shell:
        "workflow/scripts/metarib_step.sh --r1 {input.r1} --r2 {input.r2} --contigs {input.contigs} --output-dir {params.output_dir} --num-reads {params.n} --em-para '{params.EM_PARA}' --map-para '{params.MAP_PARA}' --cls-para '{params.CLS_PARA}' --ref-db {params.ref_db} --bt-idx {params.bt_idx} >> {log} 2>&1"


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
        run_iteration,
    output:
        "results/metarib/final_contigs.fasta",
    log:
        "logs/metarib/final_assembly.log",
    threads: 1
    shell:
        "touch {output} && "
        "cat {input} > {output} && "
        "echo 'Final contigs copied to {output}' >> {log} 2>&1 && "
        "echo 'All MetaRib reconstructions completed' >> {log} &&"
        "echo 'Final contigs for all samples are available in results/metarib/final_contigs/' >> {log} &&"
        "echo 'Cleaning up intermediate files...' >> {log} && "
        "rm -rf {WORK_DIR}"
