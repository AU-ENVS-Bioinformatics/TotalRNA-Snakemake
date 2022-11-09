import subprocess
from collections import ChainMap

RRNA_FILEPATH = config.get("RRNA_FILEPATH", "rrna/")
METARIB_FILEPATH = config.get("METARIB_FILEPATH", "MetaRib/")
AVAILABLE_THREADS = int(workflow.cores * 0.75)
metarib_params = dict(ChainMap(*config.get("metarib")))
private_metarib_params = dict(
    ChainMap(
        *config.get("PRIVATE_METARIB"),
    )
)


rule prepare_metarib:
    input:
        R1=expand(
            f"results/rrna/{{sample}}_fwd.fq.gz",
            sample=unique_samples,
        ),
        R2=expand(
            f"results/rrna/{{sample}}_rev.fq.gz",
            sample=unique_samples,
        ),
    output:
        R1=f"results/MetaRib/data/all.1.fq",
        R2=f"results/MetaRib/data/all.2.fq",
        sample_list=report(
            f"results/MetaRib/data/samples.list.txt",
            caption="report/prepare_metarib.rst",
        ),
    log:
        "logs/metarib/prepare_metarib.log",
    conda:
        "../envs/base_python.yaml"
    threads: AVAILABLE_THREADS
    params:
        samples_names="\n".join(unique_samples),
    shell:
        "cat {input.R1} | pigz -d -k -p{threads} > {output.R1} && "
        "echo 'Forward files were successfully concatenated' >> {log} && "
        "cat {input.R2} | pigz -d -k -p{threads} > {output.R2} && "
        "echo 'Reverse files were successfully concatenated' >> {log} && "
        "echo '{params.samples_names}' > {output.sample_list} && "
        "echo 'Copy samples names into samples_list.txt' >> {log} "


rule move_files_to_metarib:
    input:
        R1=expand(
            f"results/rrna/{{sample}}_fwd.fq.gz",
            sample=unique_samples,
        ),
        R2=expand(
            f"results/rrna/{{sample}}_rev.fq.gz",
            sample=unique_samples,
        ),
    output:
        R1=expand(
            f"results/MetaRib/data/{{sample}}.1.fq",
            sample=unique_samples,
        ),
        R2=expand(
            f"results/MetaRib/data/{{sample}}.2.fq",
            sample=unique_samples,
        ),
    conda:
        "../envs/base_python.yaml"
    log:
        "logs/metarib/move_files_to_metarib.log",
    script:
        "../scripts/cp_metarib.py"


rule config_file_metarib:
    output:
        f"results/MetaRib/MetaRib.cfg",
    params:
        PROJECT_DIR=lambda wildcards, output: output[0][:-12],
        DATA_DIR=f"results/MetaRib/data",
        SAMPLING_NUM=metarib_params.get("SAMPLING_NUM", "1000000")
        if metarib_params
        else "1000000",
        EM_PARA=metarib_params.get("EM_PARA", ""),
        EM_REF=private_metarib_params.get("EM_REF", ""),
        EM_BT=private_metarib_params.get("EM_BT", ""),
        MAP_PARA=metarib_params.get("MAP_PARA", ""),
        CLS_PARA=metarib_params.get("CLS_PARA", ""),
        MIN_COV=metarib_params.get("MIN_COV", "2"),
        MIN_PER=metarib_params.get("MIN_PER", "80"),
        BBTOOL="workflow/scripts/BBMap/sh",
    threads: AVAILABLE_THREADS
    conda:
        "../envs/metarib.yaml"
    log:
        "logs/metarib/config_file.log",
    shell:
        "echo [BASE] > {output} && "
        "echo PROJECT_DIR : $(pwd)/{params.PROJECT_DIR} >> {output} && "
        "echo DATA_DIR : $(pwd)/{params.PROJECT_DIR}/data >> {output} && "
        "echo SAMPLING_NUM : {params.SAMPLING_NUM} >> {output} && "
        "echo THREAD : {threads} >> {output} && "
        "echo [EMIRGE] >> {output} && "
        "echo  EM_PATH :  $(which emirge_amplicon.py) >> {output} && "
        "echo EM_PARA : {params.EM_PARA} >> {output} && "
        "echo EM_REF : {params.EM_REF} >> {output} && "
        "echo EM_BT : {params.EM_BT} >> {output} && "
        "echo [BBTOOL] >> {output} && "
        "echo  BB_PATH : $CONDA_PREFIX/bin >> {output} && "
        "echo MAP_PARA : {params.MAP_PARA} >> {output} && "
        "echo CLS_PARA : {params.CLS_PARA} >> {output} && "
        "echo [FILTER] >> {output} && "
        "echo MIN_COV : {params.MIN_COV} >> {output} && "
        "echo MIN_PER : {params.MIN_PER} >> {output} && "
        "cat  {output} > {log}"


rule MetaRib:
    input:
        config=f"results/MetaRib/MetaRib.cfg",
        R1=f"results/MetaRib/data/all.1.fq",
        R2=f"results/MetaRib/data/all.2.fq",
        samples=f"results/MetaRib/data/samples.list.txt",
        R1_samples=expand(
            f"results/MetaRib/data/{{sample}}.1.fq",
            sample=unique_samples,
        ),
    output:
        outdir=directory(f"results/MetaRib/MetaRib"),
        filtered_fasta=f"results/MetaRib/MetaRib/Abundance/all.dedup.filtered.fasta",
    params:
        outdir=f"results/MetaRib/Iteration/",
    log:
        "logs/metarib.log",
    conda:
        "../envs/metarib.yaml"
    shell:
        "python2 workflow/scripts/MetaRib/run_MetaRib.py -cfg {input.config}"
        ">> {log} 2>&1 "
