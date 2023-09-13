from collections import ChainMap

metarib_params = dict(ChainMap(*config.get("metarib")))
private_metarib_params = dict(
    ChainMap(
        *config.get("PRIVATE_METARIB"),
    )
)


rule decompress_rrna:
    input:
        fwd="results/rrna/{sample}_fwd.fq.gz",
        rev="results/rrna/{sample}_rev.fq.gz",
    output:
        fwd=temp("results/MetaRib/data/{sample}.1.fq"),
        rev=temp("results/MetaRib/data/{sample}.2.fq"),
    log:
        "logs/rrna/decompress_{sample}.log",
    conda:
        "../envs/pigz.yaml"
    threads: config["threads"]["pigz"]
    shell:
        "pigz -dkf -p{threads} < {input.fwd} > {output.fwd} && "
        "echo 'Forward file was successfully decompressed' >> {log} && "
        "pigz -dkf -p{threads} < {input.rev} > {output.rev} && "
        "echo 'Reverse file was successfully decompressed' >> {log} "


rule data_preparation:
    input:
        R1=expand(
            "results/MetaRib/data/{sample}.1.fq",
            sample=unique_samples,
        ),
        R2=expand(
            "results/MetaRib/data/{sample}.1.fq",
            sample=unique_samples,
        ),
    output:
        R1=temp("results/MetaRib/data/all.1.fq"),
        R2=temp("results/MetaRib/data/all.2.fq"),
        sample_list="results/MetaRib/data/samples.list.txt",
    log:
        "logs/metarib/data_preparation.log",
    conda:
        "../envs/pigz.yaml"
    threads: config["threads"]["pigz"]
    params:
        samples_names="\n".join(unique_samples),
    shell:
        "cat {input.R1} > {output.R1} && "
        "echo 'Forward files were successfully concatenated' >> {log} && "
        "cat {input.R2} > {output.R2} && "
        "echo 'Reverse files were successfully concatenated' >> {log} && "
        "echo '{params.samples_names}' > {output.sample_list} && "
        "echo 'Copy samples names into samples_list.txt' >> {log} "


rule config_file_metarib:
    output:
        f"results/MetaRib/MetaRib.cfg",
    params:
        PROJECT_DIR="results/MetaRib/",
        SAMPLING_NUM=metarib_params.get("SAMPLING_NUM", "1000000"),
        EM_PARA=metarib_params.get("EM_PARA", ""),
        EM_REF=private_metarib_params.get("EM_REF", ""),
        EM_BT=private_metarib_params.get("EM_BT", ""),
        MAP_PARA=metarib_params.get("MAP_PARA", ""),
        CLS_PARA=metarib_params.get("CLS_PARA", ""),
        MIN_COV=metarib_params.get("MIN_COV", "2"),
        MIN_PER=metarib_params.get("MIN_PER", "80"),
    threads: config["threads"]["metarib"]
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
        R1="results/MetaRib/data/all.1.fq",
        R2="results/MetaRib/data/all.2.fq",
        sample_list="results/MetaRib/data/samples.list.txt",
        config="results/MetaRib/MetaRib.cfg",
        R1_samples=expand(
            "results/MetaRib/data/{sample}.1.fq",
            sample=unique_samples,
        ),
        R2_samples=expand(
            "results/MetaRib/data/{sample}.2.fq",
            sample=unique_samples,
        ),
    output:
        outdir=directory("results/MetaRib/MetaRib"),
        filtered_fasta="results/MetaRib/all.dedup.filtered.fasta",
    log:
        "logs/metarib.log",
    conda:
        "../envs/metarib.yaml"
    threads: config["threads"]["metarib"]
    params:
        script="workflow/external_scripts/run_MetaRib.py",
    shell:
        "python2 {params.script} -cfg {input.config} "
        ">> {log} 2>&1 && "
        "mv {output.outdir}/Abundance/all.dedup.filtered.fasta {output.filtered_fasta} "
