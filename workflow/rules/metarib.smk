import subprocess
from collections import ChainMap
RRNA_FILEPATH = config.get("RRNA_FILEPATH", "rrna/")
METARIB_FILEPATH = config.get("METARIB_FILEPATH", "metarib/")
AVAILABLE_THREADS = int(workflow.cores * 0.75)

metarib_params = dict(ChainMap(*config.get('metarib')))
private_metarib_params = dict(ChainMap(*config.get('PRIVATE_METARIB'),))
rule prepare_metarib:
    input:
        R1= expand(f"{DEFAULT_DEST_FILEPATH}{RRNA_FILEPATH}{{sample}}_fwd.fq.gz", sample = sample),
        R2= expand(f"{DEFAULT_DEST_FILEPATH}{RRNA_FILEPATH}{{sample}}_rev.fq.gz", sample = sample),
    output:
        R1 = f"{DEFAULT_DEST_FILEPATH}{METARIB_FILEPATH}data/all.1.fq.gz",   
        R2 = f"{DEFAULT_DEST_FILEPATH}{METARIB_FILEPATH}data/all.2.fq.gz",
        sample_list = report(
            f"{DEFAULT_DEST_FILEPATH}{METARIB_FILEPATH}samples_list.txt",
            caption = "report/prepare_metarib.rst"
            ),
    log:
        "logs/sortmerna/concatenate_all.log",
    conda:
        "../envs/base_python.yaml"
    threads: AVAILABLE_THREADS
    params:
        samples_names = "\n".join(sorted(set(sample))),
    shell:
        "cat {input.R1} > {output.R1} && "
        "echo 'Forward files were successfully concatenated' >> {log} && "
        "cat {input.R2} > {output.R2} && "
        "echo 'Reverse files were successfully concatenated' >> {log} && "
        "echo '{params.samples_names}' > {output.sample_list} && "
        "echo 'Copy samples names into samples_list.txt' >> {log}"

rule config_file_metarib:
    output: f"{DEFAULT_DEST_FILEPATH}{METARIB_FILEPATH}MetaRib.cfg",
    params: 
        PROJECT_DIR = f"{DEFAULT_DEST_FILEPATH}{METARIB_FILEPATH}",
        DATA_DIR = f"{DEFAULT_DEST_FILEPATH}{METARIB_FILEPATH}data",
        SAMPLING_NUM = metarib_params.get('SAMPLING_NUM', "1000000") if metarib_params else "1000000",
        EM_PARA = metarib_params.get('EM_PARA', ""),
        EM_REF = private_metarib_params.get('EM_REF', ""),
        EM_BT = private_metarib_params.get('EM_BT', ""),
        MAP_PARA = metarib_params.get('MAP_PARA', ""),
        CLS_PARA = metarib_params.get('CLS_PARA', ""),
        MIN_COV = metarib_params.get('MIN_COV', "2"),
        MIN_PER = metarib_params.get('MIN_PER', "80"),
    threads: AVAILABLE_THREADS,
    conda:
        "../envs/metarib.yaml",
    shell:
        "echo [BASE] > {output} && "
        "echo PROJECT_DIR : {params.PROJECT_DIR} >> {output} && "
        "echo DATA_DIR : {params.DATA_DIR} >> {output} && "
        "echo SAMPLING_NUM : {params.SAMPLING_NUM} >> {output} && "
        "echo THREAD : {threads} >> {output} && "
        "echo [EMIRGE] >> {output} && "
        "echo  EM_PATH :  $(which emirge_amplicon.py) >> {output} && "
        "echo EM_PARA : {params.EM_PARA} >> {output} && "
        "echo EM_REF : {params.EM_REF} >> {output} && "
        "echo EM_BT : {params.EM_BT} >> {output} && "
        "echo [BBTOOL] >> {output} && "
        "echo  BB_PATH :  $(which bbmap.sh) >> {output} && "
        "echo MAP_PARA : {params.MAP_PARA} >> {output} && "
        "echo CLS_PARA : {params.CLS_PARA} >> {output} && "
        "echo [FILTER] >> {output} && "
        "echo MIN_COV : {params.MIN_COV} >> {output} && "
        "echo MIN_PER : {params.MIN_PER} >> {output} "