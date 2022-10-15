import subprocess
params_dict = snakemake.config

config_metarib_dict = {
    "BASE": {
        "PROJECT_DIR" : snakemake.params['PROJECT_DIR'],
        "DATA_DIR" : snakemake.params['DATA_DIR'],
        "SAMPLING_NUM": params_dict['metarib']['SAMPLING_NUM'],
        "THREAD" : snakemake.threads
    } ,
    "EMIRGE" : {
        "EM_PATH" : subprocess.check_output(['which', 'emirge_amplicon.py']).rstrip('\n'),
        "EM_PARA" : params_dict['metarib']['EM_PARA'],
        "EM_REF" : params_dict['PRIVATE_METARIB']['EM_REF'],
        "EM_BT" : params_dict['PRIVATE_METARIB']['EM_BT'],
    },
    "BBTOOL" : {
        "BB_PATH": subprocess.check_output(['which', 'bbmap.sh']).rstrip('\n'),
        "MAP_PARA" : params_dict['metarib']['MAP_PARA'],   
        "CLS_PARA" : params_dict['metarib']['CLS_PARA'],
    } ,
    "FILTER" : {
        "MIN_COV" : params_dict['metarib']['MIN_COV'],
        "MIN_PER" : params_dict['metarib']['MIN_PER'],

    }    
}
print(config_metarib_dict)