import os
import shutil
import re

DEFAULT_SOURCE_FILEPATH = snakemake.config.get("DEFAULT_SOURCE_FILEPATH", "reads/")
DEFAULT_DEST_FILEPATH =snakemake.config.get("DEFAULT_DEST_FILEPATH", "results/")
RENAMED_READS_FILEPATH = snakemake.config.get("RENAMED_READS_FILEPATH", "renamed_raw_reads/")

regular_expression = ".+[-|_](.+)_.+_(.+)_.+\.fastq.gz"
input_dir = DEFAULT_SOURCE_FILEPATH
output_dir = f"{DEFAULT_DEST_FILEPATH}{RENAMED_READS_FILEPATH}"

def copy_file(original: str, target: str) -> None:
    os.makedirs(os.path.dirname(target), exist_ok=True)
    shutil.copyfile(original, target)

for file in os.listdir(input_dir):
    match = re.match(regular_expression, file)
    if match:
        original = input_dir + match.string
        print(original)
        target = output_dir + match.group(1) + "_" + match.group(2) + ".fastq.gz"
        print(target)
        if os.path.exists(target):
            print('The previous file was already there.')
        else:
            copy_file (original, target)
  