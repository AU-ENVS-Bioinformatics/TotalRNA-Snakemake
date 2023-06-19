import os
from pathlib import Path
import shutil
import re
import sys

DEFAULT_SOURCE_FILEPATH = "reads/"
DEFAULT_DEST_FILEPATH = "results/"
RENAMED_READS_FILEPATH = "renamed_raw_reads/"

regular_expression = snakemake.config.get(
    "READS_REGEX", ".+[-|_](.+)_.+_(.+)_.+\.fastq.gz"
)
input_dir = DEFAULT_SOURCE_FILEPATH
output_dir = f"{DEFAULT_DEST_FILEPATH}{RENAMED_READS_FILEPATH}"


def copy_file(original: str, target: str) -> None:
    src, dest = Path(original).absolute(), Path(target).absolute()
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(src)
    print(dest)
    dest.symlink_to(src)


with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    for file in os.listdir(input_dir):
        match = re.match(regular_expression, file)
        if match:
            original = input_dir + match.string
            target = output_dir + match.group(1) + "_" + match.group(2) + ".fastq.gz"
            if os.path.exists(target):
                print("The previous file was already there.", file=sys.stderr)
            else:
                copy_file(original, target)
