from pathlib import Path
import sys

DEFAULT_SOURCE_FILEPATH = "reads/"
DEFAULT_DEST_FILEPATH = "results/"
RENAMED_READS_FILEPATH = "renamed_raw_reads/"

input_dir = DEFAULT_SOURCE_FILEPATH
output_dir = f"{DEFAULT_DEST_FILEPATH}{RENAMED_READS_FILEPATH}"


with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    input_dir = Path(input_dir).absolute()
    output_dir = Path(output_dir).absolute()
    output_dir.symlink_to(input_dir)
    print("Symlinked", input_dir, "to", output_dir)