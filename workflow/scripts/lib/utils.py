from pathlib import Path


def regex_filename(filename):
    # Get filename
    filename = Path(filename).name
    # Remove everything after "_sorted"
    filename = filename.split("_sorted")[0]
    return filename
