from pathlib import Path


def regex_filename(filename):
    # Get filename
    filename = Path(filename).name
    # Remove everything after either extension_fwd or extension_rev
    filename = filename.split("_fwd")[0]
    filename = filename.split("_rev")[0]
    return filename
