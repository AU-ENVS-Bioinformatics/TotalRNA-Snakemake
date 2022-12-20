import dask.dataframe as dd
# Read input from snakemake
infile = snakemake.input[0]
outfile = snakemake.output[0]
evalue_threshold = snakemake.params["evalue"]
logfile = snakemake.log[0]

# Set up new colnames
colnames = [
    "query_id","subject_id","identity","alignment_length","mismatches","gap_openings","query_start",
    "query_end","subject_start","subject_end","evalue","score"
    ]

with open(logfile, "w") as f:
    print("Starting parsing the file...",file=f)
    # Read tsv file ignoring non #Query and headers lines
    ddf = dd.read_csv(
        infile, sep = "\t", comment='#',
        header = None,names = colnames, assume_missing=True
        ).dropna()
    # Filter by evalue
    ddf2 = ddf.query(f'evalue <= {evalue_threshold}').compute()
    # Parse query_id to remove last 2 characters ("_5")
    ddf2["query_id_parsed"] = ddf2["query_id"].str[:-2]
    # Get most significant hit for each query sequence 
    ddf2.sort_values(
        by=['evalue', 'score'], ascending=[True, False]
        ).drop_duplicates('query_id_parsed').drop('query_id_parsed', axis = 1).sort_index().to_csv(outfile, index = False)
    print("Done!",file=f)