infile_path = snakemake.input
outfile_path = snakemake.output

TAXONOMY_COLUMN_POSITION = 31
PVC_POSITION_IN_TAXONOMY = 1
def get_taxonomy_list(x: [str]) -> [str]:
    taxonomy= x[TAXONOMY_COLUMN_POSITION] \
    .replace(" (superkingdom)","") \
    .replace(" (superphylum)","") \
    .split(";")
    return taxonomy[2:]

def add_prefix_taxonomy(x):
    y = ["k__","t__", "p__","c__","o__","f__","g__","s__"] 
    for i, v in enumerate(x):
        y[i] += v
    return y

with open(infile_path) as infile, open(outfile_path, 'w') as outfile:
        outfile.write(next(infile).replace("classification", "taxonomy")+'\n')
        for line in infile:
            cols = [col.strip() for col in line.strip().split("\t")]
            taxonomy = add_prefix_taxonomy(get_taxonomy_list(cols))
            taxonomy.pop(PVC_POSITION_IN_TAXONOMY)
            cols[TAXONOMY_COLUMN_POSITION] = "; ".join(taxonomy)
            outfile.write("\t".join(cols)+'\n')