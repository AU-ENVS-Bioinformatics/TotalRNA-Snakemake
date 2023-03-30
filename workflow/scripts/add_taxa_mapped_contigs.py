import sys
import pandas as pd

taxa = snakemake.input.taxa
otu = snakemake.input.otu

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    print("Reading taxa...")
    taxa_df = pd.read_csv(str(taxa), sep = "\t", names = ["OTU", "taxonomy"]) 
    print("Reading OTUs...")
    otu_df = pd.read_csv(str(otu), sep = "\t") 
    otu_df.rename(columns={'ContigID':'OTU'}, inplace=True)
    print("Merging both files...")
    otu_df.merge(
        taxa_df, on='OTU', how='left'
        ).to_csv(str(snakemake.output),sep='\t',index=False)
    print("Finished succesfully!")
