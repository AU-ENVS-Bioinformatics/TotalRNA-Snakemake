import sys
import sqlite3
import pandas as pd
from contextlib import closing
from ete3 import NCBITaxa
import json

# Parse Snakemake inputs and parameters
mapping_db = snakemake.input.mapping
counts_file = snakemake.input.counts
brite_json = snakemake.input.brite
diamond_results_file = snakemake.input.tsv
min_score_threshold = snakemake.params.min_score_threshold
threshold = snakemake.params.threshold
output_file = snakemake.output.annotated
output_brite = snakemake.output.brite
output_meta = snakemake.output.meta
log = snakemake.log[0]

ncbi = NCBITaxa()


def annotate_accessions(conn, accessions):
    """Fetch annotations for given accessions from the mapping database."""
    with closing(conn.cursor()) as cursor:
        # Create temporary table for efficient joining
        cursor.execute('CREATE TEMP TABLE temp_accessions (accession TEXT)')
        cursor.executemany('INSERT INTO temp_accessions (accession) VALUES (?)', 
                          [(acc,) for acc in accessions])
        
        # Query annotations with ORDER BY for deterministic results
        query = """
            SELECT Taxonomy, INTERPRO2GO, KEGG, mappings.Accession 
            FROM mappings 
            INNER JOIN temp_accessions ON mappings.Accession = temp_accessions.accession
            ORDER BY mappings.Accession, Taxonomy, INTERPRO2GO, KEGG
        """
        cursor.execute(query)
        df = pd.DataFrame(cursor.fetchall(), 
                         columns=['Taxonomy', 'INTERPRO2GO', 'KEGG', 'Accession'])
        df = df.drop_duplicates(subset=['Accession'])
        
        # Format KEGG and INTERPRO2GO identifiers
        df['KEGG'] = df['KEGG'].apply(
            lambda val: f"K{str(int(val)).zfill(5)}" if not pd.isna(val) else ""
        )
        df['INTERPRO2GO'] = df['INTERPRO2GO'].apply(
            lambda val: f"IPR{str(int(val)).zfill(6)}" if not pd.isna(val) else ""
        )
        
        return df


def lowest_common_ancestor(tax_ids):
    """Get lowest common ancestor of a list of tax ids"""
    # Remove zero elements
    tax_ids = [t for t in tax_ids if t != 0]
    if len(tax_ids) == 0:
        return "0"
    if len(tax_ids) == 1:
        return str(tax_ids[0])
    try:
        tree = ncbi.get_topology(tax_ids)
        root = tree.get_tree_root()
        return root.name
    except ValueError:
        # If error, try again without the last element, which is "worse" than the previous one
        return lowest_common_ancestor(tax_ids[:-1])


def get_lineage(taxid):
    """Get dictionary with rank as keys and name as values"""
    if taxid <= 0:
        return {}
    try:
        lineage = ncbi.get_lineage(taxid)
        ranks = ncbi.get_rank(lineage)
        names = ncbi.translate_to_names(lineage)
        return {ranks[k]: elm for k, elm in zip(lineage, names)}
    except ValueError:
        return {}


def get_taxonomy_row(taxid):
    """Parse taxa into greengenes-like format"""
    lineage = get_lineage(taxid)
    ranks = {
        'superkingdom': 'd__', 'phylum': 'p__', 'class': 'c__',
        'order': 'o__', 'family': 'f__', 'genus': 'g__', 'species': 's__'
    }
    for rank in ranks:
        ranks[rank] += lineage.get(rank, '')
    return ";".join(ranks.values())


def get_brite(kos, include=[]):
    """Extract BRITE annotations for given KO identifiers"""
    annotations = dict()
    for LevelA in brite['children']:
        if LevelA['name'] not in include:
            continue
        for LevelB in LevelA.get('children', []):
            for LevelC in LevelB.get('children', []):
                for LevelD in LevelC.get('children', []):
                    name = LevelD['name'].split()[0]
                    if name in kos:
                        description = LevelD['name'].split(';')[1].strip()
                        new_annotation = (LevelA['name'], LevelB['name'], LevelC['name'], description)
                        if name in annotations:
                            annotations[name].append(new_annotation)
                        annotations[name] = [new_annotation]
    return annotations


# Redirect output to log file
with open(log, "w") as f:
    sys.stderr = sys.stdout = f
    
    # ===== STEP 1: ANNOTATE ACCESSIONS =====
    print("Step 1: Annotating accessions...", file=sys.stderr)
    
    # Load input data
    counts = pd.read_csv(counts_file, sep="\t", header=0, index_col='ContigID')
    diamond_results = pd.read_csv(diamond_results_file, sep="\t", header=0, 
                                  usecols=['qseqid', 'sseqid', 'evalue', 'bitscore'])
    
    # Rename columns to standard names
    diamond_results.rename(columns={
        'qseqid': 'Query_id',
        'sseqid': 'Subject_id',
        'evalue': 'e-value',
        'bitscore': 'score'
    }, inplace=True)
    
    # Filter by minimum score threshold
    diamond_results = diamond_results[diamond_results['score'] > min_score_threshold]
    
    # Keep only hits within threshold of the best score for each query
    diamond_results = diamond_results[
        diamond_results.groupby('Query_id')['score'].transform(
            lambda x: x >= threshold * x.max()
        )
    ].reset_index(drop=True)
    
    
    print(f"Filtered to {len(diamond_results)} DIAMOND hits", file=sys.stderr)
    
    # Merge with read counts
    results = diamond_results.merge(counts, left_on='Query_id', right_index=True)
    
    # Get annotations from database
    conn = sqlite3.connect(mapping_db)
    annotations = annotate_accessions(conn, set(results['Subject_id']))
    conn.close()
    
    print(f"Retrieved {len(annotations)} annotations from database", file=sys.stderr)
    
    # Merge annotations with results
    results = results.merge(annotations, left_on='Subject_id', right_on='Accession', how='left')
    results.rename(columns={'Taxonomy': 'NCBI_TaxID'}, inplace=True)
    
    # Reorder columns with Query_id and Subject_id first
    column_order = ['Query_id', 'Subject_id'] + [
        col for col in results.columns if col not in ['Query_id', 'Subject_id', 'Accession']
    ]
    results = results[column_order]
    
    # ===== STEP 2: LCA ALGORITHM =====
    print("Step 2: Computing Lowest Common Ancestor...", file=sys.stderr)
    
    df = results.sort_values(by=["e-value", "score"], ascending=[True, False])
    df_lca = df.groupby('Query_id')['NCBI_TaxID'].apply(lowest_common_ancestor).rename('LCA').to_frame()
    df_lca['Taxonomy'] = [get_taxonomy_row(int(taxid)) for taxid in df_lca['LCA']]
    
    df = df.merge(df_lca, left_on='Query_id', right_index=True)
    
    # ===== STEP 3: FUNCTIONAL ANNOTATION =====
    print("Step 3: Adding functional annotations...", file=sys.stderr)
    
    print(f"There are {df.shape[0]} rows in the dataframe", file=sys.stderr)
    print(f"There are {(df['KEGG'] == '').sum()} rows with no KEGG identifier", file=sys.stderr)
    print(f"There are {(df['INTERPRO2GO'] == '').sum()} rows with no INTERPRO2GO identifier", file=sys.stderr)
    
    df = df.sort_values(by=["e-value", "score"], ascending=[True, False]).drop_duplicates("Query_id", keep="first")
    df.drop(columns=['e-value', 'score'], inplace=True)
    
    # Load BRITE hierarchy
    with open(brite_json, 'r') as json_file:
        json_data = json_file.read()
        brite = json.loads(json_data)
    
    # Add BRITE annotations
    annotations = get_brite(df.KEGG.to_list(), include='09180 Brite Hierarchies')
    df['Brite_Level1'] = [",".join(x[1] for x in annotations[ko]) if ko in annotations else '' for ko in df.KEGG.to_list()]
    df['Brite_Level2'] = [",".join(x[2] for x in annotations[ko]) if ko in annotations else '' for ko in df.KEGG.to_list()]
    df['Brite_Level3'] = [",".join(x[3] for x in annotations[ko]) if ko in annotations else '' for ko in df.KEGG.to_list()]
    df_brite_filtered = df[df['Brite_Level1'] != '']
    
    # Add Metabolism annotations
    df.index.name = 'ContigID'
    annotations = get_brite(df.KEGG.to_list(), include='09100 Metabolism')
    df['Metabolism_Level1'] = [",".join(x[1] for x in annotations[ko]) if ko in annotations else '' for ko in df.KEGG.to_list()]
    df['Metabolism_Level2'] = [",".join(x[2] for x in annotations[ko]) if ko in annotations else '' for ko in df.KEGG.to_list()]
    df['Metabolism_Level3'] = [",".join(x[3] for x in annotations[ko]) if ko in annotations else '' for ko in df.KEGG.to_list()]
    df_meta_filtered = df[df['Metabolism_Level1'] != '']
    
    # Save outputs
    df.to_csv(output_file, sep='\t')
    df_brite_filtered.to_csv(output_brite, sep='\t')
    df_meta_filtered.to_csv(output_meta, sep='\t')
    
    print("Done!", file=sys.stderr)
