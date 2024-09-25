import pandas as pd
import gzip
import sys

# Check if gather file is provided
if len(sys.argv) != 2:
    print("Usage: python sum_tax.py [gatherfile.csv]")
    sys.exit(1)

gather_file = sys.argv[1]

# Function to read gzipped CSV files
def read_gzipped_csv(filepath):
    with gzip.open(filepath, 'rt') as f:
        return pd.read_csv(f)

# Lineage database files
lineage_files = [
    '/media/vehuardo/Work_volume_10TB/Sourmash_gb_databases/genbank-2022.03-archaea.lineages.csv.gz',
    '/media/vehuardo/Work_volume_10TB/Sourmash_gb_databases/genbank-2022.03-bacteria.lineages.csv.gz',
    '/media/vehuardo/Work_volume_10TB/Sourmash_gb_databases/genbank-2022.03-fungi.lineages.csv.gz',
    '/media/vehuardo/Work_volume_10TB/Sourmash_gb_databases/genbank-2022.03-protozoa.lineages.csv.gz',
    '/media/vehuardo/Work_volume_10TB/Sourmash_gb_databases/genbank-2022.03-viral.lineages.csv.gz'
]

# Read lineage database files
lineage_df = pd.concat([read_gzipped_csv(file) for file in lineage_files])

# Read the gather file
gather_df = pd.read_csv(gather_file)

# Process gather files to extract accession numbers
gather_df['accession'] = gather_df['name'].apply(lambda x: x.split(' ')[0])

# Merge on accession number and sum f_unique_weighted
summary = gather_df.groupby('accession')['f_unique_weighted'].sum().reset_index()
result = summary.merge(lineage_df, left_on='accession', right_on='ident')

# Summarize at species level
final_summary = result.groupby(['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])['f_unique_weighted'].sum().reset_index()
final_summary = final_summary.rename(columns={'f_unique_weighted': 'abundance'})

# Output to a file
output_file = gather_file.replace('.csv', '_summary.csv')
final_summary.to_csv(output_file, index=False)

print(f"Summary written to {output_file}")