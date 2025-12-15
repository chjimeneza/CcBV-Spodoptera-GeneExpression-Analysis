import pandas as pd
import re

# Load UniProt annotation file (TSV format)
df_uniprot = pd.read_csv("../mapping/SRR_protein_count/uniprot_spodoptera_and_bracoviriform.tsv", sep="\t")

# Define a function to clean and simplify protein names
def clean_protein_name(name):
    name = re.sub(r'LOC\d+', '', name)  # Remove LOC identifiers (if present)
    name = re.sub(r'isoform.*', '', name, flags=re.IGNORECASE)  # Remove isoform-specific annotations
    return name.strip()

# Apply cleaning function to the 'Protein names' column
df_uniprot['clean_protein_name'] = df_uniprot['Protein names'].apply(clean_protein_name)

# Remove rows without gene names and split multiple gene names into separate rows
df_exploded = df_uniprot.dropna(subset=['Gene Names']).copy()
df_exploded['Gene Names'] = df_exploded['Gene Names'].str.split()  # Split space-separated names
df_exploded = df_exploded.explode('Gene Names')  # Expand into separate rows

# Group protein names by gene name, merge duplicates with "-"
# This ensures a clean mapping from gene → protein name(s)
grouped = df_exploded.groupby('Gene Names')['clean_protein_name'].apply(
    lambda x: '-'.join(sorted(set(x)))
).to_dict()

# Load virus-host correlation results
df_corr = pd.read_csv("../mapping/SRR_protein_count/correlations_filtered_with_pvalue_0.7_and_fdr.tsv", sep="\t")

# Map UniProt protein names to the 'gene_2' column (typically the host gene)
df_corr['protein_name'] = df_corr['gene_2'].map(grouped)

# Save the updated DataFrame to a TSV file including protein name annotations
df_corr.to_csv("../mapping/SRR_protein_count/correlations_with_protein_names_0.7.tsv", sep="\t", index=False)
