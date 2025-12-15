import pandas as pd

# Load the correlation results file (containing gene pairs with high correlation and low FDR)
df = pd.read_csv("../mapping/SRR_protein_count/correlations_filtered_with_pvalue_0.7_and_fdr.tsv", sep="\t")

# Extract the 'gene_2' column, which refers to host gene identifiers (e.g., LOC123456)
ncbi_gene_ids = df["gene_2"].dropna()  # Remove any missing entries

# Filter to keep only gene IDs that start with 'LOC' (NCBI gene identifier prefix)
ncbi_gene_ids = ncbi_gene_ids[ncbi_gene_ids.str.startswith("LOC")]

# Remove the 'LOC' prefix and keep only the numeric part (as expected by UniProt ID mapping tool)
ncbi_gene_ids = ncbi_gene_ids.str.replace("LOC", "", regex=False).unique()

# Save the resulting numeric gene IDs to a plain text file, one per line
# This can be directly pasted into UniProt's ID mapping tool to retrieve protein annotation info
with open("../mapping/SRR_protein_count/ncbi_gene_ids.txt", "w") as f:
    f.write("\n".join(ncbi_gene_ids))

# Print summary
print(f"{len(ncbi_gene_ids)} NCBI GeneIDs (numeric) saved in ncbi_gene_ids.txt")
