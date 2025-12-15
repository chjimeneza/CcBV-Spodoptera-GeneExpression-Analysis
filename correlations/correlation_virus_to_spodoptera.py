import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests

# Load the CPM-normalized expression matrix (loci as rows, samples as columns)
df = pd.read_csv("../mapping/SRR_protein_count/matrix_CPM.tsv", sep="\t", index_col=0)

# Separate virus and host genes based on naming convention (e.g., "CcBV_" for viral genes)
virus_genes = [gene for gene in df.index if gene.startswith("CcBV_")]
host_genes = [gene for gene in df.index if not gene.startswith("CcBV_")]

results = []

# Loop over each virus gene
for virus_gene in virus_genes:
    virus_vector = df.loc[virus_gene].values

    # Correlate with all host genes
    for host_gene in host_genes:
        host_vector = df.loc[host_gene].values

        # Only use positions (samples) where both virus and host genes have non-zero expression
        mask = (virus_vector != 0) & (host_vector != 0)
        filtered_virus = virus_vector[mask]
        filtered_host = host_vector[mask]

        # Only compute correlation if at least 3 paired non-zero samples exist
        if len(filtered_virus) >= 3:
            r, p = pearsonr(filtered_virus, filtered_host)
            results.append({
                "gene_1": virus_gene,
                "gene_2": host_gene,
                "pearson_r": r,
                "p_value": p
            })

    # Correlate with all other virus genes (including self-comparison)
    for other_virus_gene in virus_genes:
        other_vector = df.loc[other_virus_gene].values
        mask = (virus_vector != 0) & (other_vector != 0)
        filtered_1 = virus_vector[mask]
        filtered_2 = other_vector[mask]

        if len(filtered_1) >= 3:
            r, p = pearsonr(filtered_1, filtered_2)
            results.append({
                "gene_1": virus_gene,
                "gene_2": other_virus_gene,
                "pearson_r": r,
                "p_value": p
            })

# Convert the results list into a DataFrame
df_results = pd.DataFrame(results)

# Adjust p-values using Benjamini-Hochberg correction (FDR control)
_, fdr_corrected, _, _ = multipletests(df_results["p_value"], method='fdr_bh')
df_results["fdr"] = fdr_corrected

# Filter: retain gene pairs with absolute Pearson correlation > 0.7 and FDR < 0.05
df_filtered = df_results[
    (df_results["pearson_r"].abs() > 0.7) &
    (df_results["fdr"] < 0.05)
]

# Save the filtered results to a TSV file
df_filtered.to_csv("../mapping/SRR_protein_count/correlations_filtered_with_pvalue_0.7_and_fdr.tsv", sep="\t", index=False)

print(f" {len(df_filtered)} gene pairs saved with Pearson r > 0.7 and FDR correction.")
