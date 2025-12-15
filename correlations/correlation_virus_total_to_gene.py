import pandas as pd
from scipy.stats import pearsonr

# Load the CPM-normalized expression matrix (genes as rows, samples as columns)
expr_df = pd.read_csv("../mapping/SRR_protein_count/matrix_CPM.tsv", sep="\t", index_col=0)

# Identify viral genes (those starting with 'CcBV_')
viral_genes = [g for g in expr_df.index if g.startswith("CcBV_")]

# All other genes are considered host genes
host_genes = [g for g in expr_df.index if g not in viral_genes]

# Compute the total viral expression per sample by summing all viral genes
viral_sum = expr_df.loc[viral_genes].sum(axis=0)

# Initialize list to collect correlation results
results = []

# For each host gene, compute the correlation with total viral expression
for host in host_genes:
    x = expr_df.loc[host]     # Host gene expression vector (across samples)
    y = viral_sum             # Total viral expression vector (same length)

    # Remove samples where either value is zero (to avoid distortion in correlation)
    mask = (x != 0) & (y != 0)
    x_filt = x[mask]
    y_filt = y[mask]

    # Only compute correlation if at least 3 valid paired values remain
    if len(x_filt) >= 3:
        r, p = pearsonr(x_filt, y_filt)  # Pearson correlation and p-value
    else:
        r, p = None, None  # Not enough data to compute a valid correlation

    # Store the results for this host gene
    results.append({
        "host_gene": host,
        "pearson_r": r,
        "p_value": p
    })

# Convert the list of results to a DataFrame
cor_df = pd.DataFrame(results)

# Save the correlation results to a TSV file for downstream analysis
cor_df.to_csv("../network/correlation_spodoptera_vs_total_virus.tsv", sep="\t", index=False)
