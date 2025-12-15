import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Define base paths for input and output
base_dir = "../mapping/SRR_protein_count"
corr_file = os.path.join(base_dir, "correlations_inmuno_only.tsv")  # Correlation results file (filtered for immune-related genes)
expr_file = os.path.join(base_dir, "matrix_CPM.tsv")  # Normalized expression matrix
output_dir = os.path.join(base_dir, "scatter_plots_0.7_immuno")  # Directory to save plots

# Load the correlation table (should include at least columns: gene_1, gene_2, pearson_r, fdr, optionally GO_terms)
df_corr = pd.read_csv(corr_file, sep="\t")

# Load the CPM-normalized expression matrix
df_expr = pd.read_csv(expr_file, sep="\t", index_col=0)

# Filter gene pairs with absolute Pearson correlation > 0.7 and statistically significant FDR (< 0.05)
df_filtered = df_corr[(df_corr["fdr"] < 0.05) & (df_corr["pearson_r"].abs() > 0.7)]

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Iterate over each filtered gene pair
for _, row in df_filtered.iterrows():
    viral_gene = row['gene_1']      # Typically the viral gene
    host_gene = row['gene_2']       # Typically the host gene
    corr_val = row['pearson_r']     # Correlation coefficient
    fdr_val = row['fdr']            # FDR-corrected p-value

    # Retrieve GO terms if available; format nicely for the plot title
    go_terms = row.get('GO_terms', '').replace(";", "\n")

    # Check that both genes are present in the expression matrix
    if viral_gene not in df_expr.index or host_gene not in df_expr.index:
        continue  # Skip this pair if data is missing

    # Extract expression values for both genes
    x = df_expr.loc[host_gene]
    y = df_expr.loc[viral_gene]

    # Remove samples where either gene has zero expression (common in sparse expression data)
    mask = (x != 0) & (y != 0)
    x = x[mask]
    y = y[mask]

    # Require at least 3 data points to compute and visualize a meaningful correlation
    if len(x) < 3:
        continue

    # Set up the plot
    plt.figure(figsize=(7, 5))
    sns.regplot(
        x=x, y=y,
        ci=95,  # 95% confidence interval around the regression line
        scatter_kws={"s": 50},  # Size of points
        line_kws={"color": "gray"}  # Line color
    )

    # Axis labels
    plt.xlabel(f"{host_gene} expression")
    plt.ylabel(f"{viral_gene} expression")

    # Title with correlation info and optional GO terms
    plt.title(
        f"{viral_gene} vs {host_gene}\n"
        f"Pearson r = {corr_val:.2f}, FDR = {fdr_val:.2e}\n{go_terms}",
        fontsize=8
    )

    plt.tight_layout()

    # Construct file name and save the plot
    fname = os.path.join(output_dir, f"{viral_gene}__vs__{host_gene}.png")
    plt.savefig(fname, dpi=300)
    plt.close()  # Close the figure to avoid overlap issues
