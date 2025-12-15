import pandas as pd

# Load the raw count matrix from a TSV file
df = pd.read_csv("../mapping/SRR_protein_count/matrix_count.tsv", sep="\t")

# Set the 'Locus' column as the DataFrame index (gene/protein identifiers)
df.set_index("Locus", inplace=True)

# Apply CPM (Counts Per Million) normalization
# For each column (sample), divide each count by the total count of that column and multiply by 1 million
df_cpm = df.div(df.sum(axis=0), axis=1) * 1e6

# Save the CPM-normalized matrix to a new file
df_cpm.to_csv("../mapping/SRR_protein_count/matrix_CPM.tsv", sep="\t")

print("CPM matrix saved as 'matrix_CPM.tsv'")
