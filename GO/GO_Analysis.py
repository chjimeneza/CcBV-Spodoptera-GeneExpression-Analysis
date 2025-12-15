import pandas as pd
from goatools.obo_parser import GODag
from goatools.goea.go_enrichment_ns import GOEnrichmentStudy
from collections import defaultdict

# Load correlation table (with protein names and GO annotations)
cor_df = pd.read_csv("../mapping/SRR_protein_count/correlations_with_protein_names_0.7.tsv", sep="\t")

# Load the final merged GO annotation file
go_df = pd.read_csv("../mapping/final_GO_annotations.tsv", sep="\t")

# Select and rename relevant columns
go_df = go_df[["locus", "GO_terms"]].rename(columns={
    "locus": "gene_2",
    "GO_terms": "GO_terms"
})

# Ensure gene_2 column is a list if multiple gene IDs exist per cell
go_df["gene_2"] = go_df["gene_2"].astype(str)
go_df = go_df.assign(gene_2=go_df["gene_2"].str.split())

# Explode multi-gene entries so each gene appears on a separate row
go_df = go_df.explode("gene_2")

# All unique UniProt genes form the background universe
all_genes = go_df["gene_2"].unique().tolist()

# Merge GO annotations with the correlation results
merged_df = cor_df.merge(go_df, on="gene_2", how="left")

# Save the merged annotation file
merged_df.to_csv("../mapping/SRR_protein_count/correlation_with_GO_0.7_and_hmmer.tsv", sep="\t", index=False)

# Separate into host and viral gene subsets
df_spodoptera = merged_df[merged_df["gene_2"].str.startswith("LOC")].copy()
df_virus = merged_df[merged_df["gene_2"].str.startswith("CcBV")].copy()

# Save both subsets for reference
df_spodoptera.to_csv("../mapping/SRR_protein_count/correlation_spodoptera_0.7_and_hmmer.tsv", sep="\t", index=False)
df_virus.to_csv("../mapping/SRR_protein_count/correlation_virus_0.7_and_hmmer.tsv", sep="\t", index=False)

# Work only with host genes for the enrichment analysis
df = df_spodoptera

# Extract genes significantly positively correlated with viral genes
sig_genes = df[df["pearson_r"] > 0.7]["gene_2"].unique().tolist()

# Print summary
print(f"Genes in the universe: {len(all_genes)}")
print(f"Genes in the test set: {len(sig_genes)}")

# Build gene-to-GO mapping dictionary
gene2go_raw = defaultdict(set)
for _, row in go_df.dropna(subset=["GO_terms"]).iterrows():
    gene = row["gene_2"]
    terms = [term.split("[")[-1].strip("]") for term in row["GO_terms"].split(";")]
    gene2go_raw[gene].update(terms)

# Load GO ontology structure (OBO format)
go_dag = GODag(
    "../mapping/SRR_protein_count/go-basic.obo",
    optional_attrs=['relationship'],
    load_obsolete=True
)

# Filter GO terms to keep only valid ones present in the ontology DAG
gene2go_clean = defaultdict(set)
for gene, terms in gene2go_raw.items():
    valid_terms = [go for go in terms if go in go_dag]
    if valid_terms:
        gene2go_clean[gene].update(valid_terms)

print(f"Genes with valid GO terms: {len(gene2go_clean)} / {len(gene2go_raw)}")

# Initialize the GO enrichment analysis
goea = GOEnrichmentStudy(
    all_genes,            # Background gene universe
    gene2go_clean,        # Cleaned gene→GO term mapping
    go_dag,               # GO DAG
    propagate_counts=False,
    alpha=0.05,           # Significance threshold
    methods=['fdr_bh']    # Multiple testing correction: Benjamini-Hochberg
)

# Run enrichment analysis using the list of significant host genes
results = goea.run_study(sig_genes)

# Filter for significant GO terms (adjusted p-value < 0.05)
sig_results = [r for r in results if r.p_fdr_bh < 0.05]

# Print top 10 enriched GO terms
for r in sig_results[:10]:
    print(f"{r.GO} | {r.name} | p={r.p_fdr_bh:.3e} | study: {r.study_count}/{r.study_n}")

# Save the full GO enrichment results to file
goea.wr_tsv("../mapping/SRR_protein_count/go_enrichment_results_spodoptera_0.7_and_hmmer.tsv", results)
