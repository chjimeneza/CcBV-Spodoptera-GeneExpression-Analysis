import pandas as pd
import re

# Load GO annotations derived from Pfam (previously mapped to locus IDs)
pfam_go_df = pd.read_csv("../mapping/proteins_with_locus_GO.tsv", sep="\t")

# Standardize column name to 'locus' for consistency across files
pfam_go_df = pfam_go_df.rename(columns={"Protein": "locus"})

# Load UniProt annotation file for Spodoptera genes
uniprot_df = pd.read_csv("../mapping/SRR_protein_count/uniprot_spodoptera.tsv", sep="\t")

# Normalize locus identifiers: keep only the first gene name in case of multiple entries
uniprot_df["locus"] = uniprot_df["Gene Names"].str.split().str[0]

# Retain only the locus and GO annotation columns from UniProt
uniprot_go_df = uniprot_df[["locus", "Gene Ontology (GO)"]].rename(
    columns={"Gene Ontology (GO)": "GO_terms_uniprot"}
)

# Merge Pfam and UniProt GO annotations by locus ID using an outer join to retain all available annotations
merged_df = pd.merge(pfam_go_df, uniprot_go_df, on="locus", how="outer")

# Function to combine GO terms from both sources (Pfam and UniProt), extracting clean GO:000XXXX IDs
def merge_go_terms(row):
    terms = set()
    for col in ["GO_terms", "GO_terms_uniprot"]:
        if pd.notna(row.get(col)):
            parts = row[col].split(";")
            for part in parts:
                # Extract GO IDs (e.g., GO:0001234), ignoring accompanying text
                go_ids = [x.strip() for x in re.findall(r"GO:\d{7}", part)]
                terms.update(go_ids)
    return ";".join(sorted(terms))

# Apply the merging function across all rows
merged_df["GO_terms_merged"] = merged_df.apply(merge_go_terms, axis=1)

# Prepare final DataFrame with merged GO annotations per locus
final_df = merged_df[["locus", "GO_terms_merged"]].rename(columns={"GO_terms_merged": "GO_terms"})

# Save the consolidated annotation table to file
final_df.to_csv("../mapping/final_GO_annotations.tsv", sep="\t", index=False)
