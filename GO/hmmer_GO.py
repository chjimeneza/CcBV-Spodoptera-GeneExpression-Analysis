import pandas as pd
from collections import defaultdict
import re

# Parse the Pfam domtblout file to extract protein–Pfam domain relationships
pfam_hits = []
with open("../mapping/domtblout.txt") as f:
    for line in f:
        if line.startswith("#"):
            continue  # Skip comment lines
        parts = line.strip().split()
        if len(parts) < 5:
            continue  # Skip malformed lines
        pfam_acc_full = parts[1]  # Full Pfam accession, e.g., PF00102.32
        pfam_acc = pfam_acc_full.split('.')[0]  # Truncate version suffix → PF00102
        protein_id = parts[3]  # Protein identifier, e.g., YP_184757.1
        pfam_hits.append((protein_id, pfam_acc))

# Convert list of tuples into a DataFrame
pfam_df = pd.DataFrame(pfam_hits, columns=["Protein", "Pfam"])

# Load the Pfam-to-GO mapping file (Pfam2GO)
pfam2go = defaultdict(set)
with open("../mapping/pfam2go") as f:
    for line in f:
        if line.startswith("!"):
            continue  # Skip header/comment lines
        pfam_id = line.split()[0].replace("Pfam:", "")  # Extract Pfam ID
        go_ids = re.findall(r"GO:(\d{7})", line)  # Extract all GO terms (e.g., GO:0001234)
        for go_id in go_ids:
            pfam2go[pfam_id].add(f"GO:{go_id}")

# Associate GO terms with each protein via Pfam domains
protein2go = defaultdict(set)
for protein, pfam in pfam_df.values:
    go_terms = pfam2go.get(pfam, [])  # Retrieve GO terms linked to the Pfam ID
    protein2go[protein].update(go_terms)

# Build a DataFrame from the protein-to-GO associations
rows = []
for protein, go_terms in protein2go.items():
    if go_terms:
        rows.append({"Protein": protein, "GO_terms": ";".join(go_terms)})

df = pd.DataFrame(rows)

# Load the mapping file that connects raw identifiers to locus IDs
mapping_df = pd.read_csv("../mapping/SRR_protein_count/id_locus_or_gene.tsv", sep="\t", header=None, names=["raw", "locus"])

# Create a dictionary to map protein IDs (e.g., YP_...) to locus identifiers
id_to_locus = {}
for _, row in mapping_df.iterrows():
    for protein_id in df["Protein"]:
        if protein_id in row["raw"]:
            id_to_locus[protein_id] = row["locus"]

# Replace protein IDs in the GO annotation DataFrame with their corresponding locus IDs (if available)
df["Protein"] = df["Protein"].map(id_to_locus).fillna(df["Protein"])

# Group by locus and merge GO terms, ensuring uniqueness
df = df.groupby("Protein")["GO_terms"].apply(
    lambda x: ";".join(sorted(set(";".join(x).split(";"))))
).reset_index()

# Save the final GO annotation table to a TSV file
df.to_csv("../mapping/proteins_with_locus_GO.tsv", sep="\t", index=False)
