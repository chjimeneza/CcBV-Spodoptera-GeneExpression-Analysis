# Import required libraries
import os
import pandas as pd
from collections import defaultdict, Counter
import re

# Paths to input and output files
carpeta = "../mapping/SRR_protein_count"  # Directory with BLASTX output files
archivo_mapeo = "../mapping/id_locus_or_gene.tsv"  # File mapping sequence IDs to loci
archivo_salida = os.path.join(carpeta, "matrix_count.tsv")  # Output file for the count matrix

# Desired order of the columns in the final matrix
orden_columnas = [
    "SRR18458706_unparasitized_2h", "SRR18458715_unparasitized_2h", "SRR18458716_unparasitized_2h",
    "SRR18458703_unparasitized_24h", "SRR18458704_unparasitized_24h", "SRR71845705_unparasitized_24h",
    "SRR18458700_unparasitized_48h", "SRR18458701_unparasitized_48h", "SRR18458702_unparasitized_48h",
    "SRR18458699_parasitized_2h", "SRR18458713_parasitized_2h", "SRR18458714_parasitized_2h",
    "SRR18458710_parasitized_24h", "SRR18458711_parasitized_24h", "SRR18458712_parasitized_24h",
    "SRR18458707_parasitized_48h", "SRR18458708_parasitized_48h", "SRR18458709_parasitized_48h"
]

# 1. Load the sseqid to locus mapping
sseqid_to_locus = {}
with open(archivo_mapeo, "r") as f:
    for linea in f:
        partes = linea.strip().split("\t")
        if len(partes) < 2:
            continue
        clave, locus = partes[0], partes[1]
        # Extracts standardized RefSeq protein IDs (e.g., XP_123456.1)
        match = re.search(r"(XP|YP|NP)_[0-9]+\.[0-9]+", clave)
        if match:
            sseqid = match.group(0)
            sseqid_to_locus[sseqid] = locus

# 2. Search for BLASTX output files and count hits by sseqid
conteos_por_locus = defaultdict(dict)

# Select only relevant files that follow the expected naming pattern
archivos_blast = sorted([
    f for f in os.listdir(carpeta)
    if f.startswith("blastx_diamond_output_SRR") and f.endswith(".tsv")
])

# Group forward and reverse files by SRR ID
grupos = defaultdict(list)
for archivo in archivos_blast:
    match = re.search(r"(SRR\d+)", archivo)
    if match:
        id_srr = match.group(1)
        grupos[id_srr].append(os.path.join(carpeta, archivo))

# Show detected SRR groups
print("Detected SRR groups:")
for srr_id, archivos in grupos.items():
    print(f"{srr_id} -> {archivos}")

# 3. Process each SRR group and sum forward + reverse counts
for srr_id, archivos in grupos.items():
    conteos_sseqid = Counter()
    for archivo in archivos:
        print(f"Processing file: {archivo}")
        try:
            with open(archivo) as f:
                for linea in f:
                    if linea.strip() == "" or linea.startswith("#"):
                        continue
                    partes = linea.strip().split("\t")
                    if len(partes) < 2:
                        continue
                    sseqid = partes[1]
                    conteos_sseqid[sseqid] += 1
        except Exception as e:
            print(f"Error processing {archivo}: {e}")

    # Map sseqid to locus and aggregate counts
    locus_counts = defaultdict(int)
    print(f"\nDebugging {srr_id}: {len(conteos_sseqid)} sseqid found")
    for sseqid, count in conteos_sseqid.items():
        locus = sseqid_to_locus.get(sseqid)
        if locus:
            print(f"{sseqid} → {locus} → +{count}")
            locus_counts[locus] += count
        else:
            print(f"{sseqid} has no associated locus (unmapped)")

    # Store counts by SRR ID
    for locus, count in locus_counts.items():
        conteos_por_locus[locus][srr_id] = count

# 4. Create DataFrame from the aggregated dictionary
df = pd.DataFrame.from_dict(conteos_por_locus, orient="index").fillna(0).astype(int)

# 5. Rename columns using biological sample labels
mapeo_nombres = {
    "SRR18458706": "SRR18458706_unparasitized_2h",
    "SRR18458715": "SRR18458715_unparasitized_2h",
    "SRR18458716": "SRR18458716_unparasitized_2h",
    "SRR18458703": "SRR18458703_unparasitized_24h",
    "SRR18458704": "SRR18458704_unparasitized_24h",
    "SRR18458705": "SRR71845705_unparasitized_24h",  # likely a typo in SRR ID
    "SRR18458700": "SRR18458700_unparasitized_48h",
    "SRR18458701": "SRR18458701_unparasitized_48h",
    "SRR18458702": "SRR18458702_unparasitized_48h",
    "SRR18458699": "SRR18458699_parasitized_2h",
    "SRR18458713": "SRR18458713_parasitized_2h",
    "SRR18458714": "SRR18458714_parasitized_2h",
    "SRR18458710": "SRR18458710_parasitized_24h",
    "SRR18458711": "SRR18458711_parasitized_24h",
    "SRR18458712": "SRR18458712_parasitized_24h",
    "SRR18458707": "SRR18458707_parasitized_48h",
    "SRR18458708": "SRR18458708_parasitized_48h",
    "SRR18458709": "SRR18458709_parasitized_48h"
}
df.rename(columns=mapeo_nombres, inplace=True)

# 6. Reorder columns based on the desired biological experiment layout
df = df.reindex(columns=orden_columnas)

# 7. Save the final count matrix to a TSV file
df.index.name = "Locus"
df.to_csv(archivo_salida, sep="\t")

print(f"\nFile saved as {archivo_salida}")
print(f"Columns present in the final DataFrame: {df.columns.tolist()}")

# Check if any expected columns are missing
faltantes = [col for col in orden_columnas if col not in df.columns]
print(f"Missing columns (if any): {faltantes}")
