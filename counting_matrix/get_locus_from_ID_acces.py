from Bio import Entrez
from time import sleep

# Required by NCBI to identify who is making the requests
Entrez.email = "tu_email@ejemplo.com"

# List of NCBI protein accession IDs that did not have locus mappings
sseqids_sin_locus = [
    # A list of YP_ accession numbers (RefSeq protein IDs)
    "YP_184814.1", "YP_184776.1", "YP_184842.1", ...
]

# Dictionary to store the mapping of protein ID to locus
id_to_locus = {}

# Loop over each protein ID in the list
for protein_id in sseqids_sin_locus:
    try:
        # Fetch the GenBank record for the protein ID from NCBI Protein database
        handle = Entrez.efetch(db="protein", id=protein_id, rettype="gb", retmode="text")
        record = handle.read()
        handle.close()

        # Search the record line-by-line for the /locus_tag qualifier
        for line in record.split("\n"):
            if "/locus_tag=" in line:
                # Extract the locus_tag value and clean it
                locus = line.strip().split("=")[-1].replace('"', '')
                id_to_locus[protein_id] = locus
                break
        else:
            # If no locus_tag was found in the record
            id_to_locus[protein_id] = "NO_LOC_ENCONTRADO"

        # Wait to avoid overwhelming NCBI servers (polite scraping)
        sleep(0.3)

    except Exception as e:
        # If there's an error fetching or parsing, record it
        print(f"Error with {protein_id}: {e}")
        id_to_locus[protein_id] = "ERROR"

# Print the final mapping result for each protein ID
for k, v in id_to_locus.items():
    print(f"{k}\t{v}")
