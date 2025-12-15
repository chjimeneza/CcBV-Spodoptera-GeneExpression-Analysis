# Virus–Host Gene Expression and Network Analysis

This repository contains the computational analysis performed to study the impact of
Cotesia congregata bracovirus (CcBV) on gene expression and immune-related networks
in Spodoptera frugiperda larvae.

## Overview
The workflow integrates transcriptomic data processing, viral–host gene expression
correlation, functional enrichment analysis, and biological network analysis to explore
how viral genes may modulate host metabolic, immune, and epigenetic pathways.

## Data
- RNA-seq data were obtained from NCBI SRA (BioProject: PRJNA818900).
- Protein reference sequences were retrieved from GenBank and UniProt.

## Methods
Main steps of the analysis:
1. Download and preprocessing of RNA-seq data using SRA Toolkit.
2. Mapping of transcriptomic reads against host and viral proteins using DIAMOND (blastx).
3. Construction and CPM normalization of gene-level expression matrices.
4. Pearson correlation analysis between viral and host gene expression.
5. Functional enrichment analysis (GO) using UniProt, Pfam (HMMER), and GOATOOLS.
6. Projection of immune-related protein interaction networks and visualization in Gephi.

## Tools and Libraries
- Python (pandas, numpy, scipy)
- DIAMOND
- GOATOOLS
- HMMER (Pfam)
- NetworkX
- Gephi

## Results
Key results include:
- Enrichment of lipid metabolism and epigenetic-related GO terms among host genes
  correlated with viral expression.
- Identification of immune-related network modules potentially modulated during infection.

## Report
A full methodological description and discussion of results can be found in:
`docs/PAE_report.pdf`
