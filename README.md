# HIV_Microbiome_Analysis
This repository contains the workflow used to perform microbial association network analysis on 16S rRNA sequencing data from healthy controls and HIV-infected individuals. The pipeline integrates taxonomic classification, ecological network inference, and graph-theoretical analysis to investigate changes in gut microbial interactions associated with HIV infection.

Methodology

The analysis pipeline consists of the following steps:

1. Literature Search & Metadata Curation
    * Publicly available 16S rRNA datasets and associated metadata were collected and curated from multiple studies.
2. Taxonomic Classification
    * Raw sequencing reads were taxonomically classified to the species level using SPINGO.
3. Species Filtering
    * Low-prevalence and low-abundance taxa were removed to reduce sparsity and retain robust microbial species for downstream analysis.
4. Network Inference
    * Microbial association networks for control and HIV cohorts were inferred using SPIEC-EASI with the Meinshausen–Bühlmann neighborhood selection algorithm.
5. Community Detection
    * Network modules were identified using the Louvain clustering algorithm.
6. Network Visualization
    * Singleton taxa were removed prior to visualization, and networks were visualized using the qgraph package.
7. Network Analysis
    * Graph-theoretical properties, including connectivity, community structure, centrality, modularity, and other topological metrics, were computed to compare microbial interaction networks between healthy controls and HIV-infected individuals.

Repository Contents

* Codes/ – R scripts for taxonomic processing, network inference, and analysis.
* Environment/ - Contains the .RData file on which the codes were run
