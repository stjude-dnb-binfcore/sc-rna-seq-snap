# Using Analysis Modules in scRNA-Seq Snap Workflow

This repository contains a collection of analysis modules designed to process and analyze single cell and single nuclei RNA (sc/snRNA) data from 10X sequencing technology. 

Each module is self-contained and can be executed independently or as part of a larger analysis pipeline. Below is a summary of each analysis module, including whether they are required or optional. Furthermore, the analysis modules should be run in the following recommended order:

1. `fastqc-analysis` module (description="Pipeline for FastQC quality control tool for high throughput sequence data analysis.", required=True)
2. `cellranger-analysis` module (description="Pipeline for running and summarizing Cell Ranger count for single or multiple libraries.", required=True)
3. `upstream-analysis` module (description="Pipeline for estimating QC metrics and filtering low quality cells.", required=True)
4. `integrative-analysis` module (description="Pipeline for Integrative analysis.", required=False)
5. `cluster-cell-calling` module (description="Pipeline for cluster cell calling and gene marker analysis.", required=True)
6. `cell-contamination-removal-analysis` module (description="To remove clusters and repeat steps (4) and (5), e.g. for PDX experiments.", required=False)
7. `cell-types-annotation` module (description="Pipeline for annotating cell types.", required=True)
8. `rshiny-app` module (description="Pipeline for generating an R shiny app for the project.", required=False)
9. `clone-phylogeny-analysis` module (description="Pipeline for Clone phylogeny analysis tool. This is currently available for human data only.", required=False)
10. `de-go-analysis` module (description="Pipeline for Differential Expression, Volcano Plots, and Gene Ontology (GO) enrichment analysis Per Cell Type.", required=False)
11. `project-updates` module (description="Pipeline for summarizing results from all modules and generating project reports.", required=False)

## Contact

Contributions, issues, and feature requests are welcome! Please feel free to check [issues](https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap/issues).

---

*These tools and pipelines have been developed by the Bioinformatics core team at the [St. Jude Children's Research Hospital](https://www.stjude.org/). These are open access materials distributed under the terms of the [BSD 2-Clause License](https://opensource.org/license/bsd-2-clause), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
