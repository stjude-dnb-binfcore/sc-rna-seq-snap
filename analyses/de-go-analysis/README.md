# Pipeline for Differential Expression, Volcano Plots, and Gene Ontology (GO) enrichment analysis Per Cell Type for sc-RNA-Seq Analysis in 10X Genomics data


## Usage

`run-de-go-analysis.sh` is designed to be run as if it was called from this module directory even when called from outside of this directory.

Parameters according to the project and analysis strategy will need to be specified in the following scripts:
- `project_parameters.Config.yaml` located at the `root_dir`.


### Run module on an interactive session on HPC within the container

To run all of the Rscripts in this module sequentially on an interactive session on HPC, please run the following command from an interactive compute node:

```
bash run-de-go-analysis.sh
```

### Run module by using lsf on HPC within the container

There is also the option to run a lsf job on the HPC cluster by using the following command on an HPC node:

```
bsub < lsf-script.txt
```

## Folder content

This folder contains scripts designed to:

(a) compute differential expression (DE) results for each cell type in the cohort,
(b) generate volcano plots to visualize and summarize DE patterns, and
(c) perform Gene Ontology (GO) enrichment analysis using clusterProfiler to interpret the potential biological functions associated with the DE gene sets.

Together, these steps provide a complete workflow for identifying, visualizing, and functionally characterizing transcriptional changes at the cell‑type level.


## Folder structure 

The structure of this folder is as follows:

```
├── 01-de-analysis.Rmd
├── 02-volcano-plots.Rmd
├── 03-go-clusterprofiler
├── lsf_script.txt
├── plots
|   ├── 02-volcano-plots
|   ├── 03-go-clusterprofiler
|   ├── Report-de-analysis-<Sys.Date()>.html
|   ├── Report-de-analysis-<Sys.Date()>.pdf
|   ├── Report-volcano-plots-<Sys.Date()>.html
|   ├── Report-volcano-plots-<Sys.Date()>.pdf
|   ├── Report-go-clusterprofiler-<Sys.Date()>.html
|   └── Report-go-clusterprofiler-<Sys.Date()>.pdf
├── README.md
├── results
|   └── 01-de-analysis
├── run-de-go-analysis.R
├── run-de-go-analysis.sh
└── util
|   ├── function-create-volcano-plot.R
|___└── function-go-clusterprofiler-markers.R
```

