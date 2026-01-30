# Pipeline for Gene Ontology (GO) enrichment analysis using clusterProfiler per cell type for sc-RNA-Seq Analysis in 10X Genomics data


## Usage

`run-gene-ontology-go-enrichment-analysis.sh` is designed to be run as if it was called from this module directory even when called from outside of this directory.

Parameters according to the project and analysis strategy will need to be specified in the following scripts:
- `project_parameters.Config.yaml` located at the `root_dir`
- `future_globals_value` is hard-coded in each `run-gene-ontology-go-enrichment-analysis.R` script. If necessary, user can increase/decrease resources.


### Run module on an interactive session on HPC within the container

To run all of the Rscripts in this module sequentially on an interactive session on HPC, please run the following command from an interactive compute node:

```
bash run-gene-ontology-go-enrichment-analysis.sh
```

### Run module by using lsf on HPC within the container

There is also the option to run a lsf job on the HPC cluster by using the following command on an HPC node:

```
bsub < lsf-script.txt
```

## Folder content

This folder contains a script tasked to perform Gene Ontology (GO) enrichment analysis using clusterProfiler per cell type.


## Folder structure 

The structure of this folder is as follows:

```
├── 01-gene-ontology-go-enrichment-analysis.Rmd
├── lsf_script.txt
├── plots
|   ├── 01-gene-ontology-go-enrichment-analysis
|   ├── Report-gene-ontology-go-enrichment-analysis-<Sys.Date()>.html
|   └── Report-gene-ontology-go-enrichment-analysis-<Sys.Date()>.pdf
├── README.md
├── results
|   └── 01-gene-ontology-go-enrichment-analysis
├── run-gene-ontology-go-enrichment-analysis.R
├── run-gene-ontology-go-enrichment-analysis.sh
└── util
|___└── run-clusterwise-go-enrichment-markers.R
```

