# Pipeline for Integrative analysis for sc-/sn-RNA-Seq Analysis in 10X Genomics data

## Usage

`run-integrative-analysis.sh` is designed to be run as if it was called from this module directory even when called from outside of this directory.

Parameters according to the project and analysis strategy will need to be specified in the following scripts:
- `project_parameters.Config.yaml` located at the `root_dir`.

The `integrative-analysis` module is intended for use with multi-sample cohorts. It should be excluded when analyzing single-sample cohorts, as it relies on data integration across multiple samples.


### Run module on an interactive session on HPC within the container

To run all of the Rscripts in this module sequentially on an interactive session on HPC, please run the following command from an interactive compute node:

```
bash run-integrative-analysis.sh
```

### Run module by using lsf on HPC within the container

There is also the option to run a lsf job on the HPC cluster by using the following command on an HPC node:

```
bsub < lsf-script.txt
```


## Folder content

This folder contains a script tasked to integrate the samples across the project. 

The script named `custom-seurat-functions.R` is from [CellGen Programme Template Notebooks](https://github.com/cellgeni/notebooks).

## Folder structure 

The structure of this folder is as follows:

```
├── 01-integrative-analysis.Rmd
├── lsf_script.txt
├── plots
|   ├── plots_{integration_method}
|   └── Report-integrative-analysis-{integration_method}-<Sys.Date()>.html
├── README.md
├── results
|   ├── metadata_{integration_method}
|   └── seurat_obj_integrated_{integration_method}.rds
├── run-integrative-analysis.R
├── run-integrative-analysis.sh
└── util
|  ├── custom-seurat-functions.R
|__└── function-samples-integrate.R
```