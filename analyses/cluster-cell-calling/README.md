# Pipeline for cluster cell calling and gene marker analysis for sc-/sn-RNA-Seq Analysis in 10X Genomics data

## Usage

`run-cluster-cell-calling.sh` is designed to be run as if it was called from this module directory even when called from outside of this directory.

Parameters according to the project and analysis strategy will need to be specified in the following scripts:
- `project_parameters.Config.yaml` located at the `root_dir`.
- `run-cluster-cell-calling.sh`: comment in/out according to which step user wants to run, i.e., `run-cluster-cell-calling-step1.R` or `run-cluster-cell-calling-step2.R`.

### Run module on an interactive session on HPC within the container

To run all of the Rscripts in this module sequentially on an interactive session on HPC, please run the following command from an interactive compute node:

```
bash run-cluster-cell-calling.sh
```

### Run module by using lsf on HPC within the container

There is also the option to run a lsf job on the HPC cluster by using the following command on an HPC node:

```
bsub < lsf-script.txt
```

## Folder content
This folder contains a script tasked to calculate clusters and find markers for each cluster across the project.

## Analysis strategy:
We recommend the user to follow the following steps for running the current module:
- Step (1) `run-cluster-cell-calling-step1.R`: At first, the `01-cluster-cell-calling.Rmd` should be run for a set of resolutions or by default.
- Step (2) `run-cluster-cell-calling-step2.R`: After inspection of the first round of results, the single resolution that fits best the data can be provided and used to run the `02-find-markers.Rmd`.

For more information, see [Seurat clustering](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html).

In order to interpret the clustering and decide properly on the resolution, we should consider tissue type, total number of cells in the dataset, and expected known cell types in there (i.e., tissue ecosystem). That would help to determine the correct number of clusters in a biologically meaningful way, considering known cell types and avoiding splitting clusters into unstable, small ones. That can also help to explore smaller clusters that might contain unknown or disease/patient-specific clusters that are still worth considering and investigating from the research/clinical perspective. 

After inspection of the two rounds of results in the `cluster-cell-calling` module, user can remove clusters if necessary. This is relevant to projects in which we need to remove clusters from the object. This is the case, e.g., in PDX projects, there might be both human and mouse clusters identified after the `02-find-markers.Rmd` step of the the `cluster-cell-calling` module. In this case, we recommend the user to run the `contamination-remove-cells-analysis` module that allows to remove clusters, repeat normalization and integration steps. This object can then be used for cell type annotation or other type of analysis.


## Folder structure 

The structure of this folder is as follows:

```
├── 01-cluster-cell-calling.Rmd
├── 02-find-markers.Rmd
├── lsf-script.txt
├── plots
|   ├── 01_cluster_cell_calling_{resolution}
|   ├── 02_find_markers
|   ├── Report_cluster_cell_calling_{resolution}_<Sys.Date()>.html
|   ├── Report_cluster_cell_calling_{resolution}_<Sys.Date()>.pdf
|   ├── Report_find_markers_<Sys.Date()>.html
|   └── Report_find_markers_<Sys.Date()>.pdf
├── README.md
├── results
|   ├── 01_cluster_cell_calling_{resolution}
|   └── 02_find_markers
├── run-cluster-cell-calling-step1.R
├── run-cluster-cell-calling-step2.R
├── run-cluster-cell-calling.sh
└── util
|___└──function-cluster-cell-calling.R
```

