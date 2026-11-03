# Pipeline for Differential Expression, Volcano Plots, and Gene Ontology (GO) enrichment analysis Per Cell Type for sc-RNA-Seq Analysis in 10X Genomics data


## Usage

`run-de-go-analysis.sh` is designed to be run as if it was called from this module directory even when called from outside of this directory.

Parameters according to the project and analysis strategy will need to be specified in the following scripts:
- `project_parameters.Config.yaml` located at the `root_dir`.

**Condition comparison mode** (`dego_condition_comparison_mode`):

| Value | What DE computes | Report sections (01 / 02 / 03) |
|-------|------------------|--------------------------------|
| `global_interaction` | **Overall** cell types vs rest; **global interaction** (`cell_type × condition` each vs all other cells) | Overall + global interaction |
| `within_celltype` | **Overall** cell types; **within type** `ident.1` vs `ident.2` per cell type (contrast in `Res_*_Within_contrast_*.tsv`) | Overall + within cell type |
| `both` | All of the above | Full module |

**Pipeline (all modes):** `run-de-computation.R` writes TSVs → `01-de-analysis.Rmd` / volcano / GO **read** those files (script 1 does not re-run `FindAllMarkers`). Modes differ only in **which DE steps run** and **which report sections appear**, not in how the code path works.

See comments on `dego_condition_comparison_mode` in `project_parameters.Config.yaml` for full detail.

Outputs are written under mode-specific subfolders so you can change `dego_condition_comparison_mode` in YAML and re-run without overwriting a prior run:

- **Results (TSVs):** `results/results_<mode>/01-de-analysis/` (e.g. `results/results_global_interaction/01-de-analysis/`)
- **Plots and reports:** `plots/plots_<mode>/` (HTML/PDF reports and volcano/GO subfolders)

`<mode>` is `global_interaction`, `within_celltype`, or `both`.

`min_cells_per_group_dego_module` applies to **all steps and modes** (overall cell types, global interaction, and within-cell-type contrasts) in the three HTML/PDF reports. Exclusions are written to `Res_*_Low_cell_*.tsv` (overall and global interaction) and `Res_*_Skipped_celltypes_*.tsv` (within cell type), and summarized in each report.


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
|   ├── plots_global_interaction
|   |   ├── 02-volcano-plots
|   |   ├── 03-go-clusterprofiler
|   |   └── Report-*.{html,pdf}
|   ├── plots_within_celltype
|   └── plots_both
├── README.md
├── results
|   ├── results_global_interaction
|   |   └── 01-de-analysis
|   ├── results_within_celltype
|   └── results_both
├── run-de-computation.R
├── run-de-computation-by-condition.R
├── run-de-go-analysis.R
├── run-de-go-analysis.sh
└── util
|   ├── function-create-volcano-plot.R
|   ├── function-de-global-interaction.R
|   ├── function-de-within-celltype-comparison.R
|___└── function-go-clusterprofiler-markers.R
```

