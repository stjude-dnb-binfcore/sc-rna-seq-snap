# Pipeline for generating an R shiny app for the project

## Usage

`run-rshiny-app.sh` is designed to be run as if it was called from this module directory even when called from outside of this directory.

Parameters according to the project and analysis strategy will need to be specified in the following scripts:
- `project_parameters.Config.yaml` located at the `root_dir`.

Optional **de-go-analysis** integration in the Shiny app: set `include_de_go_analysis_rshiny_app_module: "YES"`.

- **DE & GO reports** (top-level tab): if HTML reports exist under `analyses/de-go-analysis/plots/plots_<dego_condition_comparison_mode>/`, they are copied to `results/shinyApp/www/de_go_reports/` and listed in a `DE & GO reports` tab (`ui.R`). Manifest: `results/de_go_shiny_manifest.tsv`.
- **Volcano (DE)** (under **Annotations**, next to Proportion and Bubbleplot / Heatmap): precomputed marker tables from `analyses/de-go-analysis/results/results_<dego_condition_comparison_mode>/01-de-analysis/` (`Res_<resolution>_Markers_*.tsv` from `run-de-computation`) are bundled into `results/shinyApp/de_go_volcano_index.rds` and a reactive tab is added to `ui.R` / `server.R` (ggplot2 scatter, brush-to-subset table, optional cell-type filter). Manifest: `results/de_go_volcano_shiny_manifest.tsv`.
- **Static volcano HTML** (optional): if `Report-*volcano*.html` exists in the de-go plots folder, copies go to `results/shinyApp/www/volcano_reports/` for offline viewing; they are not required for the interactive tab.

You need DE TSVs (and usually knitted HTML) from `run-de-go-analysis` / `run-de-computation`; then re-run this Shiny module so `ui.R` / `server.R` are patched after `makeShinyCodesMulti`.


### Run module on an interactive session on HPC within the container

To run the module on an interactive session on HPC, please run the following command from an interactive compute node:

```
bash run-rshiny-app.sh
```

### Run module by using lsf on HPC with the container

There is also the option to run a lsf job on the HPC cluster by using the following command on an HPC node:

```
bsub < lsf-script.txt
```

## Launch the R Shiny app

To launch the R Shiny app, navigate to the script located at `./results/shinyApp/server.R`. Open it in RStudio on your local machine, and click the `Run App` button in the top-right corner to start the application.


## Folder content

This folder contains scripts designed to run and generate an R shiny app for the project.


## Folder structure 

The structure of this folder is as follows:

```
├── 01-generate-rshiny-app.R
├── lsf-script.txt
├── README.md
├── run-rshiny-app.sh
└── util
    ├── makeShinyFiles_assay.R
    └── function-attach-dego-shiny-app.R  # DE/GO + Volcano tabs and www copies
```
