# Pipeline for annotating cell types for sc-/sn-RNA-Seq Analysis in 10X Genomics data

## Usage

To run the script in this module from the command line sequentially, use:

```
bash run-cell-types-annotation.sh
```

`run-cell-types-annotation.sh` is designed to be run as if it was called from this module directory even when called from outside of this directory.

Parameters according to the project and analysis strategy will need to be specified in the following scripts:
- `project_parameters.Config.yaml` located at the `root_dir`
- `run-cell-types-annotation.R`: define `root_dir`, `data_dir`, and define params for: `01-cell-types-annotation-SingleR-broad.Rmd` and `02-cell-types-annotation-SingleR-fine.Rmd` in the `rmarkdown::render code chunk`


## Folder content
This folder contains a script tasked to identify cell types annotation with SingleR by using a broad and fine resolution labeling across the project.

########----------------------------------------------------------------------------------------------------------------
## TODO based on discussion on July 1st, 2024

•	Step 1: cell type annotation based on general cell atlas reference. For now, we will be using ` HumanPrimaryCellAtlasData` and `MouseRNAseqData` from the celldex R package for human and mouse data, respectively.
1.	Option with SingleR R package: reference to be by `celldex` R package, or public or in-house 
2.	Option with Labeltransfer Seurat function: reference to be by `celldex` R package, or public or in-house 
3.	Option with `Cell type signature` function (Cody)
•	Step 2: Evaluation of results `from step 1` -> Identify and remove any cells that are due to contamination -> Save filtered object
•	Step 3: Repeat `Step 1`
•	Cell atlas reference:
o	Species/tissue/seq technology
o	Publicly available 
o	Customized
o	Create a directory with all objects in there
o	Test: mouse brain sn-RNA-seq 10x 
########----------------------------------------------------------------------------------------------------------------


## Folder structure 

The structure of this folder is as follows:

```
├── 01-cell-types-annotation-SingleR-broad.Rmd
├── 02-cell-types-annotation-SingleR-fine.Rmd
├── lsf_script.txt
├── plots
|   ├── 01-cell-types-annotation-SingleR 
|       ├── 01_annotations_broad
|       └── 02_annotations_fine
|   ├── Report_cell_types_annotation_SingleR_broad_<Sys.Date()>.html
|   ├── Report_cell_types_annotation_SingleR_broad_<Sys.Date()>.pdf
|   ├── Report_cell_types_annotation_SingleR_fine_<Sys.Date()>.html
|   └── Report_cell_types_annotation_SingleR_fine_<Sys.Date()>.pdf
├── README.md
├── results
|   ├── 01-cell-types-annotation-SingleR 
|       ├── 01_annotations_broad
|       └── 02_annotations_fine
├── run-cell-types-annotation.R
├── run-cell-types-annotation.sh
└── util
|___└── function-cell-type-fractions.R
```

