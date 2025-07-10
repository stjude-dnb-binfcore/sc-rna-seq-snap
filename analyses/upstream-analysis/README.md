# Pipeline for estimating QC metrics for sc-/sn-RNA-Seq Analysis in 10X Genomics data

## Usage

`run-upstream-analysis.sh` is designed to be run as if it was called from this module directory even when called from outside of this directory.

Parameters according to the project and analysis strategy will need to be specified in the following scripts:
- `project_parameters.Config.yaml` located at the `root_dir`.
- `future.globals.maxSize` is hardwired coded in the `02A_run_seurat_qc.Rmd`, `03_run_scDblFinder.Rmd`, and `04_run_filter_object.Rmd`. If necessary, user can increase/decrease resources.


### Run module on an interactive session on HPC within the container

To run all of the Rscripts in this module sequentially on an interactive session on HPC, please run the following command from an interactive compute node (while within the container):

```
bash run-upstream-analysis.sh
```

### Run module by using lsf on HPC with the container

There is also the option to run a lsf job on the HPC cluster by using the following command on an HPC node:

```
bsub < lsf-script.txt
```


## Folder content

This folder contains scripts tasked to:
(1) Infer QC metrics and associated plots to visually explore the quality of each library of the project.
(2) Evaluate QC metrics and set filters to remove low quality cells in 10x single-cell- and single-nuclei-RNA-sequencing libraries (without cell hashing experiment).

## QC Steps and methods

The pipeline offers flexibility for users to include or exclude methods and adjust the workflow during QC based on factors such as sequence type, expected cell number, experiment type, and genome reference.

By default, the pipeline runs all methods from steps (1-3). Step (2) is mandatory for basic QC filtering, while integration of steps (1) and/or (3) is optional. This can be configured in the `project_parameters.Config.yaml` file. However, we recommend reviewing all results, as this can provide valuable insights into the overall quality of each library.


### (1) Estimating and filtering out ambient mRNA (empty droplets)

[SoupX](https://cran.r-project.org/web/packages/SoupX/vignettes/pbmcTutorial.html) profiles “the soup”, i.e., collection of cell-free mRNAs floating in the input solution. The soup looks different for each input solution and strongly resembles the expression pattern obtained by summing all the individual cells.

SoupX calculates `Cell-specific contamination fraction` (Estimate (or manually set) the contamination fraction, the fraction of UMIs originating from the background, in each cell) and infers a `corrected expression matrix` (Correct the expression of each cell using the ambient mRNA expression profile and estimated contamination).

The CellRanger output from the `cellranger-analysis` module is used for this step.
 -  Contamination summary table and Cell-specific contamination fraction plot are generated.


Please note that if there is no cluster variabiltiy in the library, SoupX algorithm will fail when running `autoEstCont` by default:

```
sc <- autoEstCont(sc_raw, forceAccept = TRUE)
```

This means that the algorithm was not able to identify genes that are very specific to one cluster of cells. The determination of how specific is “very specific” is based on the gene’s `tf-idf` value for the cluster it is specific to. See the quickMarkers help or [this](https://constantamateur.github.io/2020-04-10-scDE/) for an explanation of what this means. The default of `tfidfMin=1` demands that genes by reasonably specific, so if the user is getting a low number of genes for estimation they can consider decreasing this value. This list is further reduced by keeping only genes that are “highly expressed” in the soup (as these give more accurate estimates of rho), where highly expressed is controlled by `soupQuantile`. The default value sounds strict, but in practice many genes with tf-idf over 1 tend to pass it.

If this is the case, the user will need to adjust as following in the `01_run_SoupX.Rmd` script:

```
sc <- autoEstCont(sc_raw, forceAccept = TRUE, tfidfMin = 0.1, soupQuantile = 0.9)
```


### (2) Seurat QC metrics

[Seurat](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) and [scooter](https://github.com/igordot/scooter) workflows are implemented to pre-process, filter and plot the RNA-sequencing data. The CellRanger output from the `cellranger-analysis` module or the corrected matrix from `step 1: SoupX` will be used for this step. User will have to define `params` as needed for their experiment. 
  - Before and after filter: Plot distribution of the number of genes, UMI, and percent mitochondrial reads per cell.
  - Summary of Cell Statistics: Percent of reads in cells, Median UMI count per cell, Median genes detected per cell, Median percent reads mitochondrial.
  - Data are normalized by using the global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result. Then, highly variable genes (HVGs) are selected to subset features that indicate high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others). Then, these HVGs are used as input to principal component analysis, and the top 30 principal components are selected. A combination of 30 dimensions and 30 neighbors are used along with the principal components to calculate the UMAP (Uniform Manifold Approximation and Projection) embeddings.
  
Here, the user can select to implement the following strategies to remove low quality cells:
- `step A` [miQC](https://bioconductor.org/packages/devel/bioc/vignettes/miQC/inst/doc/miQC.html) R package. The miQC model is based on the assumption that there are a non-trivial number of compromised cells in the dataset, which is not true in all datasets. If it is already known that the dataset is high-quality with a trivial number of compromised cells, we recommend that the user skip this step. If miQC does not identify any low quality cells, then the `step B` will automatically be used as for the filtering strategy.
- `step B` `run_QC_default` function. This is split in two filtering steps.
   - `step 1`: Filter cells with low content of genes expressed and remove mtDNA from each library (as defined in the `params`).
   - `step 2`: `Find_Outlier_Thershold` function. This is an optional step (as defined in the `params`). 

Moreover, only libraries with more than 500 cells will be kept for merging and integration purposes. This value is a commonly used threshold for many single-cell RNA-seq studies as a minimum for obtaining reliable and meaningful analysis. 
- Statistical Power: At least 500 cells are typically needed to ensure the analysis has enough statistical power to detect meaningful biological signals.
- Cell Diversity: With fewer cells, you may not capture sufficient cellular diversity, leading to incomplete or biased results.
- Clustering: Some clustering algorithms in single-cell RNA-seq require a minimum number of cells to create robust, meaningful clusters.


#### Post alignment/cell quality filtering parameters
We recommend that the user use the following parameters for initial QC, and then adjust accordingly if necessary:
- `scRNA`: min_genes = 300 (nFeature_RNA (genes detected))
           min_count = 500 (nCount_RNA (UMIs detected))
           mtDNA_pct_default = 15 (ideally 10; Percent Mitochondrial)
- `snRNA`: min_genes = 300 (nFeature_RNA (genes detected))
           min_count = 500 (nCount_RNA (UMIs detected))
           mtDNA_pct_default = 5 (ideally 1; Percent Mitochondrial)


### (3) Estimating and filtering out doublets

Popular approach of scRNAseq uses oil droplets or wells to isolate single cells along with barcoded beads. Depending on the cell density loaded, a proportion of reaction volumes (i.e. droplets or wells) will capture more than one cell, forming ‘doublets’ (or ‘multiplets’), i.e. two or more cells captured by a single reaction volume and thus sequenced as a single-cell artifact. 

The proportion of doublets is proportional to the number of cells captured. It is common in single-cell experiments to have 10-20% doublets, making accurate doublet detection critical.

Doublets are prevalent in single-cell sequencing data and can lead to artifactual findings. We will use a computational approach to calculate and remove doublets from the library. Here, we use [ScDblFinder](https://bioconductor.org/packages/devel/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html) method for identifying doublets/multiplets in single-cell data.

The `seurat_obj_raw.rds` object from step (1) is used for this step.
 -  Summary table with doublet metrics and doublets prediction plot are generated.


### (4) Merging filtered data

Next, we merge count matrices from steps (1-3) after filtering out low quality cells, ambient RNA (optional as defined in the `params`), and doublets (optional as defined in the `params`). Seurat object and metadata for the library along with UMAP embeddings are saved to be used for downstream analyses.

### (5) Final QC summary report

Lastly, we provide a final QC summary report containing graphs and summary tables across each QC step.

## Folder structure 

The structure of this folder is as follows:

```
├── 01_run_SoupX.Rmd
├── 02A_run_seurat_qc.Rmd
├── 02B_run_seurat_qc_multiple_samples.R
├── 03_run_scDblFinder.Rmd
├── 04_run_filter_object.Rmd
├── 05_run_summary_report.Rmd
├── plots
├── lsf-script.txt
├── README.md
├── results
├── run-upstream-analysis.R
├── run-upstream-analysis.sh
└── util
|   ├── function-calculate-qc-metrics.R
|   ├── function-create-UMAP.R
|   ├── function-process-Seurat.R
|___└── function-run-QC.R
```
