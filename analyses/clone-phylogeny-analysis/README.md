# Pipeline for Clone phylogeny analysis tool for sc-/sn-RNA-Seq Analysis in 10X Genomics data

## Usage

`run-clone-phylogeny-analysis.sh` is designed to be run as if it was called from this module directory even when called from outside of this directory.

Parameters according to the project and analysis strategy will need to be specified in the following scripts:
- `project_parameters.Config.yaml` located at the `root_dir`.

### Run module on an interactive session on HPC within the container

To run the script on an interactive session on HPC, please run the following command from an interactive compute node (while within the container):

```
bash run-clone-phylogeny-analysis.sh
```

### Run module by using lsf on HPC with the container

There is also the option to run a lsf job on the HPC cluster by using the following command on an HPC node:

```
bsub < lsf-script.txt
```

## Folder content

This folder contains scripts tasked to run a clone phylogeny analysis for sc-/sn-RNA-Seq Analysis in 10X Genomics data. 

## Steps and methods

The `clone phylogeny analysis` module utilizes Numbat method developed by [Gao et al., 2022](https://www.nature.com/articles/s41587-022-01468-y). Numbat is a haplotype-aware CNV caller from single-cell and spatial transcriptomics data. It integrates signals from gene expression, allelic ratio, and population-derived haplotype information to accurately infer allele-specific CNVs in single cells and reconstruct their lineage relationship. Numbat can be used to: 
(1) Detect allele-specific copy number variations from scRNA-seq and spatial transcriptomics; 
(2) Differentiate tumor versus normal cells in the tumor microenvironment; and  
(3) Infer the clonal architecture and evolutionary history of profiled tumors.

For more information, please see:

- [numbat 1.4.2](https://github.com/kharchenkolab/numbat)
- [Numbat User Guide](https://kharchenkolab.github.io/numbat/)
- [Eagle v2.4.1 User Manual](https://alkesgroup.broadinstitute.org/Eagle/)
- [Human tutorial](https://kharchenkolab.github.io/numbat/articles/numbat.html)
- [Mouse tutorial](https://kharchenkolab.github.io/numbat/articles/mouse.html)


## Folder structure 

The structure of this folder is as follows:

```
├── 01-create-count-mat.R
├── 02-prepare-files-for-pileup-and-phase.sh
├── 03-run-pileup-and-phase.sh
├── 04-create-df-allele-rda.R
├── 05-run-numbat.R
├── 06A-create-numbat-plots.Rmd
├── 06B-create-numbat-plots-multiple-samples.R
├── input
├── lsf-script.txt
├── README.md
├── plots
|   └── 05-create-numbat-plots
├── results
|   ├── 01-create-count-mat
|   ├── 02-create-pileup-and-phase
|   ├── 03-create-df-allele-rda
|   ├── 04-run-numbat
|   └── 05-create-numbat-plots
└── run-clone-phylogeny-analysis.sh
```

