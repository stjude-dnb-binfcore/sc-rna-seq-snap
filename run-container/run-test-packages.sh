#!/bin/bash

set -e
set -o pipefail

########################################################################
pwd
echo $PATH
echo $LD_LIBRARY_PATH
which python

########################################################################
# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 
pwd

########################################################################
# Read container path and file
containerdir=$(realpath "$(dirname "${BASH_SOURCE[0]}")/..")
echo "$containerdir"


################################################################################################################
# Run other dependencies
singularity exec ${containerdir}/rstudio_4.4.0_seurat_4.4.0_latest.sif python3 --version
singularity exec ${containerdir}/rstudio_4.4.0_seurat_4.4.0_latest.sif cellranger --version
singularity exec ${containerdir}/rstudio_4.4.0_seurat_4.4.0_latest.sif fastqc --version
singularity exec ${containerdir}/rstudio_4.4.0_seurat_4.4.0_latest.sif multiqc --version
singularity exec ${containerdir}/rstudio_4.4.0_seurat_4.4.0_latest.sif pandoc --version
singularity exec ${containerdir}/rstudio_4.4.0_seurat_4.4.0_latest.sif tex --version

# Run R script
singularity exec ${containerdir}/rstudio_4.4.0_seurat_4.4.0_latest.sif \
            Rscript -e "rmarkdown::render('run-test-packages.Rmd', clean = TRUE,
                              output_dir = file.path('.'),
                              output_file = c(paste('Report-', 'run-test-packages', '-', Sys.Date(), sep = '')),
                              output_format = 'all')"
