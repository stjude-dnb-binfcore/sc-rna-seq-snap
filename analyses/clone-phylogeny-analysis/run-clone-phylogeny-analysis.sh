#!/bin/bash

set -e
set -o pipefail

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

########################################################################
# Run module
Rscript --vanilla 01-create-count-mat.R
bash 02-prepare-files-for-pileup-and-phase.sh #This step will only run the first time of analysis. If files exist, this will be skipped.
bash 03-run-pileup-and-phase.sh
Rscript --vanilla 04-create-df-allele-rda.R
Rscript --vanilla 05-run-numbat.R
Rscript --vanilla 06B-create-numbat-plots-multiple-samples.R
