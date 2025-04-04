#!/bin/bash

set -e
set -o pipefail

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

################################################################################################################
# Run module
#Rscript --vanilla run-cell-types-annotation-SingleR.R
Rscript --vanilla run-cell-types-annotation-gene-markers.R
Rscript --vanilla run-cell-types-annotation-reference.R
Rscript --vanilla run-merge-cell-types-annotations-all.R
################################################################################################################
