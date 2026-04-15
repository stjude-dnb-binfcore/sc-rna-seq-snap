#!/bin/bash

set -e
set -o pipefail

# Activate virtual environment
python -m venv py39

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

################################################################################################################
# Run module
#Rscript --vanilla run-cell-contamination-removal-analysis-steps-1-2-3-skip-step2.R
Rscript --vanilla run-cell-contamination-removal-analysis-steps-1-2-3.R
Rscript --vanilla run-cell-contamination-removal-analysis-steps-4.R
