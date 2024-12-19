#!/bin/bash

set -e
set -o pipefail

# Activate virtual environment
python -m venv py39

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 
################################################################################################################
# Run module
Rscript --vanilla run-cluster-cell-calling-step1.R
# Rscript --vanilla run-cluster-cell-calling-step2.R
################################################################################################################
