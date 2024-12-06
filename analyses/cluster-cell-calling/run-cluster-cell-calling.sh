#!/bin/bash

set -e
set -o pipefail

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 
################################################################################################################
# Run module
Rscript --vanilla run-cluster-cell-calling-step1.R
# Rscript --vanilla run-cluster-cell-calling-step2.R
################################################################################################################
