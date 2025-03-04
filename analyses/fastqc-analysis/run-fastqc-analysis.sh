#!/bin/bash

set -e
set -o pipefail

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

########################################################################
# Read multiple values and assign them to variables by parsing yaml file
mapfile -t fastqc_dir < <(grep '^ *-' ../../project_parameters.Config.yaml | sed 's/^ *- *//')

##################
# Check if the variable is an array by testing if it has elements
if [ "${#fastqc_dir[@]}" -gt 0 ]; then
  echo "fastqc_dir is an array with ${#fastqc_dir[@]} elements."
else
  echo "fastqc_dir is not an array or is empty."
fi

##################

# Create results dir
mkdir -p results
mkdir -p results/01-fastqc-reports


################################################################################################################
###### STEP 1 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
################################################################################################################
# Run FastQC per library

# Loop through the array and echo each element
#for dir in "${fastqc_dir_array[@]}"; 
for dir in $(echo "${fastqc_dir[@]}" | tr ' ' '\n' | sort); 
  do
  echo "Processing directory: ${dir}"
  fastqc -o results/01-fastqc-reports ${dir}/*R2*.fastq.gz
done

################################################################################################################
###### STEP 2 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
################################################################################################################
# Run multiqc for all samples in the `results` dir
# to summarize results
 
cd results/01-fastqc-reports
multiqc . 

# rename folder
mv multiqc_data 02-multiqc-reports

# move the files related to multiqc in the main `results` dir
mv 02-multiqc-reports ../
mv multiqc_report.html ../

################################################################################################################
###### THE END ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
################################################################################################################