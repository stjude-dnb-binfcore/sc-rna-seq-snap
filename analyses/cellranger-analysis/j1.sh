#!/bin/bash

set -e
set -o pipefail

########################################################################
# Set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

########################################################################
# Read multiple values and assign them to variables by parsing yaml file
metadata_dir=$(cat ../../project_parameters.Config.yaml | grep 'metadata_dir:' | awk '{print $2}')
metadata_dir=${metadata_dir//\"/}  # Removes all double quotes
echo "$metadata_dir"  # Output: This is a string with quotes.

genome_reference_path=$(cat ../../project_parameters.Config.yaml | grep 'genome_reference_path:' | awk '{print $2}')
genome_reference_path=${genome_reference_path//\"/}  # Removes all double quotes
echo "$genome_reference_path"  # Output: This is a string with quotes.

cellranger_parameters=$(cat ../../project_parameters.Config.yaml | grep 'cellranger_parameters:' | awk '{print $2}')
cellranger_parameters=${cellranger_parameters//\"/}  # Removes all double quotes
echo "$cellranger_parameters"  # Output: This is a string with quotes.

########################################################################
# Read root path
#rootdir=$(realpath "./../..")
#echo "$rootdir"

########################################################################
# Create directories to save output files to
mkdir -p ./results/01_logs
mkdir -p ./results/02_cellranger_count
mkdir -p ./results/02_cellranger_count/${cellranger_parameters}

########################################################################
# If your `project_metadata` is not in `*.txt` file format
# use the following code line to convert it
cat "${metadata_dir}"/project_metadata.tsv | sed 's/,/\t/g' > ./input/project_metadata.txt

########################################################################
# Run CellRanger for all libraries
# let's generate the shell scripts for later submission 
singularity exec ${rootdir}/rstudio_4.4.0_seurat_4.4.0_latest.sif python3 ./util/run_cellranger.py --file=./input/project_metadata.txt \
                                 --transcriptome=${genome_reference_path}  \
                                 --create_bam=true \
                                 --output_dir=./results/02_cellranger_count/${cellranger_parameters}/
                                 # --force_cells=8000  

########################################################################
# finally launch the shell scripts
# these are the bsub jobs
for lsf_sh_file in cellranger_job1*.sh; do
   bash $lsf_sh_file
	#echo $lsf_sh_file
done
