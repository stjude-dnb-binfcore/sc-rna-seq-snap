#!/bin/bash

set -e
set -o pipefail

########################################################################
# Load modules
module load python/3.9.9
module load cellranger/8.0.1
########################################################################

# Set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

# If your `project_metadata` is not in `*.txt` file format
# use the following code line to convert it
cat ./mouse-test-dataset/project_metadata.tsv | sed 's/,/\t/g' > ./input/project_metadata.txt


########################################################################
# Run CellRanger for all libraries
python ./util/run_cellranger.py --file=./input/project_metadata.txt \
                                --transcriptome=GRCm39 \
                                --create_bam=true \
                                --output_dir=./results/02_cellranger_count/ForcedCells8000Parameters/ \
                                --force_cells=8000
                                # --output_dir=./results/02_cellranger_count/DefaultParameters/
                                