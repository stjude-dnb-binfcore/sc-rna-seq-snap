#!/bin/bash

set -e
set -o pipefail

########################################################################

# Set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

# Read root path
rootdir=$(realpath "./../..")
echo "$rootdir"

########################################################################
# Read multiple values and assign them to variables by parsing yaml file
metadata_dir=$(cat ${rootdir}/project_parameters.Config.yaml | grep 'metadata_dir:' | awk '{print $2}')
metadata_dir=${metadata_dir//\"/}  # Removes all double quotes
echo "$metadata_dir"  # Output: This is a string with quotes.

genome_reference_path=$(cat ${rootdir}/project_parameters.Config.yaml | grep 'genome_reference_path:' | awk '{print $2}')
genome_reference_path=${genome_reference_path//\"/}  # Removes all double quotes
echo "$genome_reference_path"  # Output: This is a string with quotes.

cellranger_parameters=$(cat ${rootdir}/project_parameters.Config.yaml | grep 'cellranger_parameters:' | awk '{print $2}')
cellranger_parameters=${cellranger_parameters//\"/}  # Removes all double quotes
echo "$cellranger_parameters"  # Output: This is a string with quotes.

create_bam_value=$(cat ${rootdir}/project_parameters.Config.yaml | grep 'create_bam_value:' | awk '{print $2}')
create_bam_value=${create_bam_value//\"/}  # Removes all double quotes
echo "$create_bam_value"  # Output: This is a string with quotes.

metadata_file=$(cat ${rootdir}/project_parameters.Config.yaml | grep 'metadata_file:' | awk '{print $2}')
metadata_file=${metadata_file//\"/}  # Removes all double quotes
echo "Metadata file: $metadata_file"  # Output 

########################################################################
# Create directories to save output files to
module_dir=$rootdir/analyses/cellranger-analysis
echo "$module_dir"

results_dir=results
echo "$results_dir"

mkdir -p ./$results_dir/01_logs
mkdir -p ./$results_dir/02_cellranger_count
mkdir -p ./$results_dir/02_cellranger_count/${cellranger_parameters}


logs_dir=$module_dir/$results_dir/01_logs
echo "$logs_dir"

########################################################################
# Path to the TSV file
SAMPLES_FILE="$metadata_dir/$metadata_file"

# Check if the TSV file exists
if [ ! -f "$SAMPLES_FILE" ]; then
  echo "Error: TSV file not found at $SAMPLES_FILE"
  exit 1
fi

########################################################################
# Load modules
#module load python/3.9.9
#module load cellranger/8.0.1

cellranger_version=$(cat ${rootdir}/project_parameters.Config.yaml | grep 'cellranger_version:' | awk '{print $2}')
cellranger_version=${cellranger_version//\"/}  # Removes all double quotes
echo "$cellranger_version"  # Output: This is a string with quotes.

module load $cellranger_version

########################################################################
# STEP 1 ###############################################################
########################################################################
# Loop through each line of the TSV file (skipping the header)
# Using `cut` to extract the first three columns (ID, SAMPLE, FASTQ)
# This assumes your TSV is tab-separated. If it's space-separated, adjust the delimiter accordingly.
# Sort the samples by the second column (SAMPLE) alphabetically, and process
declare -A sample_map
declare -A fastq_map

# Read the file without using a pipe (to avoid a subshell)
while IFS=$'\t' read -r ID SAMPLE FASTQ; do
  # Skip empty or malformed lines
  [ -z "$ID" ] && continue
  [ "$ID" == "ID" ] && continue  # skip header just in case

  sample_map["$ID"]+="${SAMPLE},"
  fastq_map["$ID"]+="${FASTQ},"
done < <(tail -n +2 "$SAMPLES_FILE" | cut -f1-3)

# Iterate over each unique ID
for ID in "${!sample_map[@]}"; do
  SAMPLE_str="${sample_map[$ID]%,}"     # Remove trailing comma
  FASTQ_str="${fastq_map[$ID]%,}"       # Remove trailing comma

  echo "Processing Sample ID: $ID"
  echo "  Samples: $SAMPLE_str"
  echo "  FASTQ paths: $FASTQ_str"

  # Submit job to LSF
  bsub -J "RNA_${ID}" -n 6 -M 48000 -R "rusage[mem=8000]" \
    -o "${logs_dir}/${ID}.out" -e "${logs_dir}/${ID}.err" \
  "cellranger count --id=${ID} \
      --sample=${SAMPLE_str} \
      --fastqs=${FASTQ_str} \
      --transcriptome=${genome_reference_path} \
      --create-bam=${create_bam_value} \
      --output-dir=./$results_dir/02_cellranger_count/${cellranger_parameters}/${ID} \
      --localcores=6 \
      --localmem=48 \
      --jobmode=lsf"

  if [ $? -eq 0 ]; then
    echo "Cellranger count submitted successfully for ID: $ID"
  else
    echo "Error submitting cellranger count for ID: $ID"
  fi
done

################################
