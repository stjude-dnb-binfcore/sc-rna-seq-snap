#!/bin/bash
# 
#BSUB -P project
#BSUB -J 03-create-pileup-and-phase
#BSUB -oo job.out -eo job.err
#BSUB -n 24
#BSUB -R "rusage[mem=128GB] span[hosts=1]"
#BSUB -cwd "."

set -e
set -o pipefail

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

#######################################################
# Read multiple values and assign them to variables by parsing yaml file
root_dir=$(cat ../../project_parameters.Config.yaml | grep 'root_dir:' | awk '{print $2}')
root_dir=${root_dir//\"/}  # Removes all double quotes
echo "${root_dir}"  # Output: This is a string with quotes.

module_dir=${root_dir}/analyses/clone-phylogeny-analysis
echo "${module_dir}"

data_dir=${root_dir}/analyses/cellranger-analysis/results/02_cellranger_count/DefaultParameters
echo "${data_dir}"

########################################################################
# Directories and paths to file Inputs/Outputs
mkdir -p ./results/02-prepare-files-for-pileup-and-phase
mkdir -p ./results/03-create-pileup-and-phase

input_dir=${module_dir}/results/02-prepare-files-for-pileup-and-phase
#echo "${input_dir}"

mkdir -p ${input_dir}/sample_barcode

# Define array with samples
mapfile -t sample < <(grep '^ *-' ../../project_parameters.Config.yaml | sed 's/^ *- *//')

##################
# Check if the variable is an array by testing if it has elements
if [ "${#sample[@]}" -gt 0 ]; then
  echo "sample is an array with ${#sample[@]} elements."
else
  echo "sample is not an array or is empty."
fi
echo "${sample}"  # Output: This is a string with quotes.

# Output the array elements
echo "${sample[@]}"  # Output: sample1 sample2 etc

# Run loop for each sample
for i in "${!sample[@]}"; do
  
  echo "Sample: ${sample[i]}. Processing..."
  
  sample_data_dir=${data_dir}/${sample[i]}
  echo "$sample_data_dir"
  
  mkdir -p ${module_dir}/results/03-create-pileup-and-phase/${sample[i]}
  mkdir -p ${input_dir}/sample_barcode/${sample[i]}
  
  # Step 2
  # To copy and unzip the barcodes file from the alignment folder for each sample
  cd ..
    
  cp ${sample_data_dir}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz ${input_dir}/sample_barcode/${sample[i]}
  gunzip ${input_dir}/sample_barcode/${sample[i]}/barcodes.tsv.gz
  
  # Run pileup_and_phase.R (installed in the container)
  Rscript --vanilla /numbat/inst/bin/pileup_and_phase.R \
          --label ${sample[i]} \
          --samples ${sample[i]} \
          --bams ${sample_data_dir}/outs/possorted_genome_bam.bam \
          --barcodes ${input_dir}/sample_barcode/${sample[i]}/barcodes.tsv \
          --outdir ${module_dir}/results/03-create-pileup-and-phase/${sample[i]} \
          --gmap ${input_dir}/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
          --snpvcf ${input_dir}/references/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
          --paneldir ${input_dir}/references/1000G_hg38 \
          --ncores 64 \
          --UMItag Auto
            
  echo "Sample: ${sample[i]}. Processing completed succefully."
done
