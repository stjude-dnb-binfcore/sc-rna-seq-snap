#!/bin/bash
# 
#BSUB -P project
#BSUB -J 02-create-pileup-and-phase
#BSUB -oo job.out -eo job.err
#BSUB -n 24
#BSUB -R "rusage[mem=128GB] span[hosts=1]"
#BSUB -cwd "."

set -e
set -o pipefail

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

mkdir -p ./results
mkdir -p ./results/02-prepare-files-for-pileup-and-phase
#######################################################
# Read multiple values and assign them to variables by parsing yaml file
root_dir=$(cat ../../project_parameters.Config.yaml | grep 'root_dir:' | awk '{print $2}')
root_dir=${root_dir//\"/}  # Removes all double quotes
echo "${root_dir}"  # Output: This is a string with quotes.

module_dir=${root_dir}/analyses/clone-phylogeny-analysis
echo "${module_dir}"

input_dir=${module_dir}/results/02-prepare-files-for-pileup-and-phase
echo "${input_dir}"

references_dir=${input_dir}/references
echo "${references_dir}"


genome_name=$(cat ../../project_parameters.Config.yaml | grep 'genome_name:' | awk '{print $2}')
genome_name=${genome_name//\"/}  # Removes all double quotes
echo "$genome_name"  # Output: This is a string with quotes.

########################################################################
# Check what species the dataset is from GRCh38 
if [ "$genome_name" = "GRCh38" ]; then
    echo "Dataset is GRCh38, i.e., human genome. Processing..."
    
    ##########################################################################
    # Download input files ##################################################
    # https://kharchenkolab.github.io/numbat/articles/numbat.html#installation
    ##########################################################################
    # 1000G SNP VCF
    ######################
    file_url="https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz"
    # Extract file name from URL (you can adjust this to your needs)
    file_name=$(basename "$file_url")

    # Check if the file exists
    if [ ! -f "$references_dir/$file_name" ]; then
         echo "File snpvcf does not exist. Downloading..."
         wget -P "$references_dir" "$file_url"
         gunzip ${references_dir}/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz

    else
             echo "File snpvcf already exists. Skipping download."
    fi

    ########################################################################
    # Eagle_v2.4.1
    ######################
    # https://alkesgroup.broadinstitute.org/Eagle/#x1-30002
    file_url="https://alkesgroup.broadinstitute.org/Eagle/downloads/Eagle_v2.4.1.tar.gz"
    # Extract file name from URL (you can adjust this to your needs)
    file_name=$(basename "$file_url")

    # Check if the file exists
    if [ ! -f "$input_dir/$file_name" ]; then
         echo "Folder Eagle does not exist. Downloading..."
         wget -P "$input_dir" "$file_url" 
         
         # Extract the .tar.gz file
         tar -xvzf "$input_dir/Eagle_v2.4.1.tar.gz" -C "$input_dir"
         echo "Download and unzipping complete."
    else
             echo "Folder Eagle already exists. Skipping download."
    fi

    ########################################################################
    # 1000G Reference Panel (paste link in browser to download if wget isn’t working)
    ######################
    file_url="https://pklab.org/teng/data/1000G_hg38.zip"
    # Extract file name from URL (you can adjust this to your needs)
    file_name=$(basename "$file_url")
    
    # Check if the file exists
    if [ ! -f "$references_dir/$file_name" ]; then
        echo "File 1000G Reference Panel does not exist. Downloading..."

        # Download the file to the input/eagle directory
        wget -P "$references_dir" "$file_url"
        unzip ${references_dir}/1000G_hg38.zip -d ${references_dir}
        echo "Download and unzipping complete."
    else
        echo "File 1000G Reference Panel already exists. Skipping download."
    fi
    ########################################################################

else
    # The pipeline cannot be used yet for for mouse and dual index genomes, i.e., GRCm39, GRCh38ANDGRCm39, and GRCh38_GFP_tdTomato.
    echo "Dataset is GRCm39, i.e., mouse genome. The pipeline is currently available for human species only."
    exit 1  # This will stop the script execution if the genome is GRCm39.
fi
