#!/bin/bash

#BSUB -P run_cellranger
#BSUB -J submit-multiple-jobs
#BSUB -q standard
#BSUB -n 1
#BSUB -R "rusage[mem=2500]"
#BSUB -o submit-multiple-jobs.out
#BSUB -e submit-multiple-jobs.err
#BSUB -cwd "./sc-rna-seq-snap/analyses/cellranger-analysis"

########################################################################
# Set up running directory
prefix="./sc-rna-seq-snap/analyses/cellranger-analysis"
cd "$(dirname "${BASH_SOURCE[0]}")" 

########################################################################
# Create directories to save output files to
mkdir -p ${prefix}/input
mkdir -p ${prefix}/results

########################################################################
# File in which we store the output text to verify the job execution order.
# As the jobs will run on the cluster, we pass in a path to the directory from
# which this launch script was run so everyone writes to the same file.
output_file="output.txt"

# Prints a simple LSF job file to standard output
#function print_job {
#    outfile=${1}
#    delay_seconds=${2}
#
#    # Preamble
#    echo "#!/bin/bash"
#
#    # Variables
#    echo "time=\`date\`"
#    echo "name=\${LSB_JOBNAME}"
#
#    # Describe job when it comes online, delay, then write completion text
#    echo "echo \"Online: name=\${name} time=\${time}\" >> ${outfile}"
#    echo "sleep ${delay_seconds}"
#    echo "echo \"Done: name=\${name}\" >> ${outfile}"
#}
#
## Generate the LSF job files
#for ((i=1; i<=2; i++)); do
#    jobname="j${i}"
#    print_job "${prefix}/${output_file}" 10 > "${prefix}/${jobname}.bsub"
#    echo "bash ${prefix}/${jobname}.sh" >> "${prefix}/${jobname}.bsub"
#done
#
# Submit job 1 - note no dependencies!
#bsub -P run_CellRanger -q standard -n 1 -R "rusage[mem=2GB]" -R "span[hosts=1]" -J j1 -o ${prefix}/j1.out -e ${prefix}/j1.err "bash ${prefix}/j1.bsub"

# Job 2 depend on the successful completion of job 1
#sleep 60
#bsub < ${prefix}/waiter.sh

# define the delay-seconds
delay_seconds=10

# let's do the job 1
# because the job 1 is not submitted as lsf job anymore
# so the job name is just j1.sh
outfile=${prefix}/${output_file}
time=$(date)
name="j1.sh"

# here let's run the job 1
echo "Online: name=${name} time=${time}" >> ${outfile}
bash ${prefix}/j1.sh 2> ${prefix}/j1.err 1> ${prefix}/j1.out

# Describe job when it comes online, delay, then write completion text
sleep ${delay_seconds}
time=$(date)
echo "Done: name=${name} time=${time}" >> ${outfile}

# sleep 60 seconds and check the completion of job 1
sleep 60
bash ${prefix}/waiter.sh