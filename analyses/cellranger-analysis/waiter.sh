#!/bin/bash

#BSUB -P run_cellranger
#BSUB -J waiter
#BSUB -q standard
#BSUB -n 1
#BSUB -R "rusage[mem=2500]"
#BSUB -o waiter.out
#BSUB -e waiter.err
#BSUB -cwd "."

########################################################################
# Read container path and file
rootdir=$(realpath "./../..")
echo "$rootdir" 

########################################################################
# Set up variables
prefix="${rootdir}/analyses/cellranger-analysis"
queue="standard"

########################################################################
# Function to check if there are any running jobs with the title pattern `ID.DST<some number here>`
check_jobs() {
    # Query bjobs and filter for jobs with titles matching the pattern
    bjobs_output=$(bjobs | grep "DST[0-9]*")

    # If the output is empty, there are no matching jobs
    if [ -z "$bjobs_output" ]; then
        return 0  # No matching jobs
    else
        return 1  # Matching jobs found
    fi
}
 

# Loop until no matching jobs are found
while true; do
    if check_jobs; then
        echo "$(date +"%D"): $(date +"%T"): No matching jobs found. Submitting the job..." >> waiter.log
        # Submit your job here
        bsub -P summarize_CellRanger -q ${queue} -n 1 -R "rusage[mem=2GB]" -R "span[hosts=1]" -J j2 -o ${prefix}/j2.out -e ${prefix}/j2.err ${prefix}/j2.bsub
        break
    else
        echo "$(date +"%D"): $(date +"%T"): Matching jobs found. Sleeping for 10 seconds..." >> waiter.log
        sleep 10
    fi
done