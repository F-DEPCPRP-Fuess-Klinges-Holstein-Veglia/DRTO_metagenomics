#!/bin/bash

#SBATCH -N 1            # number of nodes
#SBATCH -c 15            # number of cores 
#SBATCH -t 2-01:00:00   # time in d-hh:mm:ss
#SBATCH -p general      # partition 
#SBATCH -q public      # QOS
#SBATCH -o slurm.%j.out # file to save job's STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err # file to save job's STDERR (%j = JobId)
#SBATCH --mail-type=ALL # Send an e-mail when a job starts, stops, or fails
#SBATCH --export=NONE   # Purge the job-submitting shell environment

while IFS= read -r id; do
    echo "Processing ID: $id"
    ls "${id}"_merged_R1.fastq
    ls "${id}"_merged_R2.fastq
    fastp -i "${id}"_merged_R1.fastq -I "${id}"_merged_R2.fastq -o clean/"${id}"_clean_R1.fq -O clean/"${id}"_clean_R2.fq --detect_adapter_for_pe --average_qual 25 -q 20 -l 50 -y -Y 30 -g -x -n 2 -c --overrepresentation_analysis --thread 30
    echo "Cleaned reads for ID: $id"
done < ids_clean_all.txt
