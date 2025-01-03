#!/bin/bash

#SBATCH -N 2            # number of nodes
#SBATCH -c 60            # number of cores 
#SBATCH -t 2-18:00:00   # time in d-hh:mm:ss
#SBATCH -p general      # partition 
#SBATCH -q public      # QOS
#SBATCH -o megahit_cnat2.out # file to save job's STDOUT (%j = JobId)
#SBATCH -e megahit_cnat2.err # file to save job's STDERR (%j = JobId)
#SBATCH --mail-type=ALL # Send an e-mail when a job starts, stops, or fails
#SBATCH --export=NONE   # Purge the job-submitting shell environment

#mamba init

source activate metaG_env

while IFS= read -r id; do
    echo "Processing ID: $id"
    ls /home/jklinges/Klinges/DRTO_meta_DNA/bowtie/final_fqs/"${id}"_final_1.fastq
    ls /home/jklinges/Klinges/DRTO_meta_DNA/bowtie/final_fqs/"${id}"_final_2.fastq
    megahit -1 /home/jklinges/Klinges/DRTO_meta_DNA/bowtie/final_fqs/"${id}"_final_1.fastq -2 //home/jklinges/Klinges/DRTO_meta_DNA/bowtie/final_fqs/"${id}"_final_2.fastq -o "${id}"_megahit --presets meta-large -t 60
    echo "Assembled reads for ID: $id"
done < cnat_ids_failed.txt