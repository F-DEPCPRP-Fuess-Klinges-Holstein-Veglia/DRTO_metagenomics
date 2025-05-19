#!/bin/bash

#SBATCH -N 1            # number of nodes
#SBATCH -c 96            # number of cores 
#SBATCH -t 2-01:00:00   # time in d-hh:mm:ss
#SBATCH -p general      # partition 
#SBATCH -q public      # QOS
#SBATCH -o metaphlan_test.%j.out # file to save job's STDOUT (%j = JobId)
#SBATCH -e metaphlan_test.%j.err # file to save job's STDERR (%j = JobId)
#SBATCH --mail-type=ALL # Send an e-mail when a job starts, stops, or fails
#SBATCH --export=NONE   # Purge the job-submitting shell environment

source activate mpa

while IFS= read -r id; do
    echo "Processing ID: $id"
    ls ~/Klinges/DRTO_meta_DNA/megaHIT/done/"${id}"_megahit_final.contigs.fa
metaphlan ~/Klinges/DRTO_meta_DNA/megaHIT/done/"${id}"_megahit_final.contigs.fa --input_type fasta --read_min_len 200 --nproc 96 > "${id}"_profiles_metaphlan.txt 

done < mcav_ids.txt
