#!/bin/bash

#SBATCH -N 1            # number of nodes
#SBATCH -c 50            # number of cores 
#SBATCH -t 3-01:00:00   # time in d-hh:mm:ss
#SBATCH -p general      # partition 
#SBATCH -q public      # QOS
#SBATCH -o eggnog.%j.out # file to save job's STDOUT (%j = JobId)
#SBATCH -e eggnog.%j.err # file to save job's STDERR (%j = JobId)
#SBATCH --mail-type=ALL # Send an e-mail when a job starts, stops, or fails
#SBATCH --export=NONE   # Purge the job-submitting shell environment

source activate metaG_env

export EGGNOG_DATA_DIR=~/scratch/databases/EggNOG_V5

#tests
#emapper.py -m diamond --itype metagenome --cpu 50 --genepred prodigal -i ~/Klinges/DRTO_meta_DNA/megaHIT/done/DT25_OFAV_1H_June_2022_megahit_final.contigs.fa -o test3 --tax_scope Bacteria --excel --report_orthologs --dbmem

# while IFS= read -r id; do
#     echo "Processing ID: $id"

# emapper.py -m diamond --itype metagenome --cpu 50 --genepred prodigal -i ~/Klinges/DRTO_meta_DNA/megaHIT/done/"${id}"_megahit_final.contigs.fa -o test3 --tax_scope Bacteria --excel --report_orthologs --dbmem

# done < ofav_ids.txt
