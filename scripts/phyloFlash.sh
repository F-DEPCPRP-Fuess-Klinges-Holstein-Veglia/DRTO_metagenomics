#!/bin/bash

#SBATCH -N 1            # number of nodes
#SBATCH -c 32            # number of cores 
#SBATCH -t 3-01:00:00   # time in d-hh:mm:ss
#SBATCH -p general      # partition 
#SBATCH -q public      # QOS
#SBATCH -o phyloflash.%j.out # file to save job's STDOUT (%j = JobId)
#SBATCH -e phyloflash.%j.err # file to save job's STDERR (%j = JobId)
#SBATCH --mail-type=ALL # Send an e-mail when a job starts, stops, or fails
#SBATCH --export=NONE   # Purge the job-submitting shell environment

source activate phyloFlash

while IFS= read -r id; do
    echo "Processing ID: $id"
    mkdir pf2_"${id}"/
    cd pf2_"${id}"/
phyloFlash.pl -CPUs 32 -lib pf_"${id}" -read1 ~/Klinges/DRTO_meta_DNA/bowtie/final_fqs/"${id}"_final_1.fastq -read2 ~/Klinges/DRTO_meta_DNA/bowtie/final_fqs/"${id}"_final_2.fastq -readlength 150 -almosteverything -dbhome ~/Klinges/databases/138.1
    cd ..
done < ofav_ids.txt


#the command that worked best:
# phyloFlash.pl -CPUs 20 -lib pf4_DT28_CNAT_2H_Sept_2021 -read1 ~/Klinges/DRTO_meta_DNA/bowtie/final_fqs/DT28_CNAT_2H_Sept_2021_final_1.fastq -read2 ~/Klinges/DRTO_meta_DNA/bowtie/final_fqs/DT28_CNAT_2H_Sept_2021_final_2.fastq -readlength 150 -almosteverything -dbhome ~/Klinges/databases/138.1


# conda activate phyloFlash 

# phyloFlash.pl -CPUs 16 -lib pf3_DT28_CNAT_2H_Sept_2021 -read1 ~/Klinges/DRTO_meta_DNA/clean/CNAT/DT28_CNAT_2H_Sept_2021_clean_1.fastq -read2 ~/Klinges/DRTO_meta_DNA/clean/CNAT/DT28_CNAT_2H_Sept_2021_clean_2.fastq -readlength 150 -almosteverything -dbhome ~/Klinges/databases/138.1 -sortmerna

# phyloFlash.pl -CPUs 20 -lib pf3_DT28_CNAT_2H_Sept_2021 -read1 ~/Klinges/DRTO_meta_DNA/bowtie/final_fqs/DT28_CNAT_2H_Sept_2021_final_1.fastq -read2 ~/Klinges/DRTO_meta_DNA/bowtie/final_fqs/DT28_CNAT_2H_Sept_2021_final_2.fastq -readlength 150 -almosteverything -dbhome ~/Klinges/databases/138.1 -sortmerna


# sortmerna --ref ~/Klinges/databases/rRNA_databases_v4/smr_v4.3_default_db.fasta --reads ~/Klinges/DRTO_meta_DNA/clean/CNAT/DT25_CNAT_1H_Sept_2021_clean_1.fastq --reads ~/Klinges/DRTO_meta_DNA/clean/CNAT/DT25_CNAT_1H_Sept_2021_clean_2.fastq -sam -fastx -blast 1 -num_alignments 1 -v -aligned DT25_CNAT_1H_Sept_2021_SMR
