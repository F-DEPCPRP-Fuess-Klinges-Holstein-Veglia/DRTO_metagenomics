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


source activate metaG_env

#bowtie2-build --threads 15 Mcavernosa_July2018.fasta mcav
#bowtie2-build --threads 15 GCA_002042975.1_ofav_dov_v1_genomic.fna.fai ofav
#bowtie2-build --threads 15 GCA_003297005.1_SymA_ver_1.0_genomic.fna sfit
#bowtie2-build --threads 15 GCA_947184155.1_Cgoreaui_SCF055-01_genomic.fna cgor
#bowtie2-build --threads 15 GCA_963969995.1_Durusdinium_trenchii_SCF082_genomic.fna dur
#bowtie2-build --threads 15 GCA_000507305.1_ASM50730v1_genomic.fna brevi
#bowtie2-build --threads 15 GCA_963693305.1_jaMeaMean2.1_genomic.fna mmea
#bowtie2-build --threads 15 /home/jklinges/Klinges/genomes/GCF_000001405.40_GRCh38.p14_genomic.fna human

while IFS= read -r id; do
    echo "Processing ID: $id"
    ls "${id}"_clean_1.fastq
    ls "${id}"_clean_2.fastq
    bowtie2 -p 15 -x ofav -1 "${id}"_clean_1.fastq -2 "${id}"_clean_2.fastq --un-conc-gz --fast-local --un-conc-gz "${id}"_ofav_rem
    echo "ofav rem for: $id"
    bowtie2 -p 15 -x mcav -1 "${id}"_ofav_rem.1 -2 "${id}"_ofav_rem.2 --un-conc-gz --fast-local --un-conc-gz "${id}"_mcav_rem
    echo "mcav rem for: $id"
    bowtie2 -p 15 -x cgor -1 "${id}"_mcav_rem.1 -2 "${id}"_mcav_rem.2 --un-conc-gz --fast-local --un-conc-gz "${id}"_cgor_rem
    echo "cladicopium rem for: $id"
	bowtie2 -p 15 -x dur -1 "${id}"_cgor_rem.1 -2 "${id}"_cgor_rem.2 --un-conc-gz --fast-local --un-conc-gz "${id}"_dur_rem
    echo "durusdinium rem for: $id"
    bowtie2 -p 15 -x sfit -1 "${id}"_dur_rem.1 -2 "${id}"_dur_rem.2 --un-conc-gz --fast-local --un-conc-gz "${id}"_sfit_rem
    echo "symbiodinium rem for: $id"
    bowtie2 -p 15 -x brevi -1 "${id}"_sfit_rem.1 -2 "${id}"_sfit_rem.2 --un-conc-gz --fast-local --un-conc-gz "${id}"_symbiont_rem
    echo "breviolum rem for: $id"
    bowtie2 -p 15 -x mmea -1 "${id}"_symbiont_rem.1 -2 "${id}"_symbiont_rem.2 --un-conc-gz --fast-local --un-conc-gz "${id}"_mmea_rem
    echo "meandrina meandrites rem for: $id"
    bowtie2 -p 15 -x human -1 "${id}"_mmea_rem.1 -2 "${id}"_mmea_rem.2 --un-conc-gz --fast-local --un-conc-gz "${id}"_human_rem
    echo "human rem for: $id"
done < ofra_ids.txt #replace with IDs for other species
