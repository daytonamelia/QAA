#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --time=1-0
#SBATCH --job-name=qaa_trimmomatic
#SBATCH --output=slurm_out/slurm%j_blastp.out
#SBATCH --error=slurm_out/slurm%j_blastp.err
#SBATCH --mail-user=adayton@uoregon.edu
#SBATCH --mail-type=ALL

/usr/bin/time -v trimmomatic PE \
trimmed_2_2B_control_S2_L008_R1_001.fastq.gz \
trimmed_2_2B_control_S2_L008_R2_001.fastq.gz \
fw_paired_2_2B_control_S2_L008_R1_001.fq.gz \
fw_unpaired_2_2B_control_S2_L008_R1_001.fq.gz \
rv_paired_2_2B_control_S2_L008_R2_001.fq.gz \
rv_unpaired_2_2B_control_S2_L008_R2_001.fq.gz \
LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:35

/usr/bin/time -v trimmomatic PE \
trimmed_19_3F_fox_S14_L008_R1_001.fastq.gz \
trimmed_19_3F_fox_S14_L008_R2_001.fastq.gz \
fw_paired_19_3F_fox_S14_L008_R1_001.fq.gz \
fw_unpaired_19_3F_fox_S14_L008_R1_001.fq.gz \
rv_paired_19_3F_fox_S14_L008_R2_001.fq.gz \
rv_unpaired_19_3F_fox_S14_L008_R2_001.fq.gz \
LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:35