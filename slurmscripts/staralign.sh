#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=1-0
#SBATCH --job-name=staralign
#SBATCH --output=../slurm_out/slurm%j_staralign.out
#SBATCH --error=../slurm_out/slurm%j_staralign.err

mamba activate bgmp_star

readfile1=../2_2B_control_S2_L008/cleaned_2_2B_control_S2_L008_R1_001.fastq.gz
readfile2=../2_2B_control_S2_L008/cleaned_2_2B_control_S2_L008_R2_001.fastq.gz

readfile3=../19_3F_fox_S14_L008/cleaned_19_3F_fox_S14_L008_R1_001.fastq.gz
readfile4=../19_3F_fox_S14_L008/cleaned_19_3F_fox_S14_L008_R2_001.fastq.gz

/usr/bin/time -v STAR --runThreadN 8 --runMode alignReads \
--outFilterMultimapNmax 3 \
--outSAMunmapped Within KeepPairs \
--alignIntronMax 1000000 --alignMatesGapMax 1000000 \
--readFilesCommand zcat \
--readFilesIn $readfile1 $readfile2 \
--genomeDir ../../Mus_musculus.GRCm39.dna_ens112_STAR_27.11b \
--outFileNamePrefix ../2_2B_control_S2_L008/align_2_2B_control_S2_L008_

/usr/bin/time -v STAR --runThreadN 8 --runMode alignReads \
--outFilterMultimapNmax 3 \
--outSAMunmapped Within KeepPairs \
--alignIntronMax 1000000 --alignMatesGapMax 1000000 \
--readFilesCommand zcat \
--readFilesIn $readfile3 $readfile4 \
--genomeDir ../../Mus_musculus.GRCm39.dna_ens112_STAR_27.11b \
--outFileNamePrefix ../19_3F_fox_S14_L008/align_2_19_3F_fox_S14_L008_