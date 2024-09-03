#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=1-0
#SBATCH --job-name=makestardb
#SBATCH --output=slurm_out/slurm%j_STAR.out
#SBATCH --error=slurm_out/slurm%j_STAR.err

mamba activate bgmp_star

/usr/bin/time -v STAR \
--runThreadN 8 \
--runMode genomeGenerate \
--genomeDir ./Mus_musculus.GRCm39.dna_ens112_STAR_27.11b \
--genomeFastaFiles ./Mus_musculus.GRCm39.dna.primary_assembly.fa \
--sjdbGTFfile ./Mus_musculus.GRCm39.112.gtf
