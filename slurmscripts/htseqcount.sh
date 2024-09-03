#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=1-0
#SBATCH --job-name=htseqcount
#SBATCH --output=../slurm_out/slurm%j_htseq.out
#SBATCH --error=../slurm_out/slurm%j_htseq.err
#SBATCH --mail-user=adayton@uoregon.edu
#SBATCH --mail-type=ALL

gtf_file=/projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/mouse_data/Mus_musculus.GRCm39.112.gtf
control_file=/projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA/2_2B_control_S2_L008/align_2_2B_control_S2_L008_Aligned.out.sam
fox_file=/projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA/19_3F_fox_S14_L008/align_2_19_3F_fox_S14_L008_Aligned.out.sam

conda activate QAA

/usr/bin/time -v htseq-count \
--stranded=yes \
--stranded=reverse \
$fox_file \
$gtf_file