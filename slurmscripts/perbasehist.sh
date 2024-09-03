#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --time=1-0
#SBATCH --job-name=qaa_histograms
#SBATCH --output=slurm_out/slurm%j_blastp.out
#SBATCH --error=slurm_out/slurm%j_blastp.err
#SBATCH --mail-user=daytonamelia@gmail.com
#SBATCH --mail-type=ALL

sample1read1=/projects/bgmp/shared/2017_sequencing/demultiplexed/2_2B_control_S2_L008_R1_001.fastq.gz
sample1read2=/projects/bgmp/shared/2017_sequencing/demultiplexed/2_2B_control_S2_L008_R2_001.fastq.gz
sample2read1=/projects/bgmp/shared/2017_sequencing/demultiplexed/19_3F_fox_S14_L008_R1_001.fastq.gz
sample2read2=/projects/bgmp/shared/2017_sequencing/demultiplexed/19_3F_fox_S14_L008_R2_001.fastq.gz

/usr/bin/time -v ./perbasehist.py \
-f $sample1read1 -l 101 -u True -y 42 \
-o ./part1_plots/perbasehist.py_plots/2_2B_control_S2_L008_R1_001_perbasehist.png

/usr/bin/time -v ./perbasehist.py \
-f $sample1read2 -l 101 -u True -y 42 \
-o ./part1_plots/perbasehist.py_plots/2_2B_control_S2_L008_R2_001_perbasehist.png

/usr/bin/time -v ./perbasehist.py \
-f $sample2read1 -l 101 -u True -y 42 \
-o ./part1_plots/perbasehist.py_plots/19_3F_fox_S14_L008_R1_001_perbasehist.png

/usr/bin/time -v ./perbasehist.py \
-f $sample2read1 -l 101 -u True -y 42 \
-o ./part1_plots/perbasehist.py_plots/19_3F_fox_S14_L008_R2_001_perbasehist.png