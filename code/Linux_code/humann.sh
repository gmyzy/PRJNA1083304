#!/bin/bash
#SBATCH -o /mnt/raid6/gengmingyan/PRJNA1083304/00_slurm_humann/slurm.%x.%N.%j.out
#SBATCH -e /mnt/raid6/gengmingyan/PRJNA1083304/00_slurm_humann/slurm.%x.%N.%j.err
#SBATCH --qos=normal
#SBATCH -J humann
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --array=1-80%8
. /mnt/raid2/limin/biosoft/miniconda3/bin/activate biobakery
infile=$(cat /mnt/raid6/gengmingyan/PRJNA1083304/output.txt | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print$0}')

humann  --input  /mnt/raid6/gengmingyan/PRJNA1083304/03_cat_fq/${infile}_cat.fq.gz  --output /mnt/raid6/gengmingyan/PRJNA1083304/05_humann_result/  --threads 10 --input-format fastq.gz

