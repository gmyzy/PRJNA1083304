#!/bin/bash
#SBATCH -o /mnt/raid6/gengmingyan/PRJNA1083304/00_slurm_metaphlan/slurm.%x.%N.%j.out
#SBATCH -e /mnt/raid6/gengmingyan/PRJNA1083304/00_slurm_metaphlan/slurm.%x.%N.%j.err
#SBATCH --qos=normal
#SBATCH -J metaphlan
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --array=1-80%10
. /mnt/raid8/jiayingzhu/software/conda/miniconda3/bin/activate  biobakery

file=$(ls /mnt/raid6/gengmingyan/PRJNA1083304/03_Kneaddata_result/*_paired_1.fastq.gz | sed -n "${SLURM_ARRAY_TASK_ID}p")

sample_name=${file%_paired_1.fastq.gz}

metaphlan  --input_type fastq --nproc 8 $file,${sample_name}_paired_2.fastq.gz \
	--bowtie2db /mnt/raid6/limin/biosoft/database/metaphlan4.1_db \
	--index mpa_vJun23_CHOCOPhlAnSGB_202307  \
	--output_file ${sample_name}_metaphlan.txt \
	--bowtie2out ${sample_name}.bz2


