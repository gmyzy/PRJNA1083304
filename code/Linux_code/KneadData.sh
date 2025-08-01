#!/bin/bash
#SBATCH -o /mnt/raid6/gengmingyan/PRJNA1083304/00_slurm_KneadData/slurm.%x.%N.%j.out
#SBATCH -e /mnt/raid6/gengmingyan/PRJNA1083304/00_slurm_KneadData/slurm.%x.%N.%j.err
#SBATCH --qos=normal
#SBATCH -J knead
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --array=1-80%10

source /mnt/raid6/gengmingyan/software/miniconda3/bin/activate
conda activate kneaddata

# 路径设置
INPUT_DIR="/mnt/raid6/gengmingyan/PRJNA1083304/02_fastp_result"
OUTPUT_DIR="/mnt/raid6/gengmingyan/PRJNA1083304/03_Kneaddata_result"
DB="/mnt/raid6/gengmingyan/databases/GRCh38"  # 确认为 KneadData 格式的数据库

# 获取样本列表
SAMPLES=($(ls ${INPUT_DIR}/*_1.fq.gz | sed 's/_1.fq.gz//' | xargs -n 1 basename))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}

R1="${INPUT_DIR}/${SAMPLE}_1.fq.gz"
R2="${INPUT_DIR}/${SAMPLE}_2.fq.gz"

# 输出目录准备
mkdir -p ${OUTPUT_DIR}

# 运行 KneadData
kneaddata \
  -i1 ${R1} -i2 ${R2} \
  -o ${OUTPUT_DIR} \
  -db ${DB} \
  --threads 8 \
  --reorder \
  --remove-intermediate-output \
  --output-prefix ${SAMPLE} \
  --trimmomatic /mnt/raid6/gengmingyan/software/Trimmomatic-0.39/ \
	  --trimmomatic-options 'ILLUMINACLIP:/mnt/raid6/gengmingyan/software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:40:15 SLIDINGWINDOW:4:20 MINLEN:50' \
  --bowtie2-options "--end-to-end --very-sensitive --phred33" --bypass-trf \
  --log ${OUTPUT_DIR}/${SAMPLE}_kneaddata.log

