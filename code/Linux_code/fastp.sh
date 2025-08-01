#!/bin/bash
#SBATCH -o /mnt/raid6/gengmingyan/PRJNA1083304/00_slurm_fastp/slurm.%x.%N.%j.out
#SBATCH -e /mnt/raid6/gengmingyan/PRJNA1083304/00_slurm_fastp/slurm.%x.%N.%j.err
#SBATCH --qos=normal
#SBATCH -J fastp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --array=1-80%10
# 路径设置
INPUT_DIR="/mnt/raid6/gengmingyan/PRJNA1083304/00_rowdata"
CLEAN_DIR="/mnt/raid6/gengmingyan/PRJNA1083304/02_fastp_result"
REPORT_DIR="/mnt/raid6/gengmingyan/PRJNA1083304/01_fastp_report"
FASTP_BIN="/mnt/raid6/gengmingyan/software/fastp"

# 获取样本列表
SAMPLES=($(ls ${INPUT_DIR}/*_1.fastq.gz | sed 's/_1.fastq.gz//' | xargs -n 1 basename))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}

# 输入输出文件定义
R1="${INPUT_DIR}/${SAMPLE}_1.fastq.gz"
R2="${INPUT_DIR}/${SAMPLE}_2.fastq.gz"
OUT_R1="${CLEAN_DIR}/clean_${SAMPLE}_1.fq.gz"
OUT_R2="${CLEAN_DIR}/clean_${SAMPLE}_2.fq.gz"
REPORT="${REPORT_DIR}/${SAMPLE}_fastp.html"

# 确保输出目录存在
mkdir -p ${CLEAN_DIR}
mkdir -p ${REPORT_DIR}

# 运行 fastp
${FASTP_BIN} -i ${R1} -I ${R2} -o ${OUT_R1} -O ${OUT_R2} -z 4 -q 20 -u 30 -n 5 -h ${REPORT} -w 8
