#!/bin/bash

# 创建目标目录（如果不存在）
mkdir -p 03_cat_fq

# 循环遍历源目录中的文件
for file in 03_Kneaddata_result/*_paired_1.fastq.gz; do
  # 提取文件名
  filename=$(basename "$file")
  
  # 提取文件的标识符部分
  identifier=$(echo "$filename" | cut -d "_" -f 1,2)
  
  # 合并两个文件，并将结果保存到目标目录中
  cat "$file" "03_Kneaddata_result/${identifier}_paired_2.fastq.gz" > "03_cat_fq/${identifier}_cat.fq.gz"
done
