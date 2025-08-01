# PRJNA1083304 - 宏基因组数据处理与分析

## 项目简介

本项目对PRJNA1083304数据集进行了完整的宏基因组数据处理和分析，研究了肠道微生物群在炎症性抑郁症中的免疫调节作用。项目包含了从原始数据质控到功能分析的完整流程。

## 项目结构

```
PRJNA1083304/
├── code/                           # 分析代码
│   ├── Linux_code/                # Linux环境下的数据处理脚本
│   │   ├── fastp.sh              # 质量控制脚本
│   │   ├── KneadData.sh          # 去宿主序列脚本
│   │   ├── metaphlan.sh          # 物种注释脚本
│   │   ├── humann.sh             # 功能注释脚本
│   │   └── cat.sh                # 数据合并脚本
│   └── R_code/                    # R语言分析脚本
│       ├── alpha_diversity.R      # α多样性分析
│       ├── beta_diversity.R       # β多样性分析
│       ├── Taxonomic_composition.R # 物种组成分析
│       ├── Analysis_of_species_of_difference.R # 差异物种分析
│       └── Functional_differences_and_pathway_enrichment_analysis.R # 功能差异分析
├── fastp_report/                  # 质量控制报告
├── intermediate_result/           # 中间结果文件
├── metadata/                      # 元数据信息
├── metaphlan_result/             # 物种注释结果
├── plot/                         # 可视化结果图
├── *.rar                         # 压缩的大文件数据
└── 项目报告.pdf                  # 分析报告文档
```

## 数据说明

### 原始数据
- 数据来源：NCBI SRA数据库 (PRJNA1083304)
- 样本数量：40个样本
- 数据类型：肠道微生物组测序数据

### 压缩文件说明
由于GitHub文件大小限制，以下大文件已压缩：
- `e5.og_annotations.rar` - eggNOG注释结果
- `eggnog_relab.rar` - eggNOG相对丰度表
- `humann_genefamilies_relab.part*.rar` - 基因家族相对丰度表（分卷压缩）
- `ko_relab.rar` - KEGG直系同源基因相对丰度表
- `humann_pathabundance_relab.rar` - 代谢通路相对丰度表

## 分析流程

### 1. 数据质控 (Quality Control)
使用fastp进行原始数据质量控制，去除低质量序列和接头序列。

```bash
bash code/Linux_code/fastp.sh
```

### 2. 去宿主序列 (Host Removal)
使用KneadData去除人类宿主序列。

```bash
bash code/Linux_code/KneadData.sh
```

### 3. 物种注释 (Taxonomic Profiling)
使用MetaPhlAn进行物种组成分析。

```bash
bash code/Linux_code/metaphlan.sh
```

### 4. 功能注释 (Functional Profiling)
使用HUMAnN进行功能基因和代谢通路分析。

```bash
bash code/Linux_code/humann.sh
```

### 5. 统计分析与可视化
使用R进行多样性分析、差异分析和可视化。

## 主要结果

### 多样性分析
- α多样性：评估样本内物种丰富度和均匀度
- β多样性：评估样本间的群落结构差异

### 物种组成分析
- 门水平和种水平的物种组成
- 差异物种鉴定（LEfSe和MaAsLin分析）

### 功能分析
- KEGG通路富集分析
- eggNOG功能分类统计
- MetaCyc代谢通路差异分析

## 使用要求

### 软件依赖
- fastp v0.20.0+
- KneadData v0.10.0+
- MetaPhlAn v4.0+
- HUMAnN v3.0+
- R v4.0+ (含以下包：vegan, ggplot2, DESeq2等)
