library(phyloseq)
library(tidyverse)
library(ggpubr)
library(vegan)
library(Maaslin2)
library(ggplot2)
library(ggprism)
library(ggrepel)
library(VennDiagram)
library(grid)

# --------------------------
# 1. 数据准备 ####
# --------------------------

# 读取MetaPhlAn物种丰度表 #########
species_PRJNA1083304 <- read.table("04_metaphlan_result/06_metaphlan_species.txt", 
                                   header=T, row.names=1, sep="\t") 
colnames(species_PRJNA1083304)
# 修改列名
colnames(species_PRJNA1083304) <- gsub("clean_(SRR[0-9]+)_metaphlan", "\\1", 
                                       colnames(species_PRJNA1083304))

# 读取分组信息
meta_PRJNA1083304 <- read.table("metadata.txt", 
                                header=T, row.names=1, sep="\t")

# 确保样本顺序一致
identical(colnames(species_PRJNA1083304),rownames(meta_PRJNA1083304))


########----------Maaslin-----------########
# --------------------------
# 2. 整理Maaslin 格式所需数据 ####
# --------------------------


species_PRJNA1083304 <- species_PRJNA1083304 %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(.,var = "ID")
meta_PRJNA1083304 <- rownames_to_column(meta_PRJNA1083304,var = "sample")


# --------------------------
# 3. Maaslin 运行结果 ####
# --------------------------


maas <- Maaslin2(
  input_data = species_PRJNA1083304,
  input_metadata = meta_PRJNA1083304,
  output = "maaslin_output",
  min_abundance = 0.0, #最小的丰度阈值
  min_prevalence = 0.1, #特征在样本中出现的最小百分比
  min_variance = 0.0, #保留方差大于某个值的特征，用于筛选出具有足够变异性的特征
  normalization = "TSS", #数据标准化方法 
  transform = "LOG",#数据进行转换的方法
  analysis_method = "LM", #分析方法 
  max_significance = 0.25, #q值显著性的阈值，用于多重检验校正中的假阳性率控制
  random_effects = NULL, #随机效应变量
  fixed_effects = c("Group"),  # 🔥 只比较 Control vs Case
  correction = "BH", #q值时使用的校正方法
  standardize = TRUE, 
  cores = 1,
  plot_heatmap = TRUE,
  heatmap_first_n = 50,
  plot_scatter = TRUE
)


maas#结果
#coef: 线性模型的系数
#stderr: 系数的标准误差
#pval: p值
#qval: q值，经过多重假设校正后的p值
data <- maas$results
write.csv(data, "MaAsLin_result_sig.csv", row.names = FALSE)

head(data)

data <- data %>%
  mutate(
    logP = -log10(pval),
    regulation = case_when(
      pval < 0.05 & coef > 2 ~ "Enriched in Control",  # Green
      pval < 0.05 & coef < -2 ~ "Enriched in Case",    # Red
      TRUE ~ "Not significant"
    ),
    label = ifelse(regulation != "Not significant", feature, NA)
  )


# --------------------------
# 4. Maaslin 火山图绘制 ####
# --------------------------
ggplot(data, aes(x = coef, y = logP, colour = regulation)) +
  geom_point(alpha = 0.85, size = 2) +
  scale_color_manual(values = c("Enriched in Control" = "#4daf4a", 
                                "Enriched in Case" = "#e41a1c", 
                                "Not significant" = "grey70")) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_label_repel(aes(label = label),
                   size = 3,
                   box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.8, "lines"),
                   segment.color = "black",
                   show.legend = FALSE,
                   max.overlaps = Inf) +
  labs(
    title = "Volcano Plot of Differential Microbes (MaAsLin2)",
    x = "Effect Size (Coefficient)",
    y = "-log10(p-value)"
  ) +
  theme_prism(border = TRUE) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank()
  )



########----------lefse-----------########
# --------------------------
# 5. 整理lefse 格式所需数据 ####
# --------------------------

# 获取 OTU 表和元数据
otu <- as.matrix(species_PRJNA1083304)
meta <- meta_PRJNA1083304
taxa <- rownames(otu)

# 构建 tax_table，只包含物种信息
taxonomy <- matrix(taxa, ncol = 1)
rownames(taxonomy) <- taxa
colnames(taxonomy) <- "Species"
tax_table_obj <- tax_table(taxonomy)

# 构建 phyloseq 对象（包含 tax_table）
ps <- phyloseq(otu_table(otu, taxa_are_rows = TRUE),
               sample_data(meta),tax_table_obj)

# --------------------------
# 6. LEfSe 运行结果 ####
# --------------------------

lefse_result <- run_lefse(
  ps,
  group = "Group",               # 分组变量名称
  multigrp_strat = FALSE,        # 单因素对比
  kw_cutoff = 0.25,          # Kruskal-Wallis 检验的 p 值阈值
  wilcoxon_cutoff = 0.25,    # Wilcoxon 检验的 p 值阈值
  strict = "1",              # 启用多重比较校正（类似 q 值）
  lda_cutoff = 2                 # LDA 阈值
)



# 导出结果
lefse_df <- marker_table(lefse_result)
df <- as(lefse_df, "data.frame")  # 强制转换数据框
df$feature <- gsub("s__","", df$feature)
write.csv(df, "LEfSe_result.csv", row.names = FALSE)



data <- df %>%
  mutate(
    logP = -log10(pvalue),
    regulation = case_when(
      pvalue < 0.05 & ef_lda > 2 & enrich_group == "Control"~ "Enriched in Control",  # Green
      pvalue < 0.05 & ef_lda > 2 & enrich_group == "Case"~ "Enriched in Case",    # Red
      TRUE ~ "Not significant"
    ),
    label = ifelse(regulation != "Not significant", feature, NA)
  )

# 将 Case 组的 LDA 分数乘以 -1，实现方向调整
data <- data %>%
  mutate(ef_lda = ifelse(enrich_group == "Case", ef_lda * -1, ef_lda))



# --------------------------
# 7. lefse 火山图绘制 ####
# --------------------------
p <- ggplot(data, aes(x = ef_lda, y = logP, colour = regulation)) +
  geom_point(alpha = 0.85, size = 2) +
  scale_color_manual(values = c("Enriched in Control" = "#4daf4a", 
                                "Enriched in Case" = "#e41a1c", 
                                "Not significant" = "grey70")) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.8) +
  labs(
    title = "Volcano Plot of Differential Microbes (LEfSe)",
    x = "Effect Size (Coefficient)",
    y = "-log10(p-value)"
  ) +
  theme_prism(border = TRUE) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank()
  )

p

# 只保留非 NA 的标签用于绘图
label_data <- data %>% filter(!is.na(label))


# 设置图形设备尺寸
options(repr.plot.width = 6, repr.plot.height = 8)  

pdf("lefse_plot2.pdf", width = 6, height = 8)  

# 带标签的绘图代码
print(p + geom_label_repel(
  data = label_data,
  aes(x = ef_lda, y = logP, label = label),
  size = 3,
  box.padding = unit(0.5, "lines"),
  point.padding = unit(0.8, "lines"),
  segment.color = "black",
  show.legend = FALSE,
  max.overlaps = Inf
))

dev.off()


# --------------------------
# 8. 比较两种方法鉴定出的差异物种 ####
# --------------------------
maas_sig <- read_csv("MaAsLin_result_sig.csv")
lefse_sig <- read_csv("LEfSe_result_sig.csv")

maas_sig_feature <- maas_sig %>% filter(regulation != "Not significant") %>% pull(feature)
lefse_sig_feature <- lefse_sig %>% filter(regulation != "Not significant") %>% pull(feature)

venn.plot <- venn.diagram(
  x = list(MaAsLin2 = maas_sig_feature, 
           LEfSe = lefse_sig_feature),
  filename = NULL,
  fill = c("#a6bddb", "#fa9fb5"),
  alpha = 0.7,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = c(-20, 20),
  cat.dist = 0.05
)

grid.draw(venn.plot)

# 自动添加标签

annotate_venn <- function(maas_only, lefse_only, overlap,
                          venn_object,
                          fontsize = 7,
                          spacing = 0.04,
                          maas_pos = c(0.9, 0.6),
                          lefse_pos = c(0.24, 0.8),
                          overlap_pos = c(0.6, 0.85)) {
  
  draw_text_block <- function(items, x, y, fontsize, spacing) {
    for (i in seq_along(items)) {
      grid.text(items[i],
                x = x,
                y = y - (i - 1) * spacing,
                gp = gpar(fontsize = fontsize),
                just = "center")
    }
  }
  
  grid.draw(venn_object)
  
  draw_text_block(maas_only, x = maas_pos[1], y = maas_pos[2], fontsize, spacing)
  draw_text_block(lefse_only, x = lefse_pos[1], y = lefse_pos[2], fontsize, spacing)
  draw_text_block(overlap, x = overlap_pos[1], y = overlap_pos[2], fontsize, spacing)
}




venn.plot <- venn.diagram(
  x = list(MaAsLin2 = maas_sig_feature, 
           LEfSe = lefse_sig_feature),
  filename = NULL,
  fill = c("#a6bddb", "#fa9fb5"),
  alpha = 0.7,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = c(-20, 20),
  cat.dist = 0.05
)

annotate_venn(
  maas_only = setdiff(maas_sig_feature, lefse_sig_feature),
  lefse_only = setdiff(lefse_sig_feature, maas_sig_feature),
  overlap = intersect(maas_sig_feature, lefse_sig_feature),
  venn_object = venn.plot,
  fontsize = 6.5,
  spacing = 0.035
)








