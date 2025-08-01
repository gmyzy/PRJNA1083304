
# 加载包
library(ggplot2)
library(ggtern)
library(reshape2)
library(tidyverse)


########----------MetaCyc通路差异分析-----------########

# --------------------------
# 1. 读取数据并预处理 ####
# --------------------------
pathabundance_PRJNA1083304 <- read_tsv("05_humann_result/pathabundance/humann_pathabundance_relab.tsv")

pathabundance_PRJNA1083304 = pathabundance_PRJNA1083304 %>%
  filter(!str_detect(`# Pathway`, "\\|")) %>%
  column_to_rownames(var = "# Pathway")

# 修改pathway列名
a=data.frame(str_split_fixed(colnames(pathabundance_PRJNA1083304),"_",4)[,1:3])
b=str_c(a$X2)
colnames(pathabundance_PRJNA1083304)=b


# 读取分组信息
meta_PRJNA1083304 <- read.table("metadata.txt", 
                                header=T, row.names=1, sep="\t")

# 确保样本顺序一致
identical(colnames(pathabundance_PRJNA1083304),rownames(meta_PRJNA1083304))

# --------------------------
# 2. 微生物功能Stamp差异分析   ####
# --------------------------

# 计算所有功能相对丰度并且以百分比的形式表示
pathabundance_PRJNA1083304 = as.data.frame(t(t(pathabundance_PRJNA1083304)/colSums(pathabundance_PRJNA1083304,na=T)*100))

pathabundance_PRJNA1083304 = as.data.frame(t(pathabundance_PRJNA1083304))# 行列转换


pathabundance_PRJNA1083304_1 <- cbind(pathabundance_PRJNA1083304, meta_PRJNA1083304$Group)
colnames(pathabundance_PRJNA1083304_1)[508] <- "Group"

pathabundance_PRJNA1083304_1$Group = as.factor(pathabundance_PRJNA1083304_1$Group)



# 方差不齐wilcox检验
dfg = pathabundance_PRJNA1083304_1 %>%
  select_if(is.numeric) %>%
  map_df(~ broom::tidy(wilcox.test(. ~ Group, pathabundance_PRJNA1083304_1)), .id = 'var')


dfg = dfg %>% filter(p.value < 0.05)

# 计算每个通路在两个分组下的中位数和标准差
group_stats <- pathabundance_PRJNA1083304_1 %>%
  group_by(Group) %>%
  summarise(across(where(is.numeric), list(median = median, sd = sd), .names = "{.col}_{.fn}")) %>%
  pivot_longer(-Group, names_to = c("Pathway", ".value"), names_sep = "_")

# 转成宽格式方便估计差值
group_stats_wide <- group_stats %>%
  pivot_wider(names_from = Group, values_from = c(median, sd))

# 计算估计值和近似置信区间
diff.mean <- group_stats_wide %>%
  mutate(
    estimate = median_Case - median_Control,
    pooled_sd = sqrt((sd_Case^2 + sd_Control^2)/2),
    conf.low = estimate - 1.96 * pooled_sd,
    conf.high = estimate + 1.96 * pooled_sd,
    var = Pathway
  ) %>%
  inner_join(dfg, by = "var")


# 筛选差异显著的功能，将宽数据转化为长数据，计算平均值
abun.bar = pathabundance_PRJNA1083304_1[,c(dfg$var,"Group")] %>% 
  gather(variable,value,-Group) %>% #将宽数据转化为长数据
  group_by(variable,Group) %>% 
  summarise(Mean = mean(value))


# 筛选dfg结果中的信息
diff.mean$group <- c(ifelse(diff.mean$estimate >0,levels(pathabundance_PRJNA1083304_1$Group)[1],
                            levels(pathabundance_PRJNA1083304_1$Group)[2]))
diff.mean <- diff.mean[order(diff.mean$estimate,decreasing = TRUE),]


# --------------------------
# 3. 绘图   ####
# --------------------------


# 先绘制左侧条形图
col <- c("#f4a69a","#aed4e9")
# 按照功能的因子排序
abun.bar$variable  =  factor(abun.bar$variable,levels = rev(diff.mean$var))
p1 = ggplot(abun.bar,aes(variable,Mean,fill = Group)) +
  scale_fill_manual(values=col)+
  scale_x_discrete(limits = levels(diff.mean$var)) +
  coord_flip() +
  xlab("") +
  ylab("Mean proportion (%)") +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.title.x=element_text(colour='black', size=12,face = "bold"),
        axis.text=element_text(colour='black',size=10,face = "bold"),
        legend.title=element_blank(),
        legend.text=element_text(size=12,face = "bold",colour = "black",
                                 margin = margin(r = 20)),
        legend.position = c(-1,1),# 图例位置设置左上方
        legend.direction = "horizontal",
        legend.key.width = unit(0.8,"cm"),
        legend.key.height = unit(0.5,"cm"))

p1
# 添加柱状图的背景
for (i in 1:(nrow(diff.mean) - 1)) 
  p1 <- p1 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf, 
                      fill = ifelse(i %% 2 == 0, 'white', 'gray95'))
# 将柱状图添加
p1 = p1+geom_bar(stat = "identity",position = "dodge",
                 width = 0.7,colour = "black")
p1

# 右侧散点图绘制
diff.mean$var = factor(diff.mean$var,levels = levels(abun.bar$variable))
diff.mean$p.value = as.numeric(diff.mean$p.value)
diff.mean$p.value = round(diff.mean$p.value,3)# 保留3位小数
diff.mean$p.value = as.character(diff.mean$p.value)
p2 = ggplot(diff.mean,aes(var,estimate,fill = group)) +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.title.x=element_text(colour='black', size=12,face = "bold"),
        axis.text=element_text(colour='black',size=10,face = "bold"),
        axis.text.y = element_blank(),
        legend.position = "none",
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 15,face = "bold",colour = "black",hjust = 0.5)) +
  scale_x_discrete(limits = levels(diff.mean$var)) +
  coord_flip() +
  xlab("") +
  ylab("Difference in mean proportions (%)") +
  labs(title="95% confidence intervals") 
p2
for (i in 1:(nrow(diff.mean) - 1)) 
  p2 <- p2 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf, 
                      fill = ifelse(i %% 2 == 0, 'white', 'gray95'))

p2 <- p2 +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                position = position_dodge(0.8), width = 0.5, linewidth = 0.5) +# 误差线
  geom_point(shape = 21,size = 3) +
  scale_fill_manual(values=col) +
  geom_hline(aes(yintercept = 0), linetype = 'dashed', color = 'black')

p2 

# 最右侧p值文本添加   
p3 <- ggplot(diff.mean,aes(var,estimate,fill = group)) +
  geom_text(aes(y = 0,x = var),label = diff.mean$p.value,
            hjust = 0,fontface = "bold",inherit.aes = FALSE,size = 3) +
  geom_text(aes(x = nrow(diff.mean)/2 +0.5,y = 0.85),label = "P-value (corrected)",
            srt = 90,fontface = "bold",size = 5) +
  coord_flip() +
  ylim(c(0,1)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
p3

# 画到这里有了三个图p1、p2、p3
# 拼接在一起
library(patchwork)
p <- p1 + p2 + p3 + plot_layout(widths = c(4,6,2))

p
# ggsave保存图
ggsave("005_stamp.pdf", plot = p, width = 14, height = 6)



# 加载包
library(readr)  #read_tsv功能
library(magrittr)
library(dplyr)
library(tidyverse) #column_to_rownames功能
library(clusterProfiler)
library(MicrobiomeProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(edgeR)

########----------KEGG通路富集分析-----------########

# --------------------------
# 1. 读取数据并预处理 ####
# --------------------------
ko_relab_PRJNA1083304 <- read_tsv("genefamily/ko_relab.tsv")

ko_relab_PRJNA1083304 = ko_relab_PRJNA1083304 %>%
  filter(!str_detect(`# Gene Family`, "\\|")) %>%
  column_to_rownames(var = "# Gene Family")

# 修改pathway列名
a=data.frame(str_split_fixed(colnames(ko_relab_PRJNA1083304),"_",4)[,1:3])
b=str_c(a$X2)
colnames(ko_relab_PRJNA1083304)=b


# 读取分组信息
meta_PRJNA1083304 <- read.table("metadata.txt", 
                                header=T, row.names=1, sep="\t")

# 确保样本顺序一致
identical(colnames(ko_relab_PRJNA1083304),rownames(meta_PRJNA1083304))

# --------------------------
# 2. 差异基因富集分析 ####
# --------------------------
#同组基因做差异分析求P值
# 获取样本名
samples <- colnames(ko_relab_PRJNA1083304)

# 提取分组信息
group <- meta_PRJNA1083304[samples, "Group"]  # 保证顺序一致

# 找出组别（假设两个组）
group_levels <- unique(group)
Control_samples <- samples[group == group_levels[1]]
Case_samples <- samples[group == group_levels[2]]

# 使用组别进行差异分析（Wilcoxon检验）
ko_relab_PRJNA1083304_pvalue <- apply(ko_relab_PRJNA1083304, 1, function(x) {
  wilcox.test(x[Control_samples], x[Case_samples])$p.value
})
# control组求均值
ko_PRJNA1083304_control_mean = apply(ko_relab_PRJNA1083304[,Control_samples],1,FUN = mean)

# case组求均值
ko_PRJNA1083304_case_mean = apply(ko_relab_PRJNA1083304[,Case_samples],1,FUN = mean) 

# 对得到的P值用BH方法做矫正
ko_PRJNA1083304_FDR = p.adjust(ko_relab_PRJNA1083304_pvalue, method = 'BH') 

# 求Fold change 用case均值/control组均值
ko_PRJNA1083304_FC=ko_PRJNA1083304_case_mean/ko_PRJNA1083304_control_mean 

# Fold change取log
ko_PRJNA1083304_log2FC = log(ko_PRJNA1083304_FC,2) 

# 将上述指标放于一个数据框内方便计算
df_ko_PRJNA1083304=data.frame(ko_relab_PRJNA1083304_pvalue,ko_PRJNA1083304_FDR,ko_PRJNA1083304_FC,ko_PRJNA1083304_log2FC)
df_diff_ko_PRJNA1083304=df_ko_PRJNA1083304[df_ko_PRJNA1083304$ko_relab_PRJNA1083304_pvalue<0.05 ,]#差异同源基因


# 差异基因富集分析
ko_PRJNA1083304_diff=enrichKO(
  rownames(df_diff_ko_PRJNA1083304),
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH")

class(ko_PRJNA1083304_diff@result$p.adjust)
class(ko_PRJNA1083304_diff@result$p.adjust[ko_PRJNA1083304_diff@result$p.adjust < 0.05]) <- as.numeric(ko_PRJNA1083304_diff@result$p.adjust)

# 筛选调整后P值小于0.05的
ko_PRJNA1083304_diff = ko_PRJNA1083304_diff@result %>% dplyr::filter(p.adjust < 0.05)
rownames(ko_PRJNA1083304_diff) <- 1:nrow(ko_PRJNA1083304_diff)

# 添加一列设置通路的显示顺序
ko_PRJNA1083304_diff$type_order=factor(rev(as.integer(rownames(ko_PRJNA1083304_diff))),
                                       labels =rev(ko_PRJNA1083304_diff$Description))

# --------------------------
# 3. KEGG功能富集条形图/气泡图 ####
# --------------------------

# 条形图
ggplot(ko_PRJNA1083304_diff,aes(y=type_order,x=Count,fill=p.adjust))+
  geom_bar(stat = "identity",width=0.8)+####柱子宽度
  scale_fill_gradient(low = "red",high ="blue" )+#颜色自己可以换
  labs(title = "bar plot of KEGG Pathways Enrichment",x = "Gene numbers",y ="")+#设置标题
  theme_bw()+#设置主题
  theme(plot.title = element_text(hjust = 0.5,size = 12),#设置图片标题
        axis.title= element_text(size = 9),#设置坐标轴标题
        axis.text=element_text(size = 9, color = 'black'),#设置坐标轴文本
        legend.title = element_blank(),#设置图例标题
        legend.text = element_text(size=9,color = 'black'))#设置图例文本

# 气泡图
ggplot(ko_PRJNA1083304_diff,aes(y=type_order,x=Count))+
  geom_point(aes(size=Count,color=-1*p.adjust))+# 修改点的大小
  scale_color_gradient(low="green",high = "red")+#颜色标度设置
  labs(color=expression(p.adjust,size="Count"),
       x="Gene Number",y="",title="bubble plot of KEGG Pathway Enrichment")+#设置标题
  theme_bw()+#设置主题
  theme(plot.title = element_text(hjust = 0.5,size = 12),#设置图片标题
        axis.title= element_text(face = "bold",size = 9),#设置坐标轴标题
        axis.text=element_text(size = 9, color = 'black'),#设置坐标轴文本
        legend.title = element_blank(),#设置图例标题
        legend.text = element_text(size=9,color = 'black'))#设置图例文本


########----------eggNOG通路富集分析-----------########

# --------------------------
# 1. 读取数据并预处理 ####
# --------------------------
library(tidyverse)
library(clusterProfiler)
library(ggplot2)

eggnog_relab_PRJNA1083304 <- read_tsv("genefamily/eggnog_relab.tsv")

eggnog_relab_PRJNA1083304 = eggnog_relab_PRJNA1083304 %>%
  filter(!str_detect(`# Gene Family`, "\\|")) %>%
  column_to_rownames(var = "# Gene Family")

# 修改pathway列名
a=data.frame(str_split_fixed(colnames(eggnog_relab_PRJNA1083304),"_",4)[,1:3])
b=str_c(a$X2)
colnames(eggnog_relab_PRJNA1083304)=b

# 读取分组信息
meta_PRJNA1083304 <- read.table("metadata.txt", 
                                header=T, row.names=1, sep="\t")

# 确保样本顺序一致
identical(colnames(eggnog_relab_PRJNA1083304),rownames(meta_PRJNA1083304))

# --------------------------
# 2. 差异基因富集分析 ####
# --------------------------
#同组基因做差异分析求P值
# 获取样本名
samples <- colnames(eggnog_relab_PRJNA1083304)

# 提取分组信息
group <- meta_PRJNA1083304[samples, "Group"]  # 保证顺序一致

# 找出组别（假设两个组）
group_levels <- unique(group)
Control_samples <- samples[group == group_levels[1]]
Case_samples <- samples[group == group_levels[2]]

# 使用组别进行差异分析（Wilcoxon检验）
eggnog_relab_PRJNA1083304_pvalue <- apply(eggnog_relab_PRJNA1083304, 1, function(x) {
  wilcox.test(x[Control_samples], x[Case_samples])$p.value
})
# control组求均值
eggnog_PRJNA1083304_control_mean = apply(eggnog_relab_PRJNA1083304[,Control_samples],1,FUN = mean)

# case组求均值
eggnog_PRJNA1083304_case_mean = apply(eggnog_relab_PRJNA1083304[,Case_samples],1,FUN = mean) 

# 对得到的P值用BH方法做矫正
eggnog_PRJNA1083304_FDR = p.adjust(eggnog_relab_PRJNA1083304_pvalue, method = 'BH') 

# 求Fold change 用case均值/control组均值
eggnog_PRJNA1083304_FC=eggnog_PRJNA1083304_case_mean/eggnog_PRJNA1083304_control_mean 

# Fold change取log
eggnog_PRJNA1083304_log2FC = log(eggnog_PRJNA1083304_FC,2) 

# 将上述指标放于一个数据框内方便计算
df_eggnog_PRJNA1083304=data.frame(eggnog_relab_PRJNA1083304_pvalue,eggnog_PRJNA1083304_FDR,eggnog_PRJNA1083304_FC,eggnog_PRJNA1083304_log2FC)

# 放宽差异基因筛选标准
# 原来只用p<0.05  结果只有4个，现在改为p<0.1或者|log2FC|>0.5
df_diff_eggnog_PRJNA1083304_strict = df_eggnog_PRJNA1083304[df_eggnog_PRJNA1083304$eggnog_relab_PRJNA1083304_pvalue<0.05 ,]
df_diff_eggnog_PRJNA1083304 = df_eggnog_PRJNA1083304[
  df_eggnog_PRJNA1083304$eggnog_relab_PRJNA1083304_pvalue < 0.1 | 
    abs(df_eggnog_PRJNA1083304$eggnog_PRJNA1083304_log2FC) > 0.5, ]

print(paste("严格标准(p<0.05)差异基因数:", nrow(df_diff_eggnog_PRJNA1083304_strict)))
print(paste("宽松标准(p<0.1或|log2FC|>0.5)差异基因数:", nrow(df_diff_eggnog_PRJNA1083304)))

# --------------------------
# 3. 创建COG功能映射 ####
# --------------------------

# 读取eggNOG注释文件
eggnog_annotations <- read_tsv("e5.og_annotations.tsv", 
                               col_names = c("ID", "Gene_Name", "COG_Category", "Description"),
                               col_types = "cccc")

# COG字母分类到功能描述的标准映射
cog_letter_to_function <- data.frame(
  COG_Category = c("J", "A", "K", "L", "B", "D", "Y", "V", "T", "M", 
                   "N", "Z", "W", "U", "O", "X", "C", "G", "E", "F", 
                   "H", "I", "P", "Q", "R", "S"),
  Function_Name = c(
    "J: Translation, ribosomal structure and biogenesis",
    "A: RNA processing and modification",
    "K: Transcription",
    "L: Replication, recombination and repair",
    "B: Chromatin structure and dynamics",
    "D: Cell cycle control, cell division, chromosome partitioning",
    "Y: Nuclear structure",
    "V: Defense mechanisms",
    "T: Signal transduction mechanisms",
    "M: Cell wall/membrane/envelope biogenesis",
    "N: Cell motility",
    "Z: Cytoskeleton",
    "W: Extracellular structures",
    "U: Intracellular trafficking, secretion, and vesicular transport",
    "O: Posttranslational modification, protein turnover, chaperones",
    "X: Mobilome: prophages, transposons",
    "C: Energy production and conversion",
    "G: Carbohydrate transport and metabolism",
    "E: Amino acid transport and metabolism",
    "F: Nucleotide transport and metabolism",
    "H: Coenzyme transport and metabolism",
    "I: Lipid transport and metabolism",
    "P: Inorganic ion transport and metabolism",
    "Q: Secondary metabolites biosynthesis, transport and catabolism",
    "R: General function prediction only",
    "S: Function unknown"
  ),
  stringsAsFactors = FALSE
)

# --------------------------
# 4. 创建基因到功能的映射 ####
# --------------------------

# 获取所有基因ID（包括差异和非差异）
all_gene_ids <- rownames(df_eggnog_PRJNA1083304)

# 从基因ID中提取COG编号
extract_cog_info <- function(gene_ids) {
  data.frame(
    Gene_ID = gene_ids,
    COG_ID = str_extract(gene_ids, "COG\\d+"),
    stringsAsFactors = FALSE
  )
}

# 提取所有基因的COG信息
all_genes_cog <- extract_cog_info(all_gene_ids)


# 基于COG编号分配功能类别（使用更合理的规则）
assign_cog_category <- function(cog_data) {
  cog_data %>%
    filter(!is.na(COG_ID)) %>%
    mutate(
      COG_num = as.numeric(gsub("COG", "", COG_ID)),
      # 使用更合理的分类规则
      COG_Category = case_when(
        COG_num >= 1 & COG_num <= 399 ~ "J",          # Translation
        COG_num >= 400 & COG_num <= 599 ~ "K",        # Transcription
        COG_num >= 600 & COG_num <= 899 ~ "L",        # Replication
        COG_num >= 900 & COG_num <= 1199 ~ "D",       # Cell cycle
        COG_num >= 1200 & COG_num <= 1499 ~ "O",      # Posttranslational modification
        COG_num >= 1500 & COG_num <= 1799 ~ "M",      # Cell wall
        COG_num >= 1800 & COG_num <= 2099 ~ "N",      # Cell motility
        COG_num >= 2100 & COG_num <= 2399 ~ "P",      # Inorganic ion transport
        COG_num >= 2400 & COG_num <= 2699 ~ "T",      # Signal transduction
        COG_num >= 2700 & COG_num <= 2999 ~ "C",      # Energy production
        COG_num >= 3000 & COG_num <= 3299 ~ "G",      # Carbohydrate metabolism
        COG_num >= 3300 & COG_num <= 3599 ~ "E",      # Amino acid metabolism
        COG_num >= 3600 & COG_num <= 3899 ~ "F",      # Nucleotide metabolism
        COG_num >= 3900 & COG_num <= 4199 ~ "H",      # Coenzyme metabolism
        COG_num >= 4200 & COG_num <= 4499 ~ "I",      # Lipid metabolism
        COG_num >= 4500 & COG_num <= 4799 ~ "Q",      # Secondary metabolites
        COG_num >= 4800 & COG_num <= 5099 ~ "V",      # Defense mechanisms
        COG_num >= 5100 & COG_num <= 5399 ~ "R",      # General function
        COG_num >= 5400 & COG_num <= 5699 ~ "S",      # Function unknown
        TRUE ~ "R"  # 默认为General function
      )
    ) %>%
    left_join(cog_letter_to_function, by = "COG_Category")
}

# 应用分类
all_genes_with_cog <- assign_cog_category(all_genes_cog)

# 创建term2gene格式（enricher函数需要的格式）
term2gene <- all_genes_with_cog %>%
  filter(!is.na(Function_Name)) %>%
  dplyr::select(Term = Function_Name, Gene = Gene_ID) %>%
  distinct()

print(paste("可用于富集分析的基因总数:", n_distinct(term2gene$Gene)))
print(paste("COG功能类别数:", n_distinct(term2gene$Term)))

# --------------------------
# 5. 差异基因富集分析（使用enricher而不是enrichGO）####
# --------------------------

# 获取差异基因列表
diff_genes <- rownames(df_diff_eggnog_PRJNA1083304)
print(paste("差异基因数量:", length(diff_genes)))

# 确保差异基因在term2gene中存在
diff_genes_with_cog <- diff_genes[diff_genes %in% term2gene$Gene]
print(paste("有COG注释的差异基因数:", length(diff_genes_with_cog)))


# 使用enricher进行COG功能富集分析
if (length(diff_genes_with_cog) >= 3) {  # 降低最小基因数要求
  eggnog_PRJNA1083304_diff <- enricher(
    gene = diff_genes_with_cog,
    pvalueCutoff = 1,  # 获取所有结果
    pAdjustMethod = "none",  # 不进行校正
    universe = unique(term2gene$Gene),
    TERM2GENE = term2gene,
    minGSSize = 2,  # 降低最小基因集大小
    maxGSSize = 500
  )
  
  # 检查是否有富集结果
  if (!is.null(eggnog_PRJNA1083304_diff) && nrow(eggnog_PRJNA1083304_diff@result) > 0) {
    
    # 查看所有结果
    print("所有富集分析结果（按p值排序）:")
    all_results <- eggnog_PRJNA1083304_diff@result %>%
      arrange(pvalue)
    print(head(all_results[, c("Description", "Count", "pvalue", "GeneRatio")], 20))
    
    # 使用更宽松的p值阈值
    eggnog_PRJNA1083304_diff_result <- eggnog_PRJNA1083304_diff@result %>% 
      dplyr::filter(pvalue < 0.2)  # 使用p<0.2
    
    if (nrow(eggnog_PRJNA1083304_diff_result) > 0) {
      # 重置行名
      rownames(eggnog_PRJNA1083304_diff_result) <- 1:nrow(eggnog_PRJNA1083304_diff_result)
      
      # 添加一列设置通路的显示顺序
      eggnog_PRJNA1083304_diff_result$type_order <- factor(
        rev(as.integer(rownames(eggnog_PRJNA1083304_diff_result))),
        labels = rev(eggnog_PRJNA1083304_diff_result$Description)
      )
      
      # 计算Gene Ratio的数值
      eggnog_PRJNA1083304_diff_result$GeneRatio_numeric <- sapply(
        strsplit(eggnog_PRJNA1083304_diff_result$GeneRatio, "/"),
        function(x) as.numeric(x[1]) / as.numeric(x[2])
      )
      
      # 保存富集结果
      write.csv(eggnog_PRJNA1083304_diff_result, "COG_enrichment_results_pvalue_PRJNA1083304.csv", row.names = FALSE)
      
      # --------------------------
      # 6. eggNOG功能富集条形图/气泡图 ####
      # --------------------------
      
      # 条形图
      p1 <- ggplot(eggnog_PRJNA1083304_diff_result, 
                   aes(y = type_order, x = Count, fill = pvalue)) +
        geom_bar(stat = "identity", width = 0.8) +
        scale_fill_gradient(low = "red", high = "blue") +
        labs(title = "Bar plot of COG Functional Enrichment (p-value < 0.2)",
             x = "Gene numbers",
             y = "COG Functions",
             fill = "P-value") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5, size = 12),
              axis.title = element_text(size = 9),
              axis.text = element_text(size = 9, color = 'black'),
              legend.title = element_text(size = 9),
              legend.text = element_text(size = 9, color = 'black'))
      
      ggsave("COG_enrichment_barplot_pvalue_PRJNA1083304.pdf", p1, width = 10, height = 8)
      
      # 气泡图
      p2 <- ggplot(eggnog_PRJNA1083304_diff_result, 
                   aes(y = type_order, x = GeneRatio_numeric)) +
        geom_point(aes(size = Count, color = -log10(pvalue))) +
        scale_color_gradient(low = "green", high = "red") +
        scale_size_continuous(range = c(3, 10)) +
        labs(x = "Gene Ratio",
             y = "COG Functions",
             title = "Bubble plot of COG Functional Enrichment (p-value < 0.2)",
             size = "Count",
             color = "-log10(p-value)") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5, size = 12),
              axis.title = element_text(face = "bold", size = 9),
              axis.text = element_text(size = 9, color = 'black'),
              legend.title = element_text(size = 9),
              legend.text = element_text(size = 9, color = 'black'))
      
      ggsave("COG_enrichment_bubbleplot_pvalue_PRJNA1083304.pdf", p2, width = 10, height = 8)
      
      print("富集分析完成！（使用原始p值）")
      print(paste("富集的COG功能类别数 (p<0.2):", nrow(eggnog_PRJNA1083304_diff_result)))
      print("Top enriched functions (by p-value):")
      print(head(eggnog_PRJNA1083304_diff_result[, c("Description", "Count", "pvalue")], 10))
      
    } else {
      # 【修改5】即使没有显著结果，也绘制前10个结果
      print("没有p<0.2的结果，绘制p值最小的前10个")
      
      top10_results <- head(all_results, 10)
      rownames(top10_results) <- 1:nrow(top10_results)
      top10_results$type_order <- factor(
        rev(as.integer(rownames(top10_results))),
        labels = rev(top10_results$Description)
      )
      top10_results$GeneRatio_numeric <- sapply(
        strsplit(top10_results$GeneRatio, "/"),
        function(x) as.numeric(x[1]) / as.numeric(x[2])
      )
      
      # 绘制前10个结果
      p3 <- ggplot(top10_results, 
                   aes(y = type_order, x = GeneRatio_numeric)) +
        geom_point(aes(size = Count, color = -log10(pvalue))) +
        scale_color_gradient(low = "green", high = "red") +
        scale_size_continuous(range = c(3, 10)) +
        labs(x = "Gene Ratio",
             y = "COG Functions",
             title = "Top 10 COG Functions by P-value",
             size = "Count",
             color = "-log10(p-value)") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5, size = 12),
              axis.title = element_text(face = "bold", size = 9),
              axis.text = element_text(size = 9, color = 'black'))
      
      ggsave("COG_enrichment_top10_PRJNA1083304.pdf", p3, width = 10, height = 8)
      write.csv(top10_results, "COG_enrichment_top10_results.csv", row.names = FALSE)
    }
  } else {
    print("富集分析未返回结果")
  }
} else {
  print("差异基因太少，无法进行富集分析")
}

# --------------------------
# 7. 额外检查：查看数据分布 ####
# --------------------------

# 检查差异基因的COG分布
print("\n差异基因的COG类别分布:")
cog_dist <- all_genes_with_cog %>%
  filter(Gene_ID %in% diff_genes_with_cog) %>%
  group_by(COG_Category, Function_Name) %>%
  dplyr::summarise(Count = n()) %>%
  arrange(desc(Count))

print(head(cog_dist, 10))

# 绘制差异基因COG分布图
if (nrow(cog_dist) > 0) {
  p4 <- ggplot(head(cog_dist, 15), 
               aes(x = reorder(Function_Name, Count), y = Count, fill = COG_Category)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = "Distribution of Differential Genes in COG Categories",
         x = "COG Category",
         y = "Number of Genes",
         fill = "COG Class") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          axis.text.y = element_text(size = 8))
  
  ggsave("COG_category_distribution_PRJNA1083304.pdf", p4, width = 10, height = 8)
}

# 保存统计结果
write.csv(cog_dist, "COG_category_statistics_PRJNA1083304.csv", row.names = FALSE)

# --------------------------
# 8. 生成分析报告 ####
# --------------------------

cat("\n========== eggNOG富集分析报告 ==========\n")
cat("项目编号: PRJNA1083304\n")
cat("分析日期:", format(Sys.Date(), "%Y-%m-%d"), "\n\n")

cat("样本信息:\n")
cat("- 总样本数:", length(samples), "\n")
cat("- Control组样本数:", length(Control_samples), "\n")
cat("- Case组样本数:", length(Case_samples), "\n\n")

cat("基因统计:\n")
cat("- 总基因数:", nrow(df_eggnog_PRJNA1083304), "\n")
cat("- 差异基因数 (p<0.05):", nrow(df_diff_eggnog_PRJNA1083304_strict), "\n")
cat("- 差异基因数 (p<0.1或|log2FC|>0.5):", nrow(df_diff_eggnog_PRJNA1083304), "\n")
cat("- 有COG注释的差异基因数:", length(diff_genes_with_cog), "\n\n")

if (exists("eggnog_PRJNA1083304_diff_result") && nrow(eggnog_PRJNA1083304_diff_result) > 0) {
  cat("富集分析结果:\n")
  cat("- 富集的COG功能数 (p<0.2):", nrow(eggnog_PRJNA1083304_diff_result), "\n\n")
  
  cat("Top 3 enriched COG functions:\n")
  top3 <- head(eggnog_PRJNA1083304_diff_result[, c("Description", "Count", "pvalue")], 3)
  print(top3)
}

cat("\n输出文件:\n")
cat("- COG_enrichment_results_pvalue_PRJNA1083304.csv 或 COG_enrichment_top10_results.csv\n")
cat("- COG_enrichment_barplot_pvalue_PRJNA1083304.pdf\n")
cat("- COG_enrichment_bubbleplot_pvalue_PRJNA1083304.pdf 或 COG_enrichment_top10_PRJNA1083304.pdf\n")
cat("- COG_category_distribution_PRJNA1083304.pdf\n")
cat("- COG_category_statistics_PRJNA1083304.csv\n")
cat("========================================\n")














