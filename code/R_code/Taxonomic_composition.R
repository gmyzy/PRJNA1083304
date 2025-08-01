# --------------------------
####   物种堆积图   ####
# --------------------------


# 加载包
library(tidyverse)

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


# 读取MetaPhlAn门水平丰度表 #########
phylum_PRJNA1083304 <- read.table("04_metaphlan_result/01_metaphlan_phylum.txt", 
                                  header=T, row.names=1, sep="\t") 
colnames(phylum_PRJNA1083304)
# 修改列名
colnames(phylum_PRJNA1083304) <- gsub("clean_(SRR[0-9]+)_metaphlan", "\\1", )

                                      
identical(colnames(phylum_PRJNA1083304),rownames(meta_PRJNA1083304))
colnames(phylum_PRJNA1083304)



# 确保样本顺序一致
identical(colnames(species_PRJNA1083304),rownames(meta_PRJNA1083304))


# --------------------------
# 2. 数据处理 ####
# --------------------------


# 计算物种总和
phylum_PRJNA1083304$sum <- rowSums(phylum_PRJNA1083304)

# 按照物种总的丰度进行排序
phylum_PRJNA1083304_order=phylum_PRJNA1083304 %>% 
  rownames_to_column(var="id") %>%
  arrange(desc(sum)) %>% column_to_rownames(var="id")
phylum_PRJNA1083304$sum <- NULL
# 取总丰度排名前20的物种
phylum_20=phylum_PRJNA1083304_order[1:20, -ncol(phylum_PRJNA1083304_order)]

# 求每个物种按不同组别的平均数
mean = apply(phylum_PRJNA1083304,1,function(x){tapply(x,meta_PRJNA1083304$Group,mean)}) %>% 
  t() %>% as.data.frame() #所有的物种取平均数

# # 只计算前20的物种在两个组中的平均数
# mean1 = apply(phylum_20,1,function(x){tapply(x,meta_PRJNA1083304$Group,mean)}) %>% 
#   t() %>% as.data.frame()
# 
# # 把其余不在前20的物种当作其他放在堆积图中,那么就会有21个
# others <-apply(mean[!rownames(mean) %in% rownames(mean1),],2,sum) %>% as.data.frame(.)
# colnames(others) <- "Others"
# others <- as.data.frame(t(others))
# 
# # 合并前20的物种以及其他物种在不同组别的丰度平均数数据
# phylum <- rbind(mean1,others)
phylum <- mean
# 宽变长,-1表示不对第一列进行宽变长计算
ggdata=phylum %>% rownames_to_column(var = "phylum") %>%
  pivot_longer(-1,names_to = "group", values_to = "value")

# 为了物种顺序不乱,物种因子化
ggdata$phylum <- factor(ggdata$phylum,levels =unique(ggdata$phylum))
ggdata$group <- factor(ggdata$group,levels=c("Case","Control"))



# --------------------------
# 3. 绘制堆积图 ####
# --------------------------


pdf("./003_1_phylum_stack.pdf",width = 8,height = 8 )

ggplot(ggdata,aes(x = group,y=value,fill = phylum)) +
  geom_bar(stat = "identity",position = "fill",colour= 'black') + #堆叠柱状图
  labs(x="",y="Abundance of phylum(%)")+ #设置坐标轴标题
  scale_fill_manual(values =
                      c("#FFC1CC","#FFA500","#008000","#00FFFF","#0000FF",
                                 
                                 "#800080","#6495EE","#FF0000","#00CFD2","#A1522D",
                                 
                                 "#FFE5C5","#8A2BE3","#FB8072","#FFFF00","#7B68EF",
                                 
                                 "#cfe74f","#39cf40","#3e6926","#0a0883","#49ddd0",
                                 
                                 "#e0f8f7","#651d91","#9d3eaf","#b9a3d0","#5b5a5b","#8f8f8f"))+ #设置颜色标度其实这里也可以考虑用现成的调色板
                                   theme_classic()+#设置基本主题
  theme(plot.title = element_text(hjust = 0.5), #图片标题的位置居中
        axis.title= element_text(size=20,color="black",face="bold"),#坐标轴标题的字体大小颜色以及字体样式
        axis.text =element_text(size = 20, face = "bold", color ='black'),#坐标轴文本的字体大小颜色以及字体样式
        axis.text.x =element_text(angle = 45, vjust =0.5),#x坐标轴文本的字体方向和角度
        legend.title=element_blank(),#不要图例标题
        legend.position = "right",#图例位置为右方
        legend.text = element_text(size=12,face = "bold.italic",color ='black'),#设置图例文本
        axis.line=element_line(linetype=1,color="black",size=1),#设置坐标轴的线
        axis.ticks.length = unit(5, "pt"),#设置坐标轴刻度的长度
        axis.ticks = element_line(color = "black",size = 1))+#设置坐标轴刻度
  guides(fill = guide_legend(ncol = 1, byrow = TRUE)) #图例在一列

dev.off()



# --------------------------
# 4. 差异分析（识别显著物种） ####
# --------------------------


# 使用Wilcoxon秩和检验（Mann-Whitney U检验）进行差异分析

# 准备数据：将物种丰度表与分组信息合并
diff_data <- phylum_PRJNA1083304_order[1:20, ] %>%  # 只分析前20物种
  rownames_to_column("phylum") %>%
  pivot_longer(-phylum, names_to = "Sample", values_to = "Abundance") %>%
  left_join(meta_PRJNA1083304 %>% rownames_to_column("Sample"), by = "Sample")

# 准备数据：将物种丰度表与分组信息合并
diff_data <- phylum_PRJNA1083304_order %>%  # 只分析前20物种
  rownames_to_column("phylum") %>%
  pivot_longer(-phylum, names_to = "Sample", values_to = "Abundance") %>%
  left_join(meta_PRJNA1083304 %>% rownames_to_column("Sample"), by = "Sample")




# 对每个物种进行Wilcoxon检验
diff_results <- diff_data %>%
  group_by(phylum) %>%
  summarise(
    p_value = wilcox.test(Abundance ~ Group, exact = FALSE)$p.value
  ) %>%
  ungroup() %>%
  mutate(
    Significance = ifelse(p_value < 0.05, "*", "")  # 直接使用原始p值判断显著性
  )

# 获取显著差异物种
sig_phylum <- diff_results %>% 
  filter(p_value < 0.05) %>% 
  pull(phylum)

# 打印差异分析结果
print(diff_results)
cat("\nSignificant phylum (p < 0.05):", sig_phylum, "\n")

# --------------------------
# 5. 绘制带显著性标记的堆积图 ####
# --------------------------

# 将带显著性标记的物种添加星号
phylum_labels <- rownames(phylum)
phylum_labels <- ifelse(phylum_labels %in% 
                           sig_phylum, paste0(phylum_labels, " *"), 
                         phylum_labels)

# 将新标签应用到 ggdata$phylum 因子
ggdata$phylum <- factor(ggdata$phylum, 
                         levels = rownames(phylum), 
                         labels = phylum_labels)

pdf("./003_2_phylum_stack_sig.pdf",width = 8,height = 8 )

ggplot(ggdata,aes(x = group,y=value,fill = phylum)) +
  geom_bar(stat = "identity",position = "fill",colour= 'black') + #堆叠柱状图
  labs(x="",y="Abundance of phylum(%)")+ #设置坐标轴标题
  scale_fill_manual(values =
                      c("#FFC1CC","#FFA500","#008000","#00FFFF","#0000FF",
                                 
                                 "#800080","#6495EE","#FF0000","#00CFD2","#A1522D",
                                 
                                 "#FFE5C5","#8A2BE3","#FB8072","#FFFF00","#7B68EF",
                                 
                                 "#cfe74f","#39cf40","#3e6926","#0a0883","#49ddd0",
                                 
                                 "#e0f8f7","#651d91","#9d3eaf","#b9a3d0","#5b5a5b","#8f8f8f"))+ #设置颜色标度其实这里也可以考虑用现成的调色板
                                   theme_classic()+#设置基本主题
  theme(plot.title = element_text(hjust = 0.5), #图片标题的位置居中
        axis.title= element_text(size=20,color="black",face="bold"),#坐标轴标题的字体大小颜色以及字体样式
        axis.text =element_text(size = 20, face = "bold", color ='black'),#坐标轴文本的字体大小颜色以及字体样式
        axis.text.x =element_text(angle = 45, vjust =0.5),#x坐标轴文本的字体方向和角度
        legend.title=element_blank(),#不要图例标题
        legend.position = "right",#图例位置为右方
        legend.text = element_text(size=12,face = "bold.italic",color ='black'),#设置图例文本
        axis.line=element_line(linetype=1,color="black",size=1),#设置坐标轴的线
        axis.ticks.length = unit(5, "pt"),#设置坐标轴刻度的长度
        axis.ticks = element_line(color = "black",size = 1))+#设置坐标轴刻度
  guides(fill = guide_legend(ncol = 1, byrow = TRUE)) #图例在一列

dev.off()



