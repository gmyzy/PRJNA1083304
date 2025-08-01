# --------------------------
####   α多样性   ####
# --------------------------

rm(list = ls())
getwd()
setwd("/mnt/raid6/gengmingyan/PRJNA1083304/")


# 加载必要的R包
library(vegan)
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)

# --------------------------
# 1. 数据准备 ####
# --------------------------

# 读取MetaPhlAn丰度表
data_PRJNA1083304 <- read.table("04_metaphlan_result/00_merged_abundance.txt", 
                                header=T, row.names=1, sep="\t") 
colnames(data_PRJNA1083304)
# 修改列名
colnames(data_PRJNA1083304) <- gsub("clean_(SRR[0-9]+)_metaphlan", "\\1", 
                                    colnames(data_PRJNA1083304))

# 转置数据（样本为行，物种为列）
data_PRJNA1083304_t <- as.data.frame(t(data_PRJNA1083304))

# 读取分组信息
meta_PRJNA1083304 <- read.table("metadata.txt", 
                                header=T, row.names=1, sep="\t")

# 确保样本顺序一致
identical(rownames(data_PRJNA1083304_t),rownames(meta_PRJNA1083304))


# --------------------------
# 2. 计算α多样性指数 ####
# --------------------------

### 构建计算α多样性函数
alpha_diversity <- function(x, tree = NULL) {
  # 计算observed_species、Chao1指数、ACE（需要整数值，此处用相对丰度*10000近似）
  observed_species <- estimateR(round(x * 10000))[1, ]
  Chao1 <- estimateR(round(x * 10000))[2, ]
  ACE <- estimateR(round(x * 10000))[4, ]
  Shannon <- vegan::diversity(x, index = 'shannon',base = 2)
  Simpson <- vegan::diversity(x, index = 'simpson') 
  goods_Coverage <- 1 - rowSums(x == 1) / rowSums(x)
  # 保留四位小数
  Shannon <- sprintf("%0.4f", Shannon)
  Simpson <- sprintf("%0.4f", Simpson)
  goods_Coverage <- sprintf("%0.4f", goods_Coverage)
  result <- data.frame(observed_species, ACE,Chao1,
                       Shannon, Simpson, goods_Coverage)
}

alpha <- alpha_diversity(data_PRJNA1083304_t) %>% select(-c(2,6))
# 添加分组
alpha$group <- meta_PRJNA1083304$Group


# 把α多样性参数变成数值型
alpha[,3]=as.numeric(alpha[,3])
alpha[,4]=as.numeric(alpha[,4])

# 绘图前数据宽边长
diversity_long <- alpha %>% pivot_longer(-group,names_to="grp",values_to="val")
# 数据group变成因子型
diversity_long$group <- factor(diversity_long$group,levels =
                                     c("Control","Case"))


# --------------------------
# 3. 可视化 ####
# --------------------------

# α多样性绘图
pdf("./001_alpha_diversity.pdf",width = 8,height = 8 )
ggplot(diversity_long,aes(x=group,y=val,fill=grp)) +
  geom_violin(trim=FALSE,color="white") + 
  # 绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(以下设为背景为白色，其实表示不要轮廓线)
  # "trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  geom_boxplot(width=0.2,position=position_dodge(0.9))+
  stat_compare_means(aes(group=group),
                     label = "p.signif",label.x = 1.5)+#添加P值
  labs(title = "Alpha diversity",x="",y="")+ #图片title x轴和y轴标题
  theme_classic()+#设置主题
  theme(plot.title =element_text(size = 15,hjust = 0.5,face =
                                   "bold",color = 'black'),#标题的文本设置、标题居中、字体颜色和格式
        axis.text =element_text(size=14,color = "black"),#设置y轴文本
        legend.text = element_text(size=14),#设置图例文本字体的大小
        legend.position = "bottom",#图例位置设置
        legend.title = element_blank(),#图例标题设置
        strip.text = element_text(size=14))+ #设置分面文本的字体大小
  scale_fill_manual(values=c("#FFA500","#FFE5C5","#FB8072","#7B68EF"))+
  #颜色标度
  stat_summary(fun=mean, geom="point", size=2) + #给图添加均值
  facet_wrap(~grp,scales = "free",nrow =2)+#分面
  scale_x_discrete(limits=c("Control","Case"))#离散型变量x轴顺序
dev.off()











