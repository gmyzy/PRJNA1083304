# --------------------------
####   β多样性   ####
# --------------------------


# 加载必要的R包
library(ape)



# --------------------------
# 1. 数据准备 ####
# --------------------------

# 读取MetaPhlAn物种丰度表
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


# --------------------------
# 2. 计算距离 ####
# --------------------------


bray_dist<-vegdist(t(species_PRJNA1083304),method = "bray")

# 使用ape这个包中的pcoa()函数做PCoA分析
pcoa<-pcoa(bray_dist,correction = "cailliez")
pcoa_species <- pcoa$vectors[,c(1,2)] %>% 
  cbind(.,meta_PRJNA1083304$Group) %>%
  as.data.frame()
pcoa_species$V3 <- factor(pcoa_species$V3,levels=c("Control","Case"))
pcoa_species$Axis.1 <- as.numeric(pcoa_species$Axis.1)
pcoa_species$Axis.2 <- as.numeric(pcoa_species$Axis.2)

# 看组与组之间有没有差异：最后一列
result<- adonis2(formula = t(species_PRJNA1083304) ~ (meta_PRJNA1083304$Group),
                 data = meta_PRJNA1083304, permutations = 999, method = "bray")

# --------------------------
# 3. 可视化 ####
# --------------------------


# 提前设置绘图的颜色和坐标轴标题
color=c("#4DA2AF","#F1C8C8")
x_label<-round(pcoa$values$Rel_corr_eig[1]*100,2)
y_label<-round(pcoa$values$Rel_corr_eig[2]*100,2)


# PCoA图绘制（加P值）
pdf("./002_PCoA.pdf",width = 6,height = 6 )
ggplot(data = pcoa_species,aes(x=Axis.1,y=Axis.2,color=V3,fill=V3,group = V3))+
  geom_point(size= 4,pch = 21)+#点图
  
  stat_ellipse(type = "t", linetype = 1)+#添加置信椭圆
  
  labs(title = "PCOA Bray_Dist distance(species)",
       x=paste0("PCoA1[",x_label,"%]"),
       y=paste0("PCoA2[ ",y_label,"%]"))+ #图片title x轴和y轴标题
  geom_text(aes(x=0.5,y=0.4),label="p=0.05",size=5,color="red")+#图片文本添加
  scale_color_manual(values = c("#B2182B","#56B4E9"))+#外圈颜色标度
  scale_fill_manual(values = c("#B2182B","#56B4E9"))+#内圈颜色标度
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +   
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  theme_bw() +#设置基本主题
  theme(plot.title = element_text(size = 15,color = "black",hjust =0.5),#设置图片标题
        text=element_text(size=14,face="plain",color="black"),#设置图片文本
        axis.title=element_text(size=16,face="plain",color="black"),#设置坐标轴标题
        axis.text = element_text(size=14,face="plain",color="black"),#设置坐标轴文本
        legend.title =element_text(size=16,face="plain",color="black"),#设置图例标题
        legend.text =element_text(size=14,face="plain",color="black"),#设置图例文本
        legend.background = element_blank(),#设置图例背景
        legend.position=c(0.9,0.1))#设置图例位置
dev.off()

























