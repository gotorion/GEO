setwd("GSE32676") 
rm(list = ls())
install.packages("tidyverse")
install.packages("BiocManager")
chooseBioCmirror()
BiocManager::install('GEOquery')
library(tidyverse)
library(stringr)
#### 获取数据集 ####
gset = getGEO('GSE32676', destdir = ".", AnnotGPL = F, getGPL = F)

#通过pData函数获取分组信息
pdata <- pData(gset[[1]])
#设置参考水平
group_list <- ifelse(str_detect(pdata$source_name_ch1, "tumor"), "tumor", "normal")
group_list
group_list <- factor(group_list, levels = c("normal", "tumor"))
# 通过exprs函数获取表达矩阵并矫正
exp <- exprs(gset[[1]])
boxplot(exp, outline=FALSE, notch=T, col=group_list, las=2)
#### 矫正 ####
library(limma)
exp = normalizeBetweenArrays(exp)
#什么时候考虑批次效应
#exp = removeBatchEffect() 
boxplot(exp, outline=FALSE, notch=T, col=group_list, las=2)
#range(exp) 
#一般而言，范围在20以内的表达值基本已经经过了log对数转换
#exp <- log2(exp+1)
#### PCA图 ####
dat = as.data.frame(t(exp))
install.packages("FactoMineR")
install.packages("factoextra")
library(FactoMineR)
library(factoextra)
dat.pac <- PCA(dat, graph = FALSE)
pca_plot <- fviz_pca_ind(dat.pac,
                         geom.ind = "point",
                         col.ind = group_list,
                         palette = c("green", "red"),
                         addEllipses = FALSE,
                         legend.title = "Groups"
)
pca_plot

#### 使用R包转换id ####
index = gset[[1]]@annotation # 查看基因测序的平台
BiocManager::install("org.Hs.eg.db")
if(!require("hgu133plus2.db"))
  BiocManager::install("hgu133plus2.db")  
library(hgu133plus2.db)
ids <- toTable(hgu133plus2SYMBOL)

exp <- data.frame(exp)
exp <- exp %>% mutate(probe_id=rownames(exp))
exp <- exp %>% inner_join(ids, by="probe_id")
exp <- exp[!duplicated(exp$symbol),] #去重
rownames(exp) <- exp$symbol
exp <- exp[, -(33:34)]
write.table(exp, file = "exp.txt", sep = "\t", row.names = T, col.names = NA, quote = F)
#### 差异分析 ####
exp <- read.table("exp.txt", sep = "\t", row.names = 1, check.names = F, stringsAsFactors = F, header = T)
library(limma) #GEO都是用limma包
design = model.matrix(~group_list)
fit = lmFit(exp,design)
fit = eBayes(fit)
deg = topTable(fit, coef = 2, number = Inf, adjust.method = "BH")
write.table(deg, file = "deg_all.txt", sep = "\t", row.names = T, col.names = NA, quote = F)
# 标记上下调基因
logFC = 2
P.value = 0.05
k1 = (deg$adj.P.Val < P.value)&(deg$logFC < -logFC)
k2 = (deg$adj.P.Val < P.value)&(deg$logFC > logFC)
deg$change = ifelse(k1,"down",ifelse(k2,"up","stable"))
table(deg$change)

####火山图 ####
install.packages("ggpubr")
install.packages("ggthemes")
library(ggpubr)
library(ggthemes)
deg$logP <- -log10(deg$P.Value)
ggscatter(deg, x = "logFC", y = "logP",
          color = "change",
          palette = c("blue", "black", "red"),
          size = 1) + 
  theme_base() +
  geom_hline(yintercept = -log10(P.value), linetype = "dashed") +
  geom_vline(xintercept = c(-logFC, logFC), linetype = "dashed")
#### 热图 ####
cg = rownames(deg)[deg$change != "stable"]
diff = exp[cg, ]
install.packages("pheatmap")
library(pheatmap)
anntation_col = data.frame(group=group_list)
rownames(anntation_col) = colnames(diff)
pheatmap(diff,
         anntation_col = anntation_col,
         scale = "row",
         show_rownames = F,
         show_colnames = F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         fontsize = 10,
         fontsize_row = 3,
         fontsize_col = 3
         )

#### Go 富集分析 ####
BiocManager::install("clusterProfiler")
library(clusterProfiler)
deg <- read.table("deg_all.txt", sep = "\t", row.names = 1, check.names = F, stringsAsFactors = F, header = T)
logFC = 1
P.value = 0.05
k1 = (deg$adj.P.Val < P.value)&(deg$logFC < -logFC)
k2 = (deg$adj.P.Val < P.value)&(deg$logFC > logFC)
deg$change = ifelse(k1,"down",ifelse(k2,"up","stable"))

deg <- deg %>% filter(change != "stable")
#table(deg$change)
DEG <- deg
DEG <- DEG %>% rownames_to_column("Gene")

genelist <- bitr(DEG$Gene, fromType = "SYMBOL",
                 toType = "ENTREZID", OrgDb = 'org.Hs.eg.db')
DEG <- inner_join(DEG, genelist, by = c("Gene" = "SYMBOL"))

ego <- enrichGO(gene = DEG$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego_res <- ego@result
#### KEGG富集分析 ####
kk <- enrichKEGG(gene = DEG$ENTREZID,
                 organism = 'hsa', #Human sapiens
                 pvalueCutoff = 0.1,
                 qvalueCutoff = 0.1)
kk_res <- kk@result
