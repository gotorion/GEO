setwd("2set")
install.packages("devtools")
install_github("jmzeng1314/idmap1")
install_github("jmzeng1314/idmap2")
install_github("jmzeng1314/idmap3")
install_github("jmzeng1314/AnnoProbe")
library(tidyverse)
library(GEOquery)
library(stringr)
library(devtools)
library(idmap1)
library(idmap2)
library(idmap3)
library(AnnoProbe)
library(sva)
library(limma)
rm(list = ls())

#### 获取数据集 ####
eset1 <- getGEO("GSE28735", destdir = ".", getGPL = F)
eset2 <- getGEO("GSE62452", destdir = ".", getGPL = F)

#### 提取表达矩阵 ####
exp1 <- exprs(eset1[[1]])
exp2 <- exprs(eset2[[1]])
range(exp1)
range(exp2) #已经对数转换了
table(rownames(exp1) %in% rownames(exp2)) 
length(intersect(rownames(exp1), rownames(exp2)))#没有交集
boxplot(exp1, outline=F, notch=F, las=2, col='blue')
boxplot(exp2, outline=F, notch=F, las=2, col='red')
dev.off()
#如果数据不均一
# exp1 <- normalizeBetweenArrays(exp1)

#### 提取临床信息 ####
pd1 <- pData(eset1[[1]])
pd2 <- pData(eset2[[1]])
index1 <- eset1[[1]]@annotation
index2 <- eset2[[1]]@annotation
#### exp1 ID转换 ####
BiocManager::install("hugene10sttranscriptcluster.db")#GPL6244
library(hugene10sttranscriptcluster.db)
ids <- toTable(hugene10sttranscriptclusterSYMBOL) # ids = getIDs('gpl6244')
head(ids)

exp1 <- data.frame(exp1)
exp1 <- exp1 %>% mutate(probe_id=rownames(exp1))
exp1 <- exp1 %>% inner_join(ids, by="probe_id")
exp1 <- exp1[!duplicated(exp1$symbol),] #去重
rownames(exp1) <- exp1$symbol
exp1 <- exp1[, -(91:92)]
save(exp1, file="GSE28734.Rdata")
#### exp2 ID转换 ####
exp2 <- data.frame(exp2)
exp2 <- exp2 %>% mutate(probe_id=rownames(exp2))
exp2 <- exp2 %>% inner_join(ids, by="probe_id")
exp2 <- exp2[!duplicated(exp2$symbol),] #去重
rownames(exp2) <- exp2$symbol
exp2 <- exp2[, -(131:132)]
save(exp2, file="GSE62452.Rdata")
####  合并 ####
identical(exp1, exp2)
same_genes <- intersect(rownames(exp1), rownames(exp2))
eset_merge <- cbind(exp1[same_genes, ], exp2[same_genes, ])
col <- c(rep('red', 90), rep('blue', 130))
par(mar = c(8, 3, 1, 2) + 0.1)
boxplot(eset_merge, col=col, las=2)

#### 分组 ####
group1 = ifelse(str_detect(pd1$title, "nontumor"), "Normal", "Tumor")
group2 = ifelse(str_detect(pd2$source_name_ch1, "Non-tumor"), "Normal", "Tumor")
table(group1)
table(group2)
group_list = c(group1, group2)
table(group_list)
group_list = factor(group_list, levels = c("Normal", "Tumor"))

#### 去除批次效应 ####
GSE <- c(rep('GSE28734', 90), rep('GSE62452', 130))
GSE
table(group_list, GSE)

data <- eset_merge
batch <- c(rep('GSE28734', 90), rep('GSE62452', 130))
design <- model.matrix(~group_list)
expr_limma <- removeBatchEffect(data, batch = batch, design = design)
boxplot(expr_limma, col=col, las=2)

par(mfrow = c(1, 2))
boxplot(eset_merge, col = col, las = 2, main = "before")
boxplot(expr_limma, col = col, las = 2, main = "after")

#### PCA ####
dat = as.data.frame(t(expr_limma))
install.packages("FactoMineR")
install.packages("factoextra")
library(FactoMineR)
library(factoextra)
dat.pca <- PCA(dat, graph = FALSE)
pca_plot <- fviz_pca_ind(dat.pca,
                         geom.ind = "point",
                         col.ind = group_list,
                         palette = c("green", "red"),
                         addEllipses = TRUE,
                         legend.title = "Groups"
)
pca_plot
#### 差异分析 ####
library(limma) #GEO都是用limma包
design = model.matrix(~group_list)
fit = lmFit(expr_limma,design)
fit = eBayes(fit)
deg = topTable(fit, coef = 2, number = Inf, adjust.method = "BH")
write.table(deg, file = "deg_all.txt", sep = "\t", row.names = T, col.names = NA, quote = F)
# 标记上下调基因
logFC = 1
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
deg$logP <- -log10(deg$adj.P.Val)
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