setwd("twoSets")
install.packages("tidyverse")
install.packages("BiocManager")
chooseBioCmirror()
BiocManager::install('GEOquery')
install.packages("devtools")
install_github("jmzeng1314/idmap1")
install_github("jmzeng1314/idmap2")
install_github("jmzeng1314/idmap3")
install_github("jmzeng1314/AnnoProbe")
BiocManager::install("hugene10sttranscriptcluster.db")#GPL6244
library(hugene10sttranscriptcluster.db)
library(tidyverse)
library(GEOquery)
library(stringr)
library(devtools)
library(idmap1)
library(idmap2)
library(idmap3)
library(AnnoProbe)
library(limma)
library(stringr)
rm(list = ls())

#### 获取数据集 ####
eset1 <- getGEO("GSE28735", destdir = ".", getGPL = F)
eset2 <- getGEO("GSE62452", destdir = ".", getGPL = F)
eset3 <- getGEO('GSE15471', destdir = ".", getGPL = F)
eset4 <- getGEO('GSE16515', destdir = ".", getGPL = F)

#### 提取表达矩阵 ####
exp1 <- exprs(eset1[[1]])
exp2 <- exprs(eset2[[1]])
exp3 <- exprs(eset3[[1]])
exp4 <- exprs(eset4[[1]])

range(exp1)#对数转换
range(exp2)
range(exp3)
range(exp4)
exp4 <- log2(exp4+1)
range(exp4)

table(rownames(exp1) %in% rownames(exp2)) 
length(intersect(rownames(exp1), rownames(exp2)))
boxplot(exp1, outline=F, notch=F, las=2, col='blue')
boxplot(exp2, outline=F, notch=F, las=2, col='red')
boxplot(exp3)
dev.off()
#如果数据不均一
exp1 <- normalizeBetweenArrays(exp1)
exp2 <- normalizeBetweenArrays(exp2)
exp3 <- normalizeBetweenArrays(exp3)
exp4 <- normalizeBetweenArrays(exp4)
#### 提取临床信息 ####
pd1 <- pData(eset1[[1]])
pd2 <- pData(eset2[[1]])
pd3 <- pData(eset3[[1]])
pd4 <- pData(eset4[[1]])
index1 <- eset1[[1]]@annotation
index2 <- eset2[[1]]@annotation
index3 <- eset3[[1]]@annotation
index4 <- eset4[[1]]@annotation
#### exp1 ID转换 ####
#ids <- toTable(hugene10sttranscriptclusterSYMBOL)
ids = getIDs('gpl6244')
head(ids)

exp1 <- data.frame(exp1)
exp1 <- exp1 %>% mutate(probe_id=rownames(exp1))
exp1 <- exp1 %>% inner_join(ids, by="probe_id")
exp1 <- exp1[!duplicated(exp1$symbol),] #去重
rownames(exp1) <- exp1$symbol
exp1 <- exp1[, -(91:93)]
#### exp2 ID转换 ####
exp2 <- data.frame(exp2)
exp2 <- exp2 %>% mutate(probe_id=rownames(exp2))
exp2 <- exp2 %>% inner_join(ids, by="probe_id")
exp2 <- exp2[!duplicated(exp2$symbol),] #去重
rownames(exp2) <- exp2$symbol
exp2 <- exp2[, -(131:133)]
#### exp3 ID转换 ####
ids = getIDs('gpl570')
head(ids)

exp3 <- data.frame(exp3)
exp3 <- exp3 %>% mutate(probe_id=rownames(exp3))
exp3 <- exp3 %>% inner_join(ids, by="probe_id")
exp3 <- exp3[!duplicated(exp3$symbol),] #去重
rownames(exp3) <- exp3$symbol
exp3 <- exp3[, -(79:81)]
#### exp4 ID转换 ####
exp4 <- data.frame(exp4)
exp4 <- exp4 %>% mutate(probe_id=rownames(exp4))
exp4 <- exp4 %>% inner_join(ids, by="probe_id")
exp4 <- exp4[!duplicated(exp4$symbol),] #去重
rownames(exp4) <- exp4$symbol
exp4 <- exp4[, -(53:55)]

####  合并 ####
identical(exp1, exp2)
same_genes1 <- intersect(rownames(exp1), rownames(exp2))
eset_merge1 <- cbind(exp1[same_genes1, ], exp2[same_genes1, ])

same_genes2 <- intersect(rownames(exp3), rownames(exp4))
eset_merge2 <- cbind(exp3[same_genes2, ], exp4[same_genes2, ])

identical(eset_merge1, eset_merge2)
same_genes <- intersect(rownames(eset_merge1), rownames(eset_merge2))
eset_merge <- cbind(eset_merge1[same_genes, ], eset_merge2[same_genes, ])

col <- c(rep('red', 90), rep('blue', 130), rep('yellow', 78), rep('green', 52))
par(mar = c(8, 3, 1, 2) + 0.1)
boxplot(eset_merge, col=col, las=2)
dev.off()
#### 分组 ####
group1 = ifelse(str_detect(pd1$characteristics_ch1, "N"), "Normal", "Tumor")
group2 = ifelse(str_detect(pd2$characteristics_ch1, "adjacent"), "Normal", "Tumor")
group3 = ifelse(str_detect(pd3$characteristics_ch1.1, "normal"), "Normal", "Tumor")
group4 = ifelse(str_detect(pd4$title, "Normal"), "Normal", "Tumor")
table(group1)
table(group2)
table(group3)
table(group4)
group_list = c(group1, group2, group3, group4)
table(group_list)
group_list = factor(group_list, levels = c("Normal", "Tumor"))

#### 去除批次效应 ####
GSE <- c(rep('GSE28735', 90), rep('GSE62452', 130), rep('GSE15471', 78), rep('GSE16515', 52))
GSE
table(group_list, GSE)

data <- eset_merge
batch <- c(rep('GSE28735', 90), rep('GSE62452', 130), rep('GSE15471', 78), rep('GSE16515', 52))
design <- model.matrix(~group_list)
expr_limma <- removeBatchEffect(data, batch = batch, design = design)
boxplot(expr_limma, col=col, las=2)

par(mfrow = c(1, 2))
boxplot(eset_merge, col = col, las = 2, main = "before")
boxplot(expr_limma, col = col, las = 2, main = "after")

dev.off()

expr_limma <- normalizeBetweenArrays(expr_limma)
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
                         palette = c("blue", "red"),
                         addEllipses = FALSE,
                         legend.title = "Groups"
)
pca_plot

#### 差异分析 ####
library(limma) #GEO都是用limma包
design = model.matrix(~group_list)
fit = lmFit(expr_limma,design)
fit = eBayes(fit)
deg = topTable(fit, coef = 2, number = Inf, adjust.method = "BH")
write.table(deg, file = "deg_all.csv", sep = "\t", row.names = T, col.names = NA, quote = F)
# 标记上下调基因
logFC = 1
P.value = 0.05
k1 = (deg$adj.P.Val < P.value)&(deg$logFC < -logFC)
k2 = (deg$adj.P.Val < P.value)&(deg$logFC > logFC)
deg$change = ifelse(k1,"down",ifelse(k2,"up","stable"))
table(deg$change)
my_deg <- deg %>% filter(deg$change != "stable")
head(my_deg)
write.table(my_deg, file = "my_deg.csv", sep = "\t", row.names = T, col.names = NA, quote = F)
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
diff = expr_limma[cg, ]
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
BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
library(org.Hs.eg.db)
deg <- read.table("deg_all.csv", sep = "\t", row.names = 1, check.names = F, stringsAsFactors = F, header = T)

deg <- deg %>% filter(change != "stable")
table(deg$change)
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

barplot(ego, color = "pvalue")
dotplot(ego)

barplot(ego, drop = TRUE, showCategory = 8, split = "ONTOLOGY") +
  facet_grid(ONTOLOGY~., scales = 'free')
#### KEGG富集分析 ####
kk <- enrichKEGG(gene = DEG$ENTREZID,
                 organism = 'hsa', #Human sapiens
                 pvalueCutoff = 0.1,
                 qvalueCutoff = 0.1)
kk_res <- kk@result

dotplot(kk)
