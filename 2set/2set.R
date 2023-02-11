setwd("2set")
library(tidyverse)
library(GEOquery)
library(stringr)
install.packages("devtools")
library(devtools)
install_github("jmzeng1314/idmap1")
library(idmap1)
install_github("jmzeng1314/idmap2")
library(idmap2)
install_github("jmzeng1314/idmap3")
library(idmap3)
rm(list = ls())
# 获取数据集
eset1 <- getGEO("GSE28735", destdir = ".", getGPL = F)
eset2 <- getGEO("GSE62452", destdir = ".", getGPL = F)
# 提取表达矩阵
exp1 <- exprs(eset1[[1]])
exp2 <- exprs(eset2[[1]])
range(exp1)
range(exp2) #已经对数转换了
table(rownames(exp1) %in% rownames(exp2)) 
length(intersect(rownames(exp1), rownames(exp2)))#没有交集
boxplot(exp1)
boxplot(exp2)
dev.off()

# 提取临床信息
pd1 <- pData(eset1[[1]])
pd2 <- pData(eset2[[1]])
index1 <- eset1[[1]]@annotation
index2 <- eset2[[1]]@annotation
#### exp1 ID转换 ####
BiocManager::install("hugene10sttranscriptcluster.db")#GPL6244
library(hugene10sttranscriptcluster.db)
ids <- toTable(hugene10sttranscriptclusterSYMBOL) # ids = getIDs('gpl6244')
exp1 <- data.frame(exp1)
exp1 <- exp1 %>% mutate(probe_id=rownames(exp1))
exp1 <- exp1 %>% inner_join(ids, by="probe_id")
exp1 <- exp1[!duplicated(exp1$symbol),] #去重
rownames(exp1) <- exp1$symbol
exp <- exp[, -(33:34)]
#### exp2 ID转化 ####
options(stringsAsFactors = F)
gse <- getGEO(filename = "GSE132956_family.soft", destdir = ".")
str(gse)
length(gse)

id_probe <- gse@gpls$GPL17586@dataTable@table
dim(id_probe)
head(id_probe)
id_probe[1:4, 1:15]
probe2gene <- id_probe[, c(2, 8)]
library(stringr)
probe2gene$symbol = trimws(str_split(probe2gene$gene_assignment, '//', simplify = T)[, 2])
head(probe2gene)
ids2 <- probe2gene[, c(1, 3)]
exp2 <- data.frame(exp2)
exp2 <- exp2 %>% mutate(probe_id=rownames(exp2))
colnames(ids2)
names(ids2)[1] <- "probe_id"
exp2 <- exp2 %>% inner_join(ids2, by="probe_id")
exp2 <- exp2[!duplicated(exp2$symbol),] #去重
rownames(exp2) <- exp2$symbol
exp2 <- exp2[, -(16:17)]
####  合并 ####
exp = cbind(exp1, exp2)
