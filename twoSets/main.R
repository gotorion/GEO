setwd("twoSets")
#### install packages ####
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
#### library ####
library(BiocManager)
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
library(GOplot)
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
boxplot(exp1, outline=F, notch=F, las=2, col='blue')
boxplot(exp2)
boxplot(exp3)
dev.off()
#如果数据不均一
#exp1 <- normalizeBetweenArrays(exp1)
#exp2 <- normalizeBetweenArrays(exp2)
#exp3 <- normalizeBetweenArrays(exp3)
#exp4 <- normalizeBetweenArrays(exp4)
#### 提取临床信息 ####
pd1 <- pData(eset1[[1]])
pd2 <- pData(eset2[[1]])
pd3 <- pData(eset3[[1]])
pd4 <- pData(eset4[[1]])
#index1 <- eset1[[1]]@annotation
#index2 <- eset2[[1]]@annotation
#index3 <- eset3[[1]]@annotation
#index4 <- eset4[[1]]@annotation
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

write.table(expr_limma, file = "expr_limma.csv", sep = "\t", row.names = T, col.names = NA, quote = F)

#expr_limma <- normalizeBetweenArrays(expr_limma)
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
                         palette = c("#00AFBB", "#FC4E07"),
                         addEllipses = TRUE,
                         legend.title = "Groups",
                         title = "PCA - Biplot"
)
pca_plot
dev.off()
#### 差异分析 ####
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
head(deg)
install.packages("ggpubr")
install.packages("ggthemes")
install.packages("ggrepel")
install.packages("readxl")
library(ggpubr)
library(ggthemes)
library(ggrepel)
library(readxl)

demo <- read_excel("demo.xlsx")
this_tile <- paste0('Volcano plot',
                    '\nCutoff for LogFC is ',round(logFC,3),
                    '\nThe number of up genes is ',nrow(demo[demo$change =='up',]) ,
                    '\nThe number of down genes is ',nrow(demo[demo$change =='down',])
)
demo$label <- ifelse(demo$adj.P.Val< 0.0005& abs(demo$logFC) >= 2.5,demo$gene,"")
volcano_map <- ggplot(demo, aes(x=logFC, y=-log10(adj.P.Val),color=change)) + 
  geom_point(alpha=0.4, size=2) + 
  theme_bw(base_size = 12) + 
  xlab("Log2(Fold change)") +
  ylab("-Log10(P.adj)") +
  theme(plot.title = element_text(size=13,hjust = 0.5),
        legend.position = "bottom") + 
  scale_colour_manual(values = c('steelblue','gray','brown')) +
  geom_hline(yintercept = -log10(0.05), lty = 4) +
  geom_vline(xintercept = c(-logFC, logFC), lty = 4)+
  labs(title = this_tile)+
  geom_label_repel(data = demo, aes(label = label),
                   size = 3,box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.8, "lines"),
                   segment.color = "black",
                   show.legend = FALSE, max.overlaps = 10000)
volcano_map
dev.off()
#### Go 富集分析 ####
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("GOplot")#绘图
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
#deg <- read.table("deg_all.csv", sep = "\t", row.names = 1, check.names = F, stringsAsFactors = F, header = T)

deg_GO <- deg %>% filter(change != "stable")
table(deg_GO$change)
DEG <- deg_GO
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
#弦图
#write.table(ego_res, file = "go_result.csv", sep = "\t", row.names = T, col.names = NA, quote = F)
#write.table(DEG, file = "differential_expr", sep = "\t", row.names = T, col.names = NA, quote = F)
install.packages("RColorBrewer")
library(RColorBrewer)
df1 <- read.table("go_result.csv", header = T, sep = ",", stringsAsFactors = F)
df2 <- read.csv("differential_expr.csv", header = T, stringsAsFactors = F)
circ <- circle_dat(df1, df2)
df3 <- read.csv("genes.csv", header = T, stringsAsFactors = F)
df4 <- read.csv("process.csv", header = T, stringsAsFactors = F)
chord <- chord_dat(circ, df3, df4$BP)
GOChord(chord,   #chord对象
        title = "GO-Biological Process",
        limit = c(2, 0),
        ribbon.col = brewer.pal(8,"Set2"),
        space = 0.02,  #右侧色块之间的间距
        gene.order = 'logFC',   #基因展示顺序根据logFC来
        gene.space = 0.25,  #基因名字和色块之间的距离
        gene.size = 4  #基因名字大小
)
dev.off()

#### KEGG富集分析 ####
kk <- enrichKEGG(gene = DEG$ENTREZID,
                 organism = 'hsa', #Human sapiens
                 pvalueCutoff = 0.1,
                 qvalueCutoff = 0.1)
kk_res <- kk@result
dotplot(kk)


#### 共表达分析 ####
BiocManager::install("preprocessCore")
BiocManager::install("impute")
install.packages("WGCNA")
library("tidyverse")
library("WGCNA")            
if(!require(DESeq2))BiocManager::install('DESeq2')
library(DESeq2)
input <- expr_limma[rownames(my_deg),]
###保证样本在行，基因在列 很重要！！！！
datExpr0 = as.data.frame(t(input))
##开始WGCNA
#检查缺失值
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
###如果没有达标就需要筛选
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
# 样品聚类
# 聚类
sampleTree = hclust(dist(datExpr0), method = "average")
# 画图
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
###剪切线，是否需要剪切？
abline(h = 67, col = "red")
###删除剪切线以下的样品
clust = cutreeStatic(sampleTree, cutHeight = 67, minSize = 10)
table(clust)
keepSamples = (clust==1)
datExpr0 = datExpr0[keepSamples, ]
dev.off()
# 剪切完重新聚类
sampleTree2 = hclust(dist(datExpr0), method = "average")
plot(sampleTree2)

# 记录基因和样本数，方便后续可视化
nGenes = ncol(datExpr0)#基因数
nSamples = nrow(datExpr0)#样本数
save(datExpr0, nGenes, nSamples,file = "Step01-WGCNA_input.Rda")

# 构建网络，识别模块
# power值散点图
enableWGCNAThreads()   #多线程工作
powers = c(1:20)       #幂指数范围1:20
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)

par(mfrow = c(1,2))
cex1 = 0.9
###拟合指数与power值散点图
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red") #可以修改
###平均连通性与power值散点图
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

##6.2 邻接矩阵转换
sft #查看最佳power值
softPower =sft$powerEstimate #最佳power值
softPower = 6
adjacency = adjacency(datExpr0, power = softPower)

##6.3 TOM矩阵
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM
save(TOM,file = "TOM.Rda")

# 基因聚类
geneTree = hclust(as.dist(dissTOM), method = "average");
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

# 动态剪切模块识别
minModuleSize = 30      #模块基因数目
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# 相似模块聚类
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

#MEDissThres = 0.1 #剪切高度可修改
#abline(h=MEDissThres, col = "red")

###相似模块合并
#merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
#mergedColors = merge$colors
#mergedMEs = merge$newMEs
#plotDendroAndColors(geneTree, mergedColors,"Dynamic Tree Cut",
#                    dendroLabels = FALSE, hang = 0.03,
#                    addGuide = TRUE, guideHang = 0.05,
#                    main = "Gene dendrogram and module colors")

#moduleColors = mergedColors
##table(moduleColors)
#colorOrder = c("grey", standardColors(50))
#moduleLabels = match(moduleColors, colorOrder)-1
#MEs = mergedMEs
#dev.off()
#### 计算免疫评分 ####
install.packages("utils")
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
#读取肿瘤患者01A表达谱
expr <- read.table("expr_limma.csv",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)


#计算免疫评分
filterCommonGenes(input.f = "expr_limma.csv",   #输入文件名
                  output.f = "LIHC_fpkm_mRNA_01A.gct",   #输出文件名
                  id = "GeneSymbol")   #行名为gene symbol
estimateScore("LIHC_fpkm_mRNA_01A.gct",   #刚才的输出文件名
              "LIHC_fpkm_mRNA_01A_estimate_score.txt",   #新的输出文件名（即估计的结果文件）
              platform="affymetrix")   #默认平台

#3. 输出每个样品的打分
result <- read.table("LIHC_fpkm_mRNA_01A_estimate_score.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
result <- result[,-1]   
colnames(result) <- result[1,]   
result <- as.data.frame(t(result[-1,]))

rownames(result) <- colnames(expr)
write.table(result, file = "LIHC_fpkm_mRNA_01A_estimate_score.txt",sep = "\t",row.names = T,col.names = NA,quote = F) # 保存并覆盖得分


# 整理临床信息
clinical <- read.table("LIHC_fpkm_mRNA_01A_estimate_score.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
clinical <- clinical[rownames(datExpr0),]
identical(rownames(clinical),rownames(datExpr0))
# 查看临床信息
head(clinical)
# 对表达矩阵进行预处理
datTraits = as.data.frame(do.call(cbind,lapply(clinical, as.numeric)))
rownames(datTraits) = rownames(clinical)

# 对样本进行聚类
sampleTree2 = hclust(dist(datExpr0), method = "average")

# 将临床信息转换为颜色，白色表示低，红色表示高，灰色表示缺失
traitColors = numbers2colors(datTraits, signed = FALSE)

# 样本聚类图与样本性状热图
plotDendroAndColors(sampleTree2, 
                    traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()

#### 网络的分析
# 对模块特征矩阵进行排序
MEs=orderMEs(MEs)
#计算模型特征矩阵和样本信息矩阵的相关度。
moduleTraitCor=cor(MEs, datTraits, use="p")
write.table(file="Step04-modPhysiological.cor.xls",moduleTraitCor,sep="\t",quote=F)
moduleTraitPvalue=corPvalueStudent(moduleTraitCor, nSamples)
write.table(file="Step04-modPhysiological.p.xls",moduleTraitPvalue,sep="\t",quote=F)

#使用labeledHeatmap()将上述相关矩阵和p值可视化。
textMatrix=paste(signif(moduleTraitCor,2),"\n(",signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix)=dim(moduleTraitCor)
# 基因模块与临床信息相关性图
labeledHeatmap(Matrix=moduleTraitCor,#模块和表型的相关性矩阵，这个参数最重要，其他可以不变
               xLabels=colnames(datTraits),
               yLabels=names(MEs),
               ySymbols=names(MEs),
               colorLabels=FALSE,
               colors=blueWhiteRed(50),
               textMatrix=textMatrix,
               setStdMargins=FALSE,
               cex.text=0.7,
               cex.lab=0.7,
               zlim=c(-1,1),
               main=paste("Module-trait relationships"))
dev.off()
# 不同模块与基因性状的具体分析
##矩阵一
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
####看一下目的基因和哪个模块相关性最高
a <- geneModuleMembership
a <- a %>% rownames_to_column()

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

##矩阵二
traitNames=names(datTraits)
geneTraitSignificance = as.data.frame(cor(datExpr0, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")

##批量输出性状和模块散点图
for (trait in traitNames){
  traitColumn=match(trait,traitNames)  
  for (module in modNames){
    column = match(module, modNames)
    moduleGenes = moduleColors==module
    if (nrow(geneModuleMembership[moduleGenes,]) > 1){
      outPdf=paste(trait, "_", module,".pdf",sep="")
      pdf(file=outPdf,width=7,height=7)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = paste("Gene significance for ",trait),
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      abline(v=0.8,h=0.5,col="red")
      dev.off()
    }
  }
}

#10. 输出每个模块的基因
for (mod in 1:nrow(table(moduleColors)))
{  
  modules = names(table(moduleColors))[mod]
  probes = colnames(datExpr0)
  inModule = (moduleColors == modules)
  modGenes = probes[inModule]
  write.table(modGenes, file =paste0(modules,".txt"),sep="\t",row.names=F,col.names=F,quote=F)
}
write.table(a, file = "to_string.csv",sep="\t",row.names=F,col.names=F,quote=F)

#### ROC ####
#准备R包
install.packages("ROCR")
install.packages("rms")
library(ROCR)
library(rms)
roc_key_genes <- read.table("roc_key_genes.csv", sep = ",", row.names = 1, check.names = F, stringsAsFactors = F, header = T)
roc_attribute <- read.table("roc_attribute.csv", sep = ",", row.names = 1, check.names = F, stringsAsFactors = F, header = T)
x <- unlist(roc_key_genes$`False Positive Rate`)
y <- unlist(roc_key_genes$`True Positive Rate`)
plotdata <- data.frame(x, y)
names(plotdata) <- c("x", "y")

m <- unlist(roc_attribute$`False Positive Rate`)
n <- unlist(roc_attribute$`True Positive Rate`)
plotdata1 <- data.frame(m, n)
names(plotdata1) <- c("m", "n")

library(ggplot2)
g <- ggplot(plotdata) + 
  geom_path(aes(x = x, y = y, colour = x), linewidth = 1) +
  labs(x = "False Positive Rate", y = "Sensitivity", title = "ROC curve for key genes") +
  scale_color_gradient(name = "False Positive Rate", low = 'blue', high = 'red') +
  theme(plot.title = element_text(face = 'bold', size = 15))+
  annotate("text", x = 0.6, y = 0.8, label = "AUC = 0.943", size = 5)
  
g2 <- ggplot(plotdata1) + 
  geom_path(aes(m,n, colour = m), linewidth = 1) +
  labs(x = "False Positive Rate", y = "Sensitivity", title = "ROC curve for filter result") +
  scale_color_gradient(name = "False Positive Rate", low = 'blue', high = 'red') +
  theme(plot.title = element_text(face = 'bold', size = 15))+
  annotate("text", x = 0.6, y = 0.8, label = "AUC = 0.950", size = 5)

g
g2
dev.off()
#1.4 绘制ROC曲线

AUC <- 0.943

plot(roc_key_genes$`False Positive Rate`,roc_key_genes$`True Positive Rate`,
     col="red",   #曲线的颜色
     xlab="False positive rate", ylab="True positive rate",   #x轴和y轴的名称
     lty=1,lwd=2,
     main=paste("AUC=",AUC))
abline(0, 1, lty=2, lwd=2)   #绘制对角线
legend(0.60, 0.10,
       c("AUC = 0.943"),
       lty = 1,
       lwd = 2,
       col = 'red',
       bty = 'o')
dev.off()
