# DEG analysis
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
library(limma)
library(DESeq2)
library(edgeR)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(data.table)
library(clusterProfiler)
library(DOSE)
library(GO.db)
library(topGO)
library(GSEABase)
library(stringr)
set.seed(1)
a = read.csv('.GSE75214_Expr_ok.csv', header = T, sep = ',', quote = '', stringsAsFactors = FALSE, check.names = FALSE)
a=as.data.frame(a)
GPL = read.csv('GPL_GSE75214.csv', header = F,sep = ',', quote = '', stringsAsFactors = FALSE, check.names = FALSE)  
GPL[1:6,1:4]
GPL = GPL[,c(1,2)]
GPL$V2 = sapply(strsplit(as.character(GPL$V2),'//'), "[", 2)
head(GPL)
GPL$V2 = str_replace_all(GPL$V2, fixed(" "), "")
colnames(GPL) = c("ID_REF","Gene_Symbol")
GPL = GPL[!(GPL$Gene_Symbol %in% c(NA)), ]
a = merge(GPL, a, by="ID_REF")
a = a[,-1]
Expr_count = avereps(a[,-1],ID = a$Gene_Symbol)
Expr_count = as.data.frame(Expr_count)
Expr_count[1:6,1:4]
Expr_count = Expr_count[!(rownames(Expr_count) %in% "---"),]
pdata = read.csv('../data/GSE75214_meta.csv',row.names=1,header = T,sep = ',',quote = '', check.names = FALSE)
pdata <- as.data.frame(pdata)
pdata[1:3,]
table(pdata$group)
pp_ahead = pdata[pdata$group %in% c("controle_ileum","CD_ileum_active"),]
pp = pp_ahead[,1]
pp = as.data.frame(pp)
rownames(pp) = rownames(pp_ahead)
colnames(pp) = "type" 
Expr = Expr_count[,rownames(pp)]
group_list = as.character(pp[, 1])
table(group_list)
par(cex = 0.7)
n.sample=ncol(Expr)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
boxplot(Expr, col = cols,main="expression value",las=2)
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(Expr)
design
contrast.matrix = makeContrasts(contrasts = c('CD_ileum_active-controle_ileum'), levels = design)
fit <- lmFit(Expr,design)
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2)  
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 
write.csv(nrDEG, file = "CD_ileum_active-controle_ileum.csv", quote=F, row.names=T)
