Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
set.seed(1)
library(tidyverse)
library(pheatmap)
library(CIBERSORT)
library(limma)
a = read.csv('Symbol_matrix_TPM.csv', header = T,sep = ',', quote = '', stringsAsFactors = FALSE, check.names = FALSE)  
exp  = avereps(a[,-1],ID = a[,1]) 
exp  = as.data.frame(exp)
exp2 = exp[rowSums(exp>=1)>=10,]
Expr_count = rownames_to_column(exp2 )
f = "ciber_process.Rdata"
if(!file.exists(f)){
  lm22f = system.file("extdata", "LM22.txt", package = "CIBERSORT")
  TME.results = cibersort(lm22f, 
                          "exp2.txt", 
                          perm = 1000, 
                          QN = T)
  save(TME.results,file = f)
}

# vioplot
library(vioplot)         
load("ciber_process.Rdata")
immuedata = as.data.frame(TME.results) 
immue = immuedata[immuedata$`P-value` < 0.05,]
pdata = read.csv('CD2_meta_ok.csv',row.names=1,header = T,sep = ',',quote = '',check.names = FALSE)
pdata <- as.data.frame(pdata)
rownames(pdata) = pdata$name
pdata1 = pdata[(pdata$type2=="CD2_Str1"),]
pdata2 = pdata[(pdata$type2 == "CD2_rouyazhong"),]
pdata = rbind(pdata1,pdata2)
immue = immue[rownames(pdata),]
normal= 3
CD=  3
vioplotdata = immue[,-c((ncol(immue)-2):ncol(immue))]
vioplotdata = vioplotdata[,which(colSums(vioplotdata) > 0)] 
vioplotdata = vioplotdata[rownames(pdata),]
vioplotdata = vioplotdata+0.000001
vioplotdata[vioplotdata==0] <- NA
vioplotdata[is.na(vioplotdata)] <- min(vioplotdata,na.rm = T)*0.01
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(vioplotdata))
y=c(1:ncol(vioplotdata))
plot(x,y,
     xlim=c(0,(ncol(vioplotdata)-1)*3),关
     ylim=c(min(vioplotdata),max(vioplotdata)+0.02),
     legend=TRUE,
     main="",xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")
for(i in 1:ncol(vioplotdata)){
  normalData=vioplotdata[1:normal,i]
  CDData=vioplotdata[(normal+1):(normal+CD),i]
  vioplot(normalData,at=3*(i-1),lty=1,add = T,col = ("#5390d9"))
  vioplot(CDData,at=3*(i-1)+1,lty=1,add = T,col =  ("#aa2116"))
  wilcoxTest=wilcox.test(normalData,CDData, exact = F,alternative = c("less"))
  p=round(wilcoxTest$p.value,3)
  mx=max(c(normalData,CDData))
  text(x=3*(i-1)+0.5,y=mx+0.02,labels=ifelse(p<0.05,paste0("p<0.05"),paste0("p=",p)),cex = 0.8)
  text(seq(1,(ncol(vioplotdata)-1)*3+1,3),-0.05,xpd = NA,labels=colnames(vioplotdata),cex = 1,srt = 45,pos=2)
}
