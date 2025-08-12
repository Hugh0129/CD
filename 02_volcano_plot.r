# volcano plot
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
library(ggplot2)
library(tidyverse)
library(ggrepel)
set.seed(1)
data = read.csv('CD_ileum_active-controle_ileum.csv',  header = T,row.names=1, sep = ',',  stringsAsFactors = FALSE, check.names = FALSE)  
data = as.data.frame(data)
data = data[,c(1,4)]
colnames(data) = c("log2FC","Pvalue")
data$n = rownames(data)
m=ggplot(data=data, aes(x=log2FC, y=-log10(adj.P))) +
  geom_point(data=subset(data,data$Pvalue >= 0.05|abs(data$log2FC) <= 1),aes(size=abs(log2FC)),color="gray",alpha=0.5) +
  geom_point(data=subset(data,data$Pvalue<0.05 & data$log2FC > 1),aes(size=abs(log2FC)),color="firebrick3",alpha=0.5) +
  geom_point(data=subset(data,data$Pvalue<0.05 & data$log2FC < -1),aes(size=abs(log2FC)),color="navy",alpha=0.5) +
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+ 
  labs(x="Log2(FoldChange)",y="-Log10(Pvalue)")+
  theme(legend.position='none')
m
m+xlim(-3,3)+ ylim(0,8)#+ labs(title="Volcanoplot",x=expression(logFC),y=expression(-log10(P.Value)))
