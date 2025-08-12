Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
library(ggplot2)
library(tidyverse)
library(ggrepel)
set.seed(1)
CD13_str1 = read.csv('nrDEG_CD13_CD_imm_str3.csv', header = T,row.names=1, sep = ',', quote = '', stringsAsFactors = FALSE, check.names = FALSE)  
CD13_str1$sample = "CD_imm_vs_CD_str3"
CD13_str1$geneID = rownames(CD13_str1)
CD13_str2 = read.csv('nrDEG_CD2_CD_rou_str1.csv', header = T,row.names=1, sep = ',', quote = '', stringsAsFactors = FALSE, check.names = FALSE)  
CD13_str2$sample = "CD_rou_vs_CD_str1"
CD13_str2$geneID = rownames(CD13_str2)
df = rbind(CD13_str1,CD13_str2)
df$label <- ifelse((df$P.Value<0.05 &  abs(df$logFC)>0.75),"P_value<0.05","P_value>=0.05")
head(df)
top10sig1 <- filter(df,sample=="CD_imm_vs_CD_str3") %>% distinct(geneID,.keep_all = T) %>% top_n(10,(logFC))
top10sig2 <- filter(df,sample=="CD_rou_vs_CD_str1") %>% distinct(geneID,.keep_all = T) %>% top_n(10,(logFC))
top10sig <- rbind(top10sig1,top10sig2)
df$size <- case_when(!(df$geneID %in% top10sig$geneID)~ 1,
                     df$geneID %in% top10sig$geneID ~ 2)
dt <- filter(df,size==1)
head(dt)

p <- ggplot()+
  geom_jitter(data = dt,
              aes(x = sample, y = logFC, color = label),
              size = 0.85,
              width =0.4)
p

p <- ggplot()+
  geom_jitter(data = dt,
              aes(x = sample, y = logFC, color = label),
              size = 0.5,
              width =0.4)+
  geom_jitter(data = top10sig,
              aes(x = sample, y = logFC, color = label),
              size = 2,
              width =0.4)
p

dfbar<-data.frame(x=c("CD_imm_vs_CD_str3",
                      "CD_rou_vs_CD_str1"),
                  y=c(max(top10sig1$logFC),
                      max(top10sig2$logFC)))
dfbar1<-data.frame(x=c("CD_imm_vs_CD_str3",
                       "CD_rou_vs_CD_str1"),
                   y=c(min(CD13_str1$logFC),
                       min(CD13_str2$logFC)))

p1 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),width =0.8,
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,
           mapping = aes(x = x,y = y),width =0.8,
           fill = "#dcdcdc",alpha = 0.6)
p1

p2 <- ggplot()+
  geom_jitter(data = dt,
              aes(x = sample, y = logFC, color = label),
              size = 0.5,
              width =0.4)+
  geom_jitter(data = top10sig,
              aes(x = sample, y = logFC, color = label),
              size = 2,
              width =0.4)
p2

dfcol<-data.frame(x=c("CD_imm_vs_CD_str3",
                      "CD_rou_vs_CD_str1"),
                  y=0,
                  label=c("CD_imm_vs_CD_str3",
                          "CD_rou_vs_CD_str1"))

dfcol<-data.frame(x=c(" ",
                      " "),
                  y=0,
                  label=c(" ",
                          " "))

mycol <- c("#E64B357F","#4DBBD57F","#00A0877F","#3C54887F")
mycol <- c("#E64B357F","#4DBBD57F")
p3 <- p2 + geom_tile(data = dfcol,
                     aes(x=x,y=y),
                     height=1.5,
                     width =0.8,
                     color = "black",
                     fill = mycol,
                     alpha = 0.6,
                     show.legend = F)
p3

p4 <- p3 +
  scale_color_manual(name=NULL,
                     values = c("red","black"))
p4

p5 <- p4+
  labs(x="",y="log2FC")+
  geom_text(data=dfcol,
            aes(x=x,y=y,label=label),
            size =3.5,
            color ="white")
p5

p6 <- p5+
  theme_minimal()+
  theme(
    axis.title = element_text(size = 13,
                              color = "black",
                              face = "bold"),
    axis.line.y = element_line(color = "black",
                               size = 1.2),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1,0),
    legend.text = element_text(size = 15)
  )
p6
