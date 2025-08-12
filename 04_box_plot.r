# box plot
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
set.seed(1)
a = read.csv('mRNA_matrix_RMA.csv', header = T,row.names=1, sep = ',', quote = '', stringsAsFactors = FALSE, check.names = FALSE)  
a = as.data.frame((a))
data = a[c('CD68','SPP1'),]
data = as.data.frame(t(data))
pdata = read.csv('../data_box/GSE75214_meta.csv', header = T,row.names=1,  sep = ',', quote = '', stringsAsFactors = FALSE, check.names = FALSE)  
b = c("controle_ileum","CD_ileum_active")
pdata = pdata[pdata$group%in%b,]
data = data[rownames(data)%in%rownames(pdata),]
pdata = pdata[rownames(pdata)%in%rownames(data),]
data = data[rownames(pdata),]
all(rownames(data)==rownames(pdata))
plot_data = cbind(pdata,data)
plot_data = plot_data[,c(1,4)]
table(plot_data$group)
colnames(plot_data) = c("group","Retive_Abundance")
levels(factor(plot_data$group))
plot_data$group <- factor(plot_data$group,levels = c("controle_ileum","CD_ileum_active"))
p<- ggplot(data=plot_data)+ 
  geom_boxplot(mapping=aes(x=group,y=Retive_Abundance,colour = group ), 
               alpha = 0.5,
               size=0.75,
               width = 0.6)+ 
  geom_jitter(mapping=aes(x=group,y=Retive_Abundance,colour = group), 
              width =0.2,alpha = 1,size=1)+
  scale_color_manual(limits=c("controle_ileum","CD_ileum_active"), 
                     values=c("#85B22E","#922927","#85B22E","#5F80B4","#922927","#922927"))+ 
  geom_signif(mapping=aes(x=group,y=Retive_Abundance), 
              comparisons = list(c("controle_ileum","CD_ileum_active")),
              map_signif_level=T, 
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), 
              y_position = c(13), 
              size=0.75, 
              textsize = 4, 
              test = wilcox.test)+  theme_classic(  base_line_size = 0.75  )+
  labs(title="SPP1",x="",y="RMA normalized expression")+ 
  theme(plot.title = element_text(size = 15,
                                  colour = "black",
                                  hjust = 0.5),
        axis.title.y = element_text(size = 15, 
                                    # family = "myFont", 
                                    color = "black",
                                    face = "plain", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.title = element_text(color= "black",
                                    size=15, 
                                    face= "plain"),
        legend.text = element_text(color= "black", 
                                   size = 10, 
                                   face = "plain"),
        axis.text.x = element_text(size = 13,  
                                   color = "black", 
                                   face = "plain",
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0), 
        axis.text.y = element_text(size = 13,  
                                   color = "black",
                                   face = "plain", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0) 
  )
p= p +  ylim(0,15)
p
