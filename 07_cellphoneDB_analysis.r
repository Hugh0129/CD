# analysis in linux
cellphonedb method statistical_analysis meta_data_heal_epi_glial.txt count_heal_epi_glial.txt --database cellphone.db --counts-data=gene_name --output-path ./out  --threshold 0.01 --threads 8
cellphonedb plot heatmap_plot meta_data_heal_epi_glial.txt
cellphonedb plot dot_plot 

cellphonedb method statistical_analysis meta_data_noni_epi_glial.txt count_noni_epi_glial.txt --database cellphone.db --counts-data=gene_name --output-path ./out  --threshold 0.01 --threads 8
cellphonedb plot heatmap_plot meta_data_noni_epi_glial.txt
cellphonedb plot dot_plot 

cellphonedb method statistical_analysis meta_data_infl_epi_glial.txt count_infl_epi_glial.txt --database cellphone.db --counts-data=gene_name --output-path ./out  --threshold 0.01 --threads 8
cellphonedb plot heatmap_plot meta_data_infl_epi_glial.txt
cellphonedb plot dot_plot 

#===============================================================================================================================
#plot heatmap
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
library(tidyr)
library(reshape2) 
library(pheatmap)
library(RColorBrewer)
set.seed(1)

a1 = read.table('heal_count_network.txt', header = T, sep = '\t', quote = '',  stringsAsFactors = FALSE, check.names = FALSE)  
a1 = spread(a1,TARGET,count)
rownames(a1) = a1[,1]
a1 = a1[,-1]
a1 = a1[c("Enteroendocrine_cells_TA_cells", "Enterocytes", "Goblet_cells", "Tuft_cells", "Stem_cells","Glial_cells"),
        c("Enteroendocrine_cells_TA_cells", "Enterocytes", "Goblet_cells", "Tuft_cells", "Stem_cells","Glial_cells")]
a2 = read.table('noni_count_network.txt', header = T,sep = '\t', quote = '', stringsAsFactors = FALSE, check.names = FALSE)  
a2 = spread(a2,TARGET,count)
rownames(a2) = a2[,1]
a2 = a2[,-1]
a2 = a2[c("Enteroendocrine_cells_TA_cells", "Enterocytes", "Goblet_cells", "Tuft_cells","Paneth_cells", "Stem_cells","Glial_cells"),
        c("Enteroendocrine_cells_TA_cells", "Enterocytes", "Goblet_cells", "Tuft_cells","Paneth_cells", "Stem_cells","Glial_cells")]
a3 = read.table('infl_count_network.txt', header = T,sep = '\t', quote = '', stringsAsFactors = FALSE, check.names = FALSE)  
a3 = spread(a3,TARGET,count)
rownames(a3) = a3[,1]
a3 = a3[,-1]
a3 = a3[c("Enteroendocrine_cells_TA_cells", "Enterocytes", "Goblet_cells", "Tuft_cells","Paneth_cells", "Stem_cells","Glial_cells"),
        c("Enteroendocrine_cells_TA_cells", "Enterocytes", "Goblet_cells", "Tuft_cells","Paneth_cells", "Stem_cells","Glial_cells")]

#breaks
bk <- c(seq(1,500,by=1))
aa1 = pheatmap((a1), 
              breaks = bk,
              clustering_distance_cols = "correlation", 
              clustering_method="average",
              fontsize_row = 10,
              show_colnames = T, 
              show_rownames = T, 
              cluster_cols = F,
              cluster_row = F, 
              color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")[1:6]))(500)
)

aa2 = pheatmap((a2), 
               breaks = bk,
               clustering_distance_cols = "correlation", 
               clustering_method="average",
               fontsize_row = 10,
               show_colnames = T, 
               show_rownames = T, 
               cluster_cols = F,
               cluster_row = F, 
               color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")[1:6]))(500))

aa3 = pheatmap((a3), 
               breaks = bk,
               clustering_distance_cols = "correlation", 
               clustering_method="average",
               fontsize_row = 10,
               show_colnames = T, 
               show_rownames = T, 
               cluster_cols = F,
               cluster_row = F, 
               color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")[1:6]))(500))

#==============================================================================================================
#dotplot
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
library(ktplots)
library(Seurat)
library(SingleCellExperiment)
library(tidyverse)
library(data.table)
library(dplyr)
library(readr)
library(tidyverse)
set.seed(1)

# heal
mypvals <- read.table("./heal/out/pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
mymeans <- read.table("./heal/out/means.txt",header = T,sep = "\t",stringsAsFactors = F) 
kp = grepl(pattern = "Glial", colnames(mypvals)) & grepl(pattern = "Stem", colnames(mypvals)) 
pos = (1:ncol(mypvals))[kp] 
choose_pvalues <- mypvals[,c(c(1,5,6,8,9),pos)]
choose_means <- mymeans[,c(c(1,5,6,8,9),pos)]
logi <- apply(choose_pvalues[,6:ncol(choose_pvalues)]<0.05, 1, sum) 
choose_pvalues <- choose_pvalues[logi>=1,]
logi1 <- choose_pvalues$gene_a != ""
logi2 <- choose_pvalues$gene_b != ""
logi <- logi1 & logi2
choose_pvalues <- choose_pvalues[logi,]
choose_means <- choose_means[choose_means$id_cp_interaction %in% choose_pvalues$id_cp_interaction,]
meansdf <- choose_means %>% reshape2::melt()
meansdf <- data.frame(interacting_pair = paste0(meansdf$gene_a,"_",meansdf$gene_b),
                      CC = meansdf$variable,
                      means = meansdf$value)
pvalsdf <- choose_pvalues %>% reshape2::melt()
pvalsdf <- data.frame(interacting_pair = paste0(pvalsdf$gene_a,"_",pvalsdf$gene_b),
                      CC = pvalsdf$variable,
                      pvals = pvalsdf$value)
pvalsdf$joinlab<- paste0(pvalsdf$interacting_pair,"_",pvalsdf$CC)
meansdf$joinlab<- paste0(meansdf$interacting_pair,"_",meansdf$CC)
pldf <- merge(pvalsdf,meansdf,by = "joinlab")
# dotplot
summary((filter(pldf,means >0))$means)
head(pldf)
pcc =  pldf%>% filter(means >0.5) %>% 
  ggplot(aes(CC.x,interacting_pair.x) )+ 
  geom_point(aes(color=means,size=-log10(pvals+0.0001)) ) +
  scale_color_gradient2(high="Firebrick",mid = "Khaki1",low ="SlateBlue3",midpoint = 1  )+ 
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 0,hjust = 0.5,vjust = 0))+
  labs(x="Cell interacting",y="interacting pair")
pcc
