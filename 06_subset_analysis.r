Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(ArchR)
library(viridis)
library(DoubletFinder)
library(Rcpp)
library(harmony)
library(future)
library(Matrix)
library(Rmagic)
library(ggpubr)
set.seed(1)
analysis_parent_folder <- "str"
setwd(analysis_parent_folder)
colon = readRDS("CD_stromal_all_initial.rds")
# Define variables
sample_name <- "all_samples" 
execute_steps <- c(1,2)
vln_plot <- function(features, save_name){
	pdf(save_name, width = 20, onefile=F)
	print(VlnPlot(colon, features = features, group.by = "orig.ident", pt.size = 0, cols = c(rep("#D51F26",100)))+
	geom_boxplot(outlier.shape = NA, alpha = 0.6)+theme_ArchR()+theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))+
	scale_x_discrete(labels=paste0(data.frame(table(colon@meta.data$orig.ident))$Var1, "\n n = ",  data.frame(table(colon@meta.data$orig.ident))$Freq)))
	dev.off()
}
normalize_and_dim_reduce <- function (colon, sample_name){
	# stanfard seurat pipeline
	colon <- NormalizeData(colon, normalization.method = "LogNormalize", scale.factor = 10000)
	colon <- FindVariableFeatures(colon, selection.method = "vst", nfeatures = 2000)
	all.genes <- rownames(colon)
	colon <- ScaleData(colon, features = all.genes)
	colon <- RunPCA(colon, features = VariableFeatures(object = colon))
	colon <- FindNeighbors(colon, dims = 1:30)
	colon <- FindClusters(colon, resolution = 1)
	colon <- RunUMAP(colon, dims = 1:30)
	return(colon)
}
plotUMAPandRunHarmony <- function(colon, run_harmony, version){
	if (run_harmony){
		library(harmony)
		colon <- RunHarmony(colon, "orig.ident") # will use PCA
		colon <- RunUMAP(colon, dims = 1:30, reduction = "harmony", reduction.name = "umapharmony")
		colon <- FindNeighbors(colon, reduction = "harmony", dims = 1:20)
		colon <- FindClusters(colon, resolution = 1.0)
	}
	return(colon)
}
# 1) Normalize and scale data
if (1 %in% execute_steps){
	colon <- normalize_and_dim_reduce(colon, sample_name)
}
# 2) Plot UMAPs and run harmony
if (3 %in% execute_steps){
	colon <- plotUMAPandRunHarmony(colon, TRUE, version = "initial")
}
DimPlot(colon, reduction = "umap",label = T, cols = paletteDiscrete(values = unique(colon@meta.data$seurat_clusters), set = "stallion", reverse = FALSE))
DotPlot(colon, 
        features =c("COL1A1","COL3A1", #Fib_clusters
                     "ACTA2", "DCN", # Fib_clusters
                    "DARC", "CD36", "RGCC", "RELN", "EFNA5", "ATP2A3", "EFNB2", # Endo_clusters
                    "HIGD1B", "STEAP4","RERGL","NTRK2", # Pericytes_clusters
                    "S100B","GPM6B", "PLP1"),  # Glial_clusters
        dot.scale = 8) + RotatedAxis()

colon@meta.data$celltype_str = "NA"
Fib_clusters <- c(0:4,6,7,8,12,14,17,18,20,21,22,24)
Endo_clusters <- c(9,10,11,15,16,19,23)
Pericytes_clusters <- c(13)
Glial_clusters <- c(5)
colon@meta.data[which(colon@meta.data$seurat_clusters %in% Fib_clusters),'celltype_str'] <- "Fibroblasts"
colon@meta.data[which(colon@meta.data$seurat_clusters %in% Endo_clusters),'celltype_str'] <- "Endothelial_cells"
colon@meta.data[which(colon@meta.data$seurat_clusters %in% Pericytes_clusters),'celltype_str'] <- "Pericytes"
colon@meta.data[which(colon@meta.data$seurat_clusters %in% Glial_clusters),'celltype_str'] <- "Glial_cells"
Idents(colon) = "celltype_str"
p1 = DimPlot(colon, reduction = "umap",label = F, cols = paletteDiscrete(values = unique(colon@meta.data$celltype_str), set = "stallion", reverse = F))
p1
p2 = FeaturePlot(colon,reduction = "umap",features = c("SPP1"),sort.cell = TRUE,pt.size = 1,raster=FALSE,label = F)
p2
p3 = DotPlot(colon,  features =c("COL1A1", "COL3A1", "DCN", "DARC", "CD36", "RGCC",  "HIGD1B", "STEAP4", "RERGL", "S100B", "GPM6B", "PLP1"), dot.scale = 8) + RotatedAxis()
p3
p4 = FeaturePlot(colon, reduction = "umap", features = c("SPP1"),sort.cell = TRUE,split.by = "type", pt.size = 1, label = F)
p4 

cell.prop<-as.data.frame(prop.table(table(colon@meta.data$celltype_str,colon$type)))
colnames(cell.prop)<-c("cluster","sample","proportion")
p = ggplot(cell.prop,aes(sample,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  labs(x = '', y = 'Relative Abundance(%)') +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13)) +
  theme(legend.text = element_text(size = 11))
cell_type_cols <-paletteDiscrete(values = unique(colon@meta.data$celltype_str), set = "stallion", reverse = FALSE)
p <- p + scale_fill_manual(values = (cell_type_cols)) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank())+
  guides(color = guide_legend(ncol = 1, byrow = TRUE,reverse = T))+
  theme(axis.title.y = element_text(face = 'plain',color = 'black',size = 16),
        axis.title.x = element_text(face = 'plain',color = 'black',size = 10),
        axis.text.y = element_text(face = 'plain',color = 'black',size = 16),
        axis.text.x = element_text(face = 'plain',color = 'black',
                                   size = 16,angle = 90,vjust = 0.5, hjust=0), 
        axis.ticks.length=unit(.1,"lines"),
        axis.ticks.x = element_line(size=0.05, colour = "black"),
        #axis.ticks.margin=unit(.4,"cm"),
        panel.border = element_blank(),
        axis.line = element_line(size=0.1, colour = "black"),
        #axis.ticks.x.bottom =  = element_line(size = 0.5),
        panel.grid = element_blank(),
        #scale_fill_distiller(palette = "Spectral"),
        #scale_fill_brewer(palette = 'Paired'),
        #scale_fill_manual(values=ccc),
        legend.position = 'right',
        legend.key.width = unit(0.2,'cm'),
        legend.key.height = unit(0.2,'cm'),
        legend.text = element_text(color = 'black',size = 16))
p

Fibroblasts <- c("Fibroblasts")
Endothelial_cells <- c("Endothelial_cells")
Pericytes <- c("Pericytes")
Glial_cells <- c("Glial_cells")
# Subset to make immune, stromal, and epithelial projects
colon_Fibroblasts <- DietSeurat(subset(colon, subset = celltype_str %in% Fibroblasts))
colon_Endothelial_cells <- DietSeurat(subset(colon, subset = celltype_str %in% Endothelial_cells))
colon_Pericytes <- DietSeurat(subset(colon, subset = celltype_str %in% Pericytes))
colon_Glial_cells <- DietSeurat(subset(colon, subset = celltype_str %in% Glial_cells))
saveRDS(colon_Fibroblasts, file = "CD_str_Fibroblasts_initial.rds")
saveRDS(colon_Endothelial_cells, file = "CD_str_Endothelial_initial.rds")
saveRDS(colon_Pericytes, file = "CD_str_Pericytes_initial.rds")
saveRDS(colon_Glial_cells, file = "CD_str_Glial_cells_initial.rds")
