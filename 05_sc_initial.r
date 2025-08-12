# seurat analysis
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
set.seed(1)
setwd("/home/data/")
folders=list.files('./')
folders 
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                               project = folder )
})
sce.big <- merge(sceList[[1]], 
                 y = c(sceList[[2]],sceList[[3]]), 
                 add.cell.ids = c("epi", "imm", "str"),
                 project = "CD")
sce.big <- NormalizeData(sce.big, normalization.method = "LogNormalize", scale.factor = 10000) 
sce.big <- FindVariableFeatures(sce.big, selection.method = "vst", nfeatures = 2000) 
sce.big <- ScaleData(sce.big)
sce.big <- RunPCA(sce.big) 
sce.big <- FindNeighbors(sce.big, dims = 1:20) 
sce.big <- FindClusters(sce.big) 
sce.big <- RunUMAP(sce.big, dims = 1:20)
sce.big$cell_name = substr(colnames(sce.big),1,3)
Idents(sce.big) = "cell_name"
p1 = DimPlot(sce.big, reduction = "umap",label = F,raster=FALSE,
        cols = paletteDiscrete(values = unique(sce.big@meta.data$cell_name), set = "stallion", reverse = FALSE))
p1
p2 = FeaturePlot(sce.big,reduction = "umap",features = c("SPP1"),sort.cell = TRUE,split.by = "type", pt.size = 1, label = F )
p2
p3 = DotPlot(sce.big,  features =c("EPCAM",'KRT19', "PTPRC","CD3D",'CD79A', 'MZB1',"COL1A1","COL3A1","DCN","MYL9"),  dot.scale = 8, ) + RotatedAxis()
p3
cell.prop<-as.data.frame(prop.table(table(sce.big@meta.data$cell_name,sce.big$type)))
colnames(cell.prop)<-c("cluster","sample","proportion")
p = ggplot(cell.prop,aes(sample,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  labs(x = '', y = 'Relative Abundance(%)') +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13)) +
  theme(legend.text = element_text(size = 11))
cell_type_cols <-paletteDiscrete(values = unique(sce.big@meta.data$cell_name), set = "stallion", reverse = FALSE)
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
