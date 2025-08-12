Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE)
library(patchwork)
library(dplyr)
library(hdf5r)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(cowplot)
library(dplyr)
library(BayesSpace)
library(ArchR)
set.seed(1)
spe2 = Read10X("filtered_feature_bc_matrix/")
image2 <- Read10X_Image(image.dir = file.path("data", "spatial"), filter.matrix = TRUE)
spe2 <- CreateSeuratObject(counts = spe2, assay = "Spatial")
image2 <- image2[Cells(x = spe2)]
DefaultAssay(spe2 = image2) <- "Spatial"
spe2[["slice1"]] <- image2
for (i in colnames((spe2@images$slice1@coordinates))) {
  spe2@images$slice1@coordinates[[i]] <- as.integer(spe2@images$slice1@coordinates[[i]])
}
SpatialFeaturePlot(spe2, features = "nFeature_Spatial")
spe2 = spe2
plot1 <- VlnPlot(spe2, features = "nCount_Spatial", pt.size = 0.5) + NoLegend()
plot2 <- SpatialFeaturePlot(spe2, features = "nCount_Spatial") + theme(legend.position = "right")
plot_grid(plot1, plot2)
spe2 <- SCTransform(spe2, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
spe2 <- GroupCorrelation(spe2, group.assay = "Spatial", assay = "SCT", slot = "scale.data", do.plot = FALSE)
SpatialFeaturePlot(spe2, features = "SPP1", pt.size.factor = 2.5)
SpatialFeaturePlot(spe2, features = "PLP1", pt.size.factor = 2.5)
SpatialFeaturePlot(spe2, features = "CD68", pt.size.factor = 2.5)
SpatialFeaturePlot(spe2, features = "COL1A1", pt.size.factor = 2.5)
SpatialFeaturePlot(spe2, features = "FN1", pt.size.factor = 2.5)
