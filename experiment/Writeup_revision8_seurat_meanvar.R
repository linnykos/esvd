rm(list=ls())
load("../../raw_data/marques.RData")

dat <- marques$counts
obj <- Seurat::CreateSeuratObject(counts = t(dat), project = "marques",
                                  meta.data = NULL, min.cells = 0, min.features = 0)
obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10^4) # normalize log2(x*10^6+1)

obj <- Seurat::FindVariableFeatures(obj, selection.method = "mean.var.plot",
                                    mean.cutoff = c(0.0125, 3),
                                    dispersion.cutoff = c(1, Inf),
                                    verbose=F)
Seurat::VariableFeaturePlot(obj)


#########

obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
Seurat::VariableFeaturePlot(obj)
