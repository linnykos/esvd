rm(list=ls())
load("../../raw_data/marques.RData")

dat <- marques$counts
colnames(dat) <- gsub("_", "-", colnames(dat))
dat_count <- dat

cell_type_vec <- as.character(marques$cell.info$cell.type)
cell_type_vec <- as.factor(cell_type_vec)
dim(dat)

# # remove genes with too many 0's
zz <- apply(dat, 2, function(x){length(which(x!=0))})
idx <- which(zz > nrow(dat)/100)
dat <- dat[,idx]
dat_count <- dat_count[,idx]

dat_count <- dat
# dat <- t(apply(dat, 1, function(x){10^4 * x/sum(x)}))

dat <- t(as.matrix(obj@assays$RNA@data))
object <- as(object = as.matrix(x = t(dat)), Class = 'Matrix')
feature.mean <- Seurat:::FastExpMean(object, T)
feature.dispersion <-  Seurat:::FastLogVMR(object, T)
which(is.na(feature.dispersion))
which(is.na(feature.mean))
plot(feature.mean, feature.dispersion)

obj <- Seurat::CreateSeuratObject(counts = t(dat), project = "marques",
                                  meta.data = NULL, min.cells = 0, min.features = 0)
obj <- Seurat::FindVariableFeatures(obj, selection.method = "mean.var.plot",
                                    mean.cutoff = c(0.1, Inf),
                                    dispersion.cutoff = c(1.5, Inf))
Seurat::VariableFeaturePlot(obj)
gene_name_hvg <- Seurat::VariableFeatures(object = obj)
length(gene_name_hvg)

#################

obj <- Seurat::CreateSeuratObject(counts = t(dat), project = "marques",
                                  meta.data = NULL, min.cells = 0, min.features = 0)
obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = 1000,
                                    mean.cutoff = c(0.1, 8))
# obj <- Seurat::FindVariableFeatures(obj, selection.method = "mean.var.plot", nfeatures = 1000)
gene_name_hvg2 <- Seurat::VariableFeatures(object = obj)

length(intersect(gene_name_org, gene_name_hvg2))


