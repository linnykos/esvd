library(descend)
load("../../raw_data/marques.RData")

dat <- marques$counts
colnames(dat) <- gsub("_", "-", colnames(dat))
dat_count <- dat

cell_type_vec <- as.character(marques$cell.info$cell.type)
cell_type_vec <- as.factor(cell_type_vec)
dim(dat)

# remove genes with too many 0's
zz <- apply(dat, 2, function(x){length(which(x!=0))})
idx <- which(zz > nrow(dat)/100)
dat <- dat[,idx]
dat_count <- dat_count[,idx]

obj <- Seurat::CreateSeuratObject(counts = t(dat), project = "marques",
                                  meta.data = NULL, min.cells = 0, min.features = 0)
obj <- Seurat::NormalizeData(obj, normalization.method = "RC", scale.factor = 10^4)
obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = 500)
# Seurat::VariableFeaturePlot(obj)
dat <- t(as.matrix(obj@assays$RNA@data))
vst_hvg <- Seurat::VariableFeatures(object = obj)

print(paste0(Sys.time(), ": Finished vst screening"))
save.image(paste0("../results/step0_screening", suffix, ".RData"))

# run DESCEND
res_descend <- descend::runDescend(t(dat_count), n.cores = ncores)
descend_hvg <- descend::findHVG(res_descend, threshold = 50)$HVG.genes

idx <- which(colnames(dat) %in% c(vst_hvg, descend_hvg))
dat <- dat[,idx]
dat_count <- dat_count[,idx]

rm(list = c("idx", "zz", "obj", "marques"))
print(paste0(Sys.time(), ": Finished screening"))
save.image(paste0("../results/step0_screening", suffix, ".RData"))
