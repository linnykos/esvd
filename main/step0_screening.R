library(descend)
load("../../raw_data/marques.RData")

dat <- marques$counts
colnames(dat) <- gsub("_", "-", colnames(dat))
dat_count <- dat

cell_types <- unique(marques$cell.info$cell.type)
dim(dat)

# # remove genes with too many 0's
zz <- apply(dat, 2, function(x){length(which(x!=0))})
dat <- dat[,which(zz > 30)]
dat_count <- dat_count[,which(zz > 30)]

dat_count <- dat
dat <- t(apply(dat, 1, function(x){10^4 * x/sum(x)}))

obj <- Seurat::CreateSeuratObject(counts = t(dat), project = "marques",
                                  meta.data = NULL, min.cells = 0, min.features = 0)

obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
hvg <- Seurat::VariableFeatures(object = obj)

idx <- which(colnames(dat) %in% hvg)
dat <- dat[,hvg]
dat_count <- dat_count[,hvg]

rm(list = c("idx", "zz", "k", "lvls", "obj"))
print(paste0(Sys.time(), ": Finished screening"))
save.image(paste0("../results/step0_screening", suffix, ".RData"))
