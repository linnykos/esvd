library(descend)
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
dat <- t(apply(dat, 1, function(x){10^4 * x/sum(x)}))

obj <- Seurat::CreateSeuratObject(counts = t(dat), project = "marques",
                                  meta.data = NULL, min.cells = 0, min.features = 0)

obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = 1000)
hvg <- Seurat::VariableFeatures(object = obj)

idx <- which(colnames(dat) %in% hvg)
dat <- dat[,hvg]
dat_count <- dat_count[,hvg]

rm(list = c("idx", "zz", "k", "lvls", "obj", "marques"))
print(paste0(Sys.time(), ": Finished screening"))
save.image(paste0("../results/step0_screening", suffix, ".RData"))
