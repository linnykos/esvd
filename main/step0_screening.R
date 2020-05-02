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

# obj <- Seurat::CreateSeuratObject(counts = t(dat), project = "marques",
#                                   meta.data = NULL, min.cells = 0, min.features = 0)
# obj <- Seurat::NormalizeData(obj, normalization.method = "RC", scale.factor = 10^4)
# obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = 1000)
# # Seurat::VariableFeaturePlot(obj)
# dat <- t(as.matrix(obj@assays$RNA@data))
# vst_hvg <- Seurat::VariableFeatures(object = obj)
#
# print(paste0(Sys.time(), ": Finished vst screening"))
# save.image(paste0("../results/step0_screening", suffix, ".RData"))

# try a series of SPCAs
k <- 5
lvls <- 10
v_seq <- exp(seq(log(1), log(log(ncol(dat))), length.out = lvls))
spca_list <- vector("list", lvls)
spca_func <- function(i){
  res <- PMA::SPC(dat, sumabsv = v_seq[i], K = k, trace = F)
  print(paste0("Finished SPC for level ", i))
  res
}
spca_list <- foreach::"%dopar%"(foreach::foreach(i = 1:lvls), spca_func(i))
print(paste0(Sys.time(), ": Finished sPCA"))

spca_summary <- cbind(v_seq, t(sapply(spca_list, function(x){
  c(length(unique(sort(unlist(apply(x$v, 2, function(y){which(y != 0)}))))),
    x$prop.var.explained[5])})))
idx <- min(which(spca_summary[,2] == max(spca_summary[,2])))
target_var <- spca_summary[idx,3]
idx <- min(intersect(which(spca_summary[,2] >= 500), which(spca_summary[,3] >= 0.9*target_var)))
spca_idx <- sort(unique(unlist(apply(spca_list[[idx]]$v, 2, function(x){which(x != 0)}))))
spca_hvg <- colnames(dat)[spca_idx]
print(paste0(Sys.time(), ": Finished selecting sPCA genes"))

# # run DESCEND
# res_descend <- descend::runDescend(t(dat_count), n.cores = ncores)
# descend_hvg <- descend::findHVG(res_descend, threshold = 50)$HVG.genes
# print(paste0(Sys.time(), ": Finished selecting DESCEND genes"))

obj <- Seurat::CreateSeuratObject(counts = t(dat), project = "marques",
                                  meta.data = NULL, min.cells = 0, min.features = 0)
obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = 500)
vst_hvg <- Seurat::VariableFeatures(object = obj)

idx <- which(colnames(dat) %in% c(vst_hvg, spca_hvg))
dat <- dat[,idx]
dat_count <- dat_count[,idx]

reweight_factor <- rowSums(dat)
dat <- t(sapply(1:nrow(dat), function(i){10^4 * dat[i,]/reweight_factor[i]}))

rm(list = c("idx", "zz", "obj", "marques", "k", "lvls", "spca_func", "idx", "target_var"))
print(paste0(Sys.time(), ": Finished screening"))
save.image(paste0("../results/step0_screening", suffix, ".RData"))
print(warnings())
