library(descend)
load("../../raw_data/marques.RData")

dat <- marques$counts

cell_types <- unique(marques$cell.info$cell.type)
#cell_types <- cell_types[-which(cell_types %in% c("OPC", "PPR"))]

set.seed(10)
cell_idx <- unlist(lapply(cell_types, function(x){
  tmp <- which(marques$cell.info$cell.type == x)
  sample(tmp, round(length(tmp)/5))
}))
dat <- dat[cell_idx,]
dim(dat)

# use the predetermined set of genes
load("../data/marker_genes.Rda")
gene_vec <- as.vector(gene_mat)
gene_vec <- gene_vec[!duplicated(gene_vec)]
gene_idx <- which(colnames(dat) %in% gene_vec)
dat <- dat[,gene_idx]

column_vec <- colnames(dat)[which(colnames(dat) %in% gene_vec)]
gene_vec <- gene_vec[which(gene_vec %in% colnames(dat))]
reordered_idx <- order(column_vec)[rank(gene_vec)]
dat <- dat[,reordered_idx]
dim(dat)

# # remove genes with too many 0's
zz <- apply(dat, 2, function(x){length(which(x!=0))})
dat <- dat[,which(zz > 30)]

rm(list = c("idx", "zz", "k", "lvls", "reorder_idx", "column_vec"))
print(paste0(Sys.time(), ": Finished screening"))
save.image("../results/step0_screening.RData")
