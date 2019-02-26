rm(list=ls())
library(descend)
load("../../raw_data/marques.RData")

dat <- marques$counts

cell_types <- unique(marques$cell.info$cell.type)
#cell_types <- cell_types[-which(cell_types %in% c("OPC", "PPR"))]

cell_idx <- unlist(lapply(cell_types, function(x){
  tmp <- which(marques$cell.info$cell.type == x)
  sample(tmp, round(length(tmp)/5))
}))
dat <- dat[cell_idx,]
dim(dat)

# use the predetermined set of genes
gene_mat <- readxl::read_excel("../../raw_data/Marques_genes.xlsx", range = "B3:Y53")
gene_mat <- as.matrix(gene_mat)
gene_vec <- sort(unique(as.vector(gene_mat)))
gene_idx <- which(colnames(dat) %in% gene_vec)
dat <- dat[,gene_idx]
dim(dat)

# # remove genes with too many 0's
# zz <- apply(dat, 2, function(x){length(which(x == 0))/length(x)})
# idx <- which(zz <= .95)
# dat <- dat[,idx]
# dim(dat)
#
# # try a series of SPCAs
# k <- 5
# lvls <- 20
# v_seq <- exp(seq(log(1), log(log(ncol(dat))), length.out = lvls))
# res_list <- vector("list", lvls)
#
# spca_func <- function(i){
#   PMA::SPC(dat, sumabsv = v_seq[i], K = k, trace = F)
# }
#
# doMC::registerDoMC(cores = 10)
# res_list <- foreach::"%dopar%"(foreach::foreach(i = 1:lvls), spca_func(i))
#
# print("Finished SPC")
#
# # run DESCEND
# res_descend <- descend::runDescend(t(dat), n.cores = 10)

rm(list = c("idx", "zz", "k", "lvls", "gene_mat", "gene_vec"))
print(paste0(Sys.time(), ": Finished screening"))
save.image("../results/step0_screening.RData")
