rm(list=ls())
set.seed(10)
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

res <- VIPER::VIPER(as.data.frame(t(dat)), num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 1,
             report = FALSE, outdir = NULL, prefix = NULL)
save.image("../experiment/Week33_marques_viper.RData")
