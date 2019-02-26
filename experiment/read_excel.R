library(readxl)
gene_mat <- read_excel("../../raw_data/Marques_genes.xlsx", range = "B3:Y53")
gene_mat <- as.matrix(gene_mat)
gene_vec <- sort(unique(as.vector(gene_mat)))

load("../../raw_data/marques.RData")
gene_all <- colnames(marques$counts)
length(which(gene_vec %in% gene_all))
idx <- which(gene_all %in% gene_vec)
