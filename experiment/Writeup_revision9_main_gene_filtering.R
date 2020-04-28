rm(list=ls())
load("../results/step3_scalar_heuristic_cg_hvg_tmp.RData")
gene_name_hvg <- colnames(dat_impute)
length(gene_name_hvg)

load("../results/step5_trajectory.RData")
gene_name_org <- colnames(dat_impute)
length(gene_name_org)

length(intersect(gene_name_hvg, gene_name_org))
