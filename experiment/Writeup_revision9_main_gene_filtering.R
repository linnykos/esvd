rm(list=ls())
load("../results/step3_scalar_heuristic_cg_hvg_tmp.RData")
gene_name_hvg <- colnames(dat_impute)
length(gene_name_hvg)

load("../results/step5_trajectory.RData")
gene_name_org <- colnames(dat_impute)
length(gene_name_org)

length(intersect(colnames(dat), gene_name_org))

descend_hvg <- descend::findHVG(res_descend, threshold = 20)$HVG.genes
length(descend_hvg)
length(intersect(colnames(dat), descend_hvg))
