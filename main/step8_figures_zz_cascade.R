var <- ls()
rm(list = var[var != "suffix"])
load(paste0("../results/step8_figures", suffix, ".RData"))

# add in the zhang marker genes
gene_vec <- colnames(dat_impute)
gene_vec <- eSVD::convert_synonyms(gene_vec)
zhang_common_genes <- c('Rlbp1','Col1a2','Enpp6','Rnf122','Rras2','Itpr2','Mical3','Chn2','Elovl6','Cemip2','Shisal1','Rap2a','9630013A20Rik','Bmp4','Fyn','Cnksr3','Nfasc','Tmem163','Opalin','Mbp','Pdlim2')
zhang_genes1 <- c('Trf','Tppp3','Gsn')
zhang_genes2 <- c('Ppp1r14a','Ndrg1','Inf2','Itgb4')

manual_add_common <- which(gene_vec %in% zhang_common_genes)
manual_add_traj1 <- which(gene_vec %in% zhang_genes1)
manual_add_traj2 <- which(gene_vec %in% zhang_genes2)

# determine genes and their ordering
zz <- eSVD::order_highly_expressed_genes(segmentation_res,
                                         nrow1 = nrow(segmentation_prep$dat1),
                                         nrow2= nrow(segmentation_prep$dat2),
                                         common_n = length(segmentation_prep$cell_idx_common),
                                         threshold = 1.5, manual_add_common = manual_add_common,
                                         manual_add_traj1 = manual_add_traj1,
                                         manual_add_traj2 = manual_add_traj2)
# form the matrix
col_idx <- c(zz$common_genes, zz$traj1_genes, zz$traj2_genes)
zz_mat <- sapply(col_idx, function(j){
  tmp <- segmentation_res[[j]]$vec2_smooth
  c(segmentation_res[[j]]$vec1_smooth, tmp[(length(cell_idx_common)+1):length(tmp)])
})

