rm(list=ls())
load("../results/Writeup_revision11_continuum.RData")
load("data/zhang_genes.rda")
zhang_genes[,2] <- eSVD:::convert_synonyms(zhang_genes[,2])
gene_vec <- colnames(dat_impute)
gene_vec <- eSVD:::convert_synonyms(gene_vec)

length(which(zhang_genes[,2] %in% gene_vec))

unique(zhang_genes[,1])
chosen_idx <- which(gene_vec %in% zhang_genes[which(zhang_genes[,1] == "MO"),2])
length(chosen_idx)

for(j in chosen_idx){
  par(mfrow = c(1,2))
  max_val <- max(c(segmentation_res[[j]]$vec1_smooth), c(segmentation_res[[j]]$vec2_smooth))
  col_vec <- rep("black", length(segmentation_res[[j]]$vec1_smooth))
  col_vec[segmentation_res[[j]]$cut_1$i : segmentation_res[[j]]$cut_1$j] <- "red"
  plot(segmentation_res[[j]]$vec1_smooth, ylim = c(0, max_val), main = j, col = col_vec, pch = 16)
  lines(rep(length(cell_idx_common), 2), c(-1e4, 1e4), col = "red", lty = 2)

  col_vec <- rep("black", length(segmentation_res[[j]]$vec2_smooth))
  col_vec[segmentation_res[[j]]$cut_2$i : segmentation_res[[j]]$cut_2$j] <- "red"
  plot(segmentation_res[[j]]$vec2_smooth, ylim = c(0, max_val), main = gene_vec[j], col = col_vec, pch = 16)
  lines(rep(length(cell_idx_common), 2), c(-1e4, 1e4), col = "red", lty = 2)
}

#####################3

gene_idx <- c(791, 327, 232, 184)
gene_vec[gene_idx]
paste0("\'", paste0(gene_vec[gene_idx], collapse = "\',\'"),"\'")

# common (OPCs): 919, 809, 725, 652, 177, 173, 112, 97,
# common (NFO): 938, 866, 862, 828, 773, 761, 743, 600, 467, 339, 311, 304, 280, 68, 62, 33, 31,
# common (MO): 471 (right on boundary), 447, 296?
# traj1 (MO): 971, 940, 887, 570?, 492
# traj2 (MO): 916, 821, 791, 557?, 383, 327, 232, 184, 178, 142, 25

# deleted: 440, 808, 521, 460, 406, 360, 233
