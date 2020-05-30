load("../results/step5_trajectory_original.RData")
gene_vec <- colnames(dat_impute)
gene_vec2 <- sapply(gene_vec, function(gene){
  bool <- any(c(gene %in% aliasSymbol$alias_symbol, gene %in% aliasSymbol$symbol))
  if(!bool | gene %in% aliasSymbol$symbol) return(gene)

  idx <- which(aliasSymbol$alias_symbol %in% gene)[1]
  aliasSymbol$symbol[idx]
})
names(gene_vec2) <- NULL


load("../results/Writeup_revision11_continuum.RData")

unique(new_mat[,1])
chosen_idx <- which(gene_vec2 %in% new_mat[which(new_mat[,1] == "MO"),2])
length(chosen_idx)

for(j in chosen_idx){
  par(mfrow = c(1,2))
  max_val <- max(c(segmentation_res[[j]]$vec1_smooth), c(segmentation_res[[j]]$vec2_smooth))
  plot(segmentation_res[[j]]$vec1_smooth, ylim = c(0, max_val), main = j)
  lines(rep(length(cell_idx_common), 2), c(-1e4, 1e4), col = "red", lty = 2)

  plot(segmentation_res[[j]]$vec2_smooth, ylim = c(0, max_val), main = j)
  lines(rep(length(cell_idx_common), 2), c(-1e4, 1e4), col = "red", lty = 2)
}

# common (OPC): 919, 809, 725, 652, 177, 173, 112, 97,
# common (NFO): 938, 866, 862, 828, 773, 761, 743, 600, 467, 339, 311, 304, 280, 68, 62, 33, 31,
# common (MO): 471 (right on boundary), 447, 296?
# traj1 (MO): 971, 940, 887, 570?, 492
# traj2 (MO): 916, 821, 791, 557?, 383, 327, 232, 184, 178, 142, 25

# deleted: 440, 808, 521, 460, 406, 360, 233
