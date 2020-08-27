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

# extract relevant information from segmentation_prep for convenience
cell_idx_common <- segmentation_prep$cell_idx_common
cell_idx_traj1 <- segmentation_prep$cell_idx_traj1
cell_idx_traj2 <- segmentation_prep$cell_idx_traj2
common_length <- length(cell_idx_common)
cell_vec <- cluster_labels[c(cell_idx_common)]
dat1 <- segmentation_prep$dat1
dat2 <- segmentation_prep$dat2
n1 <- nrow(dat1); n2 <- nrow(dat2)

# determine genes and their ordering
zz <- eSVD::order_highly_expressed_genes(segmentation_res$df,
                                         nrow1 = n1, nrow2= n2,
                                         common_n = length(cell_idx_common),
                                         threshold = 1.5, manual_add_common = manual_add_common,
                                         manual_add_traj1 = manual_add_traj1,
                                         manual_add_traj2 = manual_add_traj2)

# form the matrix of smoothed signal based on segmentation_res$segmentation_fit
col_idx <- c(zz$common_genes, zz$traj1_genes, zz$traj2_genes)
zz_mat <- sapply(col_idx, function(j){
  tmp <- segmentation_res$segmentation_fit[[j]]$vec2_smooth
  c(segmentation_res$segmentation_fit[[j]]$vec1_smooth,
    tmp[(length(cell_idx_common)+1):length(tmp)])
})

# mark which genes are marker genes
# load("data/zhang_genes.rda")
# marker_genes <- unlist(lapply(c("OPCs", "NFO", "MO"), function(x){zhang_genes[which(zhang_genes$cell_type == x), "enriched_genes"]}))
marker_genes <- c(zhang_common_genes, zhang_genes1, zhang_genes2)
col_idx_idx <- which(gene_vec[col_idx] %in% marker_genes)
col_idx_marker <- col_idx[col_idx_idx]

# rescale the matrix based on the largest value in the region found
colnames(zz_mat) <- col_idx
for(i in 1:ncol(zz_mat)){
  idx <- col_idx[i]
  vec <- segmentation_res$df[idx,]
  bool <- vec$obj_1 > vec$obj_2
  if(bool){
    start <- vec$start_1; end <- vec$end_1
  } else {
    start <- vec$start_2; end <- vec$end_2
    if(start > length(cell_idx_common)) start <- start + (n1 - length(cell_idx_common))
    if(end > length(cell_idx_common)) end <- end + (n1 - length(cell_idx_common))
  }
  max_val <- max(zz_mat[start:end, i])
  min_val <- min(zz_mat[,i])
  zz_mat[,i] <- pmin((zz_mat[,i]-min_val)/(max_val - min_val), 1)^1.5
}

zz_mat <- t(zz_mat)

#########################################


clockwise90 <- function(a) { t(a[nrow(a):1,]) }
grDevices::png(filename = paste0("../../esvd_results/figure/main/gene_continuum.png"),
               height = 2500, width = 1500, res = 300,
               units = "px")
layout_matrix <- matrix(c(3,1,4,2), nrow = 2, ncol = 2)
graphics::layout(layout_matrix, heights = c(7, nrow(zz_mat)),
                 widths = c(length(c(zz$common_genes, zz$traj1_genes)), 3*length(zz$traj2_genes)))

par(mar = c(4,4,0,1))
graphics::image(clockwise90(zz_mat[,1:n1]), xaxt='n', yaxt='n',
                col = grDevices::hcl.colors(12, "ag_GrnYl"),
                xlab = "Pseudotime (Cells, rescaled)", ylab = "Genes")

num_ticks <- 4
seq_1 <- seq(0, common_length/n1, length.out = num_ticks)[-num_ticks]
seq_2 <- seq(common_length/n1, 1, length.out = num_ticks)
axis(1, at = c(seq_1, seq_2), labels = round(c(seq_1, seq_2), 2))

lines(rep(common_length/n1, 2), c(0,1), lwd = 4, col = "white")
lines(rep(common_length/n1, 2), c(0,1), lwd = 2, lty = 2)

# add dashed lines (dividing genes)
spacing <- 1/(length(col_idx)-1)
y_vec <- seq(0, 1, by = spacing)
lines(c(0,1), rep(1-y_vec[length(zz$common_genes)]-spacing/2, 2), lwd = 4, col = "white")
lines(c(0,1), rep(1-y_vec[length(zz$common_genes)+length(zz$traj1_genes)]-spacing/2, 2), lwd = 4, col = "white")
lines(c(0,1), rep(1-y_vec[length(zz$common_genes)]-spacing/2, 2), lwd = 2, lty = 2)
lines(c(0,1), rep(1-y_vec[length(zz$common_genes)+length(zz$traj1_genes)]-spacing/2, 2), lwd = 2, lty = 2)

# add the gene names
spacing <- 1/(length(col_idx)-1)
y_vec <- rev(seq(0, 1, by = spacing))
stopifnot(length(y_vec) == length(col_idx))
axis(2, at = y_vec[col_idx_idx], labels = rep(" ", length(col_idx_idx)), cex.axis = 0.5, las = 1)
axis(4, at = y_vec[col_idx_idx], labels = rep(" ", length(col_idx_idx)), cex.axis = 0.5, las = 1)

########

par(mar = c(4,0.5,0,4))
graphics::image(clockwise90(zz_mat[,(n1+1):ncol(zz_mat)]), xaxt='n', yaxt='n',
                col = grDevices::hcl.colors(12, "ag_GrnYl"))

seq_3 <- seq(common_length/n1, 1, length.out = num_ticks)
axis(1, at = seq(0, 1, length.out = num_ticks), labels = round(seq_3, 2))

# add the gene names
spacing <- 1/(length(col_idx)-1)
y_vec <- rev(seq(0, 1, by = spacing))
stopifnot(length(y_vec) == length(col_idx))
axis(4, at = y_vec[col_idx_idx], labels = gene_vec[col_idx_marker], cex.axis = 0.35, las = 1)

# add dashed lines (dividing genes)
spacing <- 1/(length(col_idx)-1)
y_vec <- seq(0, 1, by = spacing)
lines(c(0,1), rep(1-y_vec[length(zz$common_genes)]-spacing/2, 2), lwd = 4, col = "white")
lines(c(0,1), rep(1-y_vec[length(zz$common_genes)+length(zz$traj1_genes)]-spacing/2, 2), lwd = 4, col = "white")
lines(c(0,1), rep(1-y_vec[length(zz$common_genes)]-spacing/2, 2), lwd = 2, lty = 2)
lines(c(0,1), rep(1-y_vec[length(zz$common_genes)+length(zz$traj1_genes)]-spacing/2, 2), lwd = 2, lty = 2)

#######

cell_vec_1 <- c(cell_vec, cluster_labels[c(cell_idx_traj1)])
cell_vis <- rbind(cell_vec_1, cell_vec_1)

par(mar = c(1,4,1,1))
image(clockwise90(cell_vis), breaks = seq(0.5,13.5,by=1), col = col_info_esvd$col_code, axes = F)

lines(rep(common_length/n1, 2), c(0,1), lwd = 4, col = "white")
lines(rep(common_length/n1, 2), c(0,1), lwd = 2, lty = 2)

#######

cell_vec_2 <- cluster_labels[c(cell_idx_traj2)]
cell_vis <- rbind(cell_vec_2, cell_vec_2)

par(mar = c(1,0.5,1,4))
image(clockwise90(cell_vis), breaks = seq(0.5,13.5,by=1), col = col_info_esvd$col_code, axes = F)

graphics.off()



