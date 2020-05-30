rm(list=ls())
load("../results/Writeup_revision11_continuum.RData")

color_func <- function(alpha = 0.2){
  c(rgb(240/255, 228/255, 66/255, alpha), #yellow
    rgb(86/255, 180/255, 233/255, alpha), #skyblue
    rgb(0/255, 158/255, 115/255, alpha), #bluish green
    rgb(0/255, 114/255, 178/255,alpha), #blue
    rgb(230/255, 159/255, 0/255,alpha), #orange
    rgb(100/255, 100/255, 100/255, alpha)) #gray
}

color_name_vec <- c("yellow", "skyblue", "bluish green", "blue", "orange", "gray")

num_order_vec_esvd <- c(5, rep(3,2), c(6,1,1,4,4,4), rep(2,2),  rep(5,2))
col_vec_esvd <- color_func(1)[num_order_vec_esvd]
col_vec2_esvd <- color_func(0.5)[num_order_vec_esvd]
col_name_esvd <- color_name_vec[num_order_vec_esvd]
order_vec_esvd <- c(3, 5.1, 5.2, seq(6.1, 6.6, by = 0.1), 4.1, 4.2, 2, 1)
col_info_esvd <- data.frame(name = levels(cell_type_vec),
                            idx = sort(unique(cluster_labels)),
                            order = order_vec_esvd,
                            col_name = col_name_esvd,
                            col_code = col_vec_esvd)
col_info_esvd$factor_idx <- as.numeric(as.factor(col_info_esvd$col_name))
col_info_esvd[,c(5,6)] <- col_info_esvd[,c(6,5)]
colnames(col_info_esvd)[c(5,6)] <- colnames(col_info_esvd)[c(6,5)]
col_info_esvd
col_vec_short <- color_func(0.9)[c(1,4)]
plotting_order_esvd <- list(3,2,5,4,c(6,1))

###########################

res_mat <- .extract_information(segmentation_res)
zz <- order_highly_expressed_genes(res_mat, nrow(dat1), nrow(dat2), common_n = length(cell_idx_common),
                                   threshold = 1.5)

col_idx <- c(zz$common_genes, zz$traj1_genes, zz$traj2_genes)
zz_mat <- sapply(col_idx, function(j){
  tmp <- segmentation_res[[j]]$vec2_smooth
  c(segmentation_res[[j]]$vec1_smooth, tmp[(length(cell_idx_common)+1):length(tmp)])
})

stopifnot(nrow(zz_mat) == length(cell_idx_common)+length(cell_idx_traj1)+length(cell_idx_traj2))

# rescale the matrix based on the largest value in the region found
colnames(zz_mat) <- col_idx
for(i in 1:ncol(zz_mat)){
  idx <- col_idx[i]
  vec <- res_mat[idx,]
  bool <- vec$obj_1 > vec$obj_2
  if(bool){
    start <- vec$start_1; end <- vec$end_1
  } else {
    start <- vec$start_2; end <- vec$end_2
    if(start > length(cell_idx_common)) start <- start + (nrow(dat1) - length(cell_idx_common))
    if(end > length(cell_idx_common)) end <- end + (nrow(dat1) - length(cell_idx_common))
  }
  max_val <- max(zz_mat[start:end, i])
  min_val <- min(zz_mat[,i])
  zz_mat[,i] <- pmin((zz_mat[,i]-min_val)/(max_val - min_val), 1)^1.5
}

# delete some cells
cell_vec <- cluster_labels[c(cell_idx_common)]
idx <- which(cell_vec[1:1200] %in% c(2,3,4))
zz_mat <- zz_mat[-idx,]
cell_vec <- cell_vec[-idx]
common_length <- length(cell_idx_common) - length(idx)
n1 <- nrow(dat1) - length(idx)

zz_mat <- t(zz_mat)

###########3

clockwise90 <- function(a) { t(a[nrow(a):1,]) }
grDevices::png(filename = paste0("../../esvd_results/figure/experiment/Writeup_revision12_continuum.png"),
               height = 2500, width = 1500, res = 300,
               units = "px")
layout_matrix <- matrix(c(3,1,4,2), nrow = 2, ncol = 2)
graphics::layout(layout_matrix, heights = c(5, nrow(zz_mat)),
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
lines(c(0,1), rep(1-length(zz$common_genes)/nrow(zz_mat), 2), lwd = 4, col = "white")
lines(c(0,1), rep(1-(length(zz$common_genes)+length(zz$traj1_genes))/nrow(zz_mat), 2), lwd = 4, col = "white")
lines(c(0,1), rep(1-length(zz$common_genes)/nrow(zz_mat), 2), lwd = 2, lty = 2)
lines(c(0,1), rep(1-(length(zz$common_genes)+length(zz$traj1_genes))/nrow(zz_mat), 2), lwd = 2, lty = 2)

########

par(mar = c(4,2,0,4))
graphics::image(clockwise90(zz_mat[,(n1+1):ncol(zz_mat)]), xaxt='n', yaxt='n',
                col = grDevices::hcl.colors(12, "ag_GrnYl"))

seq_3 <- seq(common_length/n1, 1, length.out = num_ticks)
axis(1, at = seq(0, 1, length.out = num_ticks), labels = round(seq_3, 2))

spacing <- 1/(length(col_idx)+1)
y_vec <- seq(spacing, 1-spacing, length.out = length(col_idx))
axis(4, at = y_vec, labels = colnames(dat_impute)[col_idx], cex = 0.5)

# add dashed lines (dividing genes)
lines(c(0,1), rep(1-length(zz$common_genes)/nrow(zz_mat), 2), lwd = 4, col = "white")
lines(c(0,1), rep(1-(length(zz$common_genes)+length(zz$traj1_genes))/nrow(zz_mat), 2), lwd = 4, col = "white")
lines(c(0,1), rep(1-length(zz$common_genes)/nrow(zz_mat), 2), lwd = 2, lty = 2)
lines(c(0,1), rep(1-(length(zz$common_genes)+length(zz$traj1_genes))/nrow(zz_mat), 2), lwd = 2, lty = 2)


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

par(mar = c(1,2,1,4))
image(clockwise90(cell_vis), breaks = seq(0.5,13.5,by=1), col = col_info_esvd$col_code, axes = F)

graphics.off()
###########


