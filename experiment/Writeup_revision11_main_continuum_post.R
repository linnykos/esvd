rm(list=ls())
load("../results/Writeup_revision11_continuum.RData")

nrow_vec <- c(nrow(dat1), nrow(dat2))
start_vec_list <- vector("list", 2)
end_vec_list <- vector("list", 2)
midpoint_vec_list <- vector("list", 2)
obj_vec_list <- vector("list", 2)

par(mfrow = c(1,2))
for(k in 1:2){
  start_vec_list[[k]] <- sapply(1:length(segmentation_res), function(i){segmentation_res[[i]][[k]]$i})
  end_vec_list[[k]] <- sapply(1:length(segmentation_res), function(i){segmentation_res[[i]][[k]]$j})
  midpoint_vec_list[[k]] <- sapply(1:length(start_vec_list[[k]]),
                                   function(i){(start_vec_list[[k]][i] + end_vec_list[[k]][i])/2})
  obj_vec_list[[k]] <- sapply(1:length(segmentation_res), function(i){segmentation_res[[i]][[k]]$obj_val})

  plot(NA, ylim = range(obj_vec_list[[k]]), xlim = c(0, nrow_vec[k]),
       main = paste0(nrow_vec[k]))
  lines(rep(length(cell_idx_common), 2), c(-1e5, 1e5), col = "red", lwd = 2)
  for(i in 1:length(start_vec_list[[k]])){
    lines(x = c(start_vec_list[[k]][i], end_vec_list[[k]][i]), y = rep(obj_vec_list[[k]][i], 2), lwd = 2)
  }
}

##############################

res_mat <- .extract_information(segmentation_res)
zz <- order_highly_expressed_genes(res_mat, nrow(dat1), nrow(dat2), common_n = length(cell_idx_common),
                                   threshold = 1.5)

for(j in zz$common_genes){
  par(mfrow = c(1,2))
  vec1 <- dat1[,j]
  vec2 <- dat2[,j]
  col_vec <- rep("black", length(vec2)); col_vec[segmentation_res[[j]]$cut_2$i:segmentation_res[[j]]$cut_2$j] <- "red"
  plot(vec1, pch = 16, cex = 0.5, main = j)
  lines(rep(length(cell_idx_common), 2), c(-1e5, 1e5), col = "red", lwd = 2)
  points(segmentation_res[[j]]$vec1_smooth, col = "white", cex = 0.75, pch = 16)
  points(segmentation_res[[j]]$vec1_smooth, col = "red", cex = 0.5, pch = 16)

  plot(vec2, col = col_vec, pch = 16, cex = 0.5, main = j)
  lines(rep(length(cell_idx_common), 2), c(-1e5, 1e5), col = "red", lwd = 2)
  points(segmentation_res[[j]]$vec2_smooth, col = "white", cex = 0.75, pch = 16)
  points(segmentation_res[[j]]$vec2_smooth, col = "red", cex = 0.5, pch = 16)
}

for(j in zz$traj1_genes){
  par(mfrow = c(1,2))
  vec1 <- dat1[,j]
  col_vec <- rep("black", length(vec1)); col_vec[segmentation_res[[j]]$cut_1$i:segmentation_res[[j]]$cut_1$j] <- "red"
  vec2 <- dat2[,j]
  plot(vec1, col = col_vec, pch = 16, cex = 0.5, main = j)
  lines(rep(length(cell_idx_common), 2), c(-1e5, 1e5), col = "red", lwd = 2)
  points(segmentation_res[[j]]$vec1_smooth, col = "white", cex = 0.75, pch = 16)
  points(segmentation_res[[j]]$vec1_smooth, col = "red", cex = 0.5, pch = 16)
  plot(vec2, pch = 16, cex = 0.5, main = j)
  lines(rep(length(cell_idx_common), 2), c(-1e5, 1e5), col = "red", lwd = 2)
  points(segmentation_res[[j]]$vec2_smooth, col = "white", cex = 0.75, pch = 16)
  points(segmentation_res[[j]]$vec2_smooth, col = "red", cex = 0.5, pch = 16)
}

for(j in zz$traj2_genes){
  par(mfrow = c(1,2))
  vec1 <- dat1[,j]
  vec2 <- dat2[,j]
  col_vec <- rep("black", length(vec2)); col_vec[segmentation_res[[j]]$cut_2$i:segmentation_res[[j]]$cut_2$j] <- "red"
  plot(vec1, pch = 16, cex = 0.5, main = j)
  lines(rep(length(cell_idx_common), 2), c(-1e5, 1e5), col = "red", lwd = 2)
  points(segmentation_res[[j]]$vec1_smooth, col = "white", cex = 0.75, pch = 16)
  points(segmentation_res[[j]]$vec1_smooth, col = "red", cex = 0.5, pch = 16)

  plot(vec2, col = col_vec, pch = 16, cex = 0.5, main = j)
  lines(rep(length(cell_idx_common), 2), c(-1e5, 1e5), col = "red", lwd = 2)
  points(segmentation_res[[j]]$vec2_smooth, col = "white", cex = 0.75, pch = 16)
  points(segmentation_res[[j]]$vec2_smooth, col = "red", cex = 0.5, pch = 16)
}

# let's visualize the data in the simpliest way possible for now
length(cell_idx_common)
length(cell_idx_traj1)
length(cell_idx_traj2)
length(zz$common_genes)
length(zz$traj1_genes)
length(zz$traj2_genes)

# let's do the simplest thing possible first
## compile the matrix
zz_mat <- sapply(c(zz$common_genes, zz$traj1_genes, zz$traj2_genes), function(j){
  tmp <- segmentation_res[[j]]$vec2_smooth
  c(segmentation_res[[j]]$vec1_smooth, tmp[(length(cell_idx_common)+1):length(tmp)])
})
stopifnot(nrow(zz_mat) == length(cell_idx_common)+length(cell_idx_traj1)+length(cell_idx_traj2))
## rescale the matrix
zz_mat <- apply(zz_mat, 2, function(x){(x-min(x))/(max(x)-min(x))})

clockwise90 = function(a) { t(a[nrow(a):1,]) }
grDevices::png(filename = paste0("../../esvd_results/figure/experiment/Writeup_revision11_continuum.png"),
               height = 2500, width = 1500, res = 300,
               units = "px")
graphics::image(clockwise90(zz_mat))

# add dashed lines (dividing cells)
lines(c(0,1), rep(1-length(cell_idx_common)/nrow(zz_mat), 2), lwd = 4, col = "white")
lines(c(0,1), rep(1-(length(cell_idx_common)+length(cell_idx_traj1))/nrow(zz_mat), 2), lwd = 4, col = "white")
lines(c(0,1), rep(1-length(cell_idx_common)/nrow(zz_mat), 2), lwd = 2, lty = 2)
lines(c(0,1), rep(1-(length(cell_idx_common)+length(cell_idx_traj1))/nrow(zz_mat), 2), lwd = 2, lty = 2)

# add dashed lines (dividing genes)
lines(rep(length(zz$common_genes)/ncol(zz_mat), 2), c(0,1), lwd = 4, col = "white")
lines(rep((length(zz$common_genes)+length(zz$traj1_genes))/ncol(zz_mat), 2), c(0,1), lwd = 4, col = "white")
lines(rep(length(zz$common_genes)/ncol(zz_mat), 2), c(0,1), lwd = 2, lty = 2)
lines(rep((length(zz$common_genes)+length(zz$traj1_genes))/ncol(zz_mat), 2), c(0,1), lwd = 2, lty = 2)
graphics.off()

########################
#
# zz <- intersect(intersect(which(obj_vec_list[[1]]>1.5), which(start_vec_list[[1]] >= length(cell_idx_common))),
#                 which(end_vec_list[[1]] <= 3000))
# zz_obj <- sapply(zz, function(x){obj_vec_list[[1]][x]})
# zz <- zz[order(zz_obj, decreasing = T)]
# head(zz); length(zz)
#
# for(k in 1:min(10, length(zz))){
#   j <- zz[k]
#   par(mfrow = c(1,2))
#   vec1 <- dat1[,j]
#   col_vec <- rep("black", length(vec1)); col_vec[segmentation_res[[j]]$cut_1$i:segmentation_res[[j]]$cut_1$j] <- "red"
#   vec2 <- dat2[,j]
#   plot(vec1, col = col_vec, pch = 16, cex = 0.5, main = j)
#   lines(rep(length(cell_idx_common), 2), c(-1e5, 1e5), col = "red", lwd = 2)
#   points(segmentation_res[[j]]$vec1_smooth, col = "white", cex = 0.75, pch = 16)
#   points(segmentation_res[[j]]$vec1_smooth, col = "red", cex = 0.5, pch = 16)
#   plot(vec2, pch = 16, cex = 0.5, main = j)
#   lines(rep(length(cell_idx_common), 2), c(-1e5, 1e5), col = "red", lwd = 2)
#   points(segmentation_res[[j]]$vec2_smooth, col = "white", cex = 0.75, pch = 16)
#   points(segmentation_res[[j]]$vec2_smooth, col = "red", cex = 0.5, pch = 16)
# }
#
# #######
#
# zz <- intersect(intersect(which(obj_vec_list[[2]]>1.5), which(start_vec_list[[2]] >= length(cell_idx_common))),
#                 which(end_vec_list[[2]] <= 3250))
# zz_obj <- sapply(zz, function(x){obj_vec_list[[2]][x]})
# zz <- zz[order(zz_obj, decreasing = T)]
# head(zz); length(zz)
#
# for(k in 1:min(10, length(zz))){
#   j <- zz[k]
#   par(mfrow = c(1,2))
#   vec1 <- dat1[,j]
#   vec2 <- dat2[,j]
#   col_vec <- rep("black", length(vec2)); col_vec[segmentation_res[[j]]$cut_2$i:segmentation_res[[j]]$cut_2$j] <- "red"
#   plot(vec1, pch = 16, cex = 0.5)
#   lines(rep(length(cell_idx_common), 2), c(-1e5, 1e5), col = "red", lwd = 2)
#   points(segmentation_res[[j]]$vec1_smooth, col = "white", cex = 0.75, pch = 16)
#   points(segmentation_res[[j]]$vec1_smooth, col = "red", cex = 0.5, pch = 16)
#
#   plot(vec2, col = col_vec, pch = 16, cex = 0.5)
#   lines(rep(length(cell_idx_common), 2), c(-1e5, 1e5), col = "red", lwd = 2)
#   points(segmentation_res[[j]]$vec2_smooth, col = "white", cex = 0.75, pch = 16)
#   points(segmentation_res[[j]]$vec2_smooth, col = "red", cex = 0.5, pch = 16)
# }

#################################

col_palatte <- colorRampPalette(c("azure","darkviolet"))(100)
max_val <- max(pseudotime_df2$pseudotime)
min_val <- min(pseudotime_df2$pseudotime)
col_vec <-  col_palatte[pmin(pmax(round((pseudotime_df2$pseudotime-min_val)/(max_val - min_val) * 100), 1), 100)]

cluster_center_esvd <- .compute_cluster_center(esvd_embedding$u_mat[,1:3], .construct_cluster_matrix(cluster_labels))


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

combn_mat <- combn(3,2)

grDevices::png(filename = paste0("../../esvd_results/figure/experiment/Writeup_revision11_continuum_pseudotime.png"),
               height = 830, width = 2300, res = 300,
               units = "px")
par(mfrow = c(1,3), mar = c(4,4,4,1))
for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]

  plot(NA, xlim = range(esvd_embedding$u_mat[,i]), ylim = range(esvd_embedding$u_mat[,j]),
       asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
       main = "eSVD embedding and trajectories\n(Curved Gaussian)")

  points(x = esvd_embedding$u_mat[pseudotime_df2$cell_idx,i],
         y = esvd_embedding$u_mat[pseudotime_df2$cell_idx,j], pch = 16,
         col = col_vec)
  curves <- esvd_curves_short$curves
  for(ll in rev(1:length(curves))) {
    ord <- curves[[ll]]$ord
    graphics::lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "white", lwd = 15)
    graphics::lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = col_vec_short[ll], lwd = 7)
    graphics::lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "black",
                    lty = 3, lwd = 3)
  }

  for(ll in 1:nrow(cluster_center_esvd)){
    points(cluster_center_esvd[ll,i], cluster_center_esvd[ll,j], pch = 16, cex = 2.25, col = "black")
    points(cluster_center_esvd[ll,i], cluster_center_esvd[ll,j], pch = 16, cex = 1.5, col = col_vec_esvd[ll])
  }
}
graphics.off()
