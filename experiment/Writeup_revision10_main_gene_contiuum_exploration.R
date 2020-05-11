rm(list=ls())
load("../results/step5_trajectory.RData")

# cluster_labels <- as.numeric(cell_type_vec)
# order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
# cluster_group_list <- lapply(order_vec, function(x){
#   grep(paste0("^", x), levels(cell_type_vec))
# })
#
# upscale_factor <- 1
#
# p <- 3
# set.seed(10)
# esvd_curves <- eSVD::slingshot(esvd_embedding$u_mat[,1:p], cluster_labels, starting_cluster = cluster_group_list[[1]][1],
#                                cluster_group_list = cluster_group_list,
#                                verbose = T, upscale_factor = upscale_factor,
#                                reduction_percentage = 0.25,
#                                squared = T)
#
# esvd_curves$lineages

################################################

zz <- esvd_curves$curves[[1]]$ord
quantile(zz)
nrow(dat_impute)


length(esvd_curves$curves[[1]]$lambda)
head(esvd_curves$curves[[1]]$lambda)
length(esvd_curves$curves[[1]]$ord)
head(esvd_curves$curves[[1]]$ord)
length(esvd_curves$curves[[1]]$idx)
head(esvd_curves$curves[[1]]$idx)

length(unique(esvd_curves$idx))
length(esvd_curves$idx) # idx tells which "ghost" points belong to which cells
table(esvd_curves$curves[[1]]$W) # W tells you which "ghost" points are in which curve
length(esvd_curves$curves[[1]]$W)

plot(esvd_curves$curves[[1]]$lambda[esvd_curves$curves[[1]]$ord])

# ask if points assigned to the same cell have the same order
idx_all <- esvd_curves$idx[which(esvd_curves$curves[[1]]$W == 1)]
uniq_val <- sort(unique(idx_all))
len_uniq <- length(uniq_val)
zz = table(idx_all)
table(zz)

sd_vec <- rep(NA, len_uniq)
for(i in 1:len_uniq){
  idx <- which(idx_all == uniq_val[i]) # which ghost points belong the cell i?
  if(length(idx) == 1) {
    sd_vec[i] <- 1
  } else {
    # sd_vec[i] <- sd(esvd_curves$curves[[1]]$lambda[idx])
    sd_vec[i] <- max(diff(sort(which(esvd_curves$curves[[1]]$ord %in% idx))))
  }
}
plot(sd_vec)

# verify that the points with non-trivial sd have lambda == 0
lambda_vec <- rep(NA, len_uniq)
for(i in 1:len_uniq){
  idx <- which(idx_all == uniq_val[i]) # which ghost points belong the cell i?
  lambda_vec[i] <- max(esvd_curves$curves[[1]]$lambda[idx])
}
plot(sd_vec[which(sd_vec > 1)], lambda_vec[which(sd_vec > 1)])

###############################################################
curve_idx <- 2
idx_all <- esvd_curves$idx[esvd_curves$curves[[curve_idx]]$idx] # which ghost points are in the first lineage
lambda_all <- esvd_curves$curves[[curve_idx]]$lambda
stopifnot(length(idx_all) == length(lambda_all))
idx_cell <- sort(unique(idx_all)) # which cells do these correspond to
order_vec <- rep(NA, length(idx_cell))
lambda_vec <- rep(NA, length(idx_cell))

for(i in 1:length(idx_cell)){
  tmp_idx <- which(idx_all == idx_cell[i]) # idx of the ghost points
  order_vec[i] <- mean(which(esvd_curves$curves[[1]]$ord %in% tmp_idx))
  lambda_vec[i] <- mean(lambda_all[tmp_idx])
}

# plot(sort(order_vec), asp = T)

sd_vec <- rep(NA, length(idx_cell))
for(i in 1:length(idx_cell)){
  tmp_idx <- which(idx_all == idx_cell[i])
  if(length(tmp_idx) == 1) {
    sd_vec[i] <- 0
  } else {
    sd_vec[i] <- sd(lambda_all[tmp_idx])
  }
}
quantile(sd_vec)


nat_mat <- esvd_embedding$u_mat %*% t(esvd_embedding$v_mat)
pred_mat <- compute_mean(nat_mat = nat_mat, family = "curved_gaussian")
zz = order(apply(dat_impute, 2, max), decreasing = T)
zz = order(apply(pred_mat, 2, function(x){quantile(x, probs = 0.75)}), decreasing = T)
zz = order(apply(pred_mat, 2, function(x){mean(x[3000:4000])}), decreasing = T)
zz = order(apply(dat_impute, 2, function(x){quantile(x, probs = 0.75)}), decreasing = T)
zz[1:20]
# kk <- 501; ylim <- c(0, 800)

for(kk in zz[1:5]){
  png(filename = paste0("../../esvd_results/figure/experiment/Writeup_revision10_main_gene_continuum_", kk, ".png"),
      height = 1200, width = 2400, res = 300,
      units = "px")
  par(mfrow = c(1,2))
  plot(pred_mat[idx_cell[order(order_vec)], kk], ylim = range(c(pred_mat[,kk], dat_impute[,kk])),
       col = rgb(0,0,0,0.1), pch = 16,
       xlab = "Cell ordering (via estimated trajectory)", ylab = "Expression",
       main = paste0("Predicted gene expression\n(Gene ", colnames(dat_impute)[kk], ")"))
  plot(dat_impute[idx_cell[order(order_vec)], kk], ylim = range(c(pred_mat[,kk], dat_impute[,kk])),
       col = rgb(0,0,0,0.1), pch = 16,
       xlab = "Cell ordering (via estimated trajectory)", ylab = "Expression",
       main = paste0("Observed gene expression\n(Gene ", colnames(dat_impute)[kk], ")"))
  graphics.off()
}


#######

# verify that cell type is roughly aligned with trajectory
color_func <- function(alpha = 0.2){
  c(rgb(240/255, 228/255, 66/255, alpha), #yellow
    rgb(86/255, 180/255, 233/255, alpha), #skyblue
    rgb(0/255, 158/255, 115/255, alpha), #bluish green
    rgb(0/255, 114/255, 178/255,alpha), #blue
    rgb(230/255, 159/255, 0/255,alpha), #orange
    rgb(100/255, 100/255, 100/255, alpha)) #gray
}
num_order_vec_svd <- c(5, rep(3,2), c(1,1,1,1,1,1), rep(2,2),  rep(5,2))
col_vec_svd <- color_func(1)[num_order_vec_svd]
col_vec2_svd <- color_func(0.5)[num_order_vec_svd]
order_vec_svd <- c(3, 5.1, 5.2, seq(6.1, 6.6, by = 0.1), 4.1, 4.2, 2, 1)
row_idx <- floor(order_vec_svd)[cluster_labels]

plot(NA, xlim = c(0, max(lambda_vec)), ylim = c(0,7))
for(i in 1:6){
  tmp <- lambda_vec[which(row_idx[idx_cell] == i)]
  den_res <- density(tmp)
  y_vec <- (c(0, den_res$y, 0 , 0))
  y_vec <- y_vec/1.5
  polygon(x = c(den_res$x[1], den_res$x, den_res$x[length(den_res$x)], den_res$x[1]),
          y = y_vec + i,
          col = col_vec_svd[which(floor(order_vec_svd) == i)][1])
  # points(tmp, rep(i, length(tmp)), pch = 1, col = col_vec2_svd[which(floor(order_vec_svd) == i)][1])
}

####
# kk <- 1063
# vec <- pred_mat[idx_cell[order(order_vec)], kk]
# res <- circular_segmentation(vec, resolution = 50)
# col_vec <- rep("blue", length(vec))
# col_vec[res$i : res$j] <- "red"
# plot(vec, col = col_vec, pch = 16)

circular_list <- lapply(1:ncol(pred_mat), function(j){
  print(j)

  vec <- pred_mat[idx_cell[order(order_vec)], j]
  circular_segmentation(vec, resolution = 50)
})

# now extract the midpoints
midpoint_vec <- sapply(1:length(circular_list), function(i){
  (circular_list[[i]]$i + circular_list[[i]]$j)/2
})
start_vec <- sapply(1:length(circular_list), function(i){circular_list[[i]]$i})
end_vec <- sapply(1:length(circular_list), function(i){circular_list[[i]]$j})
plot(sort(midpoint_vec))
obj_vec <- sapply(1:length(circular_list), function(i){circular_list[[i]]$obj_val})
# plot(midpoint_vec, log(pmax(obj_vec,0)+1))
plot(NA, ylim = range(log(pmax(obj_vec,0)+1)), xlim = c(0, length(idx_cell)))
for(i in 1:length(start_vec)){
  lines(x = c(start_vec[i], end_vec[i]), y = rep(log(max(obj_vec[i],0)+1), 2), lwd = 2)
}

# select highly expressed genes
gene_idx <- which(log(pmax(obj_vec,0)+1) >= 1)
# order these genes
gene_idx <- gene_idx[order(midpoint_vec[gene_idx], decreasing = F)]
length(gene_idx)
image(pred_mat[idx_cell[order(order_vec)], gene_idx])
