rm(list=ls())
load("../results/step4_factorization_cg_spca-vst_before_rescaling_300_all.RData")

cluster_labels <- as.numeric(cell_type_vec)
order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
cluster_group_list <- lapply(order_vec, function(x){
  grep(paste0("^", x), levels(cell_type_vec))
})

upscale_factor <- 1

p <- 3
set.seed(10)
esvd_curves <- eSVD::slingshot(esvd_embedding$u_mat[,1:p], cluster_labels, starting_cluster = cluster_group_list[[1]][1],
                               cluster_group_list = cluster_group_list,
                               verbose = T, upscale_factor = upscale_factor,
                               reduction_percentage = 0.25,
                               squared = T)

esvd_curves$lineages

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
idx_all <- esvd_curves$idx[which(esvd_curves$curves[[1]]$W == 1)] # which ghost points are in the first lineage
idx_cell <- sort(unique(idx_all)) # which cells do these correspond to
order_vec <- rep(NA, length(idx_cell))

for(i in 1:length(idx_cell)){
  tmp_idx <- which(idx_all == idx_cell[i]) # idx of the ghost points
  order_vec[i] <- mean(which(esvd_curves$curves[[1]]$ord %in% tmp_idx))
}

plot(sort(order_vec), asp = T)


nat_mat <- esvd_embedding$u_mat %*% t(esvd_embedding$v_mat)
pred_mat <- compute_mean(nat_mat = nat_mat, family = "curved_gaussian")
zz = order(apply(dat_impute, 2, max), decreasing = T)
zz = order(apply(pred_mat, 2, function(x){quantile(x, probs = 0.75)}), decreasing = T)
zz = order(apply(pred_mat, 2, function(x){mean(x[3000:4000])}), decreasing = T)
zz = order(apply(dat_impute, 2, function(x){quantile(x, probs = 0.75)}), decreasing = T)
zz[1:20]
# kk <- 501; ylim <- c(0, 800)
kk <- 494
par(mfrow = c(1,2))
plot(pred_mat[idx_cell[order(order_vec)], kk], ylim = range(c(pred_mat[,kk], dat_impute[,kk])), col = rgb(0,0,0,0.1), pch = 16)
plot(dat_impute[idx_cell[order(order_vec)], kk], ylim = range(c(pred_mat[,kk], dat_impute[,kk])), col = rgb(0,0,0,0.1), pch = 16)


####
kk <- 1063
vec <- pred_mat[idx_cell[order(order_vec)], kk]
res <- circular_segmentation(vec, resolution = 50)
col_vec <- rep("blue", length(vec))
col_vec[res$i : res$j] <- "red"
plot(vec, col = col_vec, pch = 16)

circular_list <- lapply(1:ncol(pred_mat), function(j){
  print(j)

  vec <- pred_mat[idx_cell[order(order_vec)], j]
  circular_segmentation(vec, resolution = 50)
})

# now extract the midpoints
midpoint_vec <- sapply(1:length(circular_list), function(i){
  (circular_list[[i]]$i + circular_list[[i]]$j)/2
})
plot(sort(midpoint_vec))
obj_vec <- sapply(1:length(circular_list), function(i){circular_list[[i]]$obj_val})
plot(midpoint_vec, obj_vec, ylim = c(0,5))
