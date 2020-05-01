rm(list=ls())
load("../results/step4_factorization_cg_vst-spca.RData")
paramMat_esvd

zz1 <- esvd_embedding$u_mat

cluster_labels <- as.numeric(cell_type_vec)
cluster_center1 <- .compute_cluster_center(zz1, .construct_cluster_matrix(cluster_labels))

###########################

idx_choice <- 1
nat_mat_list <- lapply(1:length(esvd_missing_list[[idx_choice]]), function(i){
  esvd_missing_list[[idx_choice]][[i]]$u_mat %*% t(esvd_missing_list[[idx_choice]][[i]]$v_mat)
})

training_idx_list <- lapply(1:length(missing_idx_list), function(i){
  c(1:prod(dim(dat_impute)))[-missing_idx_list[[i]]]
})

par(mfrow = c(1,2))
plot_prediction_against_observed(dat_impute, nat_mat_list = nat_mat_list,
                                 missing_idx_list = training_idx_list,
                                 family = fitting_distr,
                                 scalar = paramMat_esvd[idx_choice, "scalar"],
                                 main = "eSVD embedding:\nMatrix-completion diagnostic\n(Training set)",
                                 max_points = 1e6)


plot_prediction_against_observed(dat_impute, nat_mat_list = nat_mat_list,
                                 missing_idx_list = missing_idx_list,
                                 family = fitting_distr,
                                 scalar = paramMat_esvd[idx_choice, "scalar"],
                                 main = "eSVD embedding:\nMatrix-completion diagnostic\n(Testing set)")

###################3
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
col_name_svd <- c("orange", rep("bluish green", 2), rep("yellow", 6), rep("skyblue", 2), rep("orange", 2))
order_vec_svd <- c(3, 5.1, 5.2, seq(6.1, 6.6, by = 0.1), 4.1, 4.2, 2, 1)
col_info_svd <- data.frame(name = levels(cell_type_vec),
                           idx = sort(unique(cluster_labels)),
                           order = order_vec_svd,
                           col_name = col_name_svd,
                           col_code = col_vec_svd)
col_info_svd$factor_idx <- as.numeric(as.factor(col_info_svd$col_name))
col_info_svd[,c(5,6)] <- col_info_svd[,c(6,5)]
colnames(col_info_svd)[c(5,6)] <- colnames(col_info_svd)[c(6,5)]
col_info_svd
plotting_order_svd <- c(2,3,1,4)

combn_mat <- combn(3,2)

#########################

# plotting the first 3 dimensions
par(mfrow = c(1,3))

for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]

  plot(x = zz1[,i], y = zz1[,j],
       asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
       main = "eSVD embedding and trajectories\n(Curved Gaussian)",
       pch = 16, col = col_info_svd$col_code[cluster_labels])

  for(ll in 1:nrow(cluster_center1)){
    points(cluster_center1[ll,i], cluster_center1[ll,j], pch = 16, cex = 2.25, col = "black")
    points(cluster_center1[ll,i], cluster_center1[ll,j], pch = 16, cex = 1.5, col = col_vec_svd[ll])
  }
}

## # plot each of the 6 mature oligos
for(k in 1:ncol(combn_mat)){
  par(mfrow = c(2,3))
  i <- combn_mat[1,k]; j <- combn_mat[2,k]

  for(ii in 1:6){
    ll <- which(col_info_svd$factor_idx == 4)[ii]
    point_idx <- which(cluster_labels == ll)

    plot(x = NA, y = NA, xlim = range(zz1[,i]), ylim = range(zz1[,j]),
         asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
         main = "eSVD embedding and trajectories\n(Curved Gaussian)")

    lines(x = rep(0, 2), y = c(-1e5,1e5), col = "red", lwd = 2, lty = 2)
    lines(x = c(-1e5,1e5), y = rep(0, 2), col = "red", lwd = 2, lty = 2)
    points(x = zz1[point_idx,i], y = zz1[point_idx,j], pch = 16, col = col_info_svd$col_code[ll])

    for(ll in 1:nrow(cluster_center1)){
      points(cluster_center1[ll,i], cluster_center1[ll,j], pch = 16, cex = 2.25, col = "black")
      points(cluster_center1[ll,i], cluster_center1[ll,j], pch = 16, cex = 1.5, col = col_vec_svd[ll])
    }
  }
}

###############

cluster_labels <- as.numeric(cell_type_vec)
order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
cluster_group_list <- lapply(order_vec, function(x){
  grep(paste0("^", x), levels(cell_type_vec))
})

upscale_factor <- 1

p <- 5
set.seed(10)
esvd_curves <- eSVD::slingshot(zz1[,1:p], cluster_labels, starting_cluster = cluster_group_list[[1]][1],
                               cluster_group_list = cluster_group_list,
                               verbose = T, upscale_factor = upscale_factor, reduction_percentage = reduction_percentage,
                               squared = T)

par(mfrow = c(1,3))

for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]

  plot(x = zz1[,i], y = zz1[,j],
       asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
       main = "eSVD embedding and trajectories\n(Curved Gaussian)",
       pch = 16, col = col_info_svd$col_code[cluster_labels])

  for(ll in 1:nrow(cluster_center1)){
    points(cluster_center1[ll,i], cluster_center1[ll,j], pch = 16, cex = 2.25, col = "black")
    points(cluster_center1[ll,i], cluster_center1[ll,j], pch = 16, cex = 1.5, col = col_vec_svd[ll])
  }

  curves <- esvd_curves$curves
  for(ll in 1:length(curves)) {
    ord <- curves[[ll]]$ord
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "white", lwd = 8)
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "black", lwd = 5)
  }
}



###############

zz_pca <- stats::prcomp(zz1, scale. = F, center = T)
plot(zz_pca$sdev)
cluster_labels <- as.numeric(cell_type_vec)
order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
cluster_group_list <- lapply(order_vec, function(x){
  grep(paste0("^", x), levels(cell_type_vec))
})
for(i in 3:15){
  p <- i
  dat <- zz1[,1:p]
  starting_cluster = cluster_group_list[[1]][1]
  verbose = T
  squared = T
  stopifnot(!is.list(cluster_group_list) || starting_cluster %in% cluster_group_list[[1]])
  stopifnot(all(cluster_labels > 0), all(cluster_labels %% 1 == 0), length(unique(cluster_labels)) == max(cluster_labels))
  if(all(!is.na(cluster_group_list))){
    tmp <- unlist(cluster_group_list)
    stopifnot(length(tmp) == length(unique(tmp)), length(tmp) == length(unique(cluster_labels)))
  }

  ### construct the distance matrix
  dist_mat <- .compute_cluster_distances(dat, cluster_labels)
  if(squared) dist_mat <- dist_mat^2

  if(all(is.na(cluster_group_list))){
    ### construct the spt (shortest path tree)
    g <- .construct_graph(dist_mat)

    ### identify lineages (paths through trees)
    lineages <- .construct_lineages(g, starting_cluster = starting_cluster)

  } else {
    lineages <- .construct_lineage_from_hierarchy(dist_mat, cluster_group_list,
                                                  starting_cluster)
  }
  print(paste0(i, " : ", length(lineages)))
}
