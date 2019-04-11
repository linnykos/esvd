rm(list=ls())
load("../results/step5_clustering_spca.RData")

# let's start with some basic investigation
our_sd_val$sd_val
naive_sd_val$sd_val

###################################

# basic 3d plots
# colors from http://mkweb.bcgsc.ca/colorblind/img/colorblindness.palettes.trivial.png
col_func <- function(alpha){
  c(rgb(240/255, 228/255, 66/255, alpha), #yellow
     rgb(86/255, 180/255, 233/255, alpha), #skyblue
     rgb(0/255, 158/255, 115/255, alpha), #bluish green
     rgb(0/255, 114/255, 178/255,alpha), #blue
     rgb(230/255, 159/255, 0/255,alpha)) #orange
}

reorder_cluster_labels <- function(cluster_labels, cluster_group_list,
                                   lineage = NA){
  cluster_labels2 <- rep(NA, length(cluster_labels))
  cluster_labels2[which(cluster_labels %in% unlist(cluster_group_list[1:3]))] <- 1
  cluster_labels2[which(cluster_labels %in% cluster_group_list[[4]])] <- 2
  cluster_labels2[which(cluster_labels %in% cluster_group_list[[5]])] <- 3
  if(all(is.na(lineage))){
    cluster_labels2[which(cluster_labels %in% cluster_group_list[[6]])] <- 4
  } else {
    idx <- lineage[[1]][which(!lineage[[1]] %in% lineage[[2]])]
    cluster_labels2[which(cluster_labels %in% idx)] <- 4
    idx <- lineage[[2]][which(!lineage[[2]] %in% lineage[[1]])]
    cluster_labels2[which(cluster_labels %in% idx)] <- 5
  }
  cluster_labels2[is.na(cluster_labels2)] <- 3

  center_labels <- rep(NA, max(cluster_labels))
  center_labels[unlist(cluster_group_list[1:3])] <- 1
  center_labels[cluster_group_list[[4]]] <- 2
  center_labels[cluster_group_list[[5]]] <- 3
  if(all(is.na(lineage))){
    center_labels[cluster_group_list[[6]]] <- 4
  } else {
    idx <- lineage[[1]][which(!lineage[[1]] %in% lineage[[2]])]
    center_labels[idx] <- 4
    idx <- lineage[[2]][which(!lineage[[2]] %in% lineage[[1]])]
    center_labels[idx] <- 5
  }
  center_labels[is.na(center_labels)] <- 3

  list(cluster_labels = cluster_labels2, center_labels = center_labels)
}

cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
cell_type_vec <- as.factor(cell_type_vec)
cluster_labels <- as.numeric(cell_type_vec)
order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
cluster_group_list <- lapply(order_vec, function(x){
  grep(paste0("^", x), levels(cell_type_vec))
})
label_res <- reorder_cluster_labels(cluster_labels, cluster_group_list,
                                          our_curves$lineages)
cluster_center <- .compute_cluster_center(res_our$u_mat[,1:3], .construct_cluster_matrix(cluster_labels))

slingshot_3dplot(res_our$u_mat[,1:3], label_res$cluster_labels,
                 bg_col_vec = col_func(0.1),
                 breaks = seq(0.5, 5.5, by = 1),
                 cluster_center = cluster_center,
                 center_col_vec = col_func(1),
                 center_labels = label_res$center_labels,
                 curves = our_curves$curves,
                 pch = 16, lwd = 2, main = "Our lineages",
                 xlab = "Latent dimension 1",
                 ylab = "Latent dimension 2",
                 zlab = "Latent dimension 3")

# try a grid of angles
paramMat <- as.matrix(expand.grid(seq(0, 360, length.out = 9)[-1],
                                  seq(0, 360, length.out = 9)[-1]))
sapply(1:nrow(paramMat), function(x){
  if(x %% round(nrow(paramMat)/10) == 0) cat('*')

  png(paste0("../figure/experiment/Writeup22_3dplots/Writeup22_our_lineage_theta",
             paramMat[x,1], "_phi", paramMat[x,2], ".png"),
      height = 2000, width = 2000, res = 300, units = "px")
  par(mar = c(0.5, 0.5, 4, 0.5))
  slingshot_3dplot(res_our$u_mat[,1:3], label_res$cluster_labels,
                   bg_col_vec = col_func(0.1),
                   breaks = seq(0.5, 5.5, by = 1),
                   cluster_center = cluster_center,
                   center_col_vec = col_func(1),
                   center_labels = label_res$center_labels,
                   curves = our_curves$curves,
                   theta = paramMat[x,1], phi = paramMat[x,2],
                   pch = 16, lwd = 2, main = "Our lineages",
                   xlab = "Latent dimension 1",
                   ylab = "Latent dimension 2",
                   zlab = "Latent dimension 3")
  graphics.off()
})

####

label_res <- reorder_cluster_labels(cluster_labels, cluster_group_list)
cluster_center <- .compute_cluster_center(naive_embedding, .construct_cluster_matrix(cluster_labels))

sapply(1:nrow(paramMat), function(x){
  if(x %% round(nrow(paramMat)/10) == 0) cat('*')

  png(paste0("../figure/experiment/Writeup22_3dplots/Writeup22_naive_lineage_theta",
             paramMat[x,1], "_phi", paramMat[x,2], ".png"),
      height = 2000, width = 2000, res = 300, units = "px")
  par(mar = c(0.5, 0.5, 4, 0.5))
  slingshot_3dplot(naive_embedding, label_res$cluster_labels,
                   bg_col_vec = col_func(0.1),
                   breaks = seq(0.5, 5.5, by = 1),
                   cluster_center = cluster_center,
                   center_col_vec = col_func(1),
                   center_labels = label_res$center_labels,
                   curves = naive_curves$curves,
                   theta = paramMat[x,1], phi = paramMat[x,2],
                   pch = 16, lwd = 2, main = "Naive lineages",
                   xlab = "Latent dimension 1",
                   ylab = "Latent dimension 2",
                   zlab = "Latent dimension 3")
  graphics.off()
})

#######################

# add the tubes
our_tube_list <- lapply(1:length(our_curves$curves), function(x){
  s_mat <- our_curves$curves[[x]]$s[our_curves$curves[[x]]$ord,]
  construct_3d_tube(s_mat, radius = our_sd_val$sd_val)
})

col_tube_vec <- c(rgb(0/255, 114/255, 178/255), #blue
                  rgb(230/255, 159/255, 0/255)) #orange
label_res <- reorder_cluster_labels(cluster_labels, cluster_group_list,
                                    our_curves$lineages)
cluster_center <- .compute_cluster_center(res_our$u_mat[,1:3], .construct_cluster_matrix(cluster_labels))

sapply(1:nrow(paramMat), function(x){
  if(x %% round(nrow(paramMat)/10) == 0) cat('*')

  png(paste0("../figure/experiment/Writeup22_3dplots/Writeup22_our_lineage_tube_theta",
             paramMat[x,1], "_phi", paramMat[x,2], ".png"),
      height = 2000, width = 2000, res = 300, units = "px")
  par(mar = c(0.5, 0.5, 4, 0.5))
  slingshot_3dplot(res_our$u_mat[,1:3], label_res$cluster_labels,
                   bg_col_vec = col_func(0.1),
                   breaks = seq(0.5, 5.5, by = 1),
                   cluster_center = cluster_center,
                   center_col_vec = col_func(1),
                   center_labels = label_res$center_labels,
                   curves = our_curves$curves,
                   theta = paramMat[x,1], phi = paramMat[x,2],
                   pch = 16, lwd = 2, main = "Our lineages",
                   xlab = "Latent dimension 1",
                   ylab = "Latent dimension 2",
                   zlab = "Latent dimension 3")
  for(i in 1:length(our_tube_list)){
    # col_mat <- sapply(1:ncol(our_tube_list[[i]]$z_mat), function(x){
    #   rep(x, nrow(our_tube_list[[i]]$z_mat))
    # })
    col_mat <- our_tube_list[[i]]$z_mat

    plot3D::surf3D(our_tube_list[[i]]$x_mat,
                   our_tube_list[[i]]$y_mat,
                   our_tube_list[[i]]$z_mat, add = T,
                   colvar = col_mat,
                   col = colorRampPalette(c("white", col_tube_vec[i]))(100),
                   breaks = seq(min(col_mat), max(col_mat), length.out = 101),
                   colkey = F)
  }
  graphics.off()
})

###

label_res <- reorder_cluster_labels(cluster_labels, cluster_group_list)
cluster_center <- .compute_cluster_center(naive_embedding, .construct_cluster_matrix(cluster_labels))

# add the tubes
naive_tube_list <- lapply(1:length(naive_curves$curves), function(x){
  s_mat <- naive_curves$curves[[x]]$s[naive_curves$curves[[x]]$ord,]
  construct_3d_tube(s_mat, radius = naive_sd_val$sd_val)
})

col_tube_col <-rgb(0/255, 114/255, 178/255)

sapply(1:nrow(paramMat), function(x){
  if(x %% round(nrow(paramMat)/10) == 0) cat('*')

  png(paste0("../figure/experiment/Writeup22_3dplots/Writeup22_naive_lineage_tube_theta",
             paramMat[x,1], "_phi", paramMat[x,2], ".png"),
      height = 2000, width = 2000, res = 300, units = "px")
  par(mar = c(0.5, 0.5, 4, 0.5))
  slingshot_3dplot(naive_embedding, label_res$cluster_labels,
                   bg_col_vec = col_func(0.1),
                   breaks = seq(0.5, 5.5, by = 1),
                   cluster_center = cluster_center,
                   center_col_vec = col_func(1),
                   center_labels = label_res$center_labels,
                   curves = naive_curves$curves,
                   theta = paramMat[x,1], phi = paramMat[x,2],
                   pch = 16, lwd = 2, main = "Naive lineages",
                   xlab = "Latent dimension 1",
                   ylab = "Latent dimension 2",
                   zlab = "Latent dimension 3")

  i <- 3
  col_mat <- naive_tube_list[[i]]$z_mat

  plot3D::surf3D(naive_tube_list[[i]]$x_mat,
                 naive_tube_list[[i]]$y_mat,
                 naive_tube_list[[i]]$z_mat, add = T,
                 colvar = col_mat,
                 col = colorRampPalette(c(col_tube_col, "white"))(100),
                 breaks = seq(min(col_mat), max(col_mat), length.out = 101),
                 colkey = F)

  graphics.off()
})



