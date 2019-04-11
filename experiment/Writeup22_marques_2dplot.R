rm(list=ls())
load("../results/step5_clustering_spca.RData")

d <- 3
dat <- res_our$u_mat[,1:d]

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
cell_type_vec2 <- as.character(sapply(cell_type_vec, function(x){substr(x,1,2)}))
cell_type_vec <- as.factor(cell_type_vec)
cluster_labels <- as.numeric(cell_type_vec)
order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
cluster_group_list <- lapply(order_vec, function(x){
  grep(paste0("^", x), levels(cell_type_vec))
})
label_res <- reorder_cluster_labels(cluster_labels, cluster_group_list,
                                    our_curves$lineages)
col_idx <- col_func(0.2)[label_res$cluster_labels]

num_cell <- length(unique(cell_type_vec))
cluster_center <- .compute_cluster_center(dat, .construct_cluster_matrix(cluster_labels))

i1 <- 1; i2 <- 2
xlim <- range(dat[,i1])
ylim <- range(dat[,i2])
order_vec <- c(6,5,1,4,2,3) #these are in alphabetical order from "CO", "MF", "MO", "NF", "OP", "PP"
name_vec <- c("Pdgfra+ (1)", "Precusor (2)", "Differentiated-commited\nprecusor (3)", "Newly formed (4)",
              "Myelin-forming (5)", "Mature (6)")

png(paste0("../figure/experiment/Writeup22_marques_latent_our_", i1, i2, ".png"), height = 1500, width = 2000, res = 300, units = "px")
par(mfrow = c(2,3), mar = c(4,4,4,0.5))
for(i in 1:6){
  idx <- which(as.numeric(as.factor(cell_type_vec2)) == order_vec[i])
  plot(dat[-idx,i1], dat[-idx,i2], asp = T, pch = 16,
       col = rgb(0.8, 0.8, 0.8), xlim = xlim, ylim = ylim,
       main = name_vec[i], xlab = paste0("Our latent dimension ", i1),
       ylab = paste0("Our latent dimension ", i2))

  lines(c(-100,100), rep(0,2), lwd = 2, lty = 2)
  lines(rep(0,2), c(-100,100), lwd = 2, lty = 2)

  # plot the points
  points(dat[idx,i1], dat[idx,i2], pch = 16,
         col = rgb(1,1,1), cex = 1.5)

  points(dat[idx,i1], dat[idx,i2], pch = 16,
         col = col_idx[idx], cex = 1.5)


  # plot the curves
  col_vec <- c(rgb(0/255, 114/255, 178/255), rgb(230/255, 159/255, 0/255))
  for(i in 1:length(our_curves$curves)){
    s_mat <- our_curves$curves[[i]]$s[our_curves$curves[[i]]$ord,]
    lines(s_mat[,i1], s_mat[,i2], col = "black", lwd = 5)
    lines(s_mat[,i1], s_mat[,i2], col = col_vec[i], lwd = 2)
  }

  # plot the centers
  cell_subtypes <- unique(as.numeric(as.factor(cell_type_vec))[idx]) #what cell subtypes
  for(i in 1:length(cell_subtypes)){
    points(cluster_center[cell_subtypes[i],i1], cluster_center[cell_subtypes[i],i2],
           pch = 16, cex = 2.5, col = "black")
    points(cluster_center[cell_subtypes[i],i1], cluster_center[cell_subtypes[i],i2],
           pch = 16, cex = 2, col = col_func(1)[label_res$center_labels[cell_subtypes[i]]])
  }


}
graphics.off()

##############

# do the same thing for naive embedding
label_res <- reorder_cluster_labels(cluster_labels, cluster_group_list)
col_idx <- col_func(0.2)[label_res$cluster_labels]

dat <- naive_embedding
cluster_center <- .compute_cluster_center(dat, .construct_cluster_matrix(cluster_labels))

i1 <- 1; i2 <- 3
xlim <- range(dat[,i1])
ylim <- range(dat[,i2])
order_vec <- c(6,5,1,4,2,3) #these are in alphabetical order from "CO", "MF", "MO", "NF", "OP", "PP"
name_vec <- c("Pdgfra+ (1)", "Precusor (2)", "Differentiated-commited\nprecusor (3)", "Newly formed (4)",
              "Myelin-forming (5)", "Mature (6)")

png(paste0("../figure/experiment/Writeup22_marques_latent_naive_", i1, i2, ".png"), height = 1500, width = 2000, res = 300, units = "px")
par(mfrow = c(2,3), mar = c(4,4,4,0.5))
for(i in 1:6){
  idx <- which(as.numeric(as.factor(cell_type_vec2)) == order_vec[i])
  plot(dat[-idx,i1], dat[-idx,i2], asp = T, pch = 16,
       col = rgb(0.8, 0.8, 0.8), xlim = xlim, ylim = ylim,
       main = name_vec[i], xlab = paste0("Naive latent dimension ", i1),
       ylab = paste0("Naive latent dimension ", i2))

  lines(c(-100,100), rep(0,2), lwd = 2, lty = 2)
  lines(rep(0,2), c(-100,100), lwd = 2, lty = 2)

  # plot the points
  points(dat[idx,i1], dat[idx,i2], pch = 16,
         col = rgb(1,1,1), cex = 1.5)

  points(dat[idx,i1], dat[idx,i2], pch = 16,
         col = col_idx[idx], cex = 1.5)


  # plot the curves
  for(i in 1:length(naive_curves$curves)){
    s_mat <- naive_curves$curves[[i]]$s[naive_curves$curves[[i]]$ord,]
    lines(s_mat[,i1], s_mat[,i2], col = "black", lwd = 2)
  }

  # plot the centers
  cell_subtypes <- unique(as.numeric(as.factor(cell_type_vec))[idx]) #what cell subtypes
  for(i in 1:length(cell_subtypes)){
    points(cluster_center[cell_subtypes[i],i1], cluster_center[cell_subtypes[i],i2],
           pch = 16, cex = 2.5, col = "black")
    points(cluster_center[cell_subtypes[i],i1], cluster_center[cell_subtypes[i],i2],
           pch = 16, cex = 2, col = col_func(1)[label_res$center_labels[cell_subtypes[i]]])
  }


}
graphics.off()
