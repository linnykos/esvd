rm(list=ls())
load("../results/step4_factorization_spca.RData")
zz <- svd(dat)
naive <- zz$u[,1:3]%*%diag(zz$d[1:3])

cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
# keep only the first two letters
cell_type_vec <- as.character(sapply(cell_type_vec, function(x){substr(x,1,2)}))
alpha_val <- 0.2
col_idx <- c(rgb(238/255,204/255,17/255,alpha_val), #goldenrod
             rgb(227/255,73/255,86/255,alpha_val), #red
             rgb(100/255,140/255,252/255,alpha_val), #blue
             rgb(129/255,199/255,124/255,alpha_val), #green
             rgb(238/255,204/255,17/255,alpha_val),
             rgb(238/255,204/255,17/255,alpha_val))[as.numeric(as.factor(cell_type_vec))]

num_cell <- length(unique(cell_type_vec))

combn_mat <- combn(3,2)
for(x in 1:ncol(combn_mat)){
  i1 <- combn_mat[1,x]; i2 <- combn_mat[2,x]
  xlim <- range(naive[,i1])
  ylim <- range(naive[,i2])
  order_vec <- c(6,5,1,4,2,3) #these are in alphabetical order from "CO", "MF", "MO", "NF", "OP", "PP"
  name_vec <- c("Pdgfra+ (1)", "Precusor (2)", "Differentiated-commited\nprecusor (3)", "Newly formed (4)",
                "Myelin-forming (5)", "Mature (6)")

  png(paste0("../figure/experiment/Week36_marques_latent_naive_", i1, i2, ".png"), height = 1500, width = 2000, res = 300, units = "px")
  par(mfrow = c(2,3), mar = c(4,4,4,0.5))
  for(i in 1:6){
    idx <- which(as.numeric(as.factor(cell_type_vec)) == order_vec[i])
    plot(naive[-idx,i1], naive[-idx,i2], asp = T, pch = 16,
         col = rgb(0.8, 0.8, 0.8), xlim = xlim, ylim = ylim,
         main = name_vec[i], xlab = paste0("Latent dimension ", i1),
         ylab = paste0("Latent dimension ", i2))

    lines(c(-100,100), rep(0,2), lwd = 2, lty = 2)
    lines(rep(0,2), c(-100,100), lwd = 2, lty = 2)

    points(naive[idx,i1], naive[idx,i2], pch = 16,
           col = rgb(1,1,1), cex = 1.5)

    points(naive[idx,i1], naive[idx,i2], pch = 16,
           col = col_idx[idx], cex = 1.5)
  }
  graphics.off()
}

#######################
cell_type_vec_raw <- as.character(marques$cell.info$cell.type[cell_idx])
combn_mat <- utils::combn(3,2)

alpha_val <- 0.5
col_vec <- c(rgb(100/255,140/255,252/255,alpha_val), #blue
             rgb(227/255,73/255,86/255,alpha_val), #red
             rgb(238/255,204/255,17/255,alpha_val), #goldenrod
             rgb(129/255,199/255,124/255,alpha_val), #green
             rgb(100/255,100/255,200/255,alpha_val), # purple
             rgb(40/255,225/255,201/255,alpha_val) #turqouise
)

png("../figure/experiment/Week36_marques_MO_1_naive.png", height = 2250, width = 2000, res = 300, units = "px")
idx <- grep("MO", cell_type_vec_raw)
len <- length(unique(cell_type_vec_raw[idx]))
par(mfrow = c(3,3))
for(j in c(2,1,6)){
  for(i in 1:ncol(combn_mat)){
    i1 <- combn_mat[1,i]; i2 <- combn_mat[2,i]
    idx_j <- which(cell_type_vec_raw[idx] == unique(cell_type_vec_raw[idx])[j])

    if(i==2) main_vec <- paste0("Mature sub-type ", j, "\n(", length(idx_j), ")") else main_vec = ""
    plot(naive[idx,i1], naive[idx,i2], asp = T, pch = 16,
         col = rgb(0.8, 0.8, 0.8),
         xlim = range(c(0, naive[idx,i1])), ylim = range(c(0, naive[idx,i2])),
         main = main_vec,
         xlab = paste0("Latent dimension ", i1),
         ylab = paste0("Latent dimension ", i2))

    lines(c(-100,100), rep(0,2), lwd = 2, lty = 2)
    lines(rep(0,2), c(-100,100), lwd = 2, lty = 2)

    points(naive[idx[idx_j],i1], naive[idx[idx_j],i2], pch = 16,
           col = rgb(1,1,1), cex = 2)

    points(naive[idx[idx_j],i1], naive[idx[idx_j],i2], pch = 16,
           col = col_vec[j], cex = 1.5)

    mean_vec <- colMeans(naive[idx[idx_j],c(i1, i2)])
    ellipse_points <- .compute_ellipse_points(mean_vec,
                                              stats::cov(naive[idx[idx_j],c(i1, i2)]),
                                              scale_factor = 2)
    points(mean_vec[1], mean_vec[2], pch = 16)
    lines(ellipse_points[,1], ellipse_points[,2])
  }
}
graphics.off()

png("../figure/experiment/Week36_marques_MO_2_naive.png", height = 2250, width = 2000, res = 300, units = "px")
idx <- grep("MO", cell_type_vec_raw)
len <- length(unique(cell_type_vec_raw[idx]))
par(mfrow = c(3,3))
for(j in c(3,4,5)){
  for(i in 1:ncol(combn_mat)){
    i1 <- combn_mat[1,i]; i2 <- combn_mat[2,i]
    idx_j <- which(cell_type_vec_raw[idx] == unique(cell_type_vec_raw[idx])[j])

    if(i==2) main_vec <- paste0("Mature sub-type ", j, "\n(", length(idx_j), ")") else main_vec = ""
    plot(naive[idx,i1], naive[idx,i2], asp = T, pch = 16,
         col = rgb(0.8, 0.8, 0.8),
         xlim = range(c(0, naive[idx,i1])), ylim = range(c(0, naive[idx,i2])),
         main = main_vec,
         xlab = paste0("Latent dimension ", i1),
         ylab = paste0("Latent dimension ", i2))

    lines(c(-100,100), rep(0,2), lwd = 2, lty = 2)
    lines(rep(0,2), c(-100,100), lwd = 2, lty = 2)

    points(naive[idx[idx_j],i1], naive[idx[idx_j],i2], pch = 16,
           col = rgb(1,1,1), cex = 2)

    points(naive[idx[idx_j],i1], naive[idx[idx_j],i2], pch = 16,
           col = col_vec[j], cex = 1.5)

    mean_vec <- colMeans(naive[idx[idx_j],c(i1, i2)])
    ellipse_points <- .compute_ellipse_points(mean_vec,
                                              stats::cov(naive[idx[idx_j],c(i1, i2)]),
                                              scale_factor = 2)
    points(mean_vec[1], mean_vec[2], pch = 16)
    lines(ellipse_points[,1], ellipse_points[,2])
  }
}
graphics.off()

###################

# try finding the lineage

cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
cell_type_vec <- as.factor(cell_type_vec)
cluster_labels <- as.numeric(cell_type_vec)
order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
cluster_group_list <- lapply(order_vec, function(x){
  grep(paste0("^", x), levels(cell_type_vec))
})

reduction_factor <- max(apply(naive, 2, function(x){diff(range(x))}))*.15
curves <- slingshot(naive/reduction_factor, cluster_labels, starting_cluster = cluster_group_list[[1]][1], cluster_group_list = cluster_group_list, verbose = T,
                    b = 1)
lineages <- .get_lineages(naive, cluster_labels, starting_cluster = cluster_group_list[[1]][1],
                          cluster_group_list = cluster_group_list)

# let's try bootstrapping a bit
trials <- 100
curve_list <- lapply(1:trials, function(x){
  print(x)
  set.seed(10*x)
  naive2 <- naive
  for(i in 1:length(unique(cluster_labels))){
    idx <- which(cluster_labels == i)
    idx2 <- sample(idx, length(idx), replace = T)
    naive2[idx,] <- naive[idx2,]
  }

  .get_lineages(naive2, cluster_labels, starting_cluster = cluster_group_list[[1]][1],
                cluster_group_list = cluster_group_list)

})

# we need a way to get an confidence tube hm...

#####################################

naive2 <- naive/reduction_factor
combn_mat <- utils::combn(3, 2)
range_mat <- apply(naive2, 2, range)
col_vec <- numeric(length(unique(cluster_labels)))
alpha_val <- 1
for(i in 1:3){
  col_vec[cluster_group_list[[i]]] <- rgb(238/255,204/255,17/255,alpha_val)
}
col_vec[cluster_group_list[[4]]] <- rgb(129/255,199/255,124/255,alpha_val)
col_vec[cluster_group_list[[5]]] <- rgb(227/255,73/255,86/255,alpha_val)
col_vec[cluster_group_list[[6]]] <- rgb(100/255,140/255,252/255,alpha_val)
cluster_center <- .compute_cluster_center(naive2, .construct_cluster_matrix(cluster_labels))

png("../figure/experiment/Week36_marques_naive_lineage.png", height = length(curves$curves)*2000/2.5, width = 2000, res = 300, units = "px")
par(mfrow = c(length(curves$curves),3), mar = c(4,4,4,0.5))
for(k in 1:length(curves$curves)){
  for(i in 1:3){
    cell_idx <- which(cluster_labels %in% curves$lineages[[k]])

    idx1 <- combn_mat[1,i]; idx2 <- combn_mat[2,i]
    plot(naive2[cell_idx,idx1], naive2[cell_idx,idx2], pch = 16, col = rgb(0.85,0.85,0.85,1),
         asp = T, cex = 1,
         xlim = range_lis[,idx1], ylim = range_lis[,idx2],
         xlab = paste0("Latent dimension ", idx1),
         ylab = paste0("Latent dimension ", idx2),
         main = ifelse(i == 2, paste0("Naive Lineage ", k), ""))

    # plot curves
    ord <- curves$curves[[k]]$ord
    lines(curves$curves[[k]]$s[ord, idx1], curves$curves[[k]]$s[ord, idx2], lwd = 3.5,
          col = "white")
    lines(curves$curves[[k]]$s[ord, idx1], curves$curves[[k]]$s[ord, idx2], lwd = 3,
          col = "black")

    # plot points
    points(cluster_center[curves$lineages[[k]], idx1],
           cluster_center[curves$lineages[[k]], idx2], pch = 16,
           col = col_vec[curves$lineages[[k]]], cex = 2)
  }
}
graphics.off()
