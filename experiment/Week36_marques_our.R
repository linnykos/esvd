rm(list=ls())
load("../results/step4_factorization_spca.RData")

d <- 3
dat <- res_our$u_mat[,1:d]

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
  xlim <- range(dat[,i1])
  ylim <- range(dat[,i2])
  order_vec <- c(6,5,1,4,2,3) #these are in alphabetical order from "CO", "MF", "MO", "NF", "OP", "PP"
  name_vec <- c("Pdgfra+ (1)", "Precusor (2)", "Differentiated-commited\nprecusor (3)", "Newly formed (4)",
                "Myelin-forming (5)", "Mature (6)")

  png(paste0("../figure/experiment/Week36_marques_latent_our_", i1, i2, ".png"), height = 1500, width = 2000, res = 300, units = "px")
  par(mfrow = c(2,3), mar = c(4,4,4,0.5))
  for(i in 1:6){
    idx <- which(as.numeric(as.factor(cell_type_vec)) == order_vec[i])
    plot(dat[-idx,i1], dat[-idx,i2], asp = T, pch = 16,
         col = rgb(0.8, 0.8, 0.8), xlim = xlim, ylim = ylim,
         main = name_vec[i], xlab = paste0("Our latent dimension ", i1),
         ylab = paste0("Our latent dimension ", i2))

    lines(c(-100,100), rep(0,2), lwd = 2, lty = 2)
    lines(rep(0,2), c(-100,100), lwd = 2, lty = 2)

    points(dat[idx,i1], dat[idx,i2], pch = 16,
           col = rgb(1,1,1), cex = 1.5)

    points(dat[idx,i1], dat[idx,i2], pch = 16,
           col = col_idx[idx], cex = 1.5)
  }
  graphics.off()
}

####################################

cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
cell_type_vec <- as.factor(cell_type_vec)
cluster_labels <- as.numeric(cell_type_vec)
order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
cluster_group_list <- lapply(order_vec, function(x){
  grep(paste0("^", x), levels(cell_type_vec))
})



d <- 3
dat <- res_our$u_mat[,1:d]
reduction_factor <- max(apply(dat, 2, function(x){diff(range(x))}))*.25
set.seed(10)
upscale_vec <- rep(NA, length(unique(cluster_labels)))
size_vec <- sapply(cluster_group_list, function(x){length(which(cluster_labels %in% x))})
for(i in 1:length(cluster_group_list)){
  upscale_vec[cluster_group_list[[i]]] <- (max(size_vec)/size_vec[i])^(1/2)
}

# run slingshot
curves <- singlecell::slingshot(dat/reduction_factor, cluster_labels, starting_cluster = cluster_group_list[[1]][1], cluster_group_list = cluster_group_list, verbose = T,
                    b = 1, upscale_vec = upscale_vec)
save.image("tmp.RData")

# load("../experiment/tmp.RData")
# dat2 <- dat/reduction_factor
# combn_mat <- utils::combn(d, 2)
# range_mat <- apply(dat2, 2, range)
# col_vec <- numeric(length(unique(cluster_labels)))
# alpha_val <- 1
# for(i in 1:3){
#   col_vec[cluster_group_list[[i]]] <- rgb(238/255,204/255,17/255,alpha_val)
# }
# col_vec[cluster_group_list[[4]]] <- rgb(129/255,199/255,124/255,alpha_val)
# col_vec[cluster_group_list[[5]]] <- rgb(227/255,73/255,86/255,alpha_val)
# col_vec[cluster_group_list[[6]]] <- rgb(100/255,140/255,252/255,alpha_val)
# cluster_center <- .compute_cluster_center(dat2, .construct_cluster_matrix(cluster_labels))
#
# png("../figure/experiment/Week36_marques_our_lineage.png", height = length(curves$curves)*2000/2.5, width = 2000, res = 300, units = "px")
# par(mfrow = c(length(curves$curves), ncol(combn_mat)), mar = c(4,4,4,0.5))
# for(k in 1:length(curves$curves)){
#   for(i in 1:ncol(combn_mat)){
#     cell_idx <- which(cluster_labels %in% curves$lineages[[k]])
#
#     idx1 <- combn_mat[1,i]; idx2 <- combn_mat[2,i]
#     plot(dat2[cell_idx,idx1], dat2[cell_idx,idx2], pch = 16, col = rgb(0.85,0.85,0.85,1),
#          asp = T, cex = 1,
#          xlim = range_mat[,idx1], ylim = range_mat[,idx2],
#          xlab = paste0("Latent dimension ", idx1),
#          ylab = paste0("Latent dimension ", idx2),
#          main = ifelse(i == 2, paste0("Our Lineage ", k), ""))
#
#     # plot curves
#     ord <- curves$curves[[k]]$ord
#     lines(curves$curves[[k]]$s[ord, idx1], curves$curves[[k]]$s[ord, idx2], lwd = 3.5,
#           col = "white")
#     lines(curves$curves[[k]]$s[ord, idx1], curves$curves[[k]]$s[ord, idx2], lwd = 3,
#           col = "black")
#
#     # plot points
#     points(cluster_center[curves$lineages[[k]], idx1],
#            cluster_center[curves$lineages[[k]], idx2], pch = 16,
#            col = col_vec[curves$lineages[[k]]], cex = 2)
#   }
# }
# graphics.off()
