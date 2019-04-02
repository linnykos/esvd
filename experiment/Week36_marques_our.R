rm(list=ls())
load("../results/step4_factorization_spca.RData")

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
  upscale_vec[cluster_group_list[[i]]] <- max(size_vec)/size_vec[i]
}

# run slingshot
curves <- singlecell::slingshot(dat/reduction_factor, cluster_labels, starting_cluster = cluster_group_list[[1]][1], cluster_group_list = cluster_group_list, verbose = T,
                    b = 1, upscale_vec = upscale_vec)
save.image("tmp.RData")
#
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
