rm(list=ls())
load("../experiment/Week36_marques_slingshot.RData")

combn_mat <- utils::combn(d, 2)
col_vec <- numeric(length(unique(cluster_labels)))
alpha_val <- 1
for(i in 1:3){
  col_vec[cluster_group_list[[i]]] <- rgb(238/255,204/255,17/255,alpha_val)
}
col_vec[cluster_group_list[[4]]] <- rgb(129/255,199/255,124/255,alpha_val)
col_vec[cluster_group_list[[5]]] <- rgb(227/255,73/255,86/255,alpha_val)
col_vec[cluster_group_list[[6]]] <- rgb(100/255,140/255,252/255,alpha_val)

plot_function <- function(dat, curves, name){
  reduction_factor <- max(apply(dat, 2, function(x){diff(range(x))}))*.25
  range_mat <- apply(dat, 2, range)
  cluster_center <- .compute_cluster_center(dat, .construct_cluster_matrix(cluster_labels))

  png(paste0("../figure/experiment/Week36_marques_", name, "_lineage.png"), height = length(curves$curves)*2000/2.5, width = 2000, res = 300, units = "px")
  par(mfrow = c(length(curves$curves), ncol(combn_mat)), mar = c(4,4,4,0.5))
  for(k in 1:length(curves$curves)){
    for(i in 1:ncol(combn_mat)){
      cell_idx <- which(cluster_labels %in% curves$lineages[[k]])

      idx1 <- combn_mat[1,i]; idx2 <- combn_mat[2,i]
      plot(dat[cell_idx,idx1], dat[cell_idx,idx2], pch = 16, col = rgb(0.85,0.85,0.85,1),
           asp = T, cex = 1,
           xlim = range_mat[,idx1], ylim = range_mat[,idx2],
           xlab = paste0("Latent dimension ", idx1),
           ylab = paste0("Latent dimension ", idx2),
           main = ifelse(i == 2, paste0(name, " Lineage ", k), ""))

      # plot curves
      ord <- curves$curves[[k]]$ord
      lines(curves$curves[[k]]$s[ord, idx1]*reduction_factor,
            curves$curves[[k]]$s[ord, idx2]*reduction_factor, lwd = 3.5,
            col = "white")
      lines(curves$curves[[k]]$s[ord, idx1]*reduction_factor,
            curves$curves[[k]]$s[ord, idx2]*reduction_factor, lwd = 3,
            col = "black")

      # plot points
      points(cluster_center[curves$lineages[[k]], idx1],
             cluster_center[curves$lineages[[k]], idx2], pch = 16,
             col = col_vec[curves$lineages[[k]]], cex = 2)
    }
  }
  graphics.off()
  invisible()
}

plot_function(naive, naive_curves, "Naive")
plot_function(u_mat, our_curves, "Our")
