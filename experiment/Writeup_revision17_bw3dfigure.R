rm(list=ls())
load("../results/submission_round2b_revision/step7_additional_analyses.RData")

# fix the orientation for reproducibility
if(mean(esvd_embedding$u_mat[,1]) > 0) {
  esvd_embedding$u_mat[,1] <- -esvd_embedding$u_mat[,1]
  esvd_embedding$v_mat[,1] <- -esvd_embedding$v_mat[,1]
}

if(mean(esvd_embedding$u_mat[,2]) < 0) {
  esvd_embedding$u_mat[,2] <- -esvd_embedding$u_mat[,2]
  esvd_embedding$v_mat[,2] <- -esvd_embedding$v_mat[,2]
}

if(mean(esvd_embedding$u_mat[,3]) < 0) {
  esvd_embedding$u_mat[,3] <- -esvd_embedding$u_mat[,3]
  esvd_embedding$v_mat[,3] <- -esvd_embedding$v_mat[,3]
}

angle_matrix <- matrix(c(90,225, 180,225, 270,225), byrow = T, ncol = 2)
spacing <- 4.5; tmp <- spacing*c(-.5,.5)
bound_matrix <- list(list(xlim = -2+tmp, ylim = 0+tmp, zlim = 1+tmp),
                     list(xlim = -2.5+tmp, ylim = 0.75+tmp, zlim = 1+tmp),
                     list(xlim = -3+tmp, ylim = 0.75+tmp, zlim = 1+tmp))

kk <- 2
png(paste0("../../esvd_results/figure/experiment/eSVD_theta", angle_matrix[kk,1], "_phi", angle_matrix[kk,2], "_bw.png"),
    height = 2000, width = 2000, res = 300, units = "px")
par(mar = c(0,0,4,0))
eSVD::slingshot_3dplot(esvd_embedding$u_mat[,1:3], cluster_labels,
                       bg_col_vec = rep(rgb(0.5, 0.5, 0.5, 0.5), max(cluster_labels)), bg_cex = 0.8,
                       curves = NA,
                       pch = 16, main = "eSVD embedding and trajectories:\nCurved Gaussian",
                       xlab = "Latent dimension 1", ylab = "Latent dimension 2",
                       zlab = "Latent dimension 3",
                       theta = angle_matrix[kk,1], phi = angle_matrix[kk,2],
                       xlim = bound_matrix[[kk]]$xlim,
                       ylim = bound_matrix[[kk]]$ylim,
                       zlim = bound_matrix[[kk]]$zlim)
graphics.off()
