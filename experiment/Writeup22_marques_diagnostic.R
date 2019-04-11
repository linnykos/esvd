rm(list=ls())
load("../results/step5_clustering_spca.RData")


png("../figure/experiment/Writeup22_diagnostic.png",
    height = 1200, width = 2000, res = 300, units = "px")
par(mfrow = c(1,2))

idx <- which.min(quality_vec)
pred_mat <- 1/(res_list[[idx]]$u_mat %*% t(res_list[[idx]]$v_mat))
xlim <- range(c(pred_mat[missing_idx], pred_naive[missing_idx]))

tmp <- cbind(pred_mat[missing_idx], dat_impute[missing_idx])
pca_res <- princomp(tmp)
plot(pred_mat[missing_idx], dat_impute[missing_idx], asp = T, pch = 16,
     col = rgb(0,0,0,0.2), main = "Our embedding",
     xlab = "Predicted value", ylab = "Observed and masked value", xlim = xlim)
lines(c(-1e10,1e10), c(-1e10,1e10), col = "red", lwd = 2)
lines(c(-1e10,1e10)*pca_res$loadings[1,1], c(-1e10,1e10)*pca_res$loadings[1,2],
      col = "blue", lwd = 2, lty = 2)


##############

tmp <- cbind(pred_naive[missing_idx], dat_impute[missing_idx])
pca_res <- princomp(tmp)
plot(pred_naive[missing_idx], dat_impute[missing_idx], asp = T, pch = 16,
          col = rgb(0,0,0,0.2), main = "Naive embedding",
     xlab = "Predicted value", ylab = "Observed and masked value",
     xlim = xlim)
lines(c(-1e10,1e10), c(-1e10,1e10), col = "red", lwd = 2)
lines(c(-1e10,1e10)*pca_res$loadings[1,1], c(-1e10,1e10)*pca_res$loadings[1,2],
      col = "blue", lwd = 2, lty = 2)
graphics.off()

