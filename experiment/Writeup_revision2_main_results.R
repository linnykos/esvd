rm(list=ls())
load("../results/step4_factorization.RData")

quality_vec

# diagnostic of eSVD
idx <- which.min(quality_vec)
pred_mat <- 1/(res_list[[idx]]$u_mat %*% t(res_list[[idx]]$v_mat))
tmp <- cbind(pred_mat[missing_idx], dat_impute[missing_idx])
pca_res <- princomp(tmp)

x_val <- seq(1, 1e5, length.out = 100)
y_val_top <- sapply(x_val, function(x){stats::qnorm(0.75, mean = x, sd = x/scalar_val)})
y_val_bottom <- sapply(x_val, function(x){stats::qnorm(0.25, mean = x, sd = x/scalar_val)})

plot(NA, asp = T, pch = 16,
     xlim = c(0, max(c(pred_mat[missing_idx], dat_impute[missing_idx]))),
     ylim = c(0, max(c(pred_mat[missing_idx], dat_impute[missing_idx]))),
     main = "eSVD embedding:\nMatrix-completion diagnostic",
     xlab = "Predicted value", ylab = "Observed but masked value",
     cex.lab = 1.25, cex.axis = 1.25)

polygon(c(x_val, rev(x_val)), c(y_val_top, rev(y_val_bottom)), col = rgb(1,0,0,0.2),
        border = NA, density = 30, angle = -45)
points(pred_mat[missing_idx], dat_impute[missing_idx], pch = 16,
       col = rgb(0,0,0,0.2))

lines(rep(0,2), c(-1e10,1e10), col = "red", lwd = 1)
lines(c(-1e10,1e10), rep(0,2), col = "red", lwd = 1)
lines(c(-1e10,1e10), c(-1e10,1e10), col = "red", lwd = 2)
lines(c(-1e10,1e10)*pca_res$loadings[1,1], c(-1e10,1e10)*pca_res$loadings[1,2],
      col = "blue", lwd = 2, lty = 2)
lines(x_val, y_val_top, col = "red", lwd = 2, lty = 2)
lines(x_val, y_val_bottom, col = "red", lwd = 2, lty = 2)

rad <- 500
ang <- as.numeric(acos(c(1,0) %*% pca_res$loadings[,1]))
radian_seq <- seq(0, ang, length.out = 100)
x_circ <- rad * cos(radian_seq)
y_circ <- rad * sin(radian_seq)
lines(x_circ, y_circ, lty = 2)
text(x = 500, y = 200, pos = 4, label = paste0(round(ang * 180/pi, 1), " degrees"))

#############################################

plot(res_our$u_mat[,1], res_our$u_mat[,2], asp = T, pch = 16, col = rgb(0.5, 0.5, 0.5, 0.5))
plot(res_our$u_mat[,1], res_our$u_mat[,3], asp = T, pch = 16, col = rgb(0.5, 0.5, 0.5, 0.5))
plot(res_our$u_mat[,2], res_our$u_mat[,3], asp = T, pch = 16, col = rgb(0.5, 0.5, 0.5, 0.5))
