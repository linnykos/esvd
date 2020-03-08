rm(list=ls())
load("../results/darmanis.RData")
label_vec <- as.factor(as.character(label_vec))
label_num <- as.numeric(label_vec)

plot(gaussian_fit$u_mat[,1], gaussian_fit$u_mat[,2], asp = T, col = c(1:5)[label_num],
     pch = 16, main = "Gaussian")
plot(poisson_fit$u_mat[,1], poisson_fit$u_mat[,2], asp = T, col = c(1:5)[label_num],
     pch = 16, main = "Poisson")
plot(exponential_fit$u_mat[,1], exponential_fit$u_mat[,2], asp = T, col = c(1:5)[label_num],
     pch = 16, main = "Exponential")
plot(neg_binom_fit$u_mat[,1], neg_binom_fit$u_mat[,2], asp = T, col = c(1:5)[label_num],
     pch = 16, main = "Negative binomial")
plot(curved_gaussian_fit$u_mat[,1], curved_gaussian_fit$u_mat[,2], asp = T, col = c(1:5)[label_num],
     pch = 16, main = "Curved Gaussian")
neg_bin_param
curved_gaussian_param

#########################

obs_val <- dat_impute[missing_idx]
pred_val <- compute_mean(gaussian_fit_missing$u_mat %*% t(gaussian_fit_missing$v_mat), family = "gaussian")
pred_val <- pred_val[missing_idx]
all_mat <- cbind(pred_val, obs_val)
pca_res <- stats::prcomp(all_mat[,2:1], center = F, scale = F)
plot(pred_val, obs_val, asp = T, pch = 16, col = rgb(0.5,0.5,0.5,0.5))
lines(c(-1e4,1e4), c(-1e4,1e4), col = "red", lwd = 2, lty = 2)
lines(c(-1e10,1e10)*pca_res$rotation[2,1], c(-1e10,1e10)*pca_res$rotation[1,1],
      col = "blue", lwd = 2, lty = 2)

obs_val <- dat_impute[missing_idx]
pred_val <- compute_mean(poisson_fit_missing$u_mat %*% t(poisson_fit_missing$v_mat), family = "poisson")
pred_val <- pred_val[missing_idx]
all_mat <- cbind(pred_val, obs_val)
pca_res <- stats::prcomp(all_mat[,2:1], center = F, scale = F)
plot(pred_val, obs_val, asp = T, pch = 16, col = rgb(0.5,0.5,0.5,0.5))
lines(c(-1e4,1e4), c(-1e4,1e4), col = "red", lwd = 2, lty = 2)
lines(c(-1e10,1e10)*pca_res$rotation[2,1], c(-1e10,1e10)*pca_res$rotation[1,1],
      col = "blue", lwd = 2, lty = 2)

obs_val <- dat_impute[missing_idx]
pred_val <- compute_mean(exponential_fit_missing$u_mat %*% t(exponential_fit_missing$v_mat), family = "exponential")
pred_val <- pred_val[missing_idx]
all_mat <- cbind(pred_val, obs_val)
pca_res <- stats::prcomp(all_mat[,2:1], center = F, scale = F)
plot(pred_val, obs_val, asp = T, pch = 16, col = rgb(0.5,0.5,0.5,0.5))
lines(c(-1e4,1e4), c(-1e4,1e4), col = "red", lwd = 2, lty = 2)
lines(c(-1e10,1e10)*pca_res$rotation[2,1], c(-1e10,1e10)*pca_res$rotation[1,1],
      col = "blue", lwd = 2, lty = 2)

obs_val <- dat_impute[missing_idx]
pred_val <- compute_mean(neg_binom_fit_missing$u_mat %*% t(neg_binom_fit_missing$v_mat), family = "neg_binom", size = neg_bin_param)
pred_val <- pred_val[missing_idx]
all_mat <- cbind(pred_val, obs_val)
pca_res <- stats::prcomp(all_mat[,2:1], center = F, scale = F)
plot(pred_val, obs_val, asp = T, pch = 16, col = rgb(0.5,0.5,0.5,0.5))
lines(c(-1e4,1e4), c(-1e4,1e4), col = "red", lwd = 2, lty = 2)
lines(c(-1e10,1e10)*pca_res$rotation[2,1], c(-1e10,1e10)*pca_res$rotation[1,1],
      col = "blue", lwd = 2, lty = 2)
