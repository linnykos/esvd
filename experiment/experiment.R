rm(list=ls())

set.seed(10)
n <- 100
u_mat <- abs(matrix(rnorm(n), nrow = n, ncol = 1))
v_mat <- abs(matrix(rnorm(n), nrow = n, ncol = 1))
pred_mat <- u_mat %*% t(v_mat)
dat <- pred_mat

for(i in 1:n){
  for(j in 1:n){
    dat[i,j] <- max(stats::rnorm(1, mean = 1/pred_mat[i,j], sd = 1/(2*pred_mat[i,j])), 1e-3)
  }
}

missing_idx <- eSVD::construct_missing_values(n = nrow(dat), p = ncol(dat))
family <- "curved_gaussian"

fit <- .tuning_fit(dat, family = family, scalar = 2, k = 1, max_val = 100)
if(verbose) print("Finished initial fit")

nat_mat <- fit$u_mat %*% t(fit$v_mat)
# mean_mat <- compute_mean(nat_mat, family = "exponential")
# x <- ((dat - mean_mat)^2)
# y <- (mean_mat^2/4)
# plot(x,y, asp = T, pch = 16)

fn <- function(x){
  mean_mat <- compute_mean(nat_mat, family = "curved_gaussian")
  tmp_mat <- cbind(as.numeric((dat-mean_mat)^2), as.numeric(mean_mat^2/x^2))

  pca_res <- stats::prcomp(tmp_mat, center = F, scale = F)
  ang <- as.numeric(acos(abs(c(0,1) %*% pca_res$rotation[,1])))
  abs(45 - ang* 180/pi)
}

res <- stats::optimize(fn, interval = c(1, 10))

mean_mat <- compute_mean(nat_mat, family = "curved_gaussian")
vec1 = as.numeric((dat-mean_mat)^2)
vec2 = as.numeric(mean_mat^2/x^2)
plot(vec1,vec2,asp=T)
tmp_mat <- cbind(vec1, vec2)
pca_res <- stats::prcomp(tmp_mat, center = F, scale = F)
ang <- as.numeric(acos(abs(c(0,1) %*% pca_res$rotation[,1])))
ang*180/pi

