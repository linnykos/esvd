rm(list=ls())
load("../experiment/tmp.RData")
paramMat <- cbind(50, 120, 0.05, 150, 2, 2, -2000)
colnames(paramMat) <- c("n_each", "d_each", "sigma", "total", "k", "scalar", "max_val")
vec <- paramMat[1,]

plot(res_our[,1], res_our[,2], asp = T, pch = 16, col = c(1:4)[rep(1:4, each = paramMat[1,"n_each"])])

dat <- res_our
cluster_labels <- rep(1:4, each = vec["n_each"])
target_percentage <- 0.8

stopifnot(ncol(dat) == 2, all(cluster_labels > 0), all(cluster_labels %% 1 == 0),
          max(cluster_labels) == length(unique(cluster_labels)))

# form the ellipses
k <- max(cluster_labels)

dat_subset <- dat[which(cluster_labels == x),]

scaling_val <- .binary_search(dat_subset, evaluation_func = .ellipse_coverage,
                              target_val = target_percentage, verbose = T)

mean_vec <- colMeans(dat)
cov_mat <- stats::cov(dat)

.ellipse_points(mean_vec, scaling_val*cov_mat)
