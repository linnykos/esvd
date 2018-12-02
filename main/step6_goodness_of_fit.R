load("../results/step4_clustering.RData")

est_sd <- stats::sd(as.numeric(dat_impute))
trials <- 50
n <- nrow(dat_impute); d <- ncol(dat_impute)
svd_naive <- svd(dat_impute)
naive_pred_mat <- svd_naive$u[,1:5] %*% diag(svd_naive$d[1:5]) %*% t(svd_naive$v[,1:5])
our_pred_mat <- 4/res$u_mat %*% t(res$v_mat)
sum(svd_naive$d[1:k])/sum(svd_naive$d)

set.seed(10)
naive_simulated_dat <- matrix(NA, n, d)
for(i in 1:n){
 for(j in 1:d){
   naive_simulated_dat[i,j] <- stats::rnorm(1, naive_pred_mat[i,j], sd = est_sd)
 }
}
set.seed(10)
our_simulated_dat <- matrix(NA, n, d)
for(i in 1:n){
  for(j in 1:d){
    our_simulated_dat[i,j] <- stats::rnorm(1, our_pred_mat[i,j], sd = our_pred_mat[i,j]/2)
  }
}

dat_impute_vec <- as.numeric(dat_impute)
dat_impute_vec <- (dat_impute_vec-min(dat_impute_vec))/diff(range(dat_impute_vec))
naive_simulated_vec <- as.numeric(naive_simulated_dat)
naive_simulated_vec <- (naive_simulated_vec-min(naive_simulated_vec))/diff(range(naive_simulated_vec))
our_simulated_vec <- as.numeric(our_simulated_dat)
our_simulated_vec <- (our_simulated_vec-min(our_simulated_vec))/diff(range(our_simulated_vec))

png("../figure/main/distribution.png", height = 1200, width = 2250, res = 300, units = "px")
par(mfrow = c(1,2))
plot(sort(dat_impute_vec), sort(naive_simulated_vec), asp = T, pch = 16,
     col = rgb(0,0,0,0.1), xlab = "Quantiles of observed data",
     ylab = "Quantiles of simulated data",
     main = "QQ plot for fixed variance model")
lines(c(0,1), c(0,1), col = "red", lwd = 2, lty = 2)
plot(sort(dat_impute_vec), sort(our_simulated_vec), asp = T, pch = 16,
     col = rgb(0,0,0,0.1), xlab = "Quantiles of observed data",
     ylab = "Quantiles of simulated data",
     main = "QQ plot for curved Gaussian model")
lines(c(0,1), c(0,1), col = "red", lwd = 2, lty = 2)
graphics.off()

# naive_vec <- sapply(1:trials, function(x){
#   if(x %% floor(trials/10) == 0) cat('*')
#   set.seed(10*x)
#   simulated_dat <- matrix(NA, n, d)
#   for(i in 1:n){
#     for(j in 1:d){
#       simulated_dat[i,j] <- stats::rnorm(1, naive_pred_mat[i,j], sd = est_sd)
#     }
#   }
#
#   # stats::cor(as.numeric(dat_impute_log), as.numeric(simulated_dat))
#   ks.test(as.numeric(dat_impute_log), as.numeric(simulated_dat))$statistic
# })
#
# our_pred_mat <- res$u_mat[,1:k] %*% t(res$v_mat[,1:k])
# our_pred_mat <- 4/pred_mat
# our_vec <- sapply(1:trials, function(x){
#   if(x %% floor(trials/10) == 0) cat('*')
#   set.seed(10*x)
#   simulated_dat <- matrix(NA, n, d)
#   for(i in 1:n){
#     for(j in 1:d){
#       simulated_dat[i,j] <- stats::rnorm(1, our_pred_mat[i,j], sd = our_pred_mat[i,j]/2)
#     }
#   }
#
#   ks.test(as.numeric(dat_impute_log), as.numeric(simulated_dat))$statistic
# })

############################

load("../results/step3_factorization_logged_tmp.RData")
dat_vec <- as.numeric(dat_impute)
dat_vec <- (dat_vec-min(dat_vec))/diff(range(dat_vec))

res_list <- res_list[which(sapply(res_list, length) > 0)]
for(k in 1:length(res_list)){
  set.seed(10)
  print(k)
  our_pred_mat <- 4/(res_list[[k]]$u_mat %*% t(res_list[[k]]$v_mat))
  n <- nrow(res_list[[k]]$u_mat)
  d <- nrow(res_list[[k]]$v_mat)

  our_simulated_dat <- matrix(NA, n, d)
  for(i in 1:n){
    for(j in 1:d){
      our_simulated_dat[i,j] <- stats::rnorm(1, our_pred_mat[i,j], sd = our_pred_mat[i,j]/scalar_vec[k])
    }
  }

  our_simulated_vec <- as.numeric(our_simulated_dat)
  our_simulated_vec <- (our_simulated_vec-min(our_simulated_vec))/diff(range(our_simulated_vec))

  png(paste0("../figure/main/distribution_", k, ".png"), height = 2000, width = 2000, res = 300, units = "px")
  plot(sort(dat_vec), sort(our_simulated_vec), asp = T, pch = 16,
       col = rgb(0,0,0,0.1), main = k)
  lines(c(0,1), c(0,1), col = "red", lwd = 2)
  graphics.off()
}

