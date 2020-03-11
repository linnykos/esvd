rm(list=ls())
load("../results/lingxue_analysis.RData")

# which dataset to investigate?
dat_num <- 7
label_vec <- preprocessing_list[[dat_num]]$label_vec
label_num <- as.numeric(label_vec)

plot(fit_all_list[[dat_num]]$gaussian$u_mat[,1], fit_all_list[[dat_num]]$gaussian$u_mat[,2], asp = T, col = c(1:5)[label_num],
     pch = 16, main = "Gaussian")
plot(fit_all_list[[dat_num]]$poisson$u_mat[,1], fit_all_list[[dat_num]]$poisson$u_mat[,2], asp = T, col = c(1:5)[label_num],
     pch = 16, main = "Poisson")
plot(fit_all_list[[dat_num]]$exponential$u_mat[,1], fit_all_list[[dat_num]]$exponential$u_mat[,2], asp = T, col = c(1:5)[label_num],
     pch = 16, main = "Exponential")
plot(fit_all_list[[dat_num]]$neg_binom$u_mat[,1], fit_all_list[[dat_num]]$neg_binom$u_mat[,2], asp = T, col = c(1:5)[label_num],
     pch = 16, main = "Negative binomial")
plot(fit_all_list[[dat_num]]$curved_gaussian$u_mat[,1], fit_all_list[[dat_num]]$curved_gaussian$u_mat[,2], asp = T, col = c(1:5)[label_num],
     pch = 16, main = "Curved Gaussian")
fit_all_list[[dat_num]]$neg_bin_param
fit_all_list[[dat_num]]$curved_gaussian_param
