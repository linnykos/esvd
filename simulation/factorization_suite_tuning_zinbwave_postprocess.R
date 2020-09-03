rm(list=ls())
load("../results/factorization_results_tuning_zinbwave.RData")

# see if there's an empirical difference
x <- 1; i <- 5
plot(res[[i]][[x]]$fit$u_mat[,1], res[[i]][[x]]$fit$u_mat[,2], asp = T, pch = 16, col = rep(1:4, each = 50))

# evaluate the quality of each fit
angle_fit <- lapply(1:trials, function(x){
  print(x)
  nat_mat_list_list <- lapply(1:nrow(paramMat), function(i){
    u_mat <- res[[i]][[x]]$fit$u_mat
    v_mat <- res[[i]][[x]]$fit$v_mat
    list(u_mat %*% t(v_mat))
  })


  eSVD::tuning_select_scalar(dat = res[[1]][[x]]$dat,
                             nat_mat_list_list = nat_mat_list_list,
                             family = "neg_binom", compute_percentage = F,
                             missing_idx_list = list(res[[1]][[x]]$missing_idx),
                             scalar_vec = rep(c(50, 100, 500, 1000), times = 3))$all_results
})

# evaluate the accuracy of the fits



