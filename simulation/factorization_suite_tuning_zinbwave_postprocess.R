rm(list=ls())
load("../results/factorization_results_tuning_zinbwave.RData")

# see if there's an empirical difference
x <- 1; i <- 8
plot(res[[i]][[x]]$fit[[1]]$u_mat[,1], res[[i]][[x]]$fit[[1]]$u_mat[,2], asp = T, pch = 16, col = rep(1:4, each = 50))

# evaluate the quality of each fit
angle_fit <- lapply(1:trials, function(x){
  print(x)
  nat_mat_list_list <- lapply(1:nrow(paramMat), function(i){
    lapply(1:3, function(j){
      u_mat <- res[[i]][[x]]$fit[[j]]$u_mat
      v_mat <- res[[i]][[x]]$fit[[j]]$v_mat
      u_mat %*% t(v_mat)
    })
  })


  eSVD::tuning_select_scalar(dat = res[[1]][[x]]$dat,
                             nat_mat_list_list = nat_mat_list_list,
                             family = "neg_binom", compute_percentage = F,
                             missing_idx_list = res[[1]][[x]]$missing_idx,
                             scalar_vec = rep(c(50, 100, 500), times = 3))$all_results
})

loglik <- lapply(1:trials, function(x){
  print(x)
  nat_mat_list_list <- lapply(1:nrow(paramMat), function(i){
    lapply(1:3, function(j){
      u_mat <- res[[i]][[x]]$fit[[j]]$u_mat
      v_mat <- res[[i]][[x]]$fit[[j]]$v_mat
      u_mat %*% t(v_mat)
    })
  })

  loglik_vec <- sapply(1:nrow(paramMat), function(i){
    tmp_vec <- sapply(1:3, function(j){
      idx <- res[[1]][[x]]$missing_idx[[j]]
      dat_vec <- res[[1]][[x]]$dat[idx]
      nat_vec <- nat_mat_list_list[[i]][[j]][idx]
      r_val <- paramMat[i,"r_val"]

      sum(-r_val*log(1-exp(nat_vec)) - dat_vec*nat_vec)/length(idx)
    })

    mean(tmp_vec)
  })
})


# evaluate the accuracy of the fits



