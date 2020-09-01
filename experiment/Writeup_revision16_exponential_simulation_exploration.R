####################################################
# FIRST SUITE: See if when k is correct, the nuisance parameter can be correct

# for curved gaussian
rm(list=ls())
load("../results/factorization_results_exponential_families_3.RData")
param_idx <- 6

# aggregate results across all 50 trials
summary_list <- lapply(1:length(res[[param_idx]]), function(trial){
  nat_mat_list_list <- lapply(1:3, function(i){
    lapply(1, function(j){
      u_mat <- res[[param_idx]][[trial]]$fit[[i]]$u_mat
      v_mat <- res[[param_idx]][[trial]]$fit[[i]]$v_mat
      u_mat %*% t(v_mat)
    })
  })

  eSVD::tuning_select_scalar(dat = res[[param_idx]][[trial]]$dat,
                             nat_mat_list_list = nat_mat_list_list,
                             family = "curved_gaussian",
                             missing_idx_list = list(res[[param_idx]][[trial]]$missing_idx),
                             scalar_vec = alpha_vec)
})
summary_list
table(sapply(summary_list, function(x){x$idx}))

# for negative binomial
rm(list=ls())
load("../results/factorization_results_exponential_families_2.RData")
param_idx <- 5

# aggregate results across all 50 trials
summary_list <- lapply(1:length(res[[param_idx]]), function(trial){
  nat_mat_list_list <- lapply(1:3, function(i){
    lapply(1, function(j){
      u_mat <- res[[param_idx]][[trial]]$fit[[i]]$u_mat
      v_mat <- res[[param_idx]][[trial]]$fit[[i]]$v_mat
      u_mat %*% t(v_mat)
    })
  })

  eSVD::tuning_select_scalar(dat = res[[param_idx]][[trial]]$dat,
                             nat_mat_list_list = nat_mat_list_list,
                             family = "neg_binom",
                             missing_idx_list = list(res[[param_idx]][[trial]]$missing_idx),
                             scalar_vec = r_vec)
})
summary_list
table(sapply(summary_list, function(x){x$idx}))

####################################################
# SECOND SUITE: See if the distribution is correct, can pick the correct k

rm(list=ls())
load("../results/factorization_results_exponential_families_3.RData")
param_vec <- c(3,6,9)

# aggregate results across all 50 trials
summary_list <- lapply(1:length(res[[param_vec[1]]]), function(trial){
  print(trial)
  tmp <- lapply(param_vec, function(param_idx){
    lapply(1:3, function(i){
      u_mat <- res[[param_idx]][[trial]]$fit[[i]]$u_mat
      v_mat <- res[[param_idx]][[trial]]$fit[[i]]$v_mat
      list(u_mat %*% t(v_mat))
    })
  })
  nat_mat_list_list <- eSVD:::.flatten_list(tmp)
  for(i in 1:length(nat_mat_list_list)){
    nat_mat_list_list[[i]] <- list(nat_mat_list_list[[i]])
  }

  eSVD::tuning_select_scalar(dat = res[[param_vec[1]]][[trial]]$dat,
                             nat_mat_list_list = nat_mat_list_list,
                             family = "curved_gaussian", compute_percentage = T,
                             missing_idx_list = list(res[[param_vec[1]]][[trial]]$missing_idx),
                             scalar_vec = rep(alpha_vec, times = 3))
})


##################################
# THIRD SUITE: See if method can pick among the correct distribution
rm(list=ls())
load("../results/factorization_results_exponential_families_3.RData")

trial <- 2

# first do curved gaussian
param_vec <- c(3,6,9)
tmp <- lapply(param_vec, function(param_idx){
  lapply(1:3, function(i){
    u_mat <- res[[param_idx]][[trial]]$fit[[i]]$u_mat
    v_mat <- res[[param_idx]][[trial]]$fit[[i]]$v_mat
    list(u_mat %*% t(v_mat))
  })
})
nat_mat_list_list <- eSVD:::.flatten_list(tmp)
for(i in 1:length(nat_mat_list_list)){
  nat_mat_list_list[[i]] <- list(nat_mat_list_list[[i]])
}

diagnostic_1 <- eSVD::tuning_select_scalar(dat = res[[param_vec[1]]][[trial]]$dat,
                           nat_mat_list_list = nat_mat_list_list,
                           family = "curved_gaussian", compute_percentage = T,
                           missing_idx_list = list(res[[param_vec[1]]][[trial]]$missing_idx),
                           scalar_vec = rep(alpha_vec, times = 3))

###

param_vec <- c(2,5,8)
tmp <- lapply(param_vec, function(param_idx){
  lapply(1:3, function(i){
    u_mat <- res[[param_idx]][[trial]]$fit[[i]]$u_mat
    v_mat <- res[[param_idx]][[trial]]$fit[[i]]$v_mat
    list(u_mat %*% t(v_mat))
  })
})
nat_mat_list_list <- eSVD:::.flatten_list(tmp)
for(i in 1:length(nat_mat_list_list)){
  nat_mat_list_list[[i]] <- list(nat_mat_list_list[[i]])
}

diagnostic_2 <- eSVD::tuning_select_scalar(dat = res[[param_vec[1]]][[trial]]$dat,
                                           nat_mat_list_list = nat_mat_list_list,
                                           family = "neg_binom", compute_percentage = T,
                                           missing_idx_list = list(res[[param_vec[1]]][[trial]]$missing_idx),
                                           scalar_vec = rep(r_vec, times = 3))

###

param_vec <- c(1,4,7)
tmp <- lapply(param_vec, function(param_idx){
  lapply(1, function(i){
    u_mat <- res[[param_idx]][[trial]]$fit[[i]]$u_mat
    v_mat <- res[[param_idx]][[trial]]$fit[[i]]$v_mat
    list(u_mat %*% t(v_mat))
  })
})
nat_mat_list_list <- eSVD:::.flatten_list(tmp)
for(i in 1:length(nat_mat_list_list)){
  nat_mat_list_list[[i]] <- list(nat_mat_list_list[[i]])
}

diagnostic_3 <- eSVD::tuning_select_scalar(dat = res[[param_vec[1]]][[trial]]$dat,
                                           nat_mat_list_list = nat_mat_list_list,
                                           family = "poisson", compute_percentage = T,
                                           missing_idx_list = list(res[[param_vec[1]]][[trial]]$missing_idx),
                                           scalar_vec = rep(NA, 3))

###

diagnostic_1
diagnostic_2
diagnostic_3


