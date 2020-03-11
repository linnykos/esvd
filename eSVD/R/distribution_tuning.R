tuning <- function(dat, family, iter_max = 5, binary_search_min = 1,
                   binary_search_max = 2000, binary_search_iter = 15,
                   binary_search_tol = 1e-3, ...){
  # set class
  if(length(class(dat)) == 1) class(dat) <- c(family, class(dat)[length(class(dat))])

  # fit initial fit
  family_initial <- ifelse(family == "neg_binom", "poisson", "exponential")
  fit <- .tuning_fit(dat, family = family_initial, scalar = NA, ...)

  # determine initial param
  scalar_vec <- rep(NA, iter_max)
  scalar_vec[1] <- .tuning_param_binary_search(dat, fit$u_mat, fit$v_mat, family = family_initial, initial = T,
                                        binary_search_min = binary_search_min,
                                        binary_search_max = binary_search_max,
                                        binary_search_iter = binary_search_iter,
                                        binary_search_tol = binary_search_tol, ...)

  # iterate between fitting and parameter estimation
  for(i in 2:iter_max){
    fit <- .tuning_fit(dat, family = family, scalar = scalar_vec[i-1], ...)
    scalar_vec[i] <- .tuning_param_binary_search(dat, fit$u_mat, fit$v_mat, family = family, initial = T,
                                                 binary_search_min = binary_search_min,
                                                 binary_search_max = binary_search_max,
                                                 binary_search_iter = binary_search_iter,
                                                 binary_search_tol = binary_search_tol, ...)
  }

  scalar_vec
}

################

.tuning_fit <- function(dat, family, scalar, ...){
  init <- initialization(dat, family = family, scalar = scalar, ...)
  fit_factorization(dat, u_mat = init$u_mat, v_mat = init$v_mat,
                          family = family, scalar = scalar, ...)
}


###########

.tuning_param_binary_search <- function(dat, u_mat, v_mat, family, initial = F, binary_search_min = 1,
                                        binary_search_max = 2000, binary_search_iter = 15,
                                        binary_search_tol = 1e-3, ...){
  stopifnot(ncol(u_mat) == ncol(v_mat), nrow(dat) == nrow(u_mat), ncol(dat) == nrow(v_mat))
  k <- ncol(u_mat); n <- nrow(dat); p <- ncol(dat)
  df_val <- n*p - (n*k + p*k)
  stopifnot(df_val > 0)

  nat_mat <- u_mat %*% t(v_mat)
  mean_mat_tmp <- compute_mean(nat_mat, family, scalar = 1)
  recompute_mean <- (family == "neg" & intial)

  lo_val <- binary_search_min
  obj_lo <- .compute_tuning_objective(dat, family, nat_mat, mean_mat_tmp, scalar = lo_val, recompute_mean = recompute_mean)

  hi_val <- binary_search_max
  obj_hi <- .compute_tuning_objective(dat, family, nat_mat, mean_mat_tmp, scalar = binary_search_max, recompute_mean = recompute_mean)

  obj_prev <- Inf; obj_current <- Inf

  while(iter <= binary_search_iter){
    mid_val <- (lo_val + hi_val)/2
    obj_mid <- .compute_tuning_objective(dat, family, nat_mat, mean_mat_tmp, scalar = mid_val, recompute_mean = recompute_mean)

    if(sign(obj_lo - df_val) != sign(obj_mid - df_val)){
      hi_val <- mid_val
    } else {
      lo_val <- mid_val
    }

    obj_current <- obj_mid
    if(abs(obj_current - obj_prev) <= binary_search_tol) break()

    iter <- iter + 1
  }

  mid_val
}

.compute_tuning_objective <- function(dat, family, nat_mat, mean_mat, scalar, recompute_mean){
  mean_mat <- .recompute_mean(nat_mat, mean_mat, scalar = scalar, recompute_mean = recompute_mean)
  var_mat <- .compute_variance(mean_mat, family, scalar = scalar)
  sum((dat-mean_mat)^2/var_mat)
}

.recompute_mean <- function(nat_mat, mean_mat = NA, scalar = NA, recompute_mean = F){
  if(!recompute_mean){
    stopifnot(all(!is.na(mean_mat)))
    return(mean_mat)
  } else {
    # the nat_mat is for the neg_binom family then
    compute_mean(nat_mat, "neg_binom", scalar = scalar)
  }
}

.compute_variance <- function(mean_mat, family, scalar){
  if(family == "neg_binom"){
    mean_mat + mean_mat^2/scalar
  } else {
    # this means it's curved gaussian
    (mean_mat/scalar)^2
  }
}

##################

# .tuning_param_initial.neg_binom <- function(dat, ...){
#
# }
#
# .tuning_param_initial.curved_gaussian <- function(dat, ...){
#
# }
#
# .tuning_parameter <- function(dat, u_mat, v_mat, family_to, family_from){
#   if(family_to == "neg_binom" & family_from == "neg_binom"){
#     .tuning_neg_binom(dat, u_mat, v_mat)
#   } else if(family_to == "neg_binom" & family_from == "poisson"){
#     .tuning_neg_binom_from_poisson(dat, u_mat, v_mat)
#   } else if(family_to == "curved_gaussian" & family_from == "curved_gaussian"){
#     .tuning_curved_gaussian(dat, u_mat, v_mat)
#   } else {
#     .tuning_curved_gaussian_from_exponential(dat, u_mat, v_mat)
#   }
# }
#
# .tuning_neg_binom <- function(dat, u_mat, v_mat){
#   p_mat <- exp(u_mat %*% t(v_mat))
#
#   r_seq <- 1:100
#   proposed_val <- t(sapply(r_seq, function(x){
#     pred_mat <- x*p_mat/(1-p_mat)
#     c(sum((dat - pred_mat)^2)/prod(dim(dat)),
#       sum((pred_mat + pred_mat^2/x))/prod(dim(dat)))
#   }))
#
#   r_seq[which.min(abs(proposed_val[,1] - proposed_val[,2]))]
# }
#
# .tuning_neg_binom_from_poisson <- function(dat, u_mat, v_mat){
#   pred_mat <- exp(u_mat %*% t(v_mat))
#
#   target_val <- sum((dat - pred_mat)^2)/prod(dim(dat))
#   r_seq <- 1:100
#   proposed_val <- t(sapply(r_seq, function(x){
#     sum((pred_mat + pred_mat^2/x))/prod(dim(dat))
#   }))
#
#   r_seq[which.min(abs(target_val - proposed_val))]
# }
#
# .tuning_curved_gaussian <- function(dat, u_mat, v_mat){
#   pred_mat <- 1/(u_mat %*% t(v_mat))
#
#   target_val <- sum((dat - pred_mat)^2)/prod(dim(dat))
#   r_seq <- seq(1, 10, length.out = 101)
#   proposed_val <- sapply(r_seq, function(x){sum((pred_mat/x)^2)/prod(dim(dat))})
#
#   r_seq[which.min(abs(target_val - proposed_val))]
# }
#
# .tuning_curved_gaussian_from_exponential <- function(dat, u_mat, v_mat){
#   pred_mat <- 1/(-u_mat %*% t(v_mat))
#
#   target_val <- sum((dat - pred_mat)^2)/prod(dim(dat))
#   r_seq <- seq(1, 10, length.out = 101)
#   proposed_val <- sapply(r_seq, function(x){sum((pred_mat/x)^2)/prod(dim(dat))})
#
#   r_seq[which.min(abs(target_val - proposed_val))]
# }
