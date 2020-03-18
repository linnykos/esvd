rm(list=ls())
library(simulation)
library(eSVD)
source("../simulation/factorization_generator.R")

paramMat <- cbind(50, 120, 10,
                  rep(rep(1:3, each = 4), times = 4), 50, 2, 50, 10,
                  rep(1:4, each = 12),
                  rep(c(1/27, 1/800, 1/250, 1/1000), each = 12),
                  rep(1:4, times = 12),
                  rep(c(1, 1, NA, NA), times = 12),
                  rep(c(3000, rep(100, 3)), times = 4))
colnames(paramMat) <- c("n_each", "d_each", "sigma",
                        "k", "true_r",  "true_scalar", "max_iter", "fitting_iter",
                        "true_distr",
                        "modifier",
                        "fitting_distr",
                        "fitting_param",
                        "max_val")
trials <- 5
ncores <- 20

################

cell_pop <- matrix(c(4,10, 25,100, 60,80, 25,100,
                     40,10, 60,80, 60,80, 100, 25),
                   nrow = 4, ncol = 4, byrow = T)
gene_pop <- matrix(c(20,90, 25,100,
                     90,20, 100,25)/20, nrow = 2, ncol = 4, byrow = T)

rule <- function(vec){
  n_each <- vec["n_each"]
  d_each <- vec["d_each"]
  sigma <- vec["sigma"]
  total <- vec["total"]
  modifier <- vec["modifier"]

  res <- generate_natural_mat(cell_pop, gene_pop, n_each, d_each, sigma, modifier)
  # plot(res$cell_mat[,1], res$cell_mat[,2], asp = T, col = rep(1:4, each = n_each), pch = 16)
  nat_mat <- res$nat_mat

  if(vec["true_distr"] == 1){
    obs_mat <- round(generator_gaussian(nat_mat))
  } else if(vec["true_distr"] == 2){
    obs_mat <- generator_esvd_poisson(nat_mat)
  } else if(vec["true_distr"] == 3 ){
    obs_mat <- generator_esvd_nb(nat_mat, scalar = vec["true_r"])
  } else {
    obs_mat <- round(generator_curved_gaussian(nat_mat, scalar = vec["true_scalar"]))
  }

  list(dat = obs_mat, truth = res$cell_mat, nat_mat = nat_mat)
}

criterion <- function(dat, vec, y){
  set.seed(10*y)

  dat_obs <- dat$dat
  missing_idx <- eSVD::construct_missing_values(n = nrow(dat_obs), p = ncol(dat_obs), num_val = 2)
  dat_NA <- dat_obs
  dat_NA[missing_idx] <- NA

  missing_val <- dat_obs[missing_idx]

  set.seed(10)
  # fixed variance gaussian
  if(vec["fitting_distr"] == 1){
    init <- eSVD::initialization(dat_NA, family = "gaussian", k = vec["k"], max_val = vec["max_val"])
    fit <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                   family = "gaussian",
                                   max_iter = vec["max_iter"], max_val = vec["max_val"],
                                   return_path = F, cores = ncores,
                                   verbose = F)

    pred_mat <- fit$u_mat %*% t(fit$v_mat)
    pred_nat <- pred_mat[missing_idx]
    true_nat <- dat$nat_mat[missing_idx]
    pred_val <- eSVD::compute_mean(pred_mat, family = "gaussian")[missing_idx]
    expected_val <- eSVD::compute_mean(dat$nat_mat, family = "gaussian")[missing_idx]
    fitting_param <- vec["fitting_param"]
    fitting_vec <- NA

  # poisson
  } else if(vec["fitting_distr"] == 2){
    init <- eSVD::initialization(dat_NA, family = "poisson", k = vec["k"], max_val = vec["max_val"])
    fit <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                   family = "poisson",
                                   max_iter = vec["max_iter"], max_val = vec["max_val"],
                                   return_path = F, cores = ncores,
                                   verbose = F)

    pred_mat <- fit$u_mat %*% t(fit$v_mat)
    pred_nat <- pred_mat[missing_idx]
    true_nat <- dat$nat_mat[missing_idx]
    pred_val <- eSVD::compute_mean(pred_mat, family = "poisson")[missing_idx]
    expected_val <- eSVD::compute_mean(dat$nat_mat, family = "gaussian")[missing_idx]
    fitting_param <- vec["fitting_param"]
    fitting_vec <- NA

  # negative binomial
  } else if(vec["fitting_distr"] == 3){
    if(is.na(vec["fitting_param"])){
      fitting_vec <- eSVD::tuning_scalar(dat_obs, family = "neg_binom",
                                         max_iter = vec["max_iter"], max_val = vec["max_val"], k = vec["k"],
                                         return_path = F, cores = ncores, iter_max = 10,
                                         search_min = 1, search_max = 2*max(dat_obs))
      fitting_param <- fitting_vec[length(fitting_vec)]
    } else {
      fitting_param <- vec["fitting_param"]
      fitting_vec <- NA
    }

    init <- eSVD::initialization(dat_NA, family = "neg_binom", k = vec["k"], max_val = vec["max_val"],
                                 scalar = fitting_param)
    fit <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                   family = "neg_binom", scalar = fitting_param,
                                   max_iter = vec["max_iter"], max_val = vec["max_val"],
                                   return_path = F, cores = ncores,
                                   verbose = F)

    pred_mat <- fit$u_mat %*% t(fit$v_mat)
    pred_nat <- pred_mat[missing_idx]
    true_nat <- -dat$nat_mat[missing_idx]
    pred_val <- eSVD::compute_mean(pred_mat, family = "neg_binom", scalar = fitting_param)[missing_idx]
    expected_val <- eSVD::compute_mean(-dat$nat_mat, family = "neg_binom", scalar = vec["true_r"])[missing_idx]

  # curved gaussian
  } else {
    if(is.na(vec["fitting_param"])){
      fitting_vec <- eSVD::tuning_scalar(dat_obs, family = "curved_gaussian",
                                         max_iter = vec["max_iter"], max_val = vec["max_val"], k = vec["k"],
                                         return_path = F, cores = ncores, iter_max = 10,
                                         search_min = 0.5, search_max = 10)
      fitting_param <- fitting_vec[length(fitting_vec)]
    } else {
      fitting_param <- vec["fitting_param"]
      fitting_vec <- NA
    }

    init <- eSVD::initialization(dat_NA, family = "curved_gaussian", k = vec["k"], max_val = vec["max_val"],
                                 scalar = fitting_param)
    fit <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                   family = "curved_gaussian", scalar = fitting_param,
                                   max_iter = vec["max_iter"], max_val = vec["max_val"],
                                   return_path = F, cores = ncores,
                                   verbose = F)

    pred_mat <- fit$u_mat %*% t(fit$v_mat)
    pred_nat <- pred_mat[missing_idx]
    true_nat <- dat$nat_mat[missing_idx]
    pred_val <- eSVD::compute_mean(pred_mat, family = "curved_gaussian", scalar = fitting_param)[missing_idx]
    expected_val <- eSVD::compute_mean(dat$nat_mat, family = "curved_gaussian", scalar = fitting_param)[missing_idx]

  }

  list(fit = fit$u_mat, truth = dat$truth,
       pred_nat = pred_nat, true_nat = true_nat,
       pred_val = pred_val, missing_val = missing_val,
       expected_val = expected_val, fitting_param = fitting_param,
       fitting_vec = fitting_vec)
}

## i <- 9; y <- 20; dat <- rule(paramMat[i,]); quantile(dat$dat); plot(dat$truth[,1], dat$truth[,2], asp = T, col = rep(1:4, each = paramMat[i,"n_each"]), pch = 16)
## i <- 31; y <- 1; zz <- criterion(rule(paramMat[i,]), paramMat[i,], y); zz
## neg_binom: i = 31; curved_gaussian: i = 48
############

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = trials,
                                        cores = NA, as_list = T,
                                        filepath = "../results/factorization_esvd_tmp.RData",
                                        verbose = T)

save.image("../results/factorization_esvd.RData")
