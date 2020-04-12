rm(list=ls())
library(simulation)
library(eSVD)
source("../simulation/factorization_generator.R")

paramMat <- cbind(50, 120, 5,
                  rep(rep(1:3, each = 4), times = 4), 50, 2, 50,
                  rep(1:4, each = 12),
                  rep(c(1/27, 1/800, 1/250, 1/1000), each = 12),
                  rep(1:4, times = 12),
                  rep(c(1, 1, NA, NA), times = 12),
                  rep(c(3000, rep(100, 3)), times = 4))
colnames(paramMat) <- c("n_each", "d_each", "sigma",
                        "k", "true_r",  "true_alpha", "max_iter",
                        "true_distr",
                        "modifier",
                        "fitting_distr",
                        "fitting_param",
                        "max_val")
correct_idx <- c(5, 18, 31, 44)
rearrange_idx <- c(correct_idx, c(1:nrow(paramMat))[-correct_idx])
paramMat <- paramMat[rearrange_idx,]

trials <- 100
ncores <- 20
r_vec <- c(5, 50, 100)
alpha_vec <- c(0.5, 2, 50)

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
    obs_mat <- round(generator_curved_gaussian(nat_mat, scalar = vec["true_alpha"]))
  }

  obs_mat <- obs_mat * 1000/max(obs_mat)

  list(dat = obs_mat, u_mat = res$cell_mat, v_mat = res$gene_mat)
}

criterion <- function(dat, vec, y){
  dat_obs <- dat$dat

  set.seed(10*y)
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

    # poisson
  } else if(vec["fitting_distr"] == 2){
    init <- eSVD::initialization(dat_NA, family = "poisson", k = vec["k"], max_val = vec["max_val"])
    fit <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                   family = "poisson",
                                   max_iter = vec["max_iter"], max_val = vec["max_val"],
                                   return_path = F, cores = ncores,
                                   verbose = F)

    # negative binomial
  } else if(vec["fitting_distr"] == 3){
    fit <- lapply(r_vec, function(r_val){
      init <- eSVD::initialization(dat_NA, family = "neg_binom", k = vec["k"], max_val = vec["max_val"],
                                   scalar = r_val)
      eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                              family = "neg_binom", scalar = r_val,
                              max_iter = vec["max_iter"], max_val = vec["max_val"],
                              return_path = F, cores = ncores,
                              verbose = F)
    })

    # curved gaussian
  } else {
    fit <- lapply(alpha_vec, function(alpha_val){
      init <- eSVD::initialization(dat_NA, family = "curved_gaussian", k = vec["k"], max_val = vec["max_val"],
                                   scalar = alpha_val)
      eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                              family = "curved_gaussian", scalar = alpha_val,
                              max_iter = vec["max_iter"], max_val = vec["max_val"],
                              return_path = F, cores = ncores,
                              verbose = F)
    })
  }

  list(fit = fit, true_u_mat = dat$u_mat, true_v_mat = dat$v_mat,
       dat = dat_obs, missing_idx = missing_idx)
}

## i <- 2; y <- 1; dat <- rule(paramMat[i,]); quantile(dat$dat); plot(dat$truth[,1], dat$truth[,2], asp = T, col = rep(1:4, each = paramMat[i,"n_each"]), pch = 16)
## i <- 2; y <- 1; set.seed(y); zz <- criterion(rule(paramMat[i,]), paramMat[i,], y); zz
## scalar_vec <- r_vec; family_val <- "neg_binom"; quality_vec <- sapply(1:length(zz$fit), function(j){ nat_mat <- zz$fit[[j]]$u_mat %*% t(zz$fit[[j]]$v_mat); plot_prediction_against_observed(dat = zz$dat, nat_mat_list = list(nat_mat), family = family_val, missing_idx_list = list(zz$missing_idx), scalar = scalar_vec[j], plot = F)}); quality_vec
## neg_binom: i = 31; curved_gaussian: i = 44
############

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = trials,
                                        cores = NA, as_list = T,
                                        filepath = "../results/factorization_exponential_families_tmp.RData",
                                        verbose = T)

save.image("../results/factorization_exponential_families.RData")
