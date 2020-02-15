rm(list=ls())
library(simulation)
library(eSVD)
source("../simulation/factorization_generator.R")

paramMat <- cbind(50, 120, 5,
                  2, 2, 50, 50,
                  rep(1:4, each = 8),
                  rep(c(1, 1/400, 1/350, 1/1000), each = 8),
                  rep(c(1,2, rep(3,3), rep(4,3)), times = 4),
                  rep(c(1,1, c(25, 50, 200), c(1,2,4)), times = 4),
                  rep(c(3000, rep(100, 7)), times = 4))
colnames(paramMat) <- c("n_each", "d_each", "sigma",
                        "k", "true_scalar", "true_r", "max_iter",
                        "true_distr",
                        "modifier",
                        "fitting_distr",
                        "fitting_param",
                        "max_val")
trials <- 20
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
    obs_mat <- generator_esvd_nb(nat_mat, size = vec["true_r"])
  } else {
    obs_mat <- round(generator_curved_gaussian(nat_mat, scalar = vec["true_scalar"]))
  }

  list(dat = obs_mat, truth = res$cell_mat, nat_mat = nat_mat)
}

criterion <- function(dat, vec, y){
  set.seed(10*y)

  dat_obs <- dat$dat
  n <- nrow(dat_obs); d <- ncol(dat_obs)
  missing_idx <- rbind(do.call(rbind, (lapply(1:n, function(x){
    cbind(x, sample(1:d, 2))
  }))), do.call(rbind, (lapply(1:d, function(x){
    cbind(sample(1:n, 2), d)
  }))))

  dat_NA <- dat_obs
  for(i in 1:nrow(missing_idx)){
    dat_NA[missing_idx[i,1], missing_idx[i,2]] <- NA
  }
  missing_idx <- which(is.na(dat_NA))
  missing_val <- dat_obs[missing_idx]

  set.seed(10)
  if(vec["fitting_distr"] == 1){
    init <- eSVD::initialization(dat_NA, family = "gaussian", k = vec["k"], max_val = vec["max_val"])
    fit <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                   family = "gaussian",
                                   max_iter = vec["max_iter"], max_val = vec["max_val"],
                                   return_path = F, cores = ncores,
                                   verbose = F)


    pred_mat <- fit$u_mat %*% t(fit$v_mat)
    pred_val <- pred_mat[missing_idx]
    expected_val <- dat$nat_mat[missing_idx]

  } else if(vec["fitting_distr"] == 2){
    init <- eSVD::initialization(dat_NA, family = "poisson", k = vec["k"], max_val = vec["max_val"])
    fit <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                   family = "poisson",
                                   max_iter = vec["max_iter"], max_val = vec["max_val"],
                                   return_path = F, cores = ncores,
                                   verbose = F)

    pred_mat <- fit$u_mat %*% t(fit$v_mat)
    pred_val <- exp(pred_mat[missing_idx])
    expected_val <- exp(dat$nat_mat[missing_idx])

  } else if(vec["fitting_distr"] == 3){
    init <- eSVD::initialization(dat_NA, family = "neg_binom", k = vec["k"], max_val = -vec["max_val"],
                                 size = vec["fitting_param"])
    fit <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                   family = "neg_binom", size = vec["fitting_param"],
                                   max_iter = vec["max_iter"], max_val = -vec["max_val"],
                                   return_path = F, cores = ncores,
                                   verbose = F)

    pred_mat <- fit$u_mat %*% t(fit$v_mat)
    pred_val <- (vec["fitting_param"]*exp(pred_mat)/(1-exp(pred_mat)))[missing_idx]
    expected_val <- (vec["true_r"]*exp(-dat$nat_mat)/(1-exp(-dat$nat_mat)))[missing_idx]

  } else {
    init <- eSVD::initialization(dat_NA, family = "curved_gaussian", k = vec["k"], max_val = vec["max_val"],
                                 scalar = vec["fitting_param"])
    fit <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                   family = "curved_gaussian", scalar = vec["fitting_param"],
                                   max_iter = vec["max_iter"], max_val = vec["max_val"],
                                   return_path = F, cores = ncores,
                                   verbose = F)

    pred_mat <- fit$u_mat %*% t(fit$v_mat)
    pred_val <- 1/pred_mat[missing_idx]
    expected_val <- 1/dat$nat_mat[missing_idx]
  }

  list(fit = fit$u_mat, truth = dat$truth, pred_val = pred_val, missing_val = missing_val,
       expected_val = expected_val)
}

## i <- 9; y <- 20; dat <- rule(paramMat[i,]); quantile(dat$dat); plot(dat$truth[,1], dat$truth[,2], asp = T, col = rep(1:4, each = paramMat[i,"n_each"]), pch = 16)
## i <- 19; y <- 1; zz <- criterion(rule(paramMat[i,]), paramMat[i,], y); zz

############

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = trials,
                                        cores = NA, as_list = T,
                                        filepath = "../results/factorization_esvd_tmp.RData",
                                        verbose = T)

save.image("../results/factorization_esvd.RData")
