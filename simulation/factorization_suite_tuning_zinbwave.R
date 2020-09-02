rm(list=ls())
library(simulation)
library(eSVD)
source("../simulation/factorization_generator.R")

session_info <- sessionInfo()
date_of_run <- Sys.time()
source_code_info <- readLines("../simulation/factorization_generator.R")
source_code_info <- c(source_code_info, readLines("../simulation/factorization_methods.R"))
source_code_info <- c(source_code_info, readLines("../simulation/factorization_suite_tuning_zinbwave.R"))

paramMat <- cbind(50, 200, 5,
                  rep(c(2,3,10), each = 4),
                  rep(c(50, 100, 500, 1000), times = 3),
                  50, 1/250, 1000,
                  80, 120, 600,
                  1/4, 1/4, 1/2)
colnames(paramMat) <- c("n_each", "d_each", "sigma",
                        "k",
                        "r_val",
                        "max_iter", "modifier", "max_val",
                        "size_1", "size_2", "size_3",
                        "prop_1", "prop_2", "prop_3")

trials <- 50
ncores <- 15
doMC::registerDoMC(cores = ncores)

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
  nat_mat <- res$nat_mat

  r_vec <- sample(c(paramMat[1,"size_1"], paramMat[1,"size_2"], paramMat[1,"size_3"]),
                  size = ncol(nat_mat),
                  prob = c(paramMat[1,"prop_1"], paramMat[1,"prop_2"], paramMat[1,"prop_3"]),
                  replace = T)

  dat <- generator_zinb_nb(nat_mat, r_vec)
  obs_mat <- round(dat$dat * 1000/max(dat$dat))

  list(dat = obs_mat, truth = res$cell_mat)
}

criterion <- function(dat, vec, y){
  dat_obs <- dat$dat

  set.seed(10*y)
  missing_idx <- eSVD::construct_missing_values(n = nrow(dat_obs), p = ncol(dat_obs), num_val = 2)
  dat_NA <- dat_obs
  dat_NA[missing_idx] <- NA

  missing_val <- dat_obs[missing_idx]
  init <- eSVD::initialization(dat_NA, family = "neg_binom", k = vec["k"], max_val = vec["max_val"],
                               scalar = vec["r_val"])
  fit <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                          family = "neg_binom", scalar = vec["r_val"],
                          max_iter = vec["max_iter"], max_val = vec["max_val"],
                          return_path = F, cores = NA,
                          verbose = F)

  list(fit = fit, true_u_mat = dat$u_mat, true_v_mat = dat$v_mat,
       dat = dat_obs, missing_idx = missing_idx)
}

## i <- 1; y <- 1; set.seed(y); zz1 <- criterion(rule(paramMat[i,]), paramMat[i,], y); head(zz1$fit$fit$u_mat); head(zz1$truth)
## i <- 2; y <- 2; set.seed(y); zz2 <- criterion(rule(paramMat[i,]), paramMat[i,], y); head(zz2$fit$fit$u_mat); head(zz2$truth)
## i <- 1; y <- 1; set.seed(y); y <- rule(paramMat[i,]); y$dat[1:5,1:5]; head(y$truth)

## i <- 2; y <- 1; set.seed(y); zz3 <- criterion(rule(paramMat[i,]), paramMat[i,], y); head(zz3$truth)
## i <- 2; y <- 2; set.seed(y); zz4 <- criterion(rule(paramMat[i,]), paramMat[i,], y); head(zz4$truth)

## i <- 2; y <- 2; set.seed(y); tmp <- rule(paramMat[i,]);tmp$dat[1:5,1:5]

############

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = trials,
                                        cores = ncores, as_list = T,
                                        filepath = "../results/factorization_results_tuning_zinbwave_tmp.RData",
                                        verbose = T)

save.image("../results/factorization_results_tuning_zinbwave.RData")
