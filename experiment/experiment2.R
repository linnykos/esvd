rm(list=ls())
library(simulation)
library(eSVD)
source("../simulation/factorization_generator.R")
source("../simulation/factorization_methods.R")

paramMat <- cbind(50, 120, 10,
                  2, 50, 1/250, 1000,
                  80, 120, 600,
                  1/4, 1/4, 1/2,
                  1:6)
colnames(paramMat) <- c("n_each", "d_each", "sigma",
                        "k", "max_iter", "modifier", "max_val",
                        "size_1", "size_2", "size_3",
                        "prop_1", "prop_2", "prop_3",
                        "method")

trials <- 100
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

set.seed(10)
vec <- paramMat[1,]
dat <- rule(vec)

dat <- dat$dat

paramMat_esvd <- matrix(c(50, 100, 500, 1000), nrow = 4, ncol = 1)
colnames(paramMat_esvd) <- c("scalar")
paramMat <- paramMat_esvd

k <- 3
set.seed(10)
missing_idx <- eSVD::construct_missing_values(n = nrow(dat), p = ncol(dat), num_val = 2)
dat_NA <- dat
dat_NA[missing_idx] <- NA

fit_list <- lapply(1:nrow(paramMat), function(i){
  set.seed(10)
  init <- eSVD::initialization(dat_NA, family = "neg_binom", k = k, max_val = 2000,
                               scalar = paramMat[i, "scalar"])
  eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                          family = "neg_binom", scalar = paramMat[i, "scalar"],
                          max_iter = 50, max_val = 2000,
                          return_path = F, cores = ncores,
                          verbose = F)
})

quality_vec <- sapply(1:nrow(paramMat), function(i){
  nat_mat <- fit_list[[i]]$u_mat %*% t(fit_list[[i]]$v_mat)
  mean_mat <- eSVD::compute_mean(nat_mat, family = "neg_binom", scalar = paramMat[i, "scalar"])
  eSVD::plot_prediction_against_observed(dat, nat_mat_list = list(nat_mat),
                                         scalar = paramMat[i, "scalar"],
                                         family = "neg_binom", missing_idx_list = list(missing_idx),
                                         plot = F)
})

idx <- which.min(abs(quality_vec - 45))

scalar <- paramMat[idx, "scalar"]

set.seed(10)
init <- eSVD::initialization(dat, family = "neg_binom", k = k, max_val = 2000,
                             scalar = scalar)
fit <- eSVD::fit_factorization(dat, u_mat = init$u_mat, v_mat = init$v_mat,
                               family = "neg_binom", scalar = scalar,
                               max_iter = 50, max_val = 2000,
                               return_path = F, cores = ncores, verbose = F)

save.image("../experiment/experiment.RData")
