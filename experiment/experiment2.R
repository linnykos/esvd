rm(list=ls())
library(simulation)
library(eSVD)
source("../simulation/factorization_generator.R")

paramMat <- cbind(50, 120, 10,
                  2, 2, 50, 50,
                  rep(1:4, each = 9),
                  rep(c(1/27, 1/800, 1/250, 1/1000), each = 9),
                  rep(c(1,2, rep(3,3), rep(4,3), 5), times = 4),
                  rep(c(1,1, c(10, 50, 200), c(1,2,4), 1), times = 4),
                  rep(c(3000, rep(100, 8)), times = 4))
colnames(paramMat) <- c("n_each", "d_each", "sigma",
                        "k", "true_scalar", "true_r", "max_iter",
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
    obs_mat <- generator_esvd_nb(nat_mat, size = vec["true_r"])
  } else {
    obs_mat <- round(generator_curved_gaussian(nat_mat, scalar = vec["true_scalar"]))
  }

  list(dat = obs_mat, truth = res$cell_mat, nat_mat = nat_mat)
}

y <- 1
set.seed(y)
vec <- paramMat[36,]
dat <- rule(vec)

set.seed(10*y)
dat_obs <- dat$dat
n <- nrow(dat_obs); d <- ncol(dat_obs)
init <- eSVD::initialization(dat_obs, family = "exponential", k = vec["k"], max_val = vec["max_val"],
                             scalar = vec["fitting_param"])
fit <- eSVD::fit_factorization(dat_obs, u_mat = init$u_mat, v_mat = init$v_mat,
                               family = "exponential", scalar = vec["fitting_param"],
                               max_iter = vec["max_iter"], max_val = vec["max_val"],
                               return_path = F, cores = ncores,
                               verbose = F)

set.seed(10)

dat_obs <- dat$dat

###########

pred_mat <- -1/(fit$u_mat %*% t(fit$v_mat))
(prod(dim(dat_obs))/(sum((dat_obs-pred_mat)^2/(pred_mat[,2])^2)))^(1/2)

#######################################

y <- 1
set.seed(y)
vec <- paramMat[20,]
dat <- rule(vec)

set.seed(10*y)
dat_obs <- dat$dat
n <- nrow(dat_obs); d <- ncol(dat_obs)
init <- eSVD::initialization(dat_obs, family = "poisson", k = vec["k"], max_val = vec["max_val"],
                             scalar = vec["fitting_param"])
fit <- eSVD::fit_factorization(dat_obs, u_mat = init$u_mat, v_mat = init$v_mat,
                               family = "poisson", scalar = vec["fitting_param"],
                               max_iter = vec["max_iter"], max_val = vec["max_val"],
                               return_path = F, cores = ncores,
                               verbose = F)

#########

pred_mat <- exp(fit$u_mat %*% t(fit$v_mat))

r_seq <- 1:200
vec <- sapply(r_seq, function(x){
  sum((dat_obs - pred_mat)^2/(pred_mat + pred_mat^2/x))/2
})

r_seq[which.min(abs(vec - prod(dim(dat_obs))))]



