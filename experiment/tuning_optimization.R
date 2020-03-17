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
ncores <- NA

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

y <- 1
set.seed(y)
vec <- paramMat[31,]
dat <- rule(vec)

set.seed(10*y)
dat_obs <- dat$dat

init <- eSVD::initialization(dat_obs, family = "poisson", k = vec["k"], max_val = vec["max_val"])
fit <- eSVD::fit_factorization(dat_obs, u_mat = init$u_mat, v_mat = init$v_mat,
                               family = "poisson",
                               max_iter = vec["max_iter"], max_val = vec["max_val"],
                               return_path = F, cores = NA, verbose = F)

dat <- dat_obs
nat_mat <- fit$u_mat %*% t(fit$v_mat)
mean_mat <- compute_mean(nat_mat, family = "poisson")
k <- ncol(fit$u_mat); n <- nrow(dat); p <- ncol(dat)
df_val <- n*p - (n*k + p*k)

f <- function(x, ...){
  sum((dat - mean_mat)^2/(mean_mat + mean_mat^2/x)) - df_val
}

stats::uniroot(f, interval = c(1,10000))



