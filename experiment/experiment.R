rm(list=ls())
library(simulation)
library(eSVD)
source("../simulation/factorization_generator.R")

paramMat <- cbind(25, 60, 0.5, 150,
                  2, 2, 100, 50,
                  rep(1:4, each = 8),
                  rep(c(1, 1/150, 1/50, 1/1000), each = 8),
                  rep(c(1,2, rep(3,3), rep(4,3)), times = 4),
                  rep(c(1,1, c(50, 100, 200), c(1,2,4)), times = 4),
                  rep(c(1000, rep(100, 7)), times = 4))
colnames(paramMat) <- c("n_each", "d_each", "sigma", "total",
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

  list(dat = obs_mat, truth = res$cell_mat, alt = res$gene_mat, nat = res$nat_mat)
}

vec <- paramMat[2,]
dat <- rule(vec)
y <- 1

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

init <- eSVD::initialization(dat_obs, family = "poisson", k = vec["k"], max_val = vec["max_val"])
fit <- eSVD::fit_factorization(dat_obs, u_mat = init$u_mat, v_mat = init$v_mat,
                               family = "poisson",
                               max_iter = vec["max_iter"], max_val = vec["max_val"],
                               return_path = F, cores = ncores,
                               verbose = F)

pred_mat <- fit$u_mat %*% t(fit$v_mat)
pred_val <- exp(pred_mat[missing_idx])


range(dat$truth)
range(dat$alt)
range(dat$nat)
range(dat_obs)
range(fit$u_mat)
range(fit$v_mat)
range(fit$u_mat %*% t(fit$v_mat))

cbind(pred_val, missing_val)
