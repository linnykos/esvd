rm(list=ls())
library(simulation)
library(eSVD)
source("../simulation/factorization_generator.R")

paramMat <- cbind(50, 120, 10,
                  2, 50, 2, 50, 10,
                  rep(1:4, each = 7),
                  rep(c(1/27, 1/800, 1/250, 1/1000), each = 7),
                  rep(c(1 ,2, rep(3,2), rep(4,2), 5), times = 4),
                  rep(c(1, 1, 50, NA, 2, NA, 1), times = 4),
                  rep(c(3000, rep(100, 6)), times = 4))
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
    obs_mat <- generator_esvd_nb(nat_mat, size = vec["true_r"])
  } else {
    obs_mat <- round(generator_curved_gaussian(nat_mat, scalar = vec["true_scalar"]))
  }

  list(dat = obs_mat, truth = res$cell_mat, nat_mat = nat_mat)
}

#################3

i <- 18; y <- 1
set.seed(y)
vec <- paramMat[i,]
dat <- rule(vec)

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

init <- eSVD::initialization(dat_obs, family = "poisson", k = vec["k"], max_val = vec["max_val"])
fit <- eSVD::fit_factorization(dat_obs, u_mat = init$u_mat, v_mat = init$v_mat,
                               family = "poisson",
                               max_iter = vec["max_iter"], max_val = vec["max_val"],
                               return_path = F, cores = ncores,
                               verbose = F)

# repetition
fitting_vec <- rep(NA, vec["fitting_iter"])
fitting_vec[1] <- eSVD::tuning(dat_obs, fit$u_mat, fit$v_mat, family = "neg_binom")

i = 2
init <- eSVD::initialization(dat_obs, family = "neg_binom", k = vec["k"], max_val = vec["max_val"],
                             size = fitting_vec[i-1])
fit <- eSVD::fit_factorization(dat_obs, u_mat = init$u_mat, v_mat = init$v_mat,
                               family = "neg_binom", size = fitting_vec[i-1],
                               max_iter = vec["max_iter"], max_val = vec["max_val"],
                               return_path = F, cores = ncores,
                               verbose = F)

fitting_vec[i] <- eSVD::tuning(dat_obs, fit$u_mat, fit$v_mat, family = "neg_binom")

########
dat <- dat_obs
u_mat <- fit$u_mat
v_mat <- fit$v_mat

p_mat <- exp(u_mat %*% t(v_mat))
p_mat[1:5,1:5]

r_seq <- 1:100
proposed_val <- t(sapply(r_seq, function(x){
  pred_mat <- x*p_mat/(1-p_mat)
  c(sum((dat - pred_mat)^2)/prod(dim(dat)),
    sum((pred_mat + pred_mat^2/x))/prod(dim(dat)))
}))

r_seq[which.min(abs(proposed_val[,1] - proposed_val[,2]))]
