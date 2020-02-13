rm(list=ls())
library(simulation)
library(eSVD)
source("../simulation/factorization_generator.R")

paramMat <- cbind(50, 120, 0.5, 150,
                  2, 2, 100, 2000, 50,
                  rep(1:4, each = 8),
                  rep(c(1, 1/150, 1/50, 1/1000), each = 8),
                  rep(c(1,2, rep(3,3), rep(4,3)), times = 4),
                  rep(c(1,1, c(50, 100, 200), c(1,2,4)), times = 4))
colnames(paramMat) <- c("n_each", "d_each", "sigma", "total",
                        "k", "true_scalar", "true_r", "max_val", "max_iter",
                        "true_distr",
                        "modifier",
                        "fitting_distr",
                        "fitting_param")
trials <- 20
ncores <- 20


cell_pop <- matrix(c(4,10, 25,100, 60,80, 25,100,
                     40,10, 60,80, 60,80, 100, 25),
                   nrow = 4, ncol = 4, byrow = T)
gene_pop <- matrix(c(20,90, 25,100,
                     90,20, 100,25)/20, nrow = 2, ncol = 4, byrow = T)

i <- 0
vec = paramMat[8*i+1,]

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
  obs_mat <- generator_esvd_nb(nat_mat, r = vec["true_r"])
} else {
  obs_mat <- round(generator_curved_gaussian(nat_mat, scalar = vec["true_scalar"]))
}

range(nat_mat)
range(obs_mat)
range(res$cell_mat)
