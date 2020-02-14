rm(list=ls())
library(simulation)
library(eSVD)
source("../simulation/factorization_generator.R")

paramMat <- cbind(50, 120, 5,
                  2, 2, 100, 50,
                  rep(1:4, each = 8),
                  rep(c(1, 1/150, 1/50, 1/1000), each = 8),
                  rep(c(1,2, rep(3,3), rep(4,3)), times = 4),
                  rep(c(1,1, c(50, 100, 200), c(1,2,4)), times = 4),
                  rep(c(1000, rep(100, 7)), times = 4))
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

vec = paramMat[1,]
n_each <- vec["n_each"]
d_each <- vec["d_each"]
sigma <- vec["sigma"]
total <- vec["total"]
modifier <- vec["modifier"]

res <- generate_natural_mat(cell_pop, gene_pop, n_each, d_each, sigma, modifier)
