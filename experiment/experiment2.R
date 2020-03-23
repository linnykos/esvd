rm(list=ls())
library(SingleCellExperiment)
source("../simulation/factorization_generator.R")

paramMat <- cbind(50, 120, 10,
                  2, 50, 1/250, 1000,
                  80, 120, 600,
                  1/4, 1/4, 1/2)
colnames(paramMat) <- c("n_each", "d_each", "sigma",
                        "k", "max_iter", "modifier", "max_val",
                        "size_1", "size_2", "size_3",
                        "prop_1", "prop_2", "prop_3")

cell_pop <- matrix(c(4,10, 25,100, 60,80, 25,100,
                     40,10, 60,80, 60,80, 100, 25),
                   nrow = 4, ncol = 4, byrow = T)
gene_pop <- matrix(c(20,90, 25,100,
                     90,20, 100,25)/20, nrow = 2, ncol = 4, byrow = T)

vec <- paramMat[1,]
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
# r_vec <- rep(600, ncol(nat_mat))

dat <- generator_zinb_nb(nat_mat, r_vec)
obs_mat <- round(dat$dat * 1000/max(dat$dat))

cluster_labels <- rep(1:4, each = vec["n_each"])


set.seed(10)
tmp <- pCMF::pCMF(obs_mat, K = 2, sparsity = F, verbose = F)

save.image("../experiment/experiment.RData")
