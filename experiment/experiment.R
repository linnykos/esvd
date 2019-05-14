rm(list=ls())
library(simulation)
library(singlecell)
source("../simulation/factorization_generator.R")

paramMat <- cbind(50, 120, 0.05, 150, 2, 2, -2000)
colnames(paramMat) <- c("n_each", "d_each", "sigma", "total", "k", "scalar", "max_val")
trials <- 50
vec <- paramMat[1,]

#####################

cell_pop <- matrix(c(4,10, 25,100, 60,80, 25,100,
                     40,10, 60,80, 60,80, 100, 25)/10,
                   nrow = 4, ncol = 4, byrow = T)
gene_pop <- matrix(c(20,90, 25,100,
                     90,20, 100,25)/100, nrow = 2, ncol = 4, byrow = T)


n_each <- vec["n_each"]
d_each <- vec["d_each"]
sigma <- vec["sigma"]
total <- vec["total"]

res <- generate_natrual_mat(cell_pop, gene_pop, n_each, d_each, sigma)
nat_mat <- res$nat_mat

set.seed(10)
obs_mat <- generator_pcmf_poisson(nat_mat, dropout_prob = 0.5)

set.seed(10)
tmp <- pCMF::pCMF(obs_mat, K = vec["k"], verbose = T, sparsity = F)
res_pcmf <- tmp$factor$U
plot(res_pcmf[,1], res_pcmf[,2], asp = T, pch = 16, col = c(1:4)[rep(1:4, each = paramMat[1,"n_each"])])

##############

set.seed(10)
init <- singlecell::initialization(obs_mat, family = "gaussian", k = vec["k"], max_val = vec["max_val"])
tmp <- singlecell::fit_factorization(obs_mat, u_mat = init$u_mat, v_mat = init$v_mat,
                                     family = "gaussian",  reparameterize = T,
                                     max_iter = 100, max_val = vec["max_val"],
                                     scalar = vec["scalar"],
                                     return_path = F, cores = 1,
                                     verbose = F)
res_our <- tmp$u_mat
plot(res_our[,1], res_our[,2], asp = T, pch = 16, col = c(1:4)[rep(1:4, each = paramMat[1,"n_each"])])

