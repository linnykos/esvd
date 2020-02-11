rm(list=ls())
library(simulation)
library(eSVD)
library(NMF)
source("../simulation/factorization_generator.R")

paramMat <- cbind(round(exp(seq(log(10), log(200), length.out = 10))),
                  round(exp(seq(log(20), log(400), length.out = 10))),
                  0.05, 150, 2, 2, -2000)
colnames(paramMat) <- c("n_each", "d_each", "sigma", "total", "k", "scalar", "max_val")
trials <- 200

save.image("/raid6/Kevin/singlecell_results/simulation/wasserstein_simulation.RData")

################

cell_pop <- matrix(c(4,10, 25,100, 60,80, 25,100,
                     40,10, 60,80, 60,80, 100, 25)/10,
                   nrow = 4, ncol = 4, byrow = T)
gene_pop <- matrix(c(20,90, 25,100,
                     90,20, 100,25)/100, nrow = 2, ncol = 4, byrow = T)

rule <- function(vec){
  n_each <- vec["n_each"]
  d_each <- vec["d_each"]
  sigma <- vec["sigma"]
  total <- vec["total"]

  res <- generate_natural_mat(cell_pop, gene_pop, n_each, d_each, sigma)
  nat_mat <- res$nat_mat

  obs_mat <- generator_curved_gaussian(nat_mat, alpha = vec["scalar"])

  list(dat = obs_mat, truth = res$cell_mat)
}

criterion <- function(dat, vec, y){
  # Our method
  set.seed(y)
  init <- eSVD::initialization(dat$dat, family = "gaussian", k = vec["k"], max_val = vec["max_val"])
  tmp <- eSVD::fit_factorization(dat$dat, u_mat = init$u_mat, v_mat = init$v_mat,
                                       family = "gaussian",  reparameterize = T,
                                       max_iter = 50, max_val = vec["max_val"],
                                       scalar = vec["scalar"], tol = NA,
                                       return_path = F, cores = NA,
                                       verbose = F)
  res_our <- tmp$u_mat

  # try all different orientations
  vec1 <- res_our[,1]; vec2 <- res_our[,2]

  l2_loss1 <- sum((cbind(vec1, vec2) - dat$truth)^2)/nrow(res_our)
  wasserstein_loss1 <- transport::wasserstein(transport::pp(cbind(vec1, vec2)),
                                             transport::pp(dat$truth), p = 1)
  l2_loss2 <- sum((cbind(-vec1, vec2) - dat$truth)^2)/nrow(res_our)
  wasserstein_loss2 <- transport::wasserstein(transport::pp(cbind(-vec1, vec2)),
                                              transport::pp(dat$truth), p = 1)
  l2_loss3 <- sum((cbind(vec1, -vec2) - dat$truth)^2)/nrow(res_our)
  wasserstein_loss3 <- transport::wasserstein(transport::pp(cbind(vec1, -vec2)),
                                              transport::pp(dat$truth), p = 1)
  l2_loss4 <- sum((cbind(-vec1, -vec2) - dat$truth)^2)/nrow(res_our)
  wasserstein_loss4 <- transport::wasserstein(transport::pp(cbind(-vec1, -vec2)),
                                              transport::pp(dat$truth), p = 1)

  list(res_our = res_our,
       l2_loss = min(l2_loss1, l2_loss2, l2_loss3, l2_loss4),
       wasserstein_loss = min(wasserstein_loss1, wasserstein_loss2, wasserstein_loss3, wasserstein_loss4))
}

# set.seed(1); zz1 <- criterion(rule(paramMat[1,]), paramMat[1,], 1)
# set.seed(2); zz2 <- criterion(rule(paramMat[1,]), paramMat[1,], 2)

#################

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = trials,
                                        cores = 15, as_list = T,
                                        filepath = "/raid6/Kevin/singlecell_results/simulation/wasserstein_tmp.RData",
                                        verbose = T)

save.image("wasserstein_simulation.RData")
save.image("/raid6/Kevin/singlecell_results/simulation/wasserstein_simulation.RData")
