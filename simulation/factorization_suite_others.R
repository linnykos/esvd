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

criterion <- function(dat, vec, y){
  if(vec["method"] == 1){ #svd
    method_svd(dat$dat)

  } else if(vec["method"] == 2){ #esvd
    paramMat_esvd <- matrix(c(50, 100, 500, 1000), nrow = 4, ncol = 1)
    colnames(paramMat_esvd) <- c("scalar")

    method_esvd(dat$dat, paramMat = paramMat_esvd, ncores = ncores)

  } else if(vec["method"] == 3){ #zinbwave
    method_zinbwave(dat$dat)

  } else if(vec["method"] == 4){ #pcmf
    method_pcmf(dat$dat)

  } else if(vec["method"] == 5){ #umap
    paramMat_umap <- as.matrix(expand.grid(c(2, 3, 5, 15, 30, 50),
                                      c(1e-5, 1e-3, 0.1, 0.3, 0.5, 0.9)))
    colnames(paramMat_umap) <- c("n_neighbors", "min_dist")

    method_umap_oracle(dat$dat, cell_truth = dat$truth, paramMat = paramMat_umap)

  } else { #tsne
    paramMat_tsne <- matrix(round(seq(2, 50, length.out = 10)), ncol = 1, nrow = 1)
    colnames(paramMat_tsne) <- c("perplexity")

    method_tsne_oracle(dat$dat, cell_truth = dat$truth, paramMat = paramMat_tsne)
  }
}

## i <- 4; y <- 1; zz <- criterion(rule(paramMat[i,]), paramMat[i,], y); zz

############

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = trials,
                                        cores = NA, as_list = T,
                                        filepath = "../results/factorization_results_others_tmp_gen4.RData",
                                        verbose = T)

save.image("../results/factorization_results_others_gen4.RData")
