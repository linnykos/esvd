rm(list=ls())
library(simulation)
library(eSVD)
source("../simulation/factorization_generator.R")
source("../simulation/factorization_methods.R")

session_info <- sessionInfo()
date_of_run <- Sys.time()
source_code_info <- readLines("../simulation/factorization_generator.R")
source_code_info <- c(source_code_info, readLines("../simulation/factorization_methods.R"))
source_code_info <- c(source_code_info, readLines("../simulation/factorization_suite_negbinom_esvd.R"))

paramMat <- cbind(50, 120, 5,
                  2, 50, 1/250, 1000,
                  50, 1:6)
colnames(paramMat) <- c("n_each", "d_each", "sigma",
                        "k", "max_iter", "modifier", "max_val",
                        "size","method")
paramMat <- paramMat[2,,drop = F]

trials <- 100
ncores <- 20
doMC::registerDoMC(cores = ncores)

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

  dat <- generator_esvd_nb(nat_mat, vec["size"])
  obs_mat <- round(dat * 1000/max(dat))

  list(dat = obs_mat, truth = res$cell_mat)
}

criterion <- function(dat, vec, y){

  if(vec["method"] == 1){ #svd
    dat_obs <- dat$dat
    tmp <- method_svd(dat_obs)
    return(list(fit = tmp, truth = dat$truth))

  } else if(vec["method"] == 2){ #esvd
    dat_obs <- dat$dat
    paramMat_esvd <- matrix(c(5, 50, 100), nrow = 3, ncol = 1)
    colnames(paramMat_esvd) <- c("scalar")
    tmp <- method_esvd(dat_obs, paramMat = paramMat_esvd, ncores = ncores, k = 2)

    return(list(fit = tmp, truth = dat$truth))

  } else if(vec["method"] == 3) { #zinbwave
    dat_obs <- dat$dat
    tmp <- method_zinbwave(dat_obs)
    return(list(fit = tmp, truth = dat$truth))

  } else if(vec["method"] == 4) { #pcmf
    dat_obs <- dat$dat
    tmp <- method_pcmf(dat_obs)
    return(list(fit = tmp, truth = dat$truth))

  } else if(vec["method"] == 5) { #umap
    dat_obs <- dat$dat
    paramMat_umap <- as.matrix(expand.grid(c(2, 3, 5, 15, 30, 50),
                                           c(1e-5, 1e-3, 0.1, 0.3, 0.5, 0.9)))
    colnames(paramMat_umap) <- c("n_neighbors", "min_dist")

    tmp <- method_umap_oracle(dat_obs, cell_truth = dat$truth, paramMat = paramMat_umap)
    return(list(fit = tmp, truth = dat$truth))

  } else { #tsne
    dat_obs <- dat$dat
    paramMat_tsne <- matrix(round(seq(2, 50, length.out = 10)), ncol = 1, nrow = 1)
    colnames(paramMat_tsne) <- c("perplexity")

    tmp <- method_tsne_oracle(dat_obs, cell_truth = dat$truth, paramMat = paramMat_tsne)
    return(list(fit = tmp, truth = dat$truth))
  }
}

## i <- 1; y <- 1; set.seed(y); zz1 <- criterion(rule(paramMat[i,]), paramMat[i,], y); head(zz1$fit$fit$u_mat); head(zz1$truth)
## i <- 2; y <- 2; set.seed(y); zz2 <- criterion(rule(paramMat[i,]), paramMat[i,], y); head(zz2$fit$fit$u_mat); head(zz2$truth)
## i <- 1; y <- 1; set.seed(y); y <- rule(paramMat[i,]); y$dat[1:5,1:5]; head(y$truth)

## i <- 2; y <- 1; set.seed(y); zz3 <- criterion(rule(paramMat[i,]), paramMat[i,], y); head(zz3$truth)
## i <- 2; y <- 2; set.seed(y); zz4 <- criterion(rule(paramMat[i,]), paramMat[i,], y); head(zz4$truth)

## i <- 2; y <- 2; set.seed(y); tmp <- rule(paramMat[i,]);tmp$dat[1:5,1:5]

############

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = trials,
                                        cores = NA, as_list = T,
                                        filepath = "../results/factorization_results_negbinom_esvd_tmp.RData",
                                        verbose = T)

save.image("../results/factorization_results_negbinom_esvd.RData")
