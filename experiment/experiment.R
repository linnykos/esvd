rm(list=ls())
library(simulation)
library(singlecell)
library(Rtsne); library(pCMF); library(zinbwave); library(SummarizedExperiment)
source("../simulation/factorization_generator.R")

paramMat <- cbind(50, 120, 0.05, 150, 2, 2, -2000)
colnames(paramMat) <- c("n_each", "d_each", "sigma", "total", "k", "scalar", "max_val")
trials <- 1

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

  obs_mat <- generator_exponential(nat_mat)
  obs_mat2 <- round(100*obs_mat)

  # quantile(obs_mat2, probs = seq(0,1,length.out=11))
  # length(which(obs_mat2 == 0))/prod(dim(obs_mat2))

  # now do something more dramatic with dropout
  obs_mat3 <- generate_dropout(obs_mat2, total = total)
  # quantile(obs_mat3, probs = seq(0,1,length.out=11))
  # length(which(obs_mat3 == 0))/prod(dim(obs_mat3))

  list(dat = obs_mat3, truth = res$cell_mat)
}

set.seed(1)
vec <- paramMat[1,]
dat <- rule(vec)
set.seed(10)
init <- singlecell::initialization(dat$dat, family = "gaussian", k = vec["k"], max_val = vec["max_val"])

################

dat = dat$dat
family = "gaussian"
k = vec["k"]
max_val = vec["max_val"]
tol = 1e-3
verbose = F

direction <- .dictate_direction(family)

# initialize
dat <- .matrix_completion(dat, k = k)
if(length(class(dat)) == 1) class(dat) <- c(family, class(dat)[length(class(dat))])

# projected gradient descent
# pred_mat <- .projected_gradient_descent(dat, k = k, max_val = max_val,
#                                         direction = direction,
#                                         max_iter = max_iter,
#                                         tol = tol)

n <- nrow(dat); d <- ncol(dat)
pred_mat <- .determine_initial_matrix(dat, class(dat)[1], k = k, max_val = max_val)

######

min_val <- min(dat[which(dat > 0)])
dat[which(dat <= 0)] <- min_val/2
pred_mat <- .mean_transformation(dat, family)
direction <- .dictate_direction(family)

class(pred_mat) <- "matrix" #bookeeping purposes
.nonnegative_matrix_factorization(pred_mat, k = k, direction = direction,
                                  max_val = max_val)

##############


