rm(list=ls())
library(simulation)
library(singlecell)

paramMat <- cbind(round(exp(seq(log(10), log(200), length.out = 10))),
                  round(exp(seq(log(20), log(400), length.out = 10))))
colnames(paramMat) <- c("n", "d")
trials <- 5

################

# setup
cell_pop <- matrix(c(4,10, 25,100, 60,80, 25,100,
                     40,10, 60,80, 60,80, 100,25)/100,
                   nrow = 4, ncol = 4, byrow = T)
gene_pop <- matrix(c(20,90, 25,100,
                     90,20, 100,25)/100, nrow = 2, ncol = 4, byrow = T)

.data_generator <- function(cell_pop, gene_pop,
                            distr_func = function(x){stats::rnorm(1, 4/x, sd = 2/x)},
                            n_each = 50, d_each = 120, sigma = 0.05){

  #construct the cell information
  h <- nrow(cell_pop)
  cell_mat <- do.call(rbind, lapply(1:h, function(x){
    pos <- stats::runif(n_each)
    cbind(pos*cell_pop[x,1] + (1-pos)*cell_pop[x,3] + stats::rnorm(n_each, sd = sigma),
          pos*cell_pop[x,2] + (1-pos)*cell_pop[x,4] + stats::rnorm(n_each, sd = sigma))
  }))
  n <- nrow(cell_mat)
  k <- ncol(cell_mat)

  # construct the gene information

  g <- nrow(gene_pop)
  gene_mat <- do.call(rbind, lapply(1:g, function(x){
    pos <- stats::runif(d_each)
    cbind(pos*gene_pop[x,1] + (1-pos)*gene_pop[x,3] + stats::rnorm(d_each, sd = sigma),
          pos*gene_pop[x,2] + (1-pos)*gene_pop[x,4] + stats::rnorm(d_each, sd = sigma))
  }))
  d <- nrow(gene_mat)

  # form observations
  gram_mat <- cell_mat %*% t(gene_mat) #natural parameter
  svd_res <- svd(gram_mat)
  cell_mat <- svd_res$u[,1:k] %*% diag(sqrt(svd_res$d[1:k]))
  gene_mat <- svd_res$v[,1:k] %*% diag(sqrt(svd_res$d[1:k]))

  res <- singlecell:::.reparameterize(cell_mat, gene_mat)
  cell_mat <- res$X; gene_mat <- res$Y

  obs_mat <- matrix(0, ncol = ncol(gram_mat), nrow = nrow(gram_mat))
  for(i in 1:n){
    for(j in 1:d){
      obs_mat[i,j] <- distr_func(max(gram_mat[i,j], 1e-4))
    }
  }

  obs_mat[obs_mat < 0] <- 0

  list(dat = obs_mat, cell_mat = cell_mat, gene_mat = gene_mat,
       gram_mat = gram_mat, n_each = n_each, d_each = d_each,
       h = h, g = g, k = k, A = res$A)
}
#############

rule <- function(vec){
  .data_generator(cell_pop, gene_pop, n_each = vec["n"], d_each = vec["d"])$dat
}

criterion <- function(dat, vec, y){
  init <- singlecell:::.initialization(dat, family = "gaussian", max_val = 10)
  res <- singlecell:::.fit_factorization(dat, init$u_mat, init$v_mat,
                            max_val = 5, family = "gaussian", verbose = T,
                            max_iter = 50, reparameterize = T,
                            return_path = F)
}

# set.seed(1); criterion(rule(paramMat[1,]), paramMat[1,], 1)

#################

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = trials,
                                        cores = 5, as_list = T,
                                        filepath = "../results/wasserstein_tmp.RData",
                                        verbose = T)

save.image("../results/wasserstein_simulation.RData")
