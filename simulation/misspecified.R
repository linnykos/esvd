rm(list=ls())
library(simulation)
library(singlecell)

paramMat <- as.matrix(expand.grid(50, 100, c(0.05, 0.0875, 0.125),
                                  c(-1, -0.5, 0, 0.5, 1)))
colnames(paramMat) <- c("n", "d", "sigma", "misspecify")
trials <- 1

#######

cell_pop <- matrix(c(4,10, 25,100, 60,80, 25,100,
                     40,10, 60,80, 60,80, 100,25)/100,
                   nrow = 4, ncol = 4, byrow = T)
gene_pop <- matrix(c(20,90, 25,100,
                     90,20, 100,25)/100, nrow = 2, ncol = 4, byrow = T)

.data_generator <- function(cell_pop, gene_pop,
                            distr_func, vec,
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
  cell_mat <- res$u_mat; gene_mat <- res$v_mat

  obs_mat <- matrix(0, ncol = ncol(gram_mat), nrow = nrow(gram_mat))
  for(i in 1:n){
    for(j in 1:d){
      obs_mat[i,j] <- distr_func(max(gram_mat[i,j], 1e-4), vec)
    }
  }

  obs_mat[obs_mat < 0] <- 0

  list(dat = obs_mat, cell_mat = cell_mat, gene_mat = gene_mat,
       gram_mat = gram_mat, n_each = n_each, d_each = d_each,
       h = h, g = g, k = k, A = res$A)
}

distr_func <- function(x, vec){
  val1 <- stats::rnorm(1, 4/x, sd = 2/x)
  if(vec["misspecify"] >= 0){
    val2 <- stats::rnorm(1, 4/x, sd = 3)
  } else {
    val2 <- stats::rexp(1, 1/x)
  }

  (1-abs(vec["misspecify"]))*val1 + abs(vec["misspecify"])*val2
}

########

rule <- function(vec){
  .data_generator(cell_pop, gene_pop, distr_func, vec,
                  n_each = vec["n"], d_each = vec["d"],
                  sigma = vec["sigma"])$dat
}

criterion <- function(dat, vec, y){
  init <- singlecell:::.initialization(dat, family = "gaussian", max_val = 10)
  res <- singlecell:::.fit_factorization(dat, init$u_mat, init$v_mat,
                                         max_val = 5, family = "gaussian", verbose = T,
                                         max_iter = 50, reparameterize = T, cores = 5,
                                         return_path = F)
}

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = trials,
                                        cores = NA, as_list = T,
                                        filepath = "../results/misspecified_tmp.RData",
                                        verbose = T)

save.image("../results/misspecified_simulation.RData")

