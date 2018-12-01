rm(list=ls())
library(simulation)
library(singlecell)

paramMat <- cbind(round(exp(seq(log(10), log(200), length.out = 10))),
                  round(exp(seq(log(20), log(400), length.out = 10))))
colnames(paramMat) <- c("n", "d")
trials <- 50

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
  cell_mat <- res$u_mat; gene_mat <- res$v_mat

  obs_mat <- matrix(0, ncol = ncol(gram_mat), nrow = nrow(gram_mat))
  for(i in 1:n){
    for(j in 1:d){
      obs_mat[i,j] <- distr_func(max(gram_mat[i,j], 1e-4))
    }
  }

  obs_mat[obs_mat < 0] <- 0

  list(dat = obs_mat, cell_mat = cell_mat, gene_mat = gene_mat,
       gram_mat = gram_mat, n_each = n_each, d_each = d_each,
       h = h, g = g, k = k)
}

set.seed(10)
obj <- .data_generator(cell_pop, gene_pop, n_each = 50, d_each = 60)

scalar_vec <- c(0.1, 0.5, 1, 1.5, 2, 2.5, 3, 5, 10, 100)
res_list <- vector("list", length(scalar_vec))
max_val <- 10

for(i in 1:length(scalar_vec)){
  init <- singlecell:::.initialization(obj$dat, family = "gaussian", scalar = scalar_vec[i], k = 2,
                                       max_val = max_val)
  res_list[[i]] <- singlecell:::.fit_factorization(obj$dat, u_mat = init$u_mat, v_mat = init$v_mat,
                                                   family = "gaussian",
                                                   max_iter = 10, max_val = max_val,
                                                   scalar = scalar_vec[i],
                                                   verbose = T)

  save.image("../results/gaussian_scalar.RData")
}

save.image("../results/gaussian_scalar.RData")

# for(j in 1:length(res_list)){
#   plot(res_list[[j]]$u_mat[,1], res_list[[j]]$u_mat[,2],
#        col = c(1:4)[rep(1:4, each = 50)], asp = T, pch = 16,
#        main = j)
# }



