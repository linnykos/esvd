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
                     40,10, 60,80, 60,80, 100,25)/50,
                   nrow = 4, ncol = 4, byrow = T)
gene_pop <- matrix(c(20,90, 25,100,
                     90,20, 100,25)/50, nrow = 2, ncol = 4, byrow = T)
distr_func = function(x){stats::rnorm(1, 4/x, sd = 2/x)}

.data_generator <- function(cell_pop, gene_pop,
                            distr_func = function(x){stats::rnorm(1, 4/x, sd = 2/x)},
                            n_each = 50, d_each = 100, sigma = 0.05, total = 150){
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
  obs_mat2 <- obs_mat

  # now do something more dramatic with dropout
  obs_mat3 <- obs_mat2
  .dropped_indices <- function(x, total){
    vec <- 1:length(x)
    samp <- sample(vec, size = total, replace = T, prob = x)
    setdiff(vec, unique(samp))
  }

  total_vec <- rep(total, nrow(obs_mat3))
  for(i in 1:nrow(obs_mat3)){
    idx <- .dropped_indices(obs_mat3[i,], total = total_vec[i])
    obs_mat3[i,idx] <- 0
  }

  list(dat = obs_mat3, dat_nodropout = obs_mat2,
       cell_mat = cell_mat, gene_mat = gene_mat,
       gram_mat = gram_mat, n_each = n_each, d_each = d_each,
       h = h, g = g, k = k)
}


##########

set.seed(10)
n_each <- 50
d_each <- 120
sigma <- 0.05
total <- 250
obj <- .data_generator(cell_pop, gene_pop, n_each = n_each, d_each = d_each,
                       sigma = sigma, total = total)
dat <- obj$dat
impute_res <- SAVER::saver(t(dat))
dat_impute <- t(impute_res$estimate)

############

# now try all the different factorization methods, either on dat or dat_impute
# SVD
print(paste0(Sys.time(), ": SVD"))
tmp <- svd(dat_impute)
res_svd <- tmp$u[,1:2] %*% diag(sqrt(tmp$d[1:2]))

save.image("../results/factorization_results.RData")

# ICA
print(paste0(Sys.time(), ": ICA"))
tmp <- ica::icafast(dat_impute, nc = 2)
res_ica <- tmp$S

save.image("../results/factorization_results.RData")

# our method
print(paste0(Sys.time(), ": Our method"))

max_val <- 10
scalar_vec <- c(0.1, 0.25, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 8, 10, 100)
extra_weight <- apply(dat_impute, 1, mean)

res_our_list <- vector("list", length(scalar_vec))
for(i in 1:length(scalar_vec)){
  init <- singlecell::initialization(dat_impute, family = "gaussian", scalar = scalar_vec[i],
                                       k = 2, max_val = max_val)
  res_our_list[[i]] <- singlecell::fit_factorization(dat_impute, u_mat = init$u_mat, v_mat = init$v_mat,
                                                   family = "gaussian",  reparameterize = T,
                                                   max_iter = 25, max_val = max_val,
                                                   scalar = scalar_vec[i], extra_weight = extra_weight,
                                                   return_path = F, cores = 15,
                                                   verbose = T)
  save.image("../results/factorization_results.RData")
}

save.image("../results/factorization_results.RData")

# # VAMF
# print(paste0(Sys.time(), ": VAMF"))
# res_vamf <- vamf:::vamf(t(dat), 2, log2trans=T)$factors
#
# save.image("../results/factorization_results_tmp.RData")
#
# # zinbwave
# print(paste0(Sys.time(), ": ZINB"))
# dat_se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = t(dat)))
# res_zinb <- zinbwave::zinbwave(dat_se, K = 2)
#
# save.image("../results/factorization_results_tmp.RData")
#
#
#
