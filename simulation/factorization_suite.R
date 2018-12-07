rm(list=ls())
library(simulation)
library(singlecell)

paramMat <- cbind(50, 120, 0.01, 150, 3, 2, 10)
colnames(paramMat) <- c("n", "d", "sigma", "total", "k", "scalar", "max_val")
trials <- 50

################

# setup
cell_pop <- matrix(c(4,10, 25,100, 60,80, 25,100,
                     40,10, 60,80, 60,80, 100, 25)/10,
                   nrow = 4, ncol = 4, byrow = T)
gene_pop <- matrix(c(20,90, 25,100,
                     90,20, 100,25)/100, nrow = 2, ncol = 4, byrow = T)

.data_generator <- function(cell_pop, gene_pop,
                            n_each = 50, d_each = 100, sigma = 0.05,
                            scalar = 2, total = 150){
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

  extra_weight <- rep(1, nrow(cell_mat))
  pred_mat <- 1/(cell_mat %*% t(gene_mat))
  pred_mat <- t(sapply(1:nrow(pred_mat), function(x){
    pred_mat[x,] * extra_weight[x]
  }))

  obs_mat <- matrix(0, ncol = ncol(gram_mat), nrow = nrow(gram_mat))
  for(i in 1:n){
    for(j in 1:d){
      obs_mat[i,j] <- rnorm(1, pred_mat[i,j], pred_mat[i,j]/scalar)
    }
  }

  obs_mat[obs_mat < 0] <- 0
  obs_mat2 <- round(exp(obs_mat*14)-1)
  # length(which(obs_mat2 > 5000))/prod(dim(obs_mat2))
  obs_mat2[obs_mat2 > 5000] <- 5000
  # quantile(obs_mat2)

  # now do something more dramatic with dropout
  obs_mat3 <- obs_mat2
  .dropped_indices <- function(x, total){
    vec <- 1:length(x)
    samp <- sample(vec, size = total, replace = T, prob = x)
    setdiff(vec, unique(samp))
  }

  total_vec <- rep(total, nrow(obs_mat3))
  for(i in 1:nrow(obs_mat3)){
    idx <- .dropped_indices(obs_mat[i,], total = total_vec[i])
    obs_mat3[i,idx] <- 0
  }

  list(dat = obs_mat3, dat_nodropout = obs_mat2,
       extra_weight = extra_weight,
       cell_mat = cell_mat, gene_mat = gene_mat,
       gram_mat = cell_mat %*% t(gene_mat), n_each = n_each, d_each = d_each,
       h = h, g = g, k = k)
}

############

rule <- function(vec){
  obj <- .data_generator(cell_pop, gene_pop, n_each = vec["n"], d_each = vec["d"],
                         sigma = vec["sigma"], scalar = vec["scalar"], total = vec["total"])
  dat <- obj$dat
  dat <- log(dat+1)
  dropout_mat <- singlecell::dropout(dat)
  zero_mat <- singlecell::find_true_zeros(dropout_mat, num_neighbors = 50)
  idx <- which(is.na(zero_mat))

  dat_impute <- singlecell::scImpute(dat, drop_idx = idx, Kcluster = 4,
                                       verbose = F, weight = 1)

  list(dat = obj$dat, dat_impute = dat_impute)
}

criterion <- function(dat, vec, y){
  tmp <- svd(dat$dat_impute)
  res_svd <- tmp$u[,1:vec["k"]] %*% diag(sqrt(tmp$d[1:vec["k"]]))

  tmp <- ica::icafast(dat$dat_impute, nc = vec["k"])
  res_ica <- tmp$S

  extra_weight <- rep(1, nrow(dat$dat_impute))

  print("Starting our factorization")
  init <- singlecell::initialization(dat$dat_impute, family = "gaussian", max_val = vec["max_val"],
                                       k = vec["k"])
  res_our <- singlecell::fit_factorization(dat$dat_impute, init$u_mat, init$v_mat,
                                         max_val = vec["max_val"],
                                         family = "gaussian", verbose = F,
                                         max_iter = 25, reparameterize = T,
                                         extra_weight = extra_weight, scalar = vec["scalar"],
                                         return_path = F, cores = 15)

  list(res_svd = res_svd, res_ica = res_ica, res_our = res_our, dat = dat$dat,
       dat_impute = dat$dat_impute)
}

# set.seed(1); criterion(rule(paramMat[1,]), paramMat[1,], 1)

############

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = trials,
                                        cores = 1, as_list = T,
                                        filepath = "../results/factorization_results_tmp.RData",
                                        verbose = T)

save.image("../results/factorization_results.RData")



# extra_weight <- apply(dat, 1, mean)
#

# # VAMF
# print(paste0(Sys.time(), ": VAMF"))
# tmp <- vamf:::vamf(t(dat), k, log2trans=T)
# res_vamf <- tmp$factors
#
# save.image("../results/factorization_results_tmp.RData")
#
# # zinbwave
# print(paste0(Sys.time(), ": ZINB"))
# dat_se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = t(dat)))
# tmp <- zinbwave::zinbwave(dat_se, K = k, maxiter.optimize = 100)
# res_zinb <- tmp@reducedDims$zinbwave
#
# save.image("../results/factorization_results_tmp.RData")
#
#
#
