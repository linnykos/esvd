rm(list=ls())
library(simulation)
library(singlecell)

paramMat <- cbind(50, 120, 0.01, 150, 3, 2, 10)
colnames(paramMat) <- c("n", "d", "sigma", "total", "k", "scalar", "max_val")
trials <- 2

################

# setup
cell_pop <- matrix(c(4,10, 25,100, 60,80, 25,100,
                     40,10, 60,80, 60,80, 140,60)/100,
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

  extra_weight <- rmutil::rlaplace(nrow(cell_mat), m = 11, s = 0.5)
  pred_mat <- 1/(cell_mat %*% t(gene_mat))
  pred_mat <- t(sapply(1:nrow(pred_mat), function(x){
    pred_mat[x,] * extra_weight[x]
  }))

  obs_mat <- matrix(0, ncol = ncol(gram_mat), nrow = nrow(gram_mat))
  for(i in 1:n){
    for(j in 1:d){
      obs_mat[i,j] <- round(rnorm(1, pred_mat[i,j], pred_mat[i,j]/scalar))
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
  # dat2 <- dat; dat2 <- t(apply(dat, 1, function(x){x/sum(x)}))
  dropout_mat <- singlecell:::.dropout(dat)
  zero_mat <- singlecell:::.find_true_zeros(dropout_mat, num_neighbors = 50)
  idx1 <- intersect(which(dat == 0), which(obj$dat_nodropout != 0))

  dat_impute <- singlecell:::.scImpute(dat, idx1, Kcluster = 4,
                                       verbose = F, weight = 1)

  dat_impute

#   idx1 <- intersect(which(dat == 0), which(obj$dat_nodropout != 0))
#   idx2 <- which(is.na(zero_mat))
#   idx3 <- which(dropout_mat == 0)
#   length(intersect(idx1, idx2))/length(unique(c(idx1, idx2)))
#
#   plot(obj$dat_nodropout[which(dat == 0)], dat_impute[which(dat == 0)], asp = T,
#        pch = 16, col = rgb(0,0,0,0.1))
#   plot(obj$dat_nodropout[which(is.na(zero_mat))], dat_impute[which(is.na(zero_mat))], asp = T,
#        pch = 16, col = rgb(0,0,0,0.1))

  # impute_res <- SAVER::saver(t(obj$dat))
  # new_dat <- t(impute_res$estimate)
  # dat <- obj$dat
  # plot(obj$dat_nodropout[which(dat == 0)], new_dat[which(dat == 0)], asp = T)
  #
  # dat[which(dat == 0)] <- new_dat[which(dat == 0)]
  # dat
}

criterion <- function(dat, vec, y){
  tmp <- svd(dat)
  res_svd <- tmp$u[,1:vec["k"]] %*% diag(sqrt(tmp$d[1:vec["k"]]))

  tmp <- ica::icafast(dat, nc = vec["k"])
  res_ica <- tmp$S

  extra_weight <- apply(dat, 1, mean)

  print("Starting our factorization")
  init <- singlecell::initialization(dat, family = "gaussian", max_val = vec["max_val"],
                                       k = vec["k"])
  res_our <- singlecell::fit_factorization(dat, init$u_mat, init$v_mat,
                                         max_val = vec["max_val"],
                                         family = "gaussian", verbose = F,
                                         max_iter = 25, reparameterize = T,
                                         extra_weight = extra_weight, scalar = vec["scalar"],
                                         return_path = F, cores = 15)

  list(res_svd = res_svd, res_ica = res_ica, res_our = res_our)
}

# set.seed(1); criterion(rule(paramMat[1,]), paramMat[1,], 1)

############

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = trials,
                                        cores = 1, as_list = T,
                                        filepath = "../results/factorization_results_tmp.RData",
                                        verbose = T)

save.image("../results/factorization_results.RData")

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
