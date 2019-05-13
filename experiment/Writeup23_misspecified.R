rm(list=ls())
library(simulation)
library(singlecell)

paramMat <- cbind(50, 120, 0.05, 150, 2, 2, -2000)
colnames(paramMat) <- c("n_each", "d_each", "sigma", "total", "k", "scalar", "max_val")
trials <- 50

################

rule <- function(vec){
  cell_pop <- matrix(c(4,10, 25,100, 60,80, 25,100,
                       40,10, 60,80, 60,80, 100, 25)/10,
                     nrow = 4, ncol = 4, byrow = T)
  gene_pop <- matrix(c(20,90, 25,100,
                       90,20, 100,25)/100, nrow = 2, ncol = 4, byrow = T)
  n_each <- vec["n_each"]
  d_each <- vec["d_each"]
  sigma <- vec["sigma"]
  total <- vec["total"]

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

  pred_mat <- cell_mat %*% t(gene_mat)

  obs_mat <- matrix(0, ncol = ncol(gram_mat), nrow = nrow(gram_mat))
  for(i in 1:n){
    for(j in 1:d){
      obs_mat[i,j] <- rnorm(1, pred_mat[i,j], 0.5)
    }
  }

  obs_mat[obs_mat < 0] <- 0
  obs_mat2 <- round(100*obs_mat)
  # quantile(obs_mat2, probs = seq(0,1,length.out=11))
  # length(which(obs_mat2 == 0))/prod(dim(obs_mat2))

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
  # quantile(obs_mat3, probs = seq(0,1,length.out=11))
  # length(which(obs_mat3 == 0))/prod(dim(obs_mat3))

  list(dat = obs_mat3, truth = cell_mat)
}

criterion <- function(dat, vec, y){
  cluster_labels <- rep(1:4, each = vec["n_each"])

  # SVD
  tmp <- svd(dat$dat)
  res_svd <- tmp$u[,1:vec["k"]] %*% diag(sqrt(tmp$d[1:vec["k"]]))
  curves_svd <- singlecell::slingshot(res_svd[,1:vec["k"]], cluster_labels,
                                      starting_cluster = 1,
                                      verbose = F)

  # ICA
  tmp <- ica::icafast(dat$dat, nc = vec["k"])
  res_ica <- tmp$S
  curves_ica <- singlecell::slingshot(res_ica[,1:vec["k"]], cluster_labels,
                                      starting_cluster = 1,
                                      verbose = F)

  # tsne
  tmp <- Rtsne::Rtsne(dat$dat, perplexity = 30)
  res_tsne <- tmp$Y
  curves_tsne <- singlecell::slingshot(res_tsne[,1:vec["k"]], cluster_labels,
                                       starting_cluster = 1,
                                       verbose = F)

  # Our method
  set.seed(10)
  init <- singlecell::initialization(dat$dat, family = "gaussian", k = vec["k"], max_val = vec["max_val"])
  tmp <- singlecell::fit_factorization(dat$dat, u_mat = init$u_mat, v_mat = init$v_mat,
                                       family = "gaussian",  reparameterize = T,
                                       max_iter = 100, max_val = vec["max_val"],
                                       scalar = vec["scalar"],
                                       return_path = F, cores = 1,
                                       verbose = F)
  res_our <- tmp$u_mat
  curves_our <- singlecell::slingshot(res_our[,1:vec["k"]], cluster_labels,
                                      starting_cluster = 1,
                                      verbose = F)

  # zinb-wave
  set.seed(10)
  dat_se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = t(dat$dat)))
  tmp <- zinbwave::zinbwave(dat_se, K = vec["k"], maxiter.optimize = 100)
  res_zinb <- tmp@reducedDims$zinbwave
  curves_zinb <- singlecell::slingshot(res_zinb[,1:vec["k"]], cluster_labels,
                                       starting_cluster = 1,
                                       verbose = F)

  # pcmf
  set.seed(10)
  tmp <- pCMF::pCMF(dat$dat, K = vec["k"], verbose = F)
  res_pcmf <- tmp$factor$U
  curves_pcmf <- singlecell::slingshot(res_pcmf[,1:2], cluster_labels,
                                       starting_cluster = 1,
                                       verbose = F)

  curves_truth <- singlecell::slingshot(dat$truth, cluster_labels,
                                        starting_cluster = 1,
                                        verbose = F)

  list(res_svd = res_svd, curves_svd = curves_svd,
       res_ica = res_ica, curves_ica = curves_ica,
       res_tsne = res_tsne, curves_tsne = curves_tsne,
       res_our = res_our, curves_our = curves_our,
       res_zinb = res_zinb, curves_zinb = curves_zinb,
       res_pcmf = res_pcmf, curves_pcmf = curves_pcmf,
       curves_truth = curves_truth,
       dat = dat)
}

