method_svd <- function(dat, k = 2){
  stopifnot(k >= 2)
  tmp <- svd(dat)
  fit <- tmp$u[,1:k] %*% diag(sqrt(tmp$d[1:k]))

  list(fit = fit)
}

method_tsne_oracle <- function(dat, cell_truth, paramMat, k = 2){
  fit_list <- lapply(1:nrow(paramMat), function(i){
    set.seed(10)
    Rtsne::Rtsne(dat, dims = k, perplexity = paramMat[i, "perplexity"])
  })

  dist_mat_truth <- as.matrix(stats::dist(cell_truth))

  quality_vec <- sapply(fit_list, function(fit){
    dist_mat_est <- as.matrix(stats::dist(fit[,1:2]))

    mean(sapply(1:nrow(dist_mat_est), function(i){
      cor(dist_mat_truth[i,], dist_mat_est[i,], method = "kendall")
    }))
  })

  idx <- which.max(quality_vec)

  fit <- fit_list[[idx]]

  list(fit = fit, perplexity = paramMat[idx, "perplexity"])
}

method_umap_oracle <- function(dat, cell_truth, paramMat, k = 2){
  fit_list <- lapply(1:nrow(paramMat), function(i){
    custom.settings <- umap::umap.defaults
    custom.settings$n_components <- k
    custom.settings$n_neighbors <- paramMat[i, "n_neighbors"]
    custom.settings$min_dist <- paramMat[i, "min_dist"]

    set.seed(10)
    tmp <- umap::umap(dat, config = custom.settings)
    tmp$layout
  })

  dist_mat_truth <- as.matrix(stats::dist(cell_truth))

  quality_vec <- sapply(fit_list, function(fit){
    dist_mat_est <- as.matrix(stats::dist(fit[,1:2]))

    mean(sapply(1:nrow(dist_mat_est), function(i){
      cor(dist_mat_truth[i,], dist_mat_est[i,], method = "kendall")
    }))
  })

  idx <- which.max(quality_vec)

  fit <- fit_list[[idx]]

  list(fit = fit, min_dist = paramMat[idx, "min_dist"], n_neighbors = paramMat[idx, "n_neighbors"])
}

method_zinbwave <- function(dat, k = 2){
  set.seed(10)
  dat_se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = t(dat)))
  tmp <- zinbwave::zinbwave(dat_se, K = 2, maxiter.optimize = 100, normalizedValues = F,
                            commondispersion = F)
  fit <- SingleCellExperiment::reducedDims(tmp)$zinbwave

  list(fit = fit)
}

method_pcmf <- function(dat, k = 2){
  set.seed(10)
  tmp <- pCMF::pCMF(obs_mat, K = 2, sparsity = F, verbose = F)
  fit <- tmp$factor$U

  list(fit = fit)
}

method_esvd <- function(dat, paramMat, k = 2, ncores = NA){
  set.seed(10)
  missing_idx <- eSVD::construct_missing_values(n = nrow(dat), p = ncol(dat), num_val = 2)
  dat_NA <- dat
  dat_NA[missing_idx] <- NA

  fit_list <- lapply(1:nrow(paramMat), function(i){
    set.seed(10)
    init <- eSVD::initialization(dat_NA, family = "neg_binom", k = paramMat[i, "k"], max_val = 2000,
                                 scalar = paramMat[i, "scalar"])
    eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                   family = "neg_binom", scalar = paramMat[i, "scalar"],
                                   max_iter = 50, max_val = 2000,
                                   return_path = F, cores = ncores,
                                   verbose = F)
  })

  quality_vec <- sapply(1:nrow(paramMat), function(i){
    nat_mat <- fit_list[[i]]$u_mat %*% t(fit_list[[i]]$v_mat)
    mean_mat <- eSVD::compute_mean(nat_mat, family = "neg_binom", scalar = paramMat[i, "scalar"])
    eSVD::plot_prediction_against_observed(dat, nat_mat_list = list(nat_mat),
                                           family = "neg_binom", missing_idx_list = list(missing_idx),
                                           plot = F)
  })

  idx <- which.min(abs(quality_vec - 45))

  k <- paramMat[idx, "k"]
  scalar <- paramMat[idx, "scalar"]

  set.seed(10)
  init <- eSVD::initialization(dat, family = "neg_binom", k = k, max_val = 2000,
                               scalar = scalar)
  fit <- eSVD::fit_factorization(dat, u_mat = init$u_mat, v_mat = init$v_mat,
                                 family = "neg_binom", scalar = scalar,
                                 max_iter = 50, max_val = 2000,
                                 return_path = F, cores = ncores, verbose = F)

  list(fit = fit, k = k, scalar = scalar)
}
