rm(list=ls())
load("../results/factorization_results_exponential_families_3.RData")

zz <- sapply(res, function(x){
  table(sapply(x, length))
})
zz

# first, determine which scalar is picked when
idx_list <- lapply(1:length(res), function(i){
  # determine what the paramMat is for this setting
  tmp <- names(res)[i]
  tmp <- strsplit(tmp, split = "-")[[1]]
  tmp <- as.numeric(tmp[which(colnames(paramMat) == "fitting_distr")])
  if(tmp %in% c(1,2)){
    if(tmp == 1){
      family_val <- "gaussian"
      # compute the sd_val
    } else {
      family_val <- "poisson"
    }

    sapply(res[[i]], function(y){
      nat_mat <- y$fit$u_mat %*% t(y$fit$v_mat)
      scalar <- stats::sd(y$dat[y$missing_idx] - nat_mat[y$missing_idx]) # compute this even for Poisson, since it doesn't effect anything
      quality_res <- eSVD::plot_prediction_against_observed(dat = y$dat, nat_mat_list = list(nat_mat),
                                                            family = family_val,
                                                            missing_idx_list = list(y$missing_idx),
                                                            scalar = scalar, plot = F)
      if(quality_res$bool){
        c(quality = quality_res$angle_val, scalar = 0)
      } else {
        c(quality = NA, scalar = NA)
      }
    })

  } else {
    if(tmp == 3){
      family_val <- "neg_binom"; scalar_vec <- r_vec
    } else {
      family_val <- "curved_gaussian"; scalar_vec <- alpha_vec
    }

    sapply(res[[i]], function(y){

      # condense nat_mat_list_list
      nat_mat_list_list <- lapply(1:length(y$fit), function(j){
        list(y$fit[[j]]$u_mat %*% t(y$fit[[j]]$v_mat))
      })

      tmp_res <- tuning_select_scalar(dat = y$dat, nat_mat_list_list = nat_mat_list_list,
                                      family = family_val,  missing_idx_list = list(y$missing_idx),
                                      scalar_vec = scalar_vec)

      unlist(tmp_res)
    })
  }
})

# summarize results
sapply(idx_list, function(x){
  c(table(x["scalar",]), sum(is.na(x["scalar",])))
})
