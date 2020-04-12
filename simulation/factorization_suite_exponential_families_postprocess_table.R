rm(list=ls())
load("../results/factorization_exponential_families.RData")

# remove all poisson things
res <- res[-which(paramMat[,"true_distr"] == 2)]
paramMat <- paramMat[-which(paramMat[,"true_distr"] == 2),]
res <- res[-which(paramMat[,"fitting_distr"] == 2)]
paramMat <- paramMat[-which(paramMat[,"fitting_distr"] == 2),]

res_mat_list <- vector("list", 3)
paramMat_list <- vector("list", 3)
tot_num_distr <- 21

for(zz in 1:3){
  kk <- c(1,3,4)[zz]
  # grab the relevant parameter settings
  tmp_idx <-  which(paramMat[,"true_distr"] == kk)
  res_tmp <- res[tmp_idx]
  paramMat_tmp <- paramMat[tmp_idx,]

  # remove any trials that failed anywhere
  tmp_idx <- sort(unique(unlist(lapply(res_tmp, function(x){
    which(sapply(x, function(y){all(is.na(y))}) == 1)
  }))))

  if(length(tmp_idx) > 0){
    for(i in 1:length(res_tmp)){
      res_tmp[[i]] <- res_tmp[[i]][-tmp_idx]
    }
  }

  # among all the remaining trials, across all 3 gaussian, 3 possion, 9 neg binom and 9 curved guassian (24 total),
  # unwrap paramMat for easier usage, where all the fitting_param from NA become 3 different rows
  len <- length(res_tmp[[1]])
  paramMat_new <- matrix(NA, nrow = tot_num_distr, ncol = ncol(paramMat_tmp))
  colnames(paramMat_new) <- colnames(paramMat_tmp)
  res_new <- vector("list", tot_num_distr)
  counter <- 1
  for(i in 1:nrow(paramMat_tmp)){
    if(!is.na(paramMat_tmp[i, "fitting_param"])){
      paramMat_new[counter,] <- paramMat_tmp[i,]
      res_new[[counter]] <- res_tmp[[i]]

      counter <- counter+1
    } else {
      # is it negbinom
      for(j in 1:3){
        tmp <- lapply(res_tmp[[i]], function(x){
          y <- x; y$fit <- x$fit[[j]]; y
        })
        res_new[[counter]] <- tmp

        if(paramMat_tmp[i, "fitting_distr"] == 3){
          paramMat_new[counter,] <- paramMat_tmp[i,]
          paramMat_new[counter,"fitting_param"] <- r_vec[j]
        } else {
          paramMat_new[counter,] <- paramMat_tmp[i,]
          paramMat_new[counter,"fitting_param"] <- alpha_vec[j]
        }

        counter <- counter + 1
      }
    }
  }

  # report a matrix with len columns and 24 rows that contains an entry NA if
  # the principal curve is not within the 10-90 quantile region, and the principal angle otherwise
  res_mat <- matrix(NA, tot_num_distr, len)
  for(i in 1:length(res_new)){
    scalar <- paramMat_new[i, "fitting_param"]
    distr_family <- c("gaussian", "poisson", "neg_binom", "curved_gaussian")[paramMat_new[i, "fitting_distr"]]

    for(j in 1:len){
      nat_mat <- res_new[[i]][[j]]$fit$u_mat %*% t(res_new[[i]][[j]]$fit$v_mat)
      diagnostic_res <- plot_prediction_against_observed(dat = res_new[[i]][[j]]$dat,
                                                         nat_mat_list = list(nat_mat),
                                                         family = distr_family,
                                                         missing_idx_list = list(res_new[[i]][[j]]$missing_idx),
                                                         scalar = scalar, plot = F)
      if(diagnostic_res$bool){
        res_mat[i,j] <- diagnostic_res$angle_val
      } else {
        res_mat[i,j] <- NA
      }

      # res_mat[i,j] <- diagnostic_res$angle_val
    }
  }

  res_mat_list[[zz]] <- res_mat
  paramMat_list[[zz]] <- paramMat_new[,c("true_distr", "fitting_distr", "k", "fitting_param")]
}

for(i in 1:3){
  print(i)
  res_mat <- res_mat_list[[i]]

  # now tabulate the results
  chosen_model_vec <- apply(res_mat, 2, function(x){
    idx <- 1:length(x)
    idx <- idx[which(!is.na(x))]
    x2 <- x[which(!is.na(x))]
    idx[which.min(abs(45 - x2))]
  })

  chosen_model_vec <- factor(as.numeric(chosen_model_vec), levels = 1:24)
  print(table(chosen_model_vec))
  print(paramMat_list[[i]])
}


