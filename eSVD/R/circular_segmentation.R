circular_segmentation <- function(vec, resolution = 50){
  stopifnot(n > 5)
  n <- length(vec)
  lim <- round(n/resolution)

  candidate_idx1_vec <- round(seq(2, n-2*lim, length.out = resolution))
  obj_outer <- sapply(candidate_idx1_vec, function(i){
    candidate_idx2_vec <- round(seq(i+lim, n, length.out = resolution))

    obj_inner <- sapply(candidate_idx2_vec, function(j){
      # mean_mid <- mean(vec[(i+1):j]); sd_mid <- stats::sd(vec[(i+1):j])
      # mean_other <- mean(vec[-c((i+1):j)]); sd_other <- stats::sd(vec[-c((i+1):j)])
      #
      # mean_diff <- mean_mid - mean_other
      # if(mean_diff < 0) return(0)
      #
      # mean_diff/(sd_mid + sd_other)

      # quant_mid <- stats::quantile(vec[(i+1):j], probs = 0.25)
      # quant_other <- stats::quantile(vec[-c((i+1):j)], probs = 0.75)
      #
      # quant_mid - quant_other

      mean_mid <- mean(vec[(i+1):j])
      mean_other <- mean(vec[-c((i+1):j)])

      mean_mid - mean_other
    })

    c(j = candidate_idx2_vec[which.max(obj_inner)], obj_val = max(obj_inner))
  })

  i <- candidate_idx1_vec[which.max(obj_outer[2,])]
  j <- obj_outer[1, which.max(obj_outer[2,])]
  obj_val <- obj_outer[2, which.max(obj_outer[2,])]

  list(i = i, j = j, obj_val = obj_val)
}
