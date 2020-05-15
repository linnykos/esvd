.find_highly_expressed_region <- function(vec1, idx_trajectory1,
                                          vec2, idx_trajectory2){

  vec1_smooth <- .np_smoother(vec1)
  vec2_smooth <- .np_smoother(vec2)

  total_vec1 <- c(vec1_smooth, vec2_smooth[idx_trajectory2])
  res1 <- .circular_segmentation(total_vec1, hard_cut = length(vec1_smooth))

  total_vec2 <- c(vec2_smooth, vec1_smooth[idx_trajectory1])
  res2 <- .circular_segmentation(total_vec2, hard_cut = length(vec2_smooth))

  list(cut_1 = res1, cut_2 = res2)
}

.np_smoother <- function(vec){
  n <- length(vec)
  dat <- data.frame(y = vec, x = 1:n)
  utils::capture.output(bw_res <- np::npregbw(y ~ x, data = dat))
  res <- np::npreg(bw_res)

  res$mean
}

.circular_segmentation <- function(vec, resolution = 1/100, max_width_percentage = 0.1,
                                   hard_cut = NA){
  n <- length(vec)
  lim <- ifelse(is.na(hard_cut), round(n*resolution), round(hard_cut*resolution))
  max_width <- round(n*max_width_percentage)

  n2 <- ifelse(is.na(hard_cut), n, hard_cut)
  stopifnot(n2 > 101)

  candidate_idx1_vec <- round(seq(2, n2-2*lim, length.out = lim))
  obj_outer <- sapply(candidate_idx1_vec, function(i){
    candidate_idx2_vec <- round(seq(i+lim, n2, length.out = lim))

    obj_inner <- sapply(candidate_idx2_vec, function(j){
      if(abs(i-j) >= max_width) return(-Inf)
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
