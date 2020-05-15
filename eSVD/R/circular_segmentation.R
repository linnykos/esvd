.find_highly_expressed_region <- function(common_vec, specific_vec1, specific_vec2, standardize = T){
  n <- length(common_vec)
  n1 <- length(specific_vec1)
  n2 <- length(specific_vec2)

  vec1_smooth <- .np_smoother(c(common_vec, specific_vec1))
  vec2_smooth <- .np_smoother(c(common_vec, specific_vec2))

  vec1_all <- c(vec1_smooth, vec2_smooth[(n+1):length(vec2_smooth)])
  vec2_all <- c(vec2_smooth, vec1_smooth[(n+1):length(vec1_smooth)])

  if(standardize){
    vec1_all <- scale(vec1_all)
    vec2_all <- scale(vec2_all)
  }

  res1 <- .circular_segmentation(vec1_all, hard_cut = n+n1)
  res2 <- .circular_segmentation(vec2_all, hard_cut = n+n2)

  list(cut_1 = res1, cut_2 = res2, vec1_smooth = vec1_smooth, vec2_smooth = vec2_smooth)
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
