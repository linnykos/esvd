#' Find highly expressed regions in a vector
#'
#' The input vector \code{vec} should already be meaningfully ordered
#'
#' @param vec numeric vector
#' @param resolution value less than \code{length(vec)}
#'
#' @return numeric vector
#' @export
find_highly_expressed_region <- function(vec, resolution = 1/50){
  # first smooth the signal
  vec_smooth <- .np_smoother(vec)

  # then apply circular binary segmentation
  res <- .circular_segmentation(vec_smooth, resolution = resolution)
  res
}

.np_smoother <- function(vec){
  n <- length(vec)
  dat <- data.frame(y = vec, x = 1:n)
  utils::capture.output(bw_res <- np::npregbw(y ~ x, data = dat, verbose = F))
  res <- np::npreg(bw_res)

  res$mean
}

.circular_segmentation <- function(vec, resolution = 1/100, max_width_percentage = 0.1){
  n <- length(vec)
  stopifnot(n > 101)
  lim <- round(n*resolution)
  max_width <- round(n*max_width_percentage)

  candidate_idx1_vec <- round(seq(2, n-2*lim, length.out = lim))
  obj_outer <- sapply(candidate_idx1_vec, function(i){
    candidate_idx2_vec <- round(seq(i+lim, n, length.out = lim))

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
