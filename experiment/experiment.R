set.seed(10)
vec <- c(stats::rnorm(100), stats::rnorm(100, mean = 5), stats::rnorm(100))
# res <- .circular_segmentation(vec, max_width_percentage = 1)

resolution = 1/100
max_width_percentage = 1
hard_cut = NA

stopifnot(is.na(hard_cut) || hard_cut <= length(vec))
stopifnot(resolution > 0, resolution < 1, max_width_percentage > 0, max_width_percentage <= 1,
          max_width_percentage > resolution)

n <- length(vec)
lim <- ifelse(is.na(hard_cut), round(n*resolution), round(hard_cut*resolution))
max_width <- round(n*max_width_percentage)

n2 <- ifelse(is.na(hard_cut), n, hard_cut)

candidate_idx1_vec <- round(seq(2, n2-2*lim, by = lim))
obj_outer <- sapply(candidate_idx1_vec, function(i){
  candidate_idx2_vec <- round(seq(i+lim, n2, by = lim))

  obj_inner <- sapply(candidate_idx2_vec, function(j){
    if(abs(i-j) >= max_width) return(-Inf)
    val_mid <- stats::quantile(vec[(i+1):j], probs = 0.25)
    val_other <- stats::quantile(vec[-c((i+1):j)], probs = 0.75)

    val_mid - val_other
  })

  c(j = candidate_idx2_vec[which.max(obj_inner)], obj_val = max(obj_inner))
})
