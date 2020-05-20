segment_genes_along_trajectories <- function(dat1, dat2, common_n, standardize = T,
                                             ncores = NA, verbose = F){
  stopifnot(ncol(dat1) == ncol(dat2))
  stopifnot(sum(abs(dat1[1:common_n,] - dat2[1:common_n,])) <= 1e-6)
  n1 <- nrow(dat1) - common_n
  n2 <- nrow(dat2) - common_n
  p <- ncol(dat1)

  if(!is.na(ncores)) doMC::registerDoMC(cores = ncores)

  func <- function(j){
    if(verbose & j %% floor(ncol(dat1)/10) == 0) cat('*')
    .find_highly_expressed_region(common_vec = dat1[1:common_n,j],
                                  specific_vec1 = dat1[(common_n+1):(common_n+n1),j],
                                  specific_vec2 = dat2[(common_n+1):(common_n+n2),j],
                                  standardize = standardize)
  }

  if(is.na(ncores)){
    segmentation_res <- lapply(1:p, func)
  } else {
    j <- 0 #debugging purposes
    segmentation_res <- foreach::"%dopar%"(foreach::foreach(j = 1:p), func(j))
  }

  segmentation_res
}

# order_highly_expressed_genes <- function(segmentation_res){
#   stopifnot(all(sapply(segmentation_res, length) == 4))
#
#   # collect all the necessary ingredients
#   nrow_vec <- max_common_idx+c(length(idx_trajectory1), length(idx_trajectory2))
#   start_vec_list <- vector("list", 2)
#   end_vec_list <- vector("list", 2)
#   midpoint_vec_list <- vector("list", 2)
#   obj_vec_list <- vector("list", 2)
#
#   for(k in 1:2){
#     start_vec_list[[k]] <- sapply(1:length(segmentation_res), function(i){segmentation_res[[i]][[k]]$i})
#     end_vec_list[[k]] <- sapply(1:length(segmentation_res), function(i){segmentation_res[[i]][[k]]$j})
#     midpoint_vec_list[[k]] <- sapply(1:length(start_vec_list[[k]]),
#                                      function(i){(start_vec_list[[k]][i] + end_vec_list[[k]][i])/2})
#     obj_vec_list[[k]] <- sapply(1:length(segmentation_res), function(i){segmentation_res[[i]][[k]]$obj_val})
#
#     plot(NA, ylim = range(obj_vec_list[[k]]), xlim = c(0, nrow_vec[k]),
#          main = paste0(nrow_vec[k]))
#     for(i in 1:length(start_vec_list[[k]])){
#       lines(x = c(start_vec_list[[k]][i], end_vec_list[[k]][i]), y = rep(obj_vec_list[[k]][i], 2), lwd = 2)
#     }
#   }
#
#
#   # first focus on the unique tails (starting with the shorter end)
#
#   ## find the unique genes in each tail, and then for each location along time, pick the top x genes in that set (if any)
#
#   # now deal with the common head
#
#   ## remove any genes that appeared in the tails
#
#   ## for each portion along the sequence, find the top x genes
# }

######################################

.find_highly_expressed_region <- function(common_vec, specific_vec1, specific_vec2, standardize = T){
  n <- length(common_vec)
  n1 <- length(specific_vec1)
  n2 <- length(specific_vec2)

  vec1_smooth <- .np_smoother(c(common_vec, specific_vec1))
  vec2_smooth <- .np_smoother(c(common_vec, specific_vec2))

  vec1_all <- c(vec1_smooth, vec2_smooth[(n+1):length(vec2_smooth)])
  vec2_all <- c(vec2_smooth, vec1_smooth[(n+1):length(vec1_smooth)])

  if(standardize){
    vec_all <- c(vec1_all, vec2_all)
    vec_all <- scale(vec_all)
    vec1_all <- vec_all[1:length(vec1_all)]
    vec2_all <- vec_all[(length(vec1_all)+1):length(vec_all)]
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
      # val_mid <- mean(vec[(i+1):j])
      # val_other <- mean(vec[-c((i+1):j)])
      val_mid <- stats::quantile(vec[(i+1):j], probs = 0.25)
      val_other <- stats::quantile(vec[-c((i+1):j)], probs = 0.75)

      val_mid - val_other
    })

    c(j = candidate_idx2_vec[which.max(obj_inner)], obj_val = max(obj_inner))
  })

  i <- candidate_idx1_vec[which.max(obj_outer[2,])]
  j <- obj_outer[1, which.max(obj_outer[2,])]
  obj_val <- obj_outer[2, which.max(obj_outer[2,])]

  list(i = i, j = j, obj_val = obj_val)
}
