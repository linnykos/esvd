#' Segment the genes along the trajectory
#'
#' Given two datasets, \code{dat1} and \code{dat2}, whose first \code{common_n} rows
#' are the same, use the \code{eSVD:::.find_highly_expressed_region} function to output
#' statistics of each gene in both datasets.
#'
#' As suggested by the interface, this function can only segment genes when exactly two
#' trajectories are estimated.
#'
#' @param dat1 \code{n1} by \code{p} dataset
#' @param dat2 \code{n2} by \code{p} dataset
#' @param common_n positive integer value denoting that the first \code{common_n} (smaller
#' than \code{nrow(dat1)} or \code{nrow(dat2)}) rows of
#' \code{dat1} and \code{dat2} are equivalent
#' @param standardize boolean
#' @param ncores number of cores
#' @param verbose boolean
#'
#' @return a data frame with 9 columns and \code{ncol(dat1)} rows, where each row
#' contains statistics
#' for each gene across both \code{dat1} and \code{dat2}
#' @export
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

  .extract_information(segmentation_res)
}

#' Ordering the highly expressed genes
#'
#' This uses the output of \code{segment_genes_along_trajectories} to find \code{number_of_genes}
#' highly informative genes along each position of the pseudotime (defined in terms of
#' indices, as \code{1} through \code{common_n + (nrow(dat1) - commom_n) + (nrow(dat2) - common_n)},
#' the arguments into \code{segment_genes_along_trajectories}). The function then order such
#' genes by the midpoint of their significance (as defined in the output of \code{segment_genes_along_trajectories}).
#'
#' @param res_mat the 9-columned data frame that is the output of \code{segment_genes_along_trajectories}
#' @param nrow1 the total number of rows in \code{dat1} when used in \code{segment_genes_along_trajectories}
#' @param nrow2 the total number of rows in \code{dat2} when used in \code{segment_genes_along_trajectories}
#' @param common_n the argument \code{common_n} as used in \code{segment_genes_along_trajectories}
#' @param threshold the threshold on what qualifies as a gene worth considering. See
#' \code{eSVD:::.find_common_genes} or \code{eSVD:::.find_trajectory_genes}
#' @param number_of_genes maximum number of genes (who's significance exceeds the threshold) that are included
#' at each position of pseudotime
#' @param manual_add_common gene indices (from 1 through \code{max(res_mat$idx)}) (possibly \code{NA})
#' that are manually included as highly expressed genes shared in common between the two trajectories
#' @param manual_add_traj1 gene indices (from 1 through \code{max(res_mat$idx)}) (possibly \code{NA})
#' that are manually included as highly expressed genes in trajectory 1
#' @param manual_add_traj2 gene indices (from 1 through \code{max(res_mat$idx)}) (possibly \code{NA})
#' that are manually included as highly expressed genes in trajectory 2
#'
#' @return
#' @export
#'
#' @examples
order_highly_expressed_genes <- function(res_mat, nrow1, nrow2, common_n,
                                         threshold, number_of_genes = 2,
                                         manual_add_common = NA,
                                         manual_add_traj1 = NA,
                                         manual_add_traj2 = NA){
  # find the unique genes in each tail, and then for each location along time, pick the top x genes in that set (if any)
  traj1_genes <- .find_trajectory_genes(res_mat, traj = 1, common_n = common_n, n = nrow1, threshold = threshold,
                                        number_of_genes = number_of_genes, manual_add = manual_add_traj1)
  traj2_genes <- .find_trajectory_genes(res_mat, traj = 2, common_n = common_n, n = nrow2, threshold = threshold,
                                        number_of_genes = number_of_genes, manual_add = manual_add_traj2)

  ## now find the common genes
  common_genes <- .find_common_genes(res_mat, c(traj1_genes, traj2_genes),
                                     common_n = common_n, threshold = threshold,
                                     number_of_genes = number_of_genes, manual_add = manual_add_common)

  list(common_genes = common_genes, traj1_genes = traj1_genes,
       traj2_genes = traj2_genes)
}

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

.extract_information <- function(segmentation_res){
  stopifnot(all(sapply(segmentation_res, length) == 4))

  info_mat <- matrix(NA, nrow = length(segmentation_res), ncol = 8)

  for(k in 1:2){
    info_mat[,(k-1)*4+1] <- sapply(1:length(segmentation_res), function(i){
      segmentation_res[[i]][[k]]$i
      })

    info_mat[,(k-1)*4+2] <- sapply(1:length(segmentation_res), function(i){
      segmentation_res[[i]][[k]]$j
    })

    info_mat[,(k-1)*4+3] <- rowMeans(info_mat[,((k-1)*4)+c(1:2)])

    info_mat[,(k-1)*4+4] <- sapply(1:length(segmentation_res), function(i){
      segmentation_res[[i]][[k]]$obj_val
    })
  }

  info_mat <- cbind(1:nrow(info_mat), info_mat)
  colnames(info_mat) <- c("idx", "start_1", "end_1", "mid_1", "obj_1",
                          "start_2", "end_2", "mid_2", "obj_2")

  data.frame(info_mat)
}

.find_trajectory_genes <- function(res_mat, traj, common_n, n, threshold,
                                   number_of_genes = 2, manual_add = NA){
  stopifnot(ncol(res_mat) == 9, n > common_n)

  start_idx <- ifelse(traj == 1, 2, 6)
  mid_idx <- ifelse(traj == 1, 4, 8)
  end_idx <- ifelse(traj == 1, 3, 7)
  end_idx_other <- ifelse(traj == 1, 7, 3)
  obj_idx <- ifelse(traj == 1, 5, 9)
  obj_idx_other <- ifelse(traj == 1, 9, 5)

  subset_mat <- res_mat[intersect(which(res_mat[,end_idx] >= common_n), which(res_mat[,obj_idx] >= threshold)),]
  rm_idx <- intersect(which(subset_mat[,end_idx_other] >= common_n), which(subset_mat[,obj_idx_other] >= threshold))
  if(length(rm_idx) > 0){
    subset_mat <- subset_mat[-rm_idx,]
  }

  tmp_pos <- sapply((common_n+1):n, function(j){
    tmp_mat <- subset_mat[intersect(which(subset_mat[,start_idx] <= j), which(subset_mat[,end_idx] >= j)),]
    if(nrow(tmp_mat) == 0) return(rep(NA, number_of_genes))
    if(nrow(tmp_mat) < number_of_genes) return(c(tmp_mat$idx, rep(NA, number_of_genes - length(tmp_mat$idx))))
    tmp_mat$idx[order(tmp_mat[,obj_idx], decreasing = T)[1:number_of_genes]]
  })

  idx <- sort(unique(as.numeric(tmp_pos)))
  if(!all(is.na(manual_add))) idx <- sort(unique(c(idx, manual_add)))

  idx[order(sapply(idx, function(x){res_mat[which(res_mat$idx == x)[1], mid_idx]}), decreasing = F)]
}

.find_common_genes <- function(res_mat, traj_genes, common_n, threshold,
                               number_of_genes = 2, manual_add = NA){
  stopifnot(ncol(res_mat) == 9)

  if(length(traj_genes) > 0 && any(res_mat$idx %in% traj_genes)){
    res_mat <- res_mat[-which(res_mat$idx %in% traj_genes),]
  }

  rm_idx <- intersect(which(res_mat$start_1 >= common_n), which(res_mat$start_2 >= common_n))
  if(length(rm_idx) > 0) res_mat <- res_mat[-rm_idx,]

  tmp_pos <- sapply(1:common_n, function(j){
    tmp_mat <- res_mat[intersect(which(res_mat$start_1 <= j), which(res_mat$end_1 >= j)),]
    if(nrow(tmp_mat) == 0) return(rep(NA, number_of_genes))
    if(nrow(tmp_mat) == 1) return(c(tmp_mat$idx, rep(NA, number_of_genes - length(tmp_mat$idx))))
    tmp_mat$idx[order(tmp_mat$obj_1, decreasing = T)[1:number_of_genes]]
  })

  idx <- sort(unique(as.numeric(tmp_pos)))
  if(!all(is.na(manual_add))) idx <- sort(unique(c(idx, manual_add)))

  idx[order(sapply(idx, function(x){res_mat$mid_1[which(res_mat$idx == x)[1]]}), decreasing = F)]
}
