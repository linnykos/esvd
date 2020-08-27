#' Prepare data from segmentation downstrem
#'
#' This function is meant to preprocess \code{dat} so
#' \code{eSVD::segment_genes_along_trajectories} can be used afterwards.
#' Here, \code{curve_list} is the \code{slingshot} object returned by
#' \code{eSVD::slingshot}. This function requires that \code{curve_list} contains exactly two trajectories.
#'
#' Here, \code{min_traj_pseudotime} is a
#' pseudotime (according to \code{curve_list}) of which cells that are specific to one of the trajectories
#' (according to \code{eSVD::construct_pseudotime_trajectory_matrix}) with a smaller pseudotime
#' \code{min_traj_pseudotime} are removed.
#'
#' Similarly, \code{cluster_removal_idx_vec} and \code{cluster_removal_time_vec} are vectors of the same
#' length, of which the cells with a cluster label of the first element in \code{cluster_removal_idx_vec}
#' with a pseudotime smaller than the first element in \code{cluster_removal_time_vec} are also
#' removed. (And then second, and so on.)
#'
#' @param dat \code{n} by \code{p} dataset
#' @param cluster_labels vector of cluster labels, where
#' the cluster labels are consecutive positive integers from 1 to
#' \code{max(cluster_labels)}
#' @param curve_list output from \code{eSVD::slingshot}
#' @param min_traj_pseudotime positive numeric or \code{NA}
#' @param cluster_removal_idx_vec vector of positive integers or a single \code{NA}
#' @param cluster_removal_time_vec vector of positive numerics or a single \code{NA}
#'
#' @return a list containing \code{dat1} and \code{dat2} (two datasets formed from a subset of rows from \code{dat})
#' as well as three vectors of indicies, \code{cell_idx_common}, \code{cell_idx_traj1} and \code{cell_idx_traj2}.
#' Here, \code{dat1} is formed from \code{dat[c(cell_idx_common,cell_idx_traj1),]}
#' and \code{dat2} is formed from \code{dat[c(cell_idx_common,cell_idx_traj2),]}.
#' @export
prepare_data_for_segmentation <- function(dat, cluster_labels, curve_list,
                                          min_traj_pseudotime = NA,
                                          cluster_removal_idx_vec = NA, cluster_removal_time_vec = NA){
  stopifnot(class(curve_list) == "slingshot", length(curve_list$lineages) == 2,
            length(cluster_removal_idx_vec) == length(cluster_removal_time_vec))

  pseudotime_df <- construct_pseudotime_trajectory_matrix(curve_list, cluster_labels)

  # remove some problematic cells
  pseudotime_df2 <- pseudotime_df
  idx <- which(!pseudotime_df2$consensus)
  if(length(idx) > 0) pseudotime_df2 <- pseudotime_df2[-idx,]

  if(!is.na(min_traj_pseudotime)){
    idx <- intersect(which(is.na(pseudotime_df2$consensus)), which(pseudotime_df2$pseudotime <= min_traj_pseudotime))
    if(length(idx) > 0) pseudotime_df2 <- pseudotime_df2[-idx,]
  }

  if(!any(is.na(cluster_removal_idx_vec)) & !any(is.na(cluster_removal_time_vec))){
    for(i in 1:length(cluster_removal_idx_vec)){
      idx <- intersect(which(pseudotime_df2$pseudotime <= cluster_removal_time_vec[i]),
                       which(pseudotime_df2$cluster_labels == cluster_removal_idx_vec[i]))
      if(length(idx) > 0) pseudotime_df2 <- pseudotime_df2[-idx,]
    }
  }

  # determine indices in each trajectory
  traj1_cluster <- sort(setdiff(curve_list$lineages[[1]], curve_list$lineages[[2]]))
  traj2_cluster <- sort(setdiff(curve_list$lineages[[2]], curve_list$lineages[[1]]))
  common_cluster <- intersect(curve_list$lineages[[1]], curve_list$lineages[[2]])

  # organize the pseudotime_df2
  pseudotime_df2 <- pseudotime_df2[order(pseudotime_df2$pseudotime),]
  tmp1 <- pseudotime_df2[which(!pseudotime_df2$cluster_labels %in% traj2_cluster),]
  tmp2 <- pseudotime_df2[which(!pseudotime_df2$cluster_labels %in% traj1_cluster),]

  # determine which cells are in which category
  pseudotime_max_common <- min(pseudotime_df2$pseudotime[pseudotime_df2$cluster_labels %in% c(traj1_cluster, traj2_cluster)])
  cell_idx_common <- pseudotime_df2$cell_idx[which(pseudotime_df2$pseudotime < pseudotime_max_common)]

  cell_idx_traj1 <- pseudotime_df2$cell_idx[intersect(which(pseudotime_df2$pseudotime >= pseudotime_max_common),
                                                      which(pseudotime_df2$cluster_labels %in% c(traj1_cluster, common_cluster)))]
  cell_idx_traj2 <- pseudotime_df2$cell_idx[intersect(which(pseudotime_df2$pseudotime >= pseudotime_max_common),
                                                      which(pseudotime_df2$cluster_labels %in% c(traj2_cluster, common_cluster)))]

  # construct the datasets
  dat1 <- dat[c(cell_idx_common,cell_idx_traj1),]
  dat2 <- dat[c(cell_idx_common,cell_idx_traj2),]

  list(dat1 = dat1, dat2 = dat2, cell_idx_common = cell_idx_common,
       cell_idx_traj1 = cell_idx_traj1, cell_idx_traj2 = cell_idx_traj2)
}


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
#' @param ... additional arguments for \code{eSVD:::.circular_segmentation}
#'
#' @return a data frame with 9 columns and \code{ncol(dat1)} rows, where each row
#' contains statistics
#' for each gene across both \code{dat1} and \code{dat2}
#' @export
segment_genes_along_trajectories <- function(dat1, dat2, common_n, standardize = T,
                                             ncores = NA, verbose = F, ...){
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
                                  standardize = standardize, ...)
  }

  if(is.na(ncores)){
    segmentation_res <- lapply(1:p, func)
  } else {
    j <- 0 #debugging purposes
    segmentation_res <- foreach::"%dopar%"(foreach::foreach(j = 1:p), func(j))
  }

  list(df = .extract_information(segmentation_res), segmentation_fit = segmentation_res)
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
#' @return a list with entries \code{common_genes}, \code{traj1_genes} and \code{traj2_genes},
#' each containing the index of genes (ranging from \code{1} to \code{max(res_mat$idx)}) in order of
#' pseudotime
#' @export
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

#' Find highly specific genes - Pipeline
#'
#' This function takes in three vectors, which is specialized in finding the
#' highly informative (i.e., highly expressive) genes between two trajectories.
#' It first uses \code{eSVD:::.np_smoother} to smooth the concatenation of \code{common_vec}
#' and \code{specific_vec1}, as well as to smooth the concatenation of \code{common_vec}
#' and \code{specific_vec2}, via \code{np::npreg} function (i.e., kernel regression).
#' Then, it standardizes these smoothed outputs jointly (so the joint vector has
#' mean 0 and standard deviation 1) if \code{standardize=T}. Finally, it separates
#' the standardedized joint vector back into two vectors and
#' \code{eSVD:::.circular_segmentation} to segment each of the vector. The only caveat is that
#' this separation keeps the smooths counterparts of \code{common_vec} and \code{specific_vec1}
#' and \code{specific_vec2} together, just in different orders. The reason for doing this
#' is because \code{eSVD:::.circular_segmentation} will have an argument \code{hard_cut},
#' which means to not search for a segmentation index past a certain index. This means
#' even when we're segmenting the first trajectory, we can retain the values from the second
#' trajectory for reference.
#'
#' @param common_vec a vector
#' @param specific_vec1 a vector
#' @param specific_vec2 a vector
#' @param standardize boolean
#' @param ... extra arguments for \code{eSVD:::.circular_segmentation}
#'
#' @return a list
.find_highly_expressed_region <- function(common_vec, specific_vec1, specific_vec2, standardize = T, ...){
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

  res1 <- .circular_segmentation(vec1_all, hard_cut = n+n1, ...)
  res2 <- .circular_segmentation(vec2_all, hard_cut = n+n2, ...)

  list(cut_1 = res1, cut_2 = res2, vec1_smooth = vec1_smooth, vec2_smooth = vec2_smooth)
}

#' Smoother function based on kernel regression
#'
#' The bandwidth selection is done automatically via \code{np::npregbw}.
#'
#' @param vec a vector that will be smoothed
#'
#' @return a vector of the same length as \code{vec}
.np_smoother <- function(vec){
  n <- length(vec)
  dat <- data.frame(y = vec, x = 1:n)
  utils::capture.output(bw_res <- np::npregbw(y ~ x, data = dat))
  res <- np::npreg(bw_res)

  res$mean
}

#' Circular binary segmentation
#'
#' Inspired by Olshen et al., 2004, this method does binary segmentation with a few modifications:
#' First, it only looks for changepoints, spaced out to be \code{resolution*length(vec)} apart. This
#' means it does not consider every index of \code{vec} to be a potential changepoint location.
#' Second, it only looks for starts and ends that are \code{max_width_percentage*length(vec)} apart.
#' This means it does not consider extremely long segments. Third, it has a parameter \code{hard_cut},
#' which is an index between \code{1} and \code{length(vec)}. Fourth, the metric used to segment
#' is not based on the difference on means. Rather, it is based on the difference between
#' the 25th quantile in the middle segment to the 75th quantile of the left and right shoulders.
#' Fifth (related to the fourth change), the change must be in a positive direction (i.e.,
#' the values in the middle segment must be higher than those in the shoulders, roughly speaking).
#'
#' @param vec the numeric vector to segment
#' @param resolution numeric value between 0 and 1 (exclusive)
#' @param max_width_percentage numeric value between 0 and 1 (exclusive), and must be larger than
#' \code{resolution}
#' @param hard_cut index between 1 and \code{length(vec)}, or \code{NA}
#'
#' @returna a list, where \code{i} and \code{j} denote the start and end of the middle segment and
#' \code{obj_val} denotes the objective value achieved
.circular_segmentation <- function(vec, resolution = 1/100, max_width_percentage = 0.1,
                                   hard_cut = NA){

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

  i <- candidate_idx1_vec[which.max(obj_outer[2,])]
  j <- obj_outer[1, which.max(obj_outer[2,])]
  obj_val <- obj_outer[2, which.max(obj_outer[2,])]

  list(i = i, j = j, obj_val = obj_val)
}

#' Extract segmentation information
#'
#' This is used in \code{segment_genes_along_trajectories} to format the results
#' nicely.
#'
#' @param segmentation_res a list of outputs from \code{.find_highly_expressed_region}
#'
#' @return a data frame
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

#' Find trajectory genes
#'
#' Find genes based on \code{res_mat}, for genes in trajectory \code{traj} that
#' are above a pseudotime index of \code{common_n} and below a pseudotime index of \code{n}.
#' These genes are first selected to be associated with an objective value larger than
#' \code{threshold} (according to \code{res_mat}) in the desired trajectory \code{traj}
#' as well as being lower than \code{threshold} in the other trajectory,
#'  and if any exist, the first \code{number_of_genes}
#' are considered at each pseudotime index.
#'
#' (Note: It's possible that this function could've coded without the need for the \code{n} argument...)
#'
#' @param res_mat the 9-columned data frame that is the output of \code{segment_genes_along_trajectories}
#' @param traj the value \code{1} or \code{2}
#' @param common_n positive integer
#' @param n positive integer, larger than \code{common_n}
#' @param threshold non-negative numeric
#' @param number_of_genes positive integer
#' @param manual_add vector of indices with values between \code{1} and \code{max(res_mat$idx)}
#' that are manually included, regardless of whether or not the gene's associated objective value
#' in \code{res_mat} exceeds \code{threshold}
#'
#' @return a vector of unique indices between \code{1} and \code{max(res_mat$idx)}, sorted
#' by pseudotime
.find_trajectory_genes <- function(res_mat, traj, common_n, n, threshold,
                                   number_of_genes = 2, manual_add = NA){
  stopifnot(ncol(res_mat) == 9, n > common_n, traj %in% c(1,2), length(traj) == 1,
            threshold >= 0)

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

#' Find common genes
#'
#' Find genes based on \code{res_mat}, for genes in trajectory \code{traj} that
#' are below a pseudotime index of \code{common_n}, according to both segmentations (i.e.,
#' the \code{start_1} and \code{start_2} columns in \code{res_mat}), and then among
#' such that the both \code{obj_1} and \code{obj_2} columns in \code{res_mat} are above
#' threshold. These genes are then sorted according to \code{start_1} and \code{end_1}
#' (i.e., based on trajectory 1, for an arbitrary choice), so that the first \code{number_of_genes}
#' are selected at each pseudotime index based on their values in \code{obj_1}.
#'
#' @param res_mat the 9-columned data frame that is the output of \code{segment_genes_along_trajectories}
#' @param traj_genes index of genes that are unique to either trajectories that should be removed, possibly left empty
#' @param common_n positive integer
#' @param threshold non-negative numeric
#' @param number_of_genes positive integer
#' @param manual_add vector of indices with values between \code{1} and \code{max(res_mat$idx)}
#' that are manually included, regardless of whether or not the gene's associated objective value
#' in \code{res_mat} exceeds \code{threshold}
#'
#' @return a vector of unique indices between \code{1} and \code{max(res_mat$idx)}, sorted
#' by pseudotime
.find_common_genes <- function(res_mat, traj_genes, common_n, threshold,
                               number_of_genes = 2, manual_add = NA){
  stopifnot(ncol(res_mat) == 9)

  if(length(traj_genes) > 0 && any(res_mat$idx %in% traj_genes)){
    res_mat <- res_mat[-which(res_mat$idx %in% traj_genes),]
  }

  rm_idx <- intersect(which(res_mat$start_1 >= common_n), which(res_mat$start_2 >= common_n))
  if(length(rm_idx) > 0) res_mat <- res_mat[-rm_idx,]
  rm_idx <- unique(c(which(res_mat$obj_1 <= threshold), which(res_mat$obj_2 <= threshold)))
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
