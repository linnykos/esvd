#' Construst the pseudotime trajectory matrix
#'
#' Given the results in \code{eSVD::slingshot} (in \code{slingshot_res}), construct a data frame that summarizes
#' the information, where each cell is its own row. This function is somewhat complicated
#' since \code{eSVD::slingshot} can reweight clusters (meaning certain cells are over-sampled),
#' so \code{eSVD:::.extract_pseudotimes} needs to re-map all the results in \code{slingshot_res}
#' to unique cells.
#'
#' This function currently works only when \code{slingshot_res} contains exactly two estimated
#' trajectories.
#'
#' @param slingshot_res the output from the function \code{eSVD::slingshot}
#' @param cluster_labels vector of cluster labels, where
#' the cluster labels are consecutive positive integers from 1 to
#' \code{max(cluster_labels)}
#'
#' @return a dataframe 6 columns, labeled
#' \code{cell_idx}, \code{pseudotime}, \code{dist_to_curve} (Euclidean distance to the curve,
#' in the space used for \code{eSVD::slingshot}), \code{consensus} (\code{NA} if the cell is specific
#' to a trajectory, or \code{TRUE} or \code{FALSE} otherwise, depending on whether or not the
#' two trajectories agree on the cell's pseudotime), and
#' \code{status} (\code{"both"} for cells shared between both trajectories, or \code{"1"}
#' or \code{"2"} depending on which trajectory the cell belongs), and \code{cluster_labels}
#' @export
construct_pseudotime_trajectory_matrix <- function(slingshot_res, cluster_labels){
  stopifnot(length(slingshot_res$curves) == 2) #only works currently for 2 branching lineages
  len <- length(slingshot_res$curves)

  # first extract each trajectories set of pseudotimes
  pseudotime_df_list <- .extract_pseudotimes(slingshot_res)

  # next, focus on cells shared between both trajectories
  shared_df <- .compile_common_cells(pseudotime_df_list)
  shared_df$status <- "both"

  # deal with other cells
  specific_df <- .compile_unique_cells(pseudotime_df_list)
  for(i in 1:length(specific_df)){
    specific_df[[i]]$consensus <- NA
    specific_df[[i]]$status <- paste0(i)
  }

  # merge the data frames
  tmp <- do.call(rbind, specific_df)
  all_df <- rbind(shared_df, tmp)

  # all cluster labels
  all_df$cluster_labels <- cluster_labels[all_df$cell_idx]
  all_df <- all_df[order(all_df$cell_idx),]

  all_df
}

#' Extract which cells are common and which cells are trajectory-specific
#'
#' Using the output of \code{eSVD::construct_pseudotime_trajectory_matrix} (called
#' \code{pseudotime_df}) as well
#' as which clusters in the \code{cluster_labels} column of \code{pseudotime_df}
#' correspond to which trajectories, construct 3 vectors of indices (representing different cells) --
#' one for which cells are common (shared among the two trajectories), and one for each
#' of the two trajectories. Within each vectors, the output indices are ordered according to the
#' corresponding cells' pseudotimes.
#'
#' @param pseudotime_df the output from the function \code{eSVD::construct_pseudotime_trajectory_matrix}
#' @param trajectory_1_clusters vector containing integers from 1 to \code{max(pseudotime_df$cluster_labels)}
#' @param trajectory_2_clusters vector containing integers from 1 to \code{max(pseudotime_df$cluster_labels)}, and
#' does not have an intersection with \code{trajectory_1_clusters}
#' @param threshold_common_time numeric or can be \code{NA}, in which case it uses the maximum pseudotime
#' for a cell not in either trajectory
#'
#' @return a list of integers
#' @export
partition_cells_using_pseudotime <- function(pseudotime_df, trajectory_1_clusters,
                                             trajectory_2_clusters, threshold_common_time = NA){
  k <- max(pseudotime_df$cluster_labels)
  stopifnot(all(trajectory_1_clusters %in% c(1:k)), all(trajectory_2_clusters %in% c(1:k)),
            length(intersect(trajectory_1_clusters, trajectory_2_clusters)) == 0)

  if(is.na(threshold_common_time)){
    threshold_common_time <- max(pseudotime_df$pseudotime[!pseudotime_df$cluster_labels %in% c(trajectory_1_clusters, trajectory_2_clusters)])
  }

  # remove trajectory-specific cells that are too young
  pseudotime_df2 <- pseudotime_df[-intersect(which(is.na(pseudotime_df$consensus)),
                                             which(pseudotime_df$pseudotime <= threshold_common_time)),]
  # remove common cells that do not have a consensus in the pseudo-time
  if(length(which(!pseudotime_df2$consensus)) > 0){
    pseudotime_df2 <- pseudotime_df2[-which(!pseudotime_df2$consensus),]
  }

  common_cluster <- setdiff(c(1:k), c(trajectory_1_clusters, trajectory_2_clusters))

  # construct two data frames, each with all the cells (common and specific) for each respecitive trajectory
  pseudotime_df2 <- pseudotime_df2[order(pseudotime_df2$pseudotime),]
  pseudotime_max_common <- min(pseudotime_df2$pseudotime[pseudotime_df2$cluster_labels %in% c(trajectory_1_clusters, trajectory_2_clusters)])

  # define the relevant cell categories

  cell_idx_common <- pseudotime_df2$cell_idx[intersect(which(pseudotime_df2$pseudotime <= pseudotime_max_common),
                                                       which(!pseudotime_df2$cluster_labels %in% c(trajectory_1_clusters, trajectory_2_clusters)))]

  cell_idx_traj1 <- pseudotime_df2$cell_idx[intersect(which(pseudotime_df2$pseudotime >= pseudotime_max_common),
                                                      which(pseudotime_df2$cluster_labels %in% c(trajectory_1_clusters, common_cluster)))]
  cell_idx_traj2 <- pseudotime_df2$cell_idx[intersect(which(pseudotime_df2$pseudotime >= pseudotime_max_common),
                                                      which(pseudotime_df2$cluster_labels %in% c(trajectory_2_clusters, common_cluster)))]

  list(cell_idx_common = cell_idx_common, cell_idx_traj1 = cell_idx_traj1, cell_idx_traj2 = cell_idx_traj2)
}

####################################

#' Extract pseudotimes
#'
#' Maps the (possibly) duplicated cells in \code{eSVD::slingshot} due to its
#' \code{upscale_factor} argument and extract the relevant information about each cell.
#' When a cell is detected to have multiple instances in \code{slingshot_res}, this functions
#' outputs its mean pseudotime as well as mean distance to curve among all its instances.
#'
#' @param slingshot_res the output from the function \code{eSVD::slingshot}
#'
#' @return a list of dataframes, one for each trajectory in \code{slingshot_res}
.extract_pseudotimes <- function(slingshot_res){
  n <- length(slingshot_res$curves)

  df_list <- lapply(1:n, function(curve_idx){
    # which ghost points are in the curve
    idx_all <- slingshot_res$idx[slingshot_res$curves[[curve_idx]]$idx]

    # all the pseudotimes
    lambda_all <- slingshot_res$curves[[curve_idx]]$lambda
    dist_all <- slingshot_res$curves[[curve_idx]]$dist_ind
    stopifnot(length(idx_all) == length(lambda_all))

    # which cells do these correspond to
    idx_cell <- sort(unique(idx_all))
    lambda_vec <- rep(NA, length(idx_cell))
    dist_vec <- rep(NA, length(idx_cell))

    for(i in 1:length(idx_cell)){
      tmp_idx <- which(idx_all == idx_cell[i]) # idx of the ghost points
      lambda_vec[i] <- mean(lambda_all[tmp_idx])
      dist_vec[i] <- mean(dist_all[tmp_idx])
    }

    data.frame(cell_idx = idx_cell, pseudotime = lambda_vec, dist_to_curve = dist_vec)
  })
}

#' Compile the common cells
#'
#' Uses the output of \code{eSVD:::.extract_pseudotimes} (called \code{df_list})
#' and compiles all the information into
#' a data frame. Importantly, the maximum pseudotime in \code{df_list} across both data frames
#' divided by 10 is the tolerance, and if a cell's pseudotimes across both trajectories
#' differs more than than this tolerance, then its \code{consensus} in the output data frame
#' is \code{FALSE}. Additionally, when determining which pseudotime and distance-to-curve
#' to extract, this function uses the trajectory that the cell is closer to.
#'
#' This function is hard-coded so \code{df_list} has to be length 2 (meaning there are two
#' estimated trajectories).
#'
#' @param df_list output of \code{eSVD:::.extract_pseudotimes}
#'
#' @return a dataframe with four columns, \code{cell_idx}, \code{pseudotime},
#' \code{dist_to_curve} and \code{consensus}
.compile_common_cells <- function(df_list){
  stopifnot(length(df_list) == 2) # hard-coded

  idx <- intersect(df_list[[1]]$cell_idx, df_list[[2]]$cell_idx)
  tol_val <- max(c(df_list[[1]]$pseudotime, df_list[[2]]$pseudotime))/10
  lambda_vec <- rep(NA, length(idx))
  bool_vec <- rep(NA, length(idx))
  dist_vec <- rep(NA, length(idx))

  for(i in 1:length(idx)){
    idx1 <- which(df_list[[1]]$cell_idx == idx[i])[1]
    idx2 <- which(df_list[[2]]$cell_idx == idx[i])[1]

    lambda1 <- df_list[[1]]$pseudotime[idx1]
    lambda2 <- df_list[[2]]$pseudotime[idx2]

    dist1 <- df_list[[1]]$dist_to_curve[idx1]
    dist2 <- df_list[[2]]$dist_to_curve[idx2]

    if(dist1 < dist2){
      lambda_vec[i] <- lambda1
      dist_vec[i] <- dist1
    } else {
      lambda_vec[i] <- lambda2
      dist_vec[i] <- dist2
    }

    bool_vec[i] <- (abs(lambda1 - lambda2) <= tol_val)
  }

  data.frame(cell_idx = idx, pseudotime = lambda_vec, dist_to_curve = dist_vec, consensus = bool_vec)
}

#' Compile the unique cells
#'
#' Compared to \code{eSVD:::.compile_common_cells}, this function is a lot
#' simpler since there is nothing to compare between -- only extraction is needed
#'
#' @param df_list output of \code{eSVD:::.extract_pseudotimes}
#'
#' @return a dataframe with three columns, \code{cell_idx}, \code{pseudotime},
#' \code{dist_to_curve}
.compile_unique_cells <- function(df_list){
  stopifnot(length(df_list) == 2) # hard-coded
  len <- length(df_list)

  lapply(1:len, function(i){
    other_branch <- ifelse(i == 1, 2, 1)
    specific_idx <- which(!df_list[[i]]$cell_idx %in%
                            df_list[[other_branch]]$cell_idx)
    df_list[[i]][specific_idx, ]
  })
}
