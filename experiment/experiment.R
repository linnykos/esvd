rm(list=ls())
set.seed(10)
cell_pop <- matrix(c(4,10, 25,100,
                     40,10, 60,80,
                     60,80, 25,100,
                     60,80, 100,25)/10, nrow = 4, ncol = 4, byrow = T)
h <- nrow(cell_pop)
n_each <- 50
dat <- do.call(rbind, lapply(1:h, function(x){
  pos <- stats::runif(n_each)
  cbind(pos*cell_pop[x,1] + (1-pos)*cell_pop[x,3] + stats::rnorm(n_each, sd = 0.1),
        pos*cell_pop[x,2] + (1-pos)*cell_pop[x,4] + stats::rnorm(n_each, sd = 0.1))
}))
cluster_labels <- rep(1:4, each = 50)
lineages <- .get_lineages(dat, cluster_labels, starting_cluster = 1)

# res <- .get_curves(dat, cluster_labels, lineages)

########
cluster_group_list = NA
use_initialization = F
reduction_percentage = 0.25
shrink = 1
thresh = 0.001
max_iter = 15
upscale_vec =  rep(0.5, 4)
verbose = F

stopifnot(shrink >= 0 & shrink <= 1)

if(!any(is.na(upscale_vec))){
  idx_all <- unlist(lapply(1:max(cluster_labels, na.rm = T), function(x){
    idx <- which(cluster_labels == x)
    sample(idx, round(upscale_vec[x]*length(idx)), replace = T)
  }))
  dat <- dat[idx_all,]
  cluster_labels <- cluster_labels[idx_all]
}

### setup
num_lineage <- length(lineages)
if(any(is.na(cluster_labels))) {
  idx <- which(is.na(cluster_labels))
  dat <- dat[-idx,]
  cluster_labels <- cluster_labels[-idx]
}
cluster_mat <- .construct_cluster_matrix(cluster_labels)
cluster_vec <- 1:ncol(cluster_mat)
centers <- .compute_cluster_center(dat, cluster_mat)

W <- .initialize_weight_matrix(cluster_mat, lineages)

### initial curves are piecewise linear paths through the tree
if(verbose) print("Starting to initialize curves")
s_list <- .initial_curve_fit(lineages, cluster_vec, centers)
res <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat)
pcurve_list <- res$pcurve_list; D <- res$D

if(length(lineages) == 1) {
  s_list <- lapply(1:num_lineage, function(lin){
    sample_idx <- .determine_idx_lineage(lineages[[lin]], cluster_mat)
    .smoother_func(pcurve_list[[lin]]$lambda, dat[sample_idx,,drop = F])
  })

  res <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat)
  pcurve_list <- res$pcurve_list; D <- res$D
  return(pcurve_list)
}

### determine curve hierarchy
avg_order <- .initialize_curve_hierarchy(lineages, cluster_vec)

### track distances between curves and data points to determine convergence
dist_new <- sum(abs(D[W>0]))
dist_old <- Inf

iter <- 1
# while (abs((dist_old - dist_new) >= thresh * dist_old) && iter < max_iter){
  if(verbose) print(paste0("On iteration ", iter))
  dist_old <- dist_new

  ### predict each dimension as a function of lambda (pseudotime)
  s_list <- lapply(1:num_lineage, function(lin){
    sample_idx <- .determine_idx_lineage(lineages[[lin]], cluster_mat)
    .smoother_func(pcurve_list[[lin]]$lambda, dat[sample_idx,,drop = F])
  })

  if(verbose) print("Refining curves")
  res <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat)
  pcurve_list <- res$pcurve_list; D <- res$D
  dist_new <- sum(D[W>0], na.rm=TRUE)

  # shrink together lineages near shared clusters
  if(verbose) print("Shrinking curves together")
  #if(shrink > 0){
    avg_curve_list <- vector("list", length(avg_order))
    names(avg_curve_list) <- paste0("Average", 1:length(avg_order))
    pct_shrink <- vector("list", length(avg_order))

    ### determine average curves and amount of shrinkage
    for (i in 1:length(avg_order)) {
      to_avg_curves <- lapply(avg_order[[i]], function(n_element){
        if(grepl('Lineage', n_element)){
          pcurve_list[[n_element]]
        } else {
          avg_curve_list[[n_element]]
        }
      })

      avg <- .construct_average_curve(to_avg_curves, dat)
      avg_curve_list[[i]] <- avg

      # find the indicies shared by all the curves
      common_ind <- which(rowMeans(sapply(to_avg_curves, function(crv){ crv$W > 0 })) == 1)
      pct_shrink[[i]] <- lapply(to_avg_curves, function(curve) {
        .percent_shrinkage(curve, common_ind)
      })

      pct_shrink[[i]] <- .check_shrinkage(pct_shrink[[i]])
    }

    ### do the shrinking in reverse order
    #for(i in rev(1:length(avg_curve_list))){
    i = 1
      avg_curve <- avg_curve_list[[i]]
      to_shrink_list <-  lapply(avg_order[[i]], function(n_element){
        if(grepl('Lineage', n_element)){
          pcurve_list[[n_element]]
        } else {
          avg_curve_list[[n_element]]
        }
      })

      #shrunk_list <- lapply(1:length(to_shrink_list),function(j){
      j=1
        pcurve <- to_shrink_list[[j]]
        .shrink_to_avg(pcurve, avg_curve,
                       pct_shrink[[i]][[j]] * shrink, dat)
     # })

#       for(j in 1:length(avg_order[[i]])){
#         ns <- avg_order[[i]][[j]]
#         if(grepl('Lineage', ns)){
#           pcurve_list[[ns]] <- shrunk_list[[j]]
#         } else {
#           avg_curve_list[[ns]] <- shrunk_list[[j]]
#         }
#       }
#     # }
#  #  }
#
#   iter <- iter + 1
# # }
