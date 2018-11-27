rm(list = ls())
set.seed(10)
cluster_labels <- sample(1:10, 200, replace = T)
dat <- MASS::mvrnorm(200, rep(0, 5), diag(5))

# res <- slingshot(dat, cluster_labels, starting_cluster = 1)

cluster_mat <- .construct_cluster_matrix(cluster_labels)

lineages <- .get_lineages(dat, cluster_labels, starting_cluster = 1,
                          knn = NA)

#########

shrink = 1; thresh = 0.001; max_iter = 15; b = 1

num_lineage <- length(lineages)
if(any(is.na(cluster_labels))) cluster_labels <- .fill_in_labels(dat, cluster_labels)
cluster_mat <- .construct_cluster_matrix(cluster_labels)
cluster_vec <- 1:ncol(cluster_mat)
centers <- .compute_cluster_center(dat, cluster_mat)

W <- .initialize_weight_matrix(cluster_mat, lineages)

### determine curve hierarchy
avg_order <- .initialize_curve_hierarchy(lineages, cluster_vec)

### initial curves are piecewise linear paths through the tree
s_list <- .initial_curve_fit(lineages, cluster_vec, centers)
res <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat)
pcurve_list <- res$pcurve_list; D <- res$D

### track distances between curves and data points to determine convergence
dist_new <- sum(abs(D[W>0]))
dist_old <- Inf

iter <- 1

dist_old <- dist_new

### predict each dimension as a function of lambda (pseudotime)
s_list <- lapply(1:num_lineage, function(lin){
  sample_idx <- .determine_idx_lineage(lineages[[lin]], cluster_mat)
  .smoother_func(pcurve_list[[lin]]$lambda, dat[sample_idx,,drop = F], b = b)
})

res <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat)
pcurve_list <- res$pcurve_list; D <- res$D
dist_new <- sum(D[W>0], na.rm=TRUE)

# shrink together lineages near shared clusters

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

i = 3
avg_curve <- avg_curve_list[[i]]
to_shrink_list <-  lapply(avg_order[[i]], function(n_element){
  if(grepl('Lineage', n_element)){
    pcurve_list[[n_element]]
  } else {
    avg_curve_list[[n_element]]
  }
})

shrunk_list <- lapply(1:length(to_shrink_list),function(j){
  pcurve <- to_shrink_list[[j]]
  .shrink_to_avg(pcurve, avg_curve,
                 pct_shrink[[i]][[j]] * shrink, dat)
})

##### # we found the problem. Here, the dist_ind are infinity
## i think the problem is caused by avg_curve_list[[1]]$idx are NA?
## this affects line 425:
##  pcurve <- princurve::project_to_curve(dat[sample_idx,,drop = F],
##    as.matrix(s[pcurve$ord,,drop = FALSE]))
## no points to project

