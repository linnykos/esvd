rm(list=ls())
load("../results/step4_factorization.RData")

k_select <- 3
u_mat <- res_our$u_mat[,1:k_select]

cluster_labels <- singlecell::dbscan(u_mat, neighbor_count = 10, upper_cutoff = 14,
                                     size_cutoff = 20)

upscale_vec <- max(table(cluster_labels))/table(cluster_labels)
curves <- singlecell::slingshot(u_mat, cluster_labels, starting_cluster = 4, b = 0.5, shrink = 1,
                                upscale_vec = upscale_vec)

###############

starting_cluster = 4
b = 0.5
shrink = 1
upscale_vec = upscale_vec
knn = NA
remove_outlier = T
percentage = 0.05
dat = u_mat
thresh = 0.001
max_iter = 15

cluster_mat <- .construct_cluster_matrix(cluster_labels)

lineages <- .get_lineages(dat, cluster_labels, starting_cluster = starting_cluster,
                          knn = knn, remove_outlier = remove_outlier,
                          percentage = percentage)

# curves <- .get_curves(dat, cluster_labels, lineages, shrink = shrink,
#                       thresh = thresh, max_iter = max_iter, b = b, upscale_vec = upscale_vec)

######################

if(!any(is.na(upscale_vec))){
  idx_all <- unlist(lapply(1:max(cluster_labels, na.rm = T), function(x){
    idx <- which(cluster_labels == x)
    sample(idx, upscale_vec[x]*length(idx), replace = T)
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

### do the shrinking in reverse order
# for(i in rev(1:length(avg_curve_list))){
#   avg_curve <- avg_curve_list[[i]]
#   to_shrink_list <-  lapply(avg_order[[i]], function(n_element){
#     if(grepl('Lineage', n_element)){
#       pcurve_list[[n_element]]
#     } else {
#       avg_curve_list[[n_element]]
#     }
#   })
#
#   shrunk_list <- lapply(1:length(to_shrink_list),function(j){
#     pcurve <- to_shrink_list[[j]]
#     .shrink_to_avg(pcurve, avg_curve,
#                    pct_shrink[[i]][[j]] * shrink, dat)
#   })
#
#   for(j in 1:length(avg_order[[i]])){
#     ns <- avg_order[[i]][[j]]
#     if(grepl('Lineage', ns)){
#       pcurve_list[[ns]] <- shrunk_list[[j]]
#     } else {
#       avg_curve_list[[ns]] <- shrunk_list[[j]]
#     }
#   }
# }
# iter <- iter + 1
#
# ###################
#
# i <- 1
# to_avg_curves <- lapply(avg_order[[i]], function(n_element){
#   if(grepl('Lineage', n_element)){
#     pcurve_list[[n_element]]
#   } else {
#     avg_curve_list[[n_element]]
#   }
# })
#
# avg <- .construct_average_curve(to_avg_curves, dat)
# avg_curve_list[[i]] <- avg
#
# # find the indicies shared by all the curves
# common_ind <- which(rowMeans(sapply(to_avg_curves, function(crv){ crv$W > 0 })) == 1)
# pct_shrink[[i]] <- lapply(to_avg_curves, function(curve) {
#   .percent_shrinkage(curve, common_ind)
# })
#
# pct_shrink[[i]] <- .check_shrinkage(pct_shrink[[i]])
#
# i <- 2
# to_avg_curves <- lapply(avg_order[[i]], function(n_element){
#   if(grepl('Lineage', n_element)){
#     pcurve_list[[n_element]]
#   } else {
#     avg_curve_list[[n_element]]
#   }
# })
#
# avg <- .construct_average_curve(to_avg_curves, dat)
# avg_curve_list[[i]] <- avg
#
# # find the indicies shared by all the curves
# common_ind <- which(rowMeans(sapply(to_avg_curves, function(crv){ crv$W > 0 })) == 1)
# # pct_shrink[[i]] <- lapply(to_avg_curves, function(curve) {
# #   .percent_shrinkage(curve, common_ind)
# # })
# #
# # pct_shrink[[i]] <- .check_shrinkage(pct_shrink[[i]])
#
# ###########################
#
# pcurve = to_avg_curves[[1]]
# lambda <- pcurve$lambda_long
# common_idx = common_ind
#
# dens <- stats::density(0, bw = 1, kernel = "cosine")
# surv <- list(x = dens$x, y = (sum(dens$y) - cumsum(dens$y))/sum(dens$y))
# box_vals <- graphics::boxplot(lambda[common_idx], plot = FALSE)$stats
# if(box_vals[1] == box_vals[5]){
#   pct_l <- rep(0, length(lambda))
# } else {
#   surv$x <- .scale_vector(surv$x, lower = box_vals[1], upper = box_vals[5])
#   pct_l <- stats::approx(surv$x, surv$y, lambda, rule = 2)$y
# }
#
#

#################
i <- 1
avg_curve <- avg_curve_list[[i]]
to_shrink_list <-  lapply(avg_order[[i]], function(n_element){
  if(grepl('Lineage', n_element)){
    pcurve_list[[n_element]]
  } else {
    avg_curve_list[[n_element]]
  }
})

# shrunk_list <- lapply(1:length(to_shrink_list),function(j){
#   print(j)
#   pcurve <- to_shrink_list[[j]]
#   .shrink_to_avg(pcurve, avg_curve,
#                  pct_shrink[[i]][[j]] * shrink, dat)
# })

j = 1
pcurve <- to_shrink_list[[j]]
.shrink_to_avg(pcurve, avg_curve, pct_shrink[[i]][[j]] * shrink, dat)
