rm(list=ls())
set.seed(10)
cluster_labels <- rep(1:5, each = 20)
dat <- MASS::mvrnorm(100, rep(0, 5), diag(5))
lineages <- .get_lineages(dat, cluster_labels, starting_cluster = 1)
cluster_mat <- .construct_cluster_matrix(cluster_labels)
k <- ncol(cluster_mat)
centers <- .compute_cluster_center(dat, cluster_mat)
W <- .initialize_weight_matrix(cluster_mat, lineages)
cluster_vec <- 1:ncol(cluster_mat)
s_list <- .initial_curve_fit(lineages, cluster_vec, centers)
pcurve_list <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat)$pcurve_list

res <- .construct_average_curve(pcurve_list, dat)

#######################

n <- nrow(pcurve_list[[1]]$s)
p <- ncol(pcurve_list[[1]]$s)
lambdas_all <- unique(unlist(lapply(pcurve_list, function(pcv){pcv$lambda})))
max_shared_lambda <- min(sapply(pcurve_list, function(pcv){max(pcv$lambda)}))
lambdas_all <- sort(lambdas_all[lambdas_all <= max_shared_lambda])

# interpolate all the curves so they're parameterized on the same points
pcurves_dense <- lapply(pcurve_list, function(pcurve){
  sapply(1:p, function(jj){
    stats::approx(pcurve$lambda, pcurve$s[,jj], xout = lambdas_all)$y
  })
})

avg <- sapply(1:p, function(j){
  dim_all <- sapply(1:length(pcurves_dense),function(i){
    pcurves_dense[[i]][,j]
  })
  rowMeans(dim_all)
})
