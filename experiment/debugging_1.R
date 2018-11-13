rm(list=ls())
load("../tmp.RData")
i = 1
tmp <- res_list[[i]]$u_mat; svd_res <- svd(tmp); tmp <- svd_res$u[,1:res$k] %*% diag(svd_res$d[1:res$k])
# curves <- slingshot(tmp, cluster_labels = rep(1:4, each = n_seq[i]),
#                     starting_cluster = 1)

###########

lineages <- .get_lineages(tmp, cluster_labels = rep(1:4, each = n_seq[i]),
                          starting_cluster = 1)

# curves <- .get_curves(tmp, cluster_labels = rep(1:4, each = n_seq[i]),
#                       lineages = lineages)

###########

dat = tmp
cluster_labels = rep(1:4, each = n_seq[i])
shrink = 1; thresh = 0.001; max_iter = 15; b = 1
stopifnot(shrink >= 0 & shrink <= 1)

### setup
num_lineage <- length(lineages)
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
while (iter < 2){
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
  if(shrink > 0){
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
    for(i in rev(1:length(avg_curve_list))){
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

      for(j in 1:length(avg_order[[i]])){
        ns <- avg_order[[i]][[j]]
        if(grepl('Lineage', ns)){
          pcurve_list[[ns]] <- shrunk_list[[j]]
        } else {
          avg_curve_list[[ns]] <- shrunk_list[[j]]
        }
      }
    }
  }

  iter <- iter + 1
}

######### # bug happens on iter = 2

# first see what the plot looks like
plot(tmp[,1], tmp[,2],
     pch = 16, col = col_vec[rep(1:4, each = n_seq[1])], asp = T,
     xlab = "X[,1]", ylab = "X[,2]")
for(i in 1:length(pcurve_list)){
  lines(pcurve_list[[i]])
}

####
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
i= 1
to_avg_curves <- lapply(avg_order[[i]], function(n_element){
  if(grepl('Lineage', n_element)){
    pcurve_list[[n_element]]
  } else {
    avg_curve_list[[n_element]]
  }
})

avg <- .construct_average_curve(to_avg_curves, dat)
avg_curve_list[[i]] <- avg

#########

# problem occurs here
# find the indicies shared by all the curves
common_ind <- which(rowMeans(sapply(to_avg_curves, function(crv){ crv$W > 0 })) == 1)
.percent_shrinkage(to_avg_curves[[1]], common_ind)

#########

pcurve = to_avg_curves[[1]]
common_idx = common_ind
lambda <- pcurve$lambda

dens <- stats::density(0, bw = 1, kernel = "cosine")
surv <- list(x = dens$x, y = (sum(dens$y) - cumsum(dens$y))/sum(dens$y))
box_vals <- graphics::boxplot(lambda[common_idx], plot = FALSE)$stats
