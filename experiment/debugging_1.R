load("../results/step3_factorization_logged.RData")

k <- 3
u_mat <- res$u_mat[,1:k]
cluster_labels <- dbscan(u_mat, neighbor_count = 10, upper_cutoff = 14,
                         size_cutoff = 19)

#######
dat = u_mat
starting_cluster = 1
b = 1
shrink = 1
knn = NA
remove_outlier = T
percentage = 0.05
thresh = 0.001
max_iter = 15

lineages <- .get_lineages(dat, cluster_labels, starting_cluster = starting_cluster,
                          knn = knn, remove_outlier = remove_outlier,
                          percentage = percentage)

#######
table(cluster_labels)
upscale_vec <- c(.5,1,1, 10,10,1, 10,10)
idx_all <- unlist(lapply(1:8, function(x){
  idx <- which(cluster_labels == x)
  sample(idx, upscale_vec[x]*length(idx), replace = T)
}))
dat <- dat[idx_all,]
cluster_labels <- cluster_labels[idx_all]

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

# do some initial plotting
combn_mat <- utils::combn(3, 2)
mid_vec <- apply(u_mat, 2, function(x){mean(range(x))})[1:3]
rg <- max(apply(u_mat, 2, function(x){diff(range(x))})[1:3])
lim_list <- lapply(1:3, function(x){mid_vec[x]+c(-1,1)*rg/2})
par(mfrow = c(1,3), mar = c(4,4,0.5,0.5))
for(i in 1:ncol(combn_mat)){
  idx1 <- combn_mat[1,i]; idx2 <- combn_mat[2,i]
  plot(u_mat[,idx1], u_mat[,idx2], pch = 16, col = rgb(0.85,0.85,0.85),
       asp = T, cex = 1.3, xlim = lim_list[[idx1]], ylim = lim_list[[idx2]],
       xlab = paste0("Latent dimension ", idx1),
       ylab = paste0("Latent dimension ", idx2))

  for(k in 1:length(pcurve_list)){
    ord <- pcurve_list[[k]]$ord
    lines(pcurve_list[[k]]$s[ord, idx1], pcurve_list[[k]]$s[ord, idx2], lwd = 3.5,
          col = "white")
    lines(pcurve_list[[k]]$s[ord, idx1], pcurve_list[[k]]$s[ord, idx2], lwd = 3,
          col = "black")
    points(centers[,idx1], centers[,idx2], pch = 16, cex = 2)
  }
}

#################################

# do one iteration
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

# do some initial plotting
combn_mat <- utils::combn(3, 2)
mid_vec <- apply(u_mat, 2, function(x){mean(range(x))})[1:3]
rg <- max(apply(u_mat, 2, function(x){diff(range(x))})[1:3])
lim_list <- lapply(1:3, function(x){mid_vec[x]+c(-1,1)*rg/2})
par(mfrow = c(1,3), mar = c(4,4,0.5,0.5))
for(i in 1:ncol(combn_mat)){
  idx1 <- combn_mat[1,i]; idx2 <- combn_mat[2,i]
  plot(u_mat[,idx1], u_mat[,idx2], pch = 16, col = rgb(0.85,0.85,0.85),
       asp = T, cex = 1.3, xlim = lim_list[[idx1]], ylim = lim_list[[idx2]],
       xlab = paste0("Latent dimension ", idx1),
       ylab = paste0("Latent dimension ", idx2))

  for(k in 1:length(pcurve_list)){
    ord <- pcurve_list[[k]]$ord
    lines(pcurve_list[[k]]$s[ord, idx1], pcurve_list[[k]]$s[ord, idx2], lwd = 3.5,
          col = "white")
    lines(pcurve_list[[k]]$s[ord, idx1], pcurve_list[[k]]$s[ord, idx2], lwd = 3,
          col = "black")
    points(centers[,idx1], centers[,idx2], pch = 16, cex = 2)
  }
}
