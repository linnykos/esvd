rm(list=ls())
load("../results/step3_factorization_logged.RData")

set.seed(10)
u_mat <- res$u_mat[,1:3]
cluster_labels <- dbscan(u_mat, neighbor_count = 10, upper_cutoff = 14,
                         size_cutoff = 19)

lineages <- .get_lineages(u_mat, cluster_labels, starting_cluster = 1,
                          knn = NA, remove_outlier = T, percentage = 0.05)

#########

shrink = 1; thresh = 0.001; max_iter = 15; b = 1
dat <- u_mat

### setup
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

## plot
col_vec <- cluster_labels
col_vec[which(is.na(col_vec))] <- rgb(0,0,0,0.1)
col_vec[col_vec == 1] <- rgb(0,0,0,0.5)
idx1 <- 1; idx2 <- 2
plot(u_mat[,idx1], u_mat[,idx2], pch = 16, col = col_vec, asp = T)
for(i in 1:length(s_list)){
  lines(s_list[[i]][,idx1], s_list[[i]][,idx2], lwd = 2)
}

######

# diving into .refine_curve_fit
curve_list = s_list
n <- nrow(dat); num_lineage <- length(lineages)
D <- matrix(NA, nrow = n, ncol = num_lineage)
cluster_vec <- 1:ncol(cluster_mat)

pcurve_list <- vector("list", num_lineage)

lin = 1
sample_idx <- .determine_idx_lineage(lineages[[lin]], cluster_mat)

pcurve <- princurve::project_to_curve(dat[sample_idx, ,drop = FALSE],
                                      s = curve_list[[lin]],
                                      stretch = 9999)
# note: princurve::project_to_curve changes the input s
pcurve <- .clean_curve(pcurve, W[, lin], sample_idx)
pcurve_list[[lin]] <- pcurve
