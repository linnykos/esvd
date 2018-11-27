rm(list = ls())
set.seed(10)
cluster_labels <- sample(1:10, 200, replace = T)
cluster_labels[sample(1:200, 50)] <- NA
dat <- MASS::mvrnorm(200, rep(0, 5), diag(5))

# res <- slingshot(dat, cluster_labels, starting_cluster = 1)
cluster_mat <- .construct_cluster_matrix(cluster_labels)

lineages <- .get_lineages(dat, cluster_labels, starting_cluster = 1,
                          knn = NA)

#########

shrink = 1; thresh = 0.001; max_iter = 15; b = 1

num_lineage <- length(lineages)
if(any(is.na(cluster_labels))) cluster_labels <- .fill_in_labels(dat, cluster_labels)
