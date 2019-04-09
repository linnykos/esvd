rm(list=ls())
set.seed(10)
cluster_labels <- sample(1:10, 200, replace = T)
dat <- MASS::mvrnorm(200, rep(0, 5), diag(5))

# res <- slingshot(dat, cluster_labels, starting_cluster = 1)

########

starting_cluster = 1
cluster_group_list = NA
reduction_percentage = 0.25
shrink = 1
thresh = 0.001
max_iter = 15
upscale_vec = NA
verbose = F

reduction_factor <- max(apply(dat, 2, function(x){diff(range(x))}))*reduction_percentage
dat2 <- dat/reduction_factor

if(verbose) print("Starting to infer lineages")
lineages <- .get_lineages(dat2, cluster_labels, starting_cluster = starting_cluster,
                          cluster_group_list = cluster_group_list)

if(verbose) print("Starting to infer curves")
curves <- .get_curves(dat2, cluster_labels, lineages, shrink = shrink,
                      thresh = thresh, max_iter = max_iter, upscale_vec = upscale_vec,
                      verbose = verbose)

for(k in 1:length(curves)){
  curves[[k]]$s <- curves[[k]]$s*reduction_factor
}
