rm(list=ls())
set.seed(10)
cell_pop <- matrix(c(4,10, 25,100,
                     60,80, 25,100,
                     40,10, 60,80,
                     60,80, 100,25)/10, nrow = 4, ncol = 4, byrow = T)
h <- nrow(cell_pop)
n_vec <- c(30,40,50,60)
dat <- do.call(rbind, lapply(1:h, function(x){
  pos <- stats::runif(n_vec[x])
  cbind(pos*cell_pop[x,1] + (1-pos)*cell_pop[x,3] + stats::rnorm(n_vec[x], sd = 0.1),
        pos*cell_pop[x,2] + (1-pos)*cell_pop[x,4] + stats::rnorm(n_vec[x], sd = 0.1))
}))
cluster_labels <- unlist(lapply(1:length(n_vec), function(i){rep(i, n_vec[i])}))

cluster_group_list <- NA
starting_cluster <- 1
lineages <- .get_lineages(dat, cluster_labels, starting_cluster = starting_cluster,
                          cluster_group_list = cluster_group_list,
                          squared = F)
upscale_factor <- 1
# res <- .resample_all(dat, cluster_labels, cluster_group_list, lineages, upscale_factor = 1)

cluster_group_list <- list(sort(unique(cluster_labels)))
cluster_intersection <- .intersect_lineages_cluster_group_list(lineages, cluster_group_list)
upscale_vec <- .compute_upscale_factor(cluster_labels, cluster_intersection, upscale_factor)
