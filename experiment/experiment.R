rm(list=ls())
set.seed(10)
cell_pop <- matrix(c(4,10, 25,100,
                     60,80, 25,100,
                     40,10, 60,80,
                     60,80, 100,25)/10, nrow = 4, ncol = 4, byrow = T)
h <- nrow(cell_pop)
n_each <- 25
dat <- do.call(rbind, lapply(1:h, function(x){
  pos <- stats::runif(n_each)
  cbind(pos*cell_pop[x,1] + (1-pos)*cell_pop[x,3] + stats::rnorm(n_each, sd = 0.1),
        pos*cell_pop[x,2] + (1-pos)*cell_pop[x,4] + stats::rnorm(n_each, sd = 0.1))
}))
cluster_labels <- rep(1:4, each = n_each)

set.seed(10)
slingshot_res <- slingshot(dat, cluster_labels, starting_cluster = 1)
lineages <- slingshot_res$lineages
set.seed(10)
trials <- 10
bootstrap_res <- bootstrap_curves(dat, cluster_labels, lineages = lineages, trials = trials)

# res <- compute_curve_sd(slingshot_res$curves, bootstrap_res)

target_curve_list = slingshot_res$curves
bootstrap_curve_list = bootstrap_res
ncores = NA
verbose = F

########

stopifnot(length(target_curve_list) == length(bootstrap_curve_list[[1]]),
          length(unique(sapply(bootstrap_curve_list, length))) == 1)
num_curves <- length(target_curve_list)

if(!is.na(ncores)) doMC::registerDoMC(cores = ncores)

# discretize all the curves
target_curve_list <- lapply(target_curve_list, function(curve){
  .discretize_curve_by_pseudotime(s_mat = curve$s, pseudotime_vec = curve$lambda)
})

for(i in 1:length(bootstrap_curve_list)){
  bootstrap_curve_list[[i]] <- lapply(bootstrap_curve_list[[i]], function(curve){
    .discretize_curve_by_pseudotime(s_mat = curve$s, pseudotime_vec = curve$lambda)
  })
}

# for each curve in the target curve, find the minimum distance
mat_list <- lapply(1:num_curves, function(i){
  if(verbose) print(paste0("Starting curve ", i))
  curve_mat <- target_curve_list[[i]]$s
  curve_mat_collection <- lapply(bootstrap_curve_list, function(curve){curve[[i]]$s})

  .compute_l2_curve(curve_mat, curve_mat_collection, ncores = ncores, verbose = verbose)
})
