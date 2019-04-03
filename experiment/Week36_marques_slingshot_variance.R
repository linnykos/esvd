rm(list=ls())
load("../experiment/Week36_slingshot_bootstrap.RData")
load("../experiment/Week36_marques_slingshot.RData")

# first do a basic analysis counting
naive_lineages_all <- rep(NA, 50*10)
idx <- 1
for(i in 1:length(res_list_naive)){
  for(j in 1:length(res_list_naive[[i]]$lineages)){
    naive_lineages_all[idx] <- paste0(res_list_naive[[i]]$lineages[[j]], collapse = "-")
    idx <- idx+1
  }
}
naive_lineages_all <- naive_lineages_all[!is.na(naive_lineages_all)]
table(naive_lineages_all)
naive_curves$lineages

# first do a basic analysis counting
our_lineages_all <- rep(NA, 50*10)
idx <- 1
for(i in 1:length(res_list_our)){
  for(j in 1:length(res_list_our[[i]]$lineages)){
    our_lineages_all[idx] <- paste0(res_list_our[[i]]$lineages[[j]], collapse = "-")
    idx <- idx+1
  }
}
our_lineages_all <- our_lineages_all[!is.na(our_lineages_all)]
table(our_lineages_all)
our_curves$lineages

###################

# try to compute the "radius" of uncertainty by conditioning on the same curves
# first store all the relevant curves
capture_curves <- function(string, curve_list){
  res <- lapply(1:length(curve_list), function(i){
    string_vec <- sapply(curve_list[[i]]$lineages, function(x){paste0(x, collapse="-")})
    idx <- which(string_vec == string)
    if(length(idx) == 0) return(NA)
    curve_list[[i]]$curves[[idx]]$s[curve_list[[i]]$curves[[idx]]$ord,]
  })

  res[which(sapply(res, length) > 1)]
}

# our_mat_collection <- capture_curves(paste0(our_curves$lineages[[1]], collapse = "-"), res_list_our)
# our_mat <- our_curves$curves[[1]]$s[our_curves$curves[[1]]$ord,]

# for every point in our_mat, find its l2 distance to its closest neighbor in all curves in our_mat_collection
compute_l2_curve <- function(mat, mat_collection){
  n <- nrow(mat); k <- length(mat_collection)
  sapply(1:n, function(x){
    if(x %% floor(n/10) == 0) cat('*')

    vec <- mat[x,]
    sapply(1:k, function(y){
      dist_vec <- apply(mat_collection[[k]], 1, function(z){
        .l2norm(z-vec)
      })
      min(dist_vec)
    })
  })
}

# l2_dist_mat <- compute_l2_curve(our_mat, our_mat_collection)
# median(apply(l2_dist_mat, 2, max))

compute_curve_sd <- function(target_curve_list, bootstrap_curve_list){
  num_curves <- length(target_curve_list$lineages)

  mat_list <- lapply(1:num_curves, function(i){
    print(paste0("Starting curve ", i))
    curve_mat <- target_curve_list$curves[[i]]$s[target_curve_list$curves[[i]]$ord,]
    curve_mat_collection <- capture_curves(paste0(target_curve_list$lineages[[i]], collapse = "-"), bootstrap_curve_list)

    compute_l2_curve(curve_mat, curve_mat_collection)
  })

  sd_vec <- sapply(mat_list, function(x){median(apply(x, 2, max))})

  list(sd_vec = sd_vec, mat_list = mat_list)
}

our_sd <- compute_curve_sd(our_curves, res_list_our)
naive_sd <- compute_curve_sd(naive_curves, res_list_naive)

our_sd_final <- max(sapply(our_sd$mat_list, function(x){
  quantile(apply(x, 2, max), 1)
}))
naive_sd_final <- max(sapply(naive_sd$mat_list, function(x){
  quantile(apply(x, 2, max), 1)
}))

#######

