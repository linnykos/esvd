rm(list=ls())
load("../results/step5_trajectory.RData")

target_curve_list = esvd_curves
bootstrap_curve_list = esvd_bootstrap_list
cores = NA
verbose = F

i = 1
curve_mat <- target_curve_list$curves[[i]]$s[target_curve_list$curves[[i]]$ord,]
curve_mat_collection <- .capture_curves(paste0(target_curve_list$lineages[[i]], collapse = "-"), bootstrap_curve_list)

# .compute_l2_curve(curve_mat, curve_mat_collection, cores = cores)

mat = curve_mat
mat_collection = curve_mat_collection

n <- nrow(mat); k <- length(mat_collection)

func <- function(y){
  dist_vec <- apply(mat_collection[[y]], 1, function(z){
    .l2norm(z - vec)
  })

  min(dist_vec)
}

x = 1
vec <- mat[x,]

if(is.na(cores)){
  sapply(1:k, func)
} else {
  y <- 0 #bookkeeping purposes
  unlist(foreach::"%dopar%"(foreach::foreach(y = 1:k), func(y)))
}
