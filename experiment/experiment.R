rm(list=ls())
load("../results/step5_clustering_spca.RData")
x <- 1
zz <- construct_3d_tube(our_curves$curves[[x]]$s, radius = our_sd_val$sd_val)

########

x <- 1
dat <- our_curves$curves[[x]]$s
radius <- our_sd_val$sd_val
len = 20

dat <- .remove_duplicate_rows(dat)

n <- nrow(dat)
circle_list <- lapply(1:nrow(dat), function(x){
  direction <- .find_adjacent_directions(dat, x)
  res <- .find_basis_vectors(direction)
  basis_vec1 <- res$vec1; basis_vec2 <- res$vec2

  .construct_3d_circle(dat[x,], radius, basis_vec1, basis_vec2, len = len)
})

which(sapply(circle_list, function(x){any(is.nan(x))}))

x = 45
direction <- .find_adjacent_directions(dat, x)
res <- .find_basis_vectors(direction)
basis_vec1 <- res$vec1; basis_vec2 <- res$vec2

.construct_3d_circle(dat[x,], radius, basis_vec1, basis_vec2, len = len)
