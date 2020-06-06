context("Test slingshot curves")

## .intersect_lineages_cluster_group_list is correct

test_that(".intersect_lineages_cluster_group_list works", {
  lineages <- list(c(1:5), c(1:3,6:8))
  cluster_group_list <- list(1, 2:4, 5:6, 7:8)

  res <- .intersect_lineages_cluster_group_list(lineages, cluster_group_list)
  expect_true(length(res) == 6)
  expect_true(any(sapply(res, function(x){all(sort(x) == 1)})))
  expect_true(any(sapply(res, function(x){all(sort(x) == c(2:3))})))
  expect_true(any(sapply(res, function(x){all(sort(x) == c(4))})))
  expect_true(any(sapply(res, function(x){all(sort(x) == c(5))})))
  expect_true(any(sapply(res, function(x){all(sort(x) == c(6))})))
  expect_true(any(sapply(res, function(x){all(sort(x) == c(7:8))})))
})

test_that(".intersect_lineages_cluster_group_list with a trivial cluster_group_list", {
  lineages <- list(c(1:5), c(1:3,6:8))
  cluster_group_list <- list(1:8)

  res <- .intersect_lineages_cluster_group_list(lineages, cluster_group_list)
  expect_true(length(res) == 3)
  expect_true(any(sapply(res, function(x){all(x %in% c(1:3)) & all(c(1:3) %in% x)})))
  expect_true(any(sapply(res, function(x){all(x %in%c(4:5)) & all(c(4:5) %in% x)})))
  expect_true(any(sapply(res, function(x){all(x %in% c(6:8)) & all(c(6:8) %in% x)})))
})

###########################################

## .flatten_list is correct

test_that(".flatten_list works", {
  lis <- list(list(c(1:5), c(2:4)), list(c(1:2), c(3:6), c(1)), list(4))
  res <- .flatten_list(lis)

  expect_true(length(res) == 6)
  expect_equal(res, list(c(1:5), c(2:4), c(1:2), c(3:6), c(1), 4))
})

test_that(".flatten_list must take in list of lists", {
  lis <- list(list(c(1:5), c(2:4)), list(c(1:2), c(3:6), c(1)), "asdf")
  expect_error(.flatten_list(lis))
})

###########################################

## .resample_all is correct

test_that(".resample_all works", {
  set.seed(10)
  dat <- matrix(1:500, nrow = 100, ncol = 5)
  cluster_labels <- rep(1:10, each = 10)
  cluster_group_list <- list(1, 2:4, 5:8, 9:10)
  lineages <- list(1:6, c(1:4,7:10))

  res <- .resample_all(dat, cluster_labels, cluster_group_list, lineages, upscale_factor = 1)

  expect_true(length(res) == 3)
  expect_true(all(sort(names(res)) == c("cluster_labels", "dat", "idx_all")))
})

test_that(".resample_all works with cluster_group_list set to NA", {
  set.seed(10)
  dat <- matrix(1:500, nrow = 100, ncol = 5)
  cluster_labels <- rep(1:10, each = 10)
  cluster_group_list <- list(1, 2:4, 5:8, 9:10)
  lineages <- list(1:6, c(1:4,7:10))

  res <- .resample_all(dat, cluster_labels, cluster_group_list = NA, lineages, upscale_factor = 1)

  expect_true(length(res) == 3)
  expect_true(all(sort(names(res)) == c("cluster_labels", "dat", "idx_all")))
})

###########################################

## .construct_resample_idx is correct

test_that(".construct_resample_idx works", {
  set.seed(10)
  cluster_labels <- rep(1:6, each = 20)
  upscale_vec <- c(1,1,1.5,1.5,2,2)

  res <- .construct_resample_idx(cluster_labels, upscale_vec)

  expect_true(length(res) >= length(cluster_labels))
  expect_true(all(res <= length(cluster_labels)))
  expect_true(max(table(res)) > 1)
})

test_that(".construct_resample_idx actually respects the upscale_vec", {
  set.seed(10)
  cluster_labels <- rep(1:6, each = 20)
  upscale_vec <- c(1,1,1.5,1.5,2,2)

  res <- .construct_resample_idx(cluster_labels, upscale_vec)
  res_tab <- table(cluster_labels[res])

  expect_true(all(res_tab == c(20, 20, 30, 30, 40, 40)))
})

test_that(".construct_resample_idx retains all the indices", {
  trials <- 10

  bool_vec <- sapply(1:trials, function(i){
    set.seed(i)
    cluster_labels <- sample(1:6, 100, replace = T)
    upscale_vec <- runif(6, min = 1, max = 2)

    res <- .construct_resample_idx(cluster_labels, upscale_vec)

    all(1:100 %in% res)
  })

  expect_true(all(bool_vec))
})

###########################################

## .compute_upscale_factor is correct

test_that(".compute_upscale_factor works", {
  cluster_labels <- rep(1:6, each = 20)
  cluster_intersection <- list(1, 2:4, 5:6)
  res <- .compute_upscale_factor(cluster_labels, cluster_intersection, upscale_factor = 0.5)

  expect_true(all(is.numeric(res)))
  expect_true(length(res) == 6)
  expect_true(all(res >= 0))
})

test_that(".compute_upscale_factor is 1 for the largest group size", {
  set.seed(10)
  cluster_labels <- sample(1:6, 100, replace = T)
  cluster_intersection <- list(1, 2:4, 5:6)
  res <- .compute_upscale_factor(cluster_labels, cluster_intersection, upscale_factor = 0.5)

  len_vec <- sapply(cluster_intersection, function(x){length(which(cluster_labels %in% x))})
  expect_true(all(res[cluster_intersection[[which.max(len_vec)]]] == 1))
})

test_that(".compute_upscale_factor yields the same value within a cluster", {
  set.seed(10)
  cluster_labels <- sample(1:6, 100, replace = T)
  cluster_intersection <- list(c(1,6), c(2,5), c(3,4))
  res <- .compute_upscale_factor(cluster_labels, cluster_intersection, upscale_factor = 0.5)

  expect_true(abs(res[1] - res[6]) <= 1e-6)
  expect_true(abs(res[2] - res[5]) <= 1e-6)
  expect_true(abs(res[3] - res[4]) <= 1e-6)
})

test_that(".compute_upscale_factor with upscale factor of 1 yields equal group sizes", {
  cluster_labels <- rep(1:6, each = 20)
  cluster_intersection <- list(1, 2:4, 5:6)
  res <- .compute_upscale_factor(cluster_labels, cluster_intersection, upscale_factor = 1)

  effective_size <- rep(NA, length(cluster_intersection))
  for(i in 1:length(cluster_intersection)){
    effective_size[i] <- unique(res[cluster_intersection[[i]]])[1]*length(which(cluster_labels %in% cluster_intersection[[i]]))
  }

  expect_true(length(unique(effective_size)) == 1)
})

###########################################

## .initialize_weight_matrix is correct

test_that(".initialize_weight_matrix works", {
  set.seed(10)
  cluster_labels <- sample(1:10, 200, replace = T)
  dat <- MASS::mvrnorm(200, rep(0, 5), diag(5))
  lineages <- .get_lineages(dat, cluster_labels, starting_cluster = 1)
  cluster_mat <- .construct_cluster_matrix(cluster_labels)

  res <- .initialize_weight_matrix(cluster_mat, lineages)

  expect_true(is.matrix(res))
  expect_true(all(dim(res) == c(200, length(lineages))))
  expect_true(all(sort(unique(as.numeric(res)))== c(0,1)))
})

###############

## .initialize_curve_hierarchy is correct

test_that(".initialize_curve_hierarchy works", {
  lineages <- list(Lineage1 = c(1,2,3), Lineage2 = c(1,2,4))
  cluster_vec <- 1:4

  res <- .initialize_curve_hierarchy(lineages, cluster_vec)

  expect_true(is.list(res))
})

test_that(".initialize_curve_hierarchy works for a more complicated case", {
  set.seed(10)
  cluster_labels <- sample(1:10, 200, replace = T)
  dat <- MASS::mvrnorm(200, rep(0, 5), diag(5))
  lineages <- .get_lineages(dat, cluster_labels, starting_cluster = 1)
  cluster_vec <- 1:10

  res <- .initialize_curve_hierarchy(lineages, cluster_vec)

  expect_true(is.list(res))
})

#################

## .initial_curve_fit is correct

test_that(".initial_curve_fit works", {
  set.seed(10)
  cluster_labels <- rep(1:5, each = 20)
  dat <- MASS::mvrnorm(100, rep(0, 5), diag(5))
  lineages <- .get_lineages(dat, cluster_labels, starting_cluster = 1)
  cluster_mat <- .construct_cluster_matrix(cluster_labels)
  k <- ncol(cluster_mat)
  centers <- .compute_cluster_center(dat, cluster_mat)
  W <- .initialize_weight_matrix(cluster_mat, lineages)
  cluster_vec <- 1:ncol(cluster_mat)
  s_list <- .initial_curve_fit(lineages, cluster_vec, centers)

  res <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat)

  expect_true(is.list(res))
  expect_true(all(names(res) == c("pcurve_list", "D")))
  expect_true(all(sapply(res$pcurve_list, class) == "principal_curve"))
  expect_true(all(dim(res$D) == c(100, length(lineages))))
  expect_true(all(names(res$pcurve_list) == paste0("Lineage", 1:length(res$pcurve_list))))
})

########################

## .smoother_func is correct

test_that(".smoother_func works", {
  set.seed(10)
  cluster_labels <- rep(1:5, each = 20)
  dat <- MASS::mvrnorm(100, rep(0, 5), diag(5))
  lineages <- .get_lineages(dat, cluster_labels, starting_cluster = 1)
  cluster_mat <- .construct_cluster_matrix(cluster_labels)
  k <- ncol(cluster_mat)
  centers <- .compute_cluster_center(dat, cluster_mat)
  W <- .initialize_weight_matrix(cluster_mat, lineages)
  cluster_vec <- 1:ncol(cluster_mat)
  s_list <- .initial_curve_fit(lineages, cluster_vec, centers)
  pcurve_list <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat)$pcurve_list

  sample_idx <- .determine_idx_lineage(lineages[[1]], cluster_mat)
  res <- .smoother_func(pcurve_list[[1]]$lambda, dat[sample_idx,])

  expect_true(is.matrix(res))
  expect_true(all(dim(res) == dim(dat[sample_idx,])))
})

########################

## .construct_average_curve is correct

test_that(".construct_average_curve works", {
  set.seed(10)
  cluster_labels <- rep(1:5, each = 20)
  dat <- MASS::mvrnorm(100, rep(0, 5), diag(5))
  lineages <- .get_lineages(dat, cluster_labels, starting_cluster = 1)
  cluster_mat <- .construct_cluster_matrix(cluster_labels)
  k <- ncol(cluster_mat)
  centers <- .compute_cluster_center(dat, cluster_mat)
  W <- .initialize_weight_matrix(cluster_mat, lineages)
  cluster_vec <- 1:ncol(cluster_mat)
  s_list <- .initial_curve_fit(lineages, cluster_vec, centers)
  pcurve_list <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat)$pcurve_list

  res <- .construct_average_curve(pcurve_list, dat)

  expect_true(class(res) == "principal_curve")
  expect_true(all(!is.na(res$idx)))
})

###############################

## .refine_curve_fit is correct

test_that(".refine_curve_fit works", {
  set.seed(10)
  cluster_labels <- rep(1:5, each = 20)
  dat <- MASS::mvrnorm(100, rep(0, 5), diag(5))
  lineages <- .get_lineages(dat, cluster_labels, starting_cluster = 1)
  cluster_mat <- .construct_cluster_matrix(cluster_labels)
  k <- ncol(cluster_mat)
  centers <- .compute_cluster_center(dat, cluster_mat)
  W <- .initialize_weight_matrix(cluster_mat, lineages)
  cluster_vec <- 1:ncol(cluster_mat)
  s_list <- .initial_curve_fit(lineages, cluster_vec, centers)
  res <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat)

  expect_true(is.list(res))
  expect_true(length(res) == 2)
  expect_true(all(names(res) == c("pcurve_list", "D")))
  expect_true(all(dim(res$D) == c(100, length(lineages))))
  expect_true(is.list(res$pcurve_list))
  expect_true(all(sapply(res$pcurve_list, class) == "principal_curve"))
})

##############################

## .percent_shrinkage is correct

test_that(".percent_shrinkage works", {
  set.seed(10)
  cluster_labels <- rep(1:5, each = 20)
  dat <- MASS::mvrnorm(100, rep(0, 5), diag(5))
  lineages <- .get_lineages(dat, cluster_labels, starting_cluster = 1)
  cluster_mat <- .construct_cluster_matrix(cluster_labels)
  k <- ncol(cluster_mat)
  centers <- .compute_cluster_center(dat, cluster_mat)
  W <- .initialize_weight_matrix(cluster_mat, lineages)
  cluster_vec <- 1:ncol(cluster_mat)
  s_list <- .initial_curve_fit(lineages, cluster_vec, centers)
  pcurve <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat)$pcurve_list[[1]]
  common_idx <- which(W[,1] == 1)

  res <- .percent_shrinkage(pcurve, common_idx)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == length(pcurve$lambda))
})


test_that(".percent_shrinkage computes a shrinkage for all indices, even outliers ", {
  set.seed(10)
  cluster_labels <- rep(1:5, each = 20)
  dat <- MASS::mvrnorm(100, rep(0, 5), diag(5))
  lineages <- .get_lineages(dat, cluster_labels, starting_cluster = 1)
  cluster_mat <- .construct_cluster_matrix(cluster_labels)
  k <- ncol(cluster_mat)
  centers <- .compute_cluster_center(dat, cluster_mat)
  W <- .initialize_weight_matrix(cluster_mat, lineages)
  cluster_vec <- 1:ncol(cluster_mat)
  s_list <- .initial_curve_fit(lineages, cluster_vec, centers)
  pcurve <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat)$pcurve_list[[1]]
  pcurve$lambda[1:5] <- 10*max(pcurve$lambda)
  common_idx <- which(W[,1] == 1)

  res <- .percent_shrinkage(pcurve, common_idx)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == length(pcurve$lambda))
  expect_true(all(res[1:5] == 0))
})

test_that(".percent_shrinkage respects the order of lambda", {
  set.seed(20)
  cluster_labels <- rep(1:5, each = 20)
  dat <- MASS::mvrnorm(100, rep(0, 5), diag(5))
  lineages <- .get_lineages(dat, cluster_labels, starting_cluster = 1)
  cluster_mat <- .construct_cluster_matrix(cluster_labels)
  k <- ncol(cluster_mat)
  centers <- .compute_cluster_center(dat, cluster_mat)
  W <- .initialize_weight_matrix(cluster_mat, lineages)
  cluster_vec <- 1:ncol(cluster_mat)
  s_list <- .initial_curve_fit(lineages, cluster_vec, centers)
  pcurve <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat)$pcurve_list[[1]]
  common_idx <- which(W[,1] == 1)

  res <- .percent_shrinkage(pcurve, common_idx)

  res <- res[order(pcurve$lambda, decreasing = F)]
  vec <- diff(res)
  expect_true(all(vec <= 0))
})

###########################

## .shrink_to_avg is correct

test_that(".shrink_to_avg works", {
  set.seed(20)
  cluster_labels <- rep(1:5, each = 20)
  dat <- MASS::mvrnorm(100, rep(0, 5), diag(5))
  lineages <- .get_lineages(dat, cluster_labels, starting_cluster = 1)
  cluster_mat <- .construct_cluster_matrix(cluster_labels)
  k <- ncol(cluster_mat)
  centers <- .compute_cluster_center(dat, cluster_mat)
  W <- .initialize_weight_matrix(cluster_mat, lineages)
  cluster_vec <- 1:ncol(cluster_mat)
  s_list <- .initial_curve_fit(lineages, cluster_vec, centers)


  pcurve_list <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat)$pcurve_list
  avg_curve <- .construct_average_curve(pcurve_list, dat)

  pcurve <- pcurve_list[[1]]
  common_idx <- which(W[,1] == 1)
  pct_shrink <- .percent_shrinkage(pcurve, common_idx)

  res <- .shrink_to_avg(pcurve, avg_curve, pct_shrink, dat)

  expect_true(class(res) == "principal_curve")
  expect_true(all(res$W == pcurve$W))
})

######################

## .get_curves is correct

test_that(".get_curves works", {
  set.seed(10)
  cell_pop <- matrix(c(4,10, 25,100,
                       40,10, 60,80,
                       60,80, 25,100,
                       60,80, 100,25)/10, nrow = 4, ncol = 4, byrow = T)
  h <- nrow(cell_pop)
  n_each <- 50
  dat <- do.call(rbind, lapply(1:h, function(x){
    pos <- stats::runif(n_each)
    cbind(pos*cell_pop[x,1] + (1-pos)*cell_pop[x,3] + stats::rnorm(n_each, sd = 0.1),
          pos*cell_pop[x,2] + (1-pos)*cell_pop[x,4] + stats::rnorm(n_each, sd = 0.1))
  }))
  cluster_labels <- rep(1:4, each = 50)
  lineages <- .get_lineages(dat, cluster_labels, starting_cluster = 1)

  res <- .get_curves(dat, cluster_labels, lineages)

  expect_true(length(res$pcurve_list) == length(lineages))
  expect_true(all(sapply(res$pcurve_list, class) == "principal_curve"))

  #plot(dat[,1], dat[,2], asp = T)
  #lines(res$pcurve_list$Curve1, col = "red", lwd = 2); lines(res$pcurve_list$Curve2, col = "blue", lwd = 2)
})

test_that(".get_curves works for a harder example", {
  set.seed(20)
  cluster_labels <- rep(1:5, each = 20)
  dat <- MASS::mvrnorm(100, rep(0, 2), diag(2))
  lineages <- .get_lineages(dat, cluster_labels, starting_cluster = 1)

  res <- .get_curves(dat, cluster_labels, lineages)

  expect_true(length(res$pcurve_list) == length(lineages))
  expect_true(all(sapply(res$pcurve_list, class) == "principal_curve"))

  #plot(dat[,1], dat[,2], asp = T)
  #for(i in 1:length(res)){lines(res$pcurve_list[[i]], col = i, lwd = 2)}
})

test_that(".get_curves finds reasonable curves", {
  set.seed(10)
  cell_pop <- matrix(c(4,10, 25,100,
                       60,80, 25,100,
                       40,10, 60,80,
                       60,80, 100,25)/10, nrow = 4, ncol = 4, byrow = T)
  h <- nrow(cell_pop)
  n_each <- 50
  dat <- do.call(rbind, lapply(1:h, function(x){
    pos <- stats::runif(n_each)
    cbind(pos*cell_pop[x,1] + (1-pos)*cell_pop[x,3] + stats::rnorm(n_each, sd = 0.1),
          pos*cell_pop[x,2] + (1-pos)*cell_pop[x,4] + stats::rnorm(n_each, sd = 0.1))
  }))
  cluster_labels <- rep(1:4, each = 50)
  lineages <- .get_lineages(dat, cluster_labels, starting_cluster = 1)

  res <- .get_curves(dat, cluster_labels, lineages)

  bool_vec <- sapply(1:length(res), function(x){
    all(sapply(1:(length(lineages[[x]])-1), function(y){
      mean(res$pcurve_list[[x]]$lambda[((y-1)*50+1):(y*50)]) < mean(res$pcurve_list[[x]]$lambda[(y*50+1):((y+1)*50)])
    }))
  })

  expect_true(all(bool_vec))
})

##############################

## slingshot works

test_that("slingshot works", {
  set.seed(10)
  cell_pop <- matrix(c(4,10, 25,100,
                       60,80, 25,100,
                       40,10, 60,80,
                       60,80, 100,25)/10, nrow = 4, ncol = 4, byrow = T)
  h <- nrow(cell_pop)
  n_each <- 50
  dat <- do.call(rbind, lapply(1:h, function(x){
    pos <- stats::runif(n_each)
    cbind(pos*cell_pop[x,1] + (1-pos)*cell_pop[x,3] + stats::rnorm(n_each, sd = 0.1),
          pos*cell_pop[x,2] + (1-pos)*cell_pop[x,4] + stats::rnorm(n_each, sd = 0.1))
  }))
  cluster_labels <- rep(1:4, each = 50)
  res <- slingshot(dat, cluster_labels, starting_cluster = 1)

  expect_true(is.list(res))
  expect_true(length(res) == 3)
  expect_true(all(names(res) == c("lineages", "curves", "idx")))

  #plot(dat[,1], dat[,2], asp = T)
  #for(i in 1:length(res$curves)){lines(res$curves[[i]], col = i+1, lwd = 2)}
})


test_that("slingshot can give sensible lambdas", {
  set.seed(10)
  cell_pop <- matrix(c(4,10, 25,100,
                       60,80, 25,100,
                       40,10, 60,80,
                       60,80, 100,25)/10, nrow = 4, ncol = 4, byrow = T)
  h <- nrow(cell_pop)
  n_each <- 50
  dat <- do.call(rbind, lapply(1:h, function(x){
    pos <- stats::runif(n_each)
    cbind(pos*cell_pop[x,1] + (1-pos)*cell_pop[x,3] + stats::rnorm(n_each, sd = 0.1),
          pos*cell_pop[x,2] + (1-pos)*cell_pop[x,4] + stats::rnorm(n_each, sd = 0.1))
  }))
  cluster_labels <- rep(1:4, each = 50)
  res <- slingshot(dat, cluster_labels, starting_cluster = 1)

  bool_vec <- sapply(1:length(res$lineages), function(x){
    lin <- res$lineages[[x]]
    all(sapply(1:(length(lin)-1), function(y){
      idx1 <- which(cluster_labels == lin[y])
      idx2 <- which(cluster_labels == lin[y+1])

      mean(res$curves[[x]]$lambda[which(res$curves[[x]]$idx %in% idx1)]) <
        mean(res$curves[[x]]$lambda[which(res$curves[[x]]$idx %in% idx2)])
    }))
  })

  expect_true(all(bool_vec))
})

test_that("slingshot does not crash when there is only one lineage, 2D data", {
  set.seed(10)
  dat <- matrix(rep(1:100, times = 2), ncol = 2)
  cluster_labels <- rep(1:5, each = 20)
  res <- slingshot(dat, cluster_labels, starting_cluster = 1)

  expect_true(is.list(res))
  expect_true(class(res) == "slingshot")
  expect_true(all(res$lineages[[1]] == 1:5))
  expect_true(length(res$lineages) == 1)
  expect_true(length(res$curves) == 1)
  expect_true(all(res$curves[[1]]$W == 1))
})


test_that("slingshot works with an artifical example", {
  set.seed(10)
  cluster_labels <- sample(1:10, 200, replace = T)
  dat <- MASS::mvrnorm(200, rep(0, 5), diag(5))

  res <- slingshot(dat, cluster_labels, starting_cluster = 1)

  expect_true(is.list(res))
  expect_true(length(res) == 3)
  expect_true(all(names(res) == c("lineages", "curves", "idx")))
})

############

## .construct_cluster_matrix is correct

test_that(".construct_cluster_matrix works", {
  set.seed(10)
  res <- .construct_cluster_matrix(sample(1:10, 200, replace = T))

  expect_true(all(dim(res) == c(200,10)))
  expect_true(all(sort(unique(as.numeric(res))) == c(0,1)))
})

###########

## .compute_cluster_center is correct

test_that(".compute_cluster_center works", {
  set.seed(10)
  cluster_mat <- .construct_cluster_matrix(sample(1:10, 200, replace = T))
  dat <- MASS::mvrnorm(200, rep(0, 5), diag(5))

  res <- .compute_cluster_center(dat, cluster_mat)

  expect_true(all(dim(res) == c(10,5)))
})




