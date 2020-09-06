context("Test informative genes")

## generator functions needed for the tests to come

generate_natural_mat <- function(cell_pop, gene_pop, n_each, d_each, sigma, modifier, tol = 1){
  h <- nrow(cell_pop)
  cell_mat <- do.call(rbind, lapply(1:h, function(x){
    pos <- stats::runif(n_each)
    cbind(pos*cell_pop[x,1] + (1-pos)*cell_pop[x,3] + stats::rnorm(n_each, sd = sigma),
          pos*cell_pop[x,2] + (1-pos)*cell_pop[x,4] + stats::rnorm(n_each, sd = sigma))
  }))
  cell_mat <- pmax(cell_mat, tol)
  n <- nrow(cell_mat)
  k <- ncol(cell_mat)

  # construct the gene information
  g <- nrow(gene_pop)
  gene_mat <- do.call(rbind, lapply(1:g, function(x){
    pos <- stats::runif(d_each)
    cbind(pos*gene_pop[x,1] + (1-pos)*gene_pop[x,3] + stats::rnorm(d_each, sd = sigma),
          pos*gene_pop[x,2] + (1-pos)*gene_pop[x,4] + stats::rnorm(d_each, sd = sigma))
  }))
  gene_mat <- pmax(gene_mat, tol)
  d <- nrow(gene_mat)

  res <- eSVD:::.reparameterize(cell_mat*sqrt(modifier), gene_mat*sqrt(modifier))

  list(nat_mat = res$u_mat %*% t(res$v_mat), cell_mat =  res$u_mat, gene_mat = res$v_mat)
}

generator_esvd_nb <- function(nat_mat, scalar = 100,  ...){
  n <- nrow(nat_mat); d <- ncol(nat_mat)

  obs_mat <- matrix(0, ncol = d, nrow = n)
  for(i in 1:n){
    for(j in 1:d){
      p <- 1-exp(-nat_mat[i,j]) # remember the p param in rnbinom is "flipped" in R
      obs_mat[i,j] <- stats::rnbinom(1, size = scalar, prob = p)
    }
  }

  obs_mat
}

## .np_smoother is correct

test_that(".np_smoother works", {
  set.seed(10)
  vec <- rnorm(100)
  res <- .np_smoother(vec)

  expect_true(is.numeric(res))
  expect_true(length(res) == length(vec))
})

#############################

## .circular_segmentation is correct

test_that(".circular_segmentation works", {
  set.seed(10)
  vec <- c(stats::rnorm(100), stats::rnorm(100, mean = 5), stats::rnorm(100))
  res <- .circular_segmentation(vec)

  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("i", "j", "obj_val"))))
})

test_that(".circular_segmentation outputs reasonable answer", {
  set.seed(10)
  vec1 <- c(stats::rnorm(100), stats::rnorm(100, mean = 5), stats::rnorm(100))
  res1 <- .circular_segmentation(vec1, max_width_percentage = 1)

  vec2 <- c(stats::rnorm(100), stats::rnorm(100, mean = 5), stats::rnorm(100))
  res2 <- .circular_segmentation(vec2, max_width_percentage = 1)

  vec3 <- c(stats::rnorm(300))
  res3 <- .circular_segmentation(vec3, max_width_percentage = 1)

  expect_true(abs(res1$i - res2$i) < min(abs(res1$i - res3$i), abs(res2$i - res3$i)))
  expect_true(abs(res1$j - res2$j) < min(abs(res1$j - res3$j), abs(res2$j - res3$j)))
})

test_that(".circular_segmentation is faster with a smaller resolution", {
  set.seed(10)

  vec <- c(stats::rnorm(100), stats::rnorm(100, mean = 5), stats::rnorm(100))
  res <- microbenchmark::microbenchmark(.circular_segmentation(vec, max_width_percentage = 0.2),
                                        .circular_segmentation(vec, resolution = 1/10, max_width_percentage = 0.2),
                                        times = 20)
  sum_res <- summary(res)

  expect_true(sum_res$mean[2] < sum_res$mean[1])
})

test_that(".circular_segmentation is faster with a smaller max_width_percentage", {
  set.seed(10)

  vec <- c(stats::rnorm(100), stats::rnorm(100, mean = 5), stats::rnorm(100))
  res <- microbenchmark::microbenchmark(.circular_segmentation(vec, max_width_percentage = 0.5),
                                        .circular_segmentation(vec, max_width_percentage = 0.1),
                                        times = 20)
  sum_res <- summary(res)

  expect_true(sum_res$mean[2] < sum_res$mean[1])
})

test_that(".circular_segmentation respects hard_cut", {
  trials <- 10

  bool_vec <- sapply(1:trials, function(i){
    set.seed(i)
    vec <- c(stats::rnorm(100), stats::rnorm(100, mean = 5), stats::rnorm(100))

    res <- .circular_segmentation(vec, hard_cut = 250)
    res$i <= 250 & res$j <= 250
  })

  expect_true(all(bool_vec))
})

########################

## .find_highly_expressed_region is correct

test_that(".find_highly_expressed_region works", {
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
  slingshot_res <- slingshot(dat, cluster_labels, starting_cluster = 1, upscale_factor = 1)
  cell_partition <- prepare_data_for_segmentation(dat, cluster_labels = cluster_labels,
                                                  curve_list = slingshot_res)

  gene_idx <- 1
  res <- .find_highly_expressed_region(common_vec = dat[cell_partition$cell_idx_common, gene_idx],
                                       specific_vec1 = dat[cell_partition$cell_idx_traj1, gene_idx],
                                       specific_vec2 = dat[cell_partition$cell_idx_traj2, gene_idx])

  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("cut_1", "cut_2", "vec1_smooth", "vec2_smooth"))))
  expect_true(length(res$vec1_smooth) == length(cell_partition$cell_idx_common) + length(cell_partition$cell_idx_traj1))
  expect_true(length(res$vec2_smooth) == length(cell_partition$cell_idx_common) + length(cell_partition$cell_idx_traj2))
})

####################################3

## segment_genes_along_trajectories is correct

test_that("segment_genes_along_trajectories works", {
  set.seed(10)

  n_each <- 50
  d_each <- 10
  sigma <- 1
  modifier <- 1/250
  size <- 50
  cell_pop <- matrix(c(4,10, 25,100, 60,80, 25,100,
                       40,10, 60,80, 60,80, 100, 25),
                     nrow = 4, ncol = 4, byrow = T)
  gene_pop <- matrix(c(20,90, 25,100,
                       90,20, 100,25)/20, nrow = 2, ncol = 4, byrow = T)

  res <- generate_natural_mat(cell_pop, gene_pop, n_each, d_each, sigma, modifier)
  nat_mat <- res$nat_mat

  dat <- nat_mat
  dat <- round(dat * 1000/max(dat))

  cluster_labels <- rep(1:4, each = n_each)

  slingshot_res <- slingshot(dat, cluster_labels, starting_cluster = 1, upscale_factor = 1)
  cell_partition <- prepare_data_for_segmentation(dat, cluster_labels = cluster_labels,
                                                  curve_list = slingshot_res)

  res <- segment_genes_along_trajectories(cell_partition$dat1, cell_partition$dat2, common_n = length(cell_partition$cell_idx_common),
                                          resolution = 1/10, max_width_percentage = 1/5)

  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("df", "segmentation_fit"))))

  expect_true(is.data.frame(res$df))
  expect_true(nrow(res$df) == ncol(dat))
  expect_true(all(sort(colnames(res$df)) == sort(c("idx", "start_1", "end_1", "mid_1", "obj_1", "start_2",
                                                "end_2", "mid_2", "obj_2"))))

  expect_true(is.list(res$segmentation_fit))
  expect_true(length(res$segmentation_fit) == ncol(dat))
  expect_true(length(unique(sapply(res$segmentation_fit, length))) == 1)
  expect_true(all(sort(names(res$segmentation_fit[[1]])) == sort(c("cut_1", "cut_2", "vec1_smooth", "vec2_smooth"))))
  expect_true(all(sort(names(res$segmentation_fit[[1]]$cut_1)) == sort(c("i", "j", "obj_val"))))
  expect_true(all(sort(names(res$segmentation_fit[[1]]$cut_2)) == sort(c("i", "j", "obj_val"))))
  expect_true(res$segmentation_fit[[1]]$cut_1$i %% 1 == 0)
  expect_true(res$segmentation_fit[[1]]$cut_2$i %% 1 == 0)
  expect_true(min(sapply(res$segmentation_fit, function(x){x$cut_1$i})) > 0)
  expect_true(min(sapply(res$segmentation_fit, function(x){x$cut_2$i})) > 0)
  expect_true(max(sapply(res$segmentation_fit, function(x){x$cut_1$i})) <= nrow(cell_partition$dat1))
  expect_true(max(sapply(res$segmentation_fit, function(x){x$cut_2$i})) <= nrow(cell_partition$dat2))
  expect_true(min(sapply(res$segmentation_fit, function(x){x$cut_1$j})) > 0)
  expect_true(min(sapply(res$segmentation_fit, function(x){x$cut_2$j})) > 0)
  expect_true(max(sapply(res$segmentation_fit, function(x){x$cut_1$j})) <= nrow(cell_partition$dat1))
  expect_true(max(sapply(res$segmentation_fit, function(x){x$cut_2$j})) <= nrow(cell_partition$dat2))
  expect_true(all(sapply(res$segmentation_fit, function(x){x$cut_1$j - x$cut_1$i > 0})))
  expect_true(all(sapply(res$segmentation_fit, function(x){x$cut_2$j - x$cut_2$i > 0})))
  expect_true(all(sapply(res$segmentation_fit, function(x){length(x$vec1_smooth) == nrow(cell_partition$dat1)})))
  expect_true(all(sapply(res$segmentation_fit, function(x){length(x$vec2_smooth) == nrow(cell_partition$dat2)})))
})

#################################

## .extract_information is correct

test_that(".extract_information works", {
  set.seed(10)

  n_each <- 50
  d_each <- 10
  sigma <- 1
  modifier <- 1/250
  size <- 50
  cell_pop <- matrix(c(4,10, 25,100, 60,80, 25,100,
                       40,10, 60,80, 60,80, 100, 25),
                     nrow = 4, ncol = 4, byrow = T)
  gene_pop <- matrix(c(20,90, 25,100,
                       90,20, 100,25)/20, nrow = 2, ncol = 4, byrow = T)

  res <- generate_natural_mat(cell_pop, gene_pop, n_each, d_each, sigma, modifier)
  nat_mat <- res$nat_mat

  dat <- nat_mat
  dat <- round(dat * 1000/max(dat))

  cluster_labels <- rep(1:4, each = n_each)

  slingshot_res <- slingshot(res$cell_mat, cluster_labels, starting_cluster = 1, upscale_factor = 1)
  pseudotime_df <- construct_pseudotime_trajectory_matrix(slingshot_res, cluster_labels)
  cell_partition <- prepare_data_for_segmentation(dat, cluster_labels = cluster_labels,
                                                  curve_list = slingshot_res)

  dat1 <- cell_partition$dat1; dat2 <- cell_partition$dat2

  common_n = length(cell_partition$cell_idx_common)
  resolution = 1/10
  max_width_percentage = 1/5

  n1 <- nrow(dat1) - common_n
  n2 <- nrow(dat2) - common_n
  p <- ncol(dat1)

  func <- function(j){
    .find_highly_expressed_region(common_vec = dat1[1:common_n,j],
                                  specific_vec1 = dat1[(common_n+1):(common_n+n1),j],
                                  specific_vec2 = dat2[(common_n+1):(common_n+n2),j],
                                  standardize = T, resolution = resolution, max_width_percentage = max_width_percentage)
  }

  segmentation_res <- lapply(1:p, func)

  res <- .extract_information(segmentation_res)

  expect_true(is.data.frame(res))
  expect_true(nrow(res) == ncol(dat))
  expect_true(all(sort(colnames(res)) == sort(c("idx", "start_1", "end_1", "mid_1", "obj_1", "start_2",
                                                "end_2", "mid_2", "obj_2"))))
})

#################################

## order_highly_expressed_genes is correct

test_that("order_highly_expressed_genes works", {
  set.seed(10)

  n_each <- 50
  d_each <- 10
  sigma <- 1
  modifier <- 1/250
  size <- 50
  cell_pop <- matrix(c(4,10, 25,100, 60,80, 25,100,
                       40,10, 60,80, 60,80, 100, 25),
                     nrow = 4, ncol = 4, byrow = T)
  gene_pop <- matrix(c(20,90, 25,100,
                       90,20, 100,25)/20, nrow = 2, ncol = 4, byrow = T)

  res <- generate_natural_mat(cell_pop, gene_pop, n_each, d_each, sigma, modifier)
  nat_mat <- res$nat_mat

  dat <- nat_mat
  dat <- round(dat * 1000/max(dat))

  cluster_labels <- rep(1:4, each = n_each)

  slingshot_res <- slingshot(res$cell_mat, cluster_labels, starting_cluster = 1, upscale_factor = 1)
  pseudotime_df <- construct_pseudotime_trajectory_matrix(slingshot_res, cluster_labels)
  cell_partition <- prepare_data_for_segmentation(dat, cluster_labels = cluster_labels,
                                                  curve_list = slingshot_res)

  dat1 <- cell_partition$dat1; dat2 <- cell_partition$dat2

  res_mat <- segment_genes_along_trajectories(dat1, dat2, common_n = length(cell_partition$cell_idx_common),
                                          resolution = 1/10, max_width_percentage = 1/5)

  res <- order_highly_expressed_genes(res_mat$df, nrow1 = nrow(dat1), nrow2 = nrow(dat2),
                                      common_n = length(cell_partition$cell_idx_common), threshold = 0.1)

  expect_true(is.list(res))
  expect_true(length(res) == 3)
  expect_true(all(sort(names(res)) == sort(c("common_genes", "traj1_genes", "traj2_genes"))))
})

#####################################

## prepare_data_for_segmentation is correct

test_that("prepare_data_for_segmentation works", {
  set.seed(10)

  n_each <- 50
  d_each <- 10
  sigma <- 1
  modifier <- 1/250
  size <- 50
  cell_pop <- matrix(c(4,10, 25,100, 60,80, 25,100,
                       40,10, 60,80, 60,80, 100, 25),
                     nrow = 4, ncol = 4, byrow = T)
  gene_pop <- matrix(c(20,90, 25,100,
                       90,20, 100,25)/20, nrow = 2, ncol = 4, byrow = T)

  res <- generate_natural_mat(cell_pop, gene_pop, n_each, d_each, sigma, modifier)
  nat_mat <- res$nat_mat

  dat <- nat_mat
  dat <- round(dat * 1000/max(dat))
  cluster_labels <- rep(1:4, each = n_each)
  slingshot_res <- slingshot(res$cell_mat, cluster_labels, starting_cluster = 1, upscale_factor = 1)

  res <- prepare_data_for_segmentation(dat, cluster_labels, curve_list = slingshot_res)

  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("dat1", "dat2", "cell_idx_common", "cell_idx_traj1", "cell_idx_traj2"))))
  expect_true(length(res$cell_idx_common) <= nrow(dat))
  expect_true(length(res$cell_idx_traj1) <= nrow(dat))
  expect_true(length(res$cell_idx_traj2) <= nrow(dat))
  expect_true(ncol(res$dat1) == ncol(dat))
  expect_true(ncol(res$dat2) == ncol(dat))
  expect_true(nrow(res$dat1) == length(res$cell_idx_common) + length(res$cell_idx_traj1))
  expect_true(nrow(res$dat2) == length(res$cell_idx_common) + length(res$cell_idx_traj2))
  expect_true(length(intersect(res$cell_idx_common, res$cell_idx_traj1)) == 0)
  expect_true(length(intersect(res$cell_idx_common, res$cell_idx_traj2)) == 0)
})
