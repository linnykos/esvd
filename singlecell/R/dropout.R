#' Create the dropout function
#'
#' @param a intercept of logistic function
#' @param b slope of logistic function
#'
#' @return function
#' @export
dropout_create <- function(a = 0.7, b = 1){
  function(x){1/(1+exp(-(a+b*x)))}
}

#' Estimate dropout function
#'
#' The dropout function is a logistic function, where this function estimates
#' the intercept and slope.
#'
#' @param dat numeric matrix with 0's
#' @param k rank of prediction for \code{.funk_svd_prediction}
#' @param lambda tuning parameter for \code{.funk_svd_prediction}
#' @param threshold_quant_degree parameter for \code{.largest_average_degree_subgraph}
#' @param threshold_quant_logistic parameter for \code{.logistic_regression}
#'
#' @return numeric vector
#' @export
estimate_dropout <- function(dat, k = 5, lambda = 0.01,
                             threshold_quant_degree = 0.75,
                             threshold_quant_logistic = 0.5){
  stopifnot(is.matrix(dat))
  n <- nrow(dat); d <- ncol(dat)

  adj_mat <- .form_adj_nonzero(dat)

  g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected")
  igraph::V(g)$name <- as.character(1:(n+d))

  idx <- as.numeric(.largest_average_degree_subgraph(g, threshold_quant = threshold_quant_degree))

  idx1 <- idx[which(idx <= n)]
  idx2 <- idx[which(idx > n)] - n

  pred_dat <- .funk_svd_prediction(dat[idx1, idx2], k = k, lambda = lambda)
  if(any(is.nan(pred_dat))) stop("Lambda is too large, causing Funk SVD to crash")
  .logistic_regression(dat[idx1, idx2], pred_dat, threshold_quant = threshold_quant_logistic)
}

####

#' Convert data to adjacency matrix
#'
#' If \code{dat} is a \code{n} by \code{d} matrix, this returns a
#' \code{n+d} by \code{n+d} block symmetric matrix. The first \code{n} by \code{n}
#' block (entries \code{[1:n, 1:n]}) and the last \code{d} by \code{d} block
#' (entries \code{[(n+1:d):(n+1:d)]}) are all 0's, but the off-diagonal blocks represent
#' the sparsity pattern of the \code{dat}.
#'
#' @param dat a matrix
#'
#' @return a binary matrix
.form_adj_nonzero <- function(dat){
  n <- nrow(dat); d <- ncol(dat)

  adj_mat <- matrix(0, n+d, n+d)
  for(i in 1:n){
    for(j in 1:d){
      adj_mat[i, n+j] <- ifelse(dat[i,j] == 0, 0, 1)
      adj_mat[n+j, i] <- adj_mat[i, n+j]
    }
  }

  adj_mat
}

#' Find largest average degree subgraph
#'
#' We use a greedy algorithm to find a reasonably dense subgraph of the input
#' igraph \code{g}. This algorithm produces a sequence of nested graphs by pruning
#' the vertex with least degree sequentially. Along this path, each subgraph has
#' an average degree (defined as the number of edges over the number of vertices).
#' \code{threshold_quant} defines the quantile of this average degree that we want, and
#' this function returns the smallest graph with an average degree at least this quantile.
#' Hence, \code{threshold_quant = 1} gives the approximate subgraph with the highest
#' average degree.
#'
#' @param g igraph object
#' @param threshold_quant numeric
#'
#' @return a vector of vertex names (characters) in \code{g}
#' @references Charikar, Moses. "Greedy approximation algorithms for finding
#' dense components in a graph." International Workshop on Approximation
#' Algorithms for Combinatorial Optimization. Springer, Berlin, Heidelberg, 2000.
.largest_average_degree_subgraph <- function(g, threshold_quant = 0.75){
  stopifnot(class(g) == "igraph")
  if(all(is.null(igraph::V(g)$name))) igraph::V(g)$name <- as.character(1:igraph::vcount(g))

  vertex_list <- vector("list", igraph::vcount(g))
  avg_deg <- rep(0, igraph::vcount(g)-1)

  g_tmp <- g

  for(i in 1:(length(vertex_list)-1)){
    name_vec <- igraph::V(g_tmp)$name
    vec <- igraph::degree(g_tmp)
    min_val <- min(vec)
    idx <- which(vec == min_val)
    if(length(idx) > 1) idx <- idx[sample(length(idx), 1)]

    name_vec <- name_vec[-idx]
    g_tmp <- igraph::induced_subgraph(g_tmp, name_vec)

    vertex_list[[i]] <- name_vec
    avg_deg[i] <- igraph::ecount(g_tmp)/igraph::vcount(g_tmp)
  }

  val <- stats::quantile(avg_deg, probs = threshold_quant)
  idx <- max(which(avg_deg >= val))
  vertex_list[[idx]]
}

#' Perform SVD with missing values
#'
#' This uses the Funk SVD in the \code{recommenderlab} package. In this function,
#' 0's in \code{dat} are interpreted as missing values. This function returns
#' the predicted matrix of \code{dat}.
#'
#' @param dat numeric matrix with 0's.
#' @param k rank of the SVD decomposition
#' @param lambda tuning parameter (learning rate) for Funk SVD
#'
#' @return a numeric matrix
#' @importClassesFrom recommenderlab realRatingMatrix
#' @export
.funk_svd_prediction <- function(dat, k = 5, lambda = 0.01){
  dat[dat == 0] <- NA
  dat <- methods::new("realRatingMatrix", data = recommenderlab::dropNA(dat))

  funk_svd <- recommenderlab::funkSVD(dat, k = k, lambda = lambda)
  if(any(is.nan(funk_svd$U)) || any(is.nan(funk_svd$V))) stop("Lambda is too large, causing Funk SVD to crash")
  base::tcrossprod(funk_svd$U, funk_svd$V)
}

#' Logistic regression
#'
#' This function performs a logistic regression of \code{observed_dat} (where the
#' non-zeros are labeled as \code{TRUE} and the zeros are labeled as \code{FALSE})
#' onto \code{pred_dat}, both of which are matrices but are converted into vectors
#' within this function.
#'
#' Pruning is performed based on \code{threshold_quant}. Values in \code{pred_dat}
#' (and their corresponding value in \code{observed_dat}) that are lower than the
#' \code{threshold_quant} quantile of non-zero's are removed prior to the logistic
#' regression.
#'
#' @param observed_dat a numeric matrix with 0's
#' @param pred_dat a numeric matrix
#' @param threshold_quant a numeric
#'
#' @return a vector of 2 numerics
.logistic_regression <- function(observed_dat, pred_dat, threshold_quant = 0.5){
  stopifnot(all(dim(observed_dat) == dim(pred_dat)))

  zero_factor <- as.factor(as.numeric(observed_dat != 0))
  vec <- as.numeric(pred_dat)

  tmp_mat <- data.frame(vec, zero_factor)
  zz <- tmp_mat[(tmp_mat[,1] > 0), 1]

  cutoff_val <- stats::quantile(zz, probs = threshold_quant)
  tmp_mat <- tmp_mat[which(tmp_mat[,1] > cutoff_val),]

  fit <- stats::glm(zero_factor ~ vec, family=stats::binomial(link='logit'), data = tmp_mat)
  coef_vec <- stats::coef(fit)
  names(coef_vec) <- c("Intercept", "Slope")
  coef_vec
}
