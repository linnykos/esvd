dropout_create <- function(a = 0.7, b = 1){
  function(x){1/(1+exp(-(a+b*x)))}
}

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
  .logistic_regression(dat[idx1, idx2], pred_dat, threshold_quant = threshold_quant_logistic)
}

####

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

.largest_average_degree_subgraph <- function(g, threshold_quant = 0.75){
  stopifnot(class(g) == "igraph")

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

  val <- quantile(avg_deg, probs = threshold_quant)
  idx <- max(which(avg_deg >= val))
  vertex_list[[idx]]
}

.funk_svd_prediction <- function(dat, k = 5, lambda = 0.01){
  dat[dat == 0] <- NA
  dat <- as(dat, "realRatingMatrix")

  funk_svd <- recommenderlab::funkSVD(dat, k = k, lambda = lambda)
  base::tcrossprod(funk_svd$U, funk_svd$V)
}

.logistic_regression <- function(observed_dat, pred_dat, threshold_quant = 0.5){
  stopifnot(all(dim(observed_dat) == dim(pred_dat)))

  zero_factor <- as.factor(as.numeric(observed_dat != 0))
  vec <- as.numeric(pred_dat)

  tmp_mat <- data.frame(vec, zero_factor)
  zz <- tmp_mat[(tmp_mat[,1] > 0), 1]

  cutoff_val <- quantile(zz, probs = threshold_quant)
  tmp_mat <- tmp_mat[which(tmp_mat[,1] > cutoff_val),]

  fit <- stats::glm(zero_factor ~ vec, family=binomial(link='logit'), data = tmp_mat)
  coef_vec <- stats::coef(fit)
  names(coef_vec) <- c("Intercept", "Slope")
  coef_vec
}
