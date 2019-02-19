#' Finding the true zeros
#'
#' @param dropout_mat a \code{n} by \code{d} matrix
#' @param num_neighbors number of neighbors
#' @param verbose boolean
#'
#' @return a 0-1-NA matrix of size \code{n} by \code{d}
#' @export
find_true_zeros <- function(dropout_mat, num_neighbors = NA, cores = NA,
                            verbose = F){
  if(is.na(num_neighbors)) num_neighbors <- ceiling(nrow(dropout_mat)/10)

  if(!is.na(cores)){
    doMC::registerDoMC(cores = cores)
    tmp_func <- function(x){
      if(verbose & x %% floor(nrow(dropout_mat)/10) == 0) cat('*')
      singlecell:::.find_neighbor(dropout_mat, x, num_neighbors)
    }

    neighbor_list <- foreach::"%dopar%"(foreach::foreach(i = 1:nrow(dropout_mat)),
                                        tmp_func(i))
  } else {
    neighbor_list <- lapply(1:nrow(dropout_mat), function(x){
      if(verbose & x %% floor(nrow(dropout_mat)/10) == 0) cat('*')
      .find_neighbor(dropout_mat, x, num_neighbors)
    })
  }

  zero_mat <- dropout_mat
  for(i in 1:nrow(dropout_mat)){
    idx <- which(dropout_mat[i,] == 0)
    vec <- rep(0, length(idx))
    zz <- colSums(dropout_mat[neighbor_list[[i]], idx])
    vec[which(zz != 0)] <- NA
    zero_mat[i, idx] <- vec
  }

  zero_mat
}


.find_neighbor <- function(mat, i, num_neighbors){
  row_vec <- c(1:nrow(mat))[-i]
  idx_me <- which(mat[i,] == 0)
  vec <- sapply(row_vec, function(x){
    idx_other <- which(mat[x,] == 0)
    length(intersect(idx_me, idx_other))/length(unique(c(idx_me, idx_other)))
  })

  row_vec[order(vec, decreasing = T)[1:num_neighbors]]
}

