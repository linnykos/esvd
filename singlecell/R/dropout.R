.dropout <- function(dat){
  idx <- which(colSums(dat) != 0)

  dropout_res <- lapply(idx, function(i){
    .em_mixture(dat[,i])
  })

  dropout_idx <- lapply(1:length(idx), function(i){
    set.seed(10)
    x <- dat[,idx[i]]
    x2 <- .jitter_zeros(x)
    which(.compute_dropout(dropout_res[[i]], x2) == 1)
  })

  dropout_mat <- matrix(1, ncol = ncol(dat), nrow = nrow(dat))

  if(length(idx) != ncol(dat)){
    dropout_mat[,-idx] <- 0
  }

  for(i in 1:length(idx)){
    dropout_mat[dropout_idx[[i]], idx[i]] <- 0
  }

  dropout_mat
}
