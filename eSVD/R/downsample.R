#' Downsample a dataset
#'
#' @param dat dataset where the \code{n} rows represent cells and \code{d} columns represent genes
#' @param downsample_rate numeric of the rate of downsample (higher means less cells are kept)
#' @param dropoff_rate numeric of the rate of dropoff (higher means less cells are kept)
#'
#' @return list containing \code{dat}, \code{downsample_idx}, and \code{dropout_idx}
#' @export
downsample <- function(dat, downsample_rate = 0.5, dropoff_rate = 0.2){
  d <- ncol(dat)

  downsample_dat <- t(apply(dat, 1, function(x){
    count <- sum(x)
    tmp <- sample(1:d, round((1-downsample_rate)*count), prob = x, replace = T)
    tmp <- table(tmp)

    vec <- rep(0, d)
    vec[as.numeric(names(tmp))] <- tmp
    vec
  }))
  downsample_idx <- intersect(which(downsample_dat == 0), which(dat != 0))

  idx <- which(downsample_dat != 0)
  dropout_idx <- sample(idx, round(dropoff_rate*length(idx)))
  dropout_dat <- downsample_dat
  dropout_dat[dropout_idx] <- 0

  list(dat = dropout_dat, downsample_idx = downsample_idx, dropout_idx = dropout_idx)
}
