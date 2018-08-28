rm(list=ls())
load("../../SOUP/data/zeisel.rda")

dat <- zeisel$counts

dat <- t(apply(dat, 1, function(x){x/sum(x)}))
dat <- log10(dat + 1.01)

idx <- which(colnames(dat) %in% zeisel$select.genes)
dat <- dat[,idx]
dim(dat)

table(zeisel$cell.info$cell.type)
idx <- grep("oligodendrocytes", zeisel$cell.info$cell.type)
dat_subset <- dat[idx,]
dim(dat_subset)

diff_vec <- sapply(1:ncol(dat_subset), function(i){
  if(i %% floor(ncol(dat_subset)/10) == 0) cat('*')
  x <- dat_subset[,i]
  zz <- hist(x, breaks = 100, plot = F)
  # sum(zz$counts[-1])/sum(zz$counts)
  len <- length(which(x == log10(1.01)))

  if(len > 0){
    zz$counts[2] + zz$counts[1] - len
  } else {
    zz$counts[2]
  }
})

hist_augment <- function(x, breaks = 50, min_val = log10(1.01), ...){
  break_vec <- seq(min_val, max(x), length.out = breaks)

  min_nonzero <- min(x[x != min_val])
  interval <- diff(break_vec)[1]

  if(interval + min_val > min_nonzero){
    interval <- min_nonzero - min_val
  }

  break_vec <- break_vec + interval/2
  break_vec <- c(break_vec[1] - diff(break_vec)[1], break_vec)

  zz <- hist(x, breaks = break_vec, ...)
  len <- length(which(x == min_val))
  if(len != 0){
    rect(zz$breaks[1], 0, zz$breaks[2], len, col = rgb(0.803, 0.156, 0.211))
  }

  invisible()
}

quantile(diff_vec, na.rm = T)
range(diff_vec, na.rm = T)

idx <- order(diff_vec, decreasing = T)

hist_augment(dat_subset[,idx[1000]], breaks = 100, main = "")
hist_augment(dat_subset[,idx[1]], breaks = 100, main = "")

vec1 <- dat_subset[,idx[1000]]
vec2 <- dat_subset[,idx[1]]

############

source("../experiment/em_gamma_normal_switch.R")

compute_dropout <- function(vec, param){
  estimated <- rep(2, length(vec))

  if(class(param) == "Gamma-Normal" | class(param) == "Gamma-TNormal"){
    gamma_dist <- param[1] * sapply(vec, stats::dgamma, shape = param[2], rate = param[3])
    norm_dist <- (1-param[1]) * sapply(vec, stats::dnorm, mean = param[4], sd = param[5])
    if(class(param) == "Gamma-TNormal") norm_dist <- norm_dist / (1 - stats::pnorm(0, mean = param[4], sd = param[5]))

    ratio <- gamma_dist/(gamma_dist + norm_dist)
    estimated[ratio > 0.5] <- 1
  } else {
    estimated[which(vec == log10(1.01))] <- 1
  }

  estimated
}

draw_curve <- function(vec_list, ymax = NA, max_val = 15,
                       min_val = log10(1.01)/10, ...){
  #compute distribution
  x_seq <- seq(min_val, max_val, length.out = 1000)

  density_list <- lapply(vec_list, function(vec){
    if(class(vec) == "Gamma-Normal" | class(vec) == "Gamma-TNormal"){
      gamma_dist <- vec[1] * sapply(x_seq, stats::dgamma, shape = vec[2], rate = vec[3])
      norm_dist <- (1-vec[1]) * sapply(x_seq, stats::dnorm, mean = vec[4], sd = vec[5])
      if(class(vec) == "Gamma-TNormal") norm_dist <- norm_dist / (1 - stats::pnorm(0, mean = vec[4], sd = vec[5]))

      list(gamma_dist, norm_dist)
    } else {
      norm_dist <- (1-vec[1]) * sapply(x_seq, stats::dnorm, mean = vec[3], sd = vec[4])
      list(NA, norm_dist)
    }
  })

  if(is.na(ymax)) ymax <- max(unlist(density_list), na.rm = T)

  plot(NA, xlim = c(min_val, max_val), ylim = c(0, ymax), ...)
  for(i in 1:length(vec_list)){
    if(!any(is.na(density_list[[i]][[1]]))) {
      lines(x_seq, density_list[[i]][[1]], col = rgb(0.584, 0.858, 0.564), lty = i, ... )
    } else {
      lines(rep(log10(1.01), 2), c(-1e6, 1e6), col = rgb(0.584, 0.858, 0.564), lty = i, ... )
    }
    lines(x_seq, density_list[[i]][[2]], col = rgb(0.803, 0.156, 0.211), lty = i, ...)
  }

  invisible()
}

res_original <- get_mix_switch(vec1)
res_trunc <- get_mix_switch(vec1, truncated = T)
# res_trunc_mom <- .get_mix(vec1, prop_init = .1, MoM = T)

class1 <- compute_dropout(vec1, res_original, F)
class2 <- compute_dropout(vec1, res_trunc, T)

table(class1, class2)

# let's also draw the density
draw_curve(list(res_original), max_val = max(vec1),
           min_val = min(vec1), ymax = 100, lwd = 2)
lines(x = rep(log10(1.01), 2), y = c(-1e6, 1e6), lwd = 2, lty = 2)
draw_curve(list(res_trunc), max_val = 2*max(vec1), ymax = 100, lwd = 2)
lines(x = rep(log10(1.01), 2), y = c(-1e6, 1e6), lwd = 2, lty = 2)

##########

res_original <- get_mix_switch(vec2)

draw_curve(list(res_original), max_val = max(vec2),
           min_val = min(vec2), ymax = 20, lwd = 2)

class1 <- compute_dropout(vec2, res_original)

table(class1, class2)

###########################

#enumerate over all genes and see if there is EVER a instance where
# log10(1.01) is classified as "true"

cutoff_vec <- sapply(1:ncol(dat_subset), function(x){
  if(x %% floor(ncol(dat_subset)/10) == 0) cat('*')

  vec <- dat_subset[,x]
  if(length(vec[vec > log10(1.01)]) < 5) return(NA)

  res <- get_mix_switch(vec)
  class_vec <- compute_dropout(vec, res)
  if(length(which(class_vec == 2)) > 0){
    min(vec[class_vec == 2])
  } else {
    Inf
  }
})
length(which(cutoff_vec == log10(1.01)))

class_vec <- sapply(1:ncol(dat_subset), function(x){
  if(x %% floor(ncol(dat_subset)/10) == 0) cat('*')

  vec <- dat_subset[,x]
  if(length(vec[vec > log10(1.01)]) < 5) return(NA)

  res <- get_mix_switch(vec)
  class(res)
})
table(class_vec)

count_vec <- sapply(1:ncol(dat_subset), function(x){
  vec <- dat_subset[,x]
  length(vec[vec != log10(1.01)])
})

cutoff_vec2 <- cutoff_vec[intersect(which(!is.na(cutoff_vec)),
                                    which(!is.infinite(cutoff_vec)))]




