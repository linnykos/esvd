rm(list=ls())
load("../../SOUP/data/camp.rda")

dat <- camp$counts

idx <- which(colnames(dat) %in% camp$select.genes)
dat <- dat[,idx]

dat <- t(apply(dat, 1, function(x){x/sum(x)}))
dat <- log10(dat + 1.01)

dim(dat)

####

# let's look at only the N's (since there are the most samples)
table(camp$cell.info$cell.type)
idx <- grep("N", camp$cell.info$cell.type)
dat_subset <- dat[idx,]
dim(dat_subset)

#search for differentiation using scImpute's code
source("../experiment/em_gamma_normal.R")
# diff_vec <- sapply(1:ncol(dat_subset), function(i){
#   x <- dat_subset[,i]
#   tryCatch({
#     res <- get_mix(x, prop_init = 0.3)
#     res[4] - res[2]/res[3]
#   }, error = function(e) {
#     NA
#   }
#   )
# })

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
hist_augment(dat_subset[,idx[1]], breaks = 50)

