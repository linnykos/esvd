rm(list=ls())
load("../../SOUP/data/zeisel.rda")

dat = zeisel$counts
dim(dat)

sparsity_vec <- apply(dat, 2, function(x){
  length(which(x!=0))
})
sum_vec <- colSums(dat)

plot(sparsity_vec, sum_vec)

quantile_sum <- quantile(sum_vec, prob = 0.95)
plot(sparsity_vec, sum_vec, ylim = c(0, quantile_sum))

idx <- which(sum_vec <= quantile_sum)
dat <- dat[,idx]

anova_func <- function(x){
  uniq <- unique(zeisel$cell.info[,2])
  k <- length(uniq)
  n <- length(x)

  between <- sum(sapply(uniq, function(k){
    idx <- which(zeisel$cell.info[,2] == k)
    length(idx)*(mean(x[idx]) - mean(x))^2
  }))/(k-1)

  within <- sum(sapply(unique(zeisel$cell.info[,2]), function(k){
    idx <- which(zeisel$cell.info[,2] == k)
    sum((x[idx]-mean(x[idx]))^2)
  }))/(n-k)

  between/within
}

anova_vec <- apply(dat, 2, anova_func)
plot(sort(anova_vec))

cell_idx <- as.numeric(as.factor(zeisel$cell.info[,2]))
idx <- which.max(anova_vec)
plot(jitter(sort(dat[,idx])), col = c(1:7)[cell_idx[order(dat[,idx])]])

idx <- which.min(anova_vec)
plot(jitter(sort(dat[,idx])), col = c(1:7)[cell_idx[order(dat[,idx])]])
