rm(list=ls())
load("../../SOUP/data/camp.rda")

dat <- camp$counts

idx <- which(colnames(dat) %in% camp$select.genes)
dat <- dat[,idx]

dat <- t(apply(dat, 1, function(x){x/sum(x)}))
dat <- 10*log(dat + 1)

dim(dat)
dat <- dat[order(camp$cell.info[,2]),]
camp$cell.info <- camp$cell.info[order(camp$cell.info[,2]),]

length(which(dat == 0))/prod(dim(dat))

##################################

par(mfrow = c(1,3))
zz <- apply(dat, 1, function(x){length(which(x != 0))/length(x)})
plot(sort(zz))
zz <- apply(dat, 2, function(x){length(which(x != 0))/length(x)})
plot(sort(zz))
plot(sort(dat[dat != 0]))
