rm(list=ls())
load("../../SOUP/data/camp.rda")

dat <- camp$counts

idx <- which(colnames(dat) %in% camp$select.genes)
dat <- dat[,idx]

dat <- t(apply(dat, 1, function(x){x/sum(x)}))
dat <- 10*log(dat + 1)

dim(dat)

##################

dropout_res <- vector("list", length = ncol(dat))
for(i in 1:ncol(dat)){
  set.seed(10)
  if(i %% floor(ncol(dat)/10) == 0) cat('*')
  x <- dat[,i]
  if(all(x == 0)) {
    dropout_res[[i]] <- NA
    next()
  }
  if(length(x[x > 0]) <= 10) {
    dropout_res[[i]] <- NA
    next()
  }

  dropout_res[[i]] <- tryCatch({
    .em_mixture(x)
  }, error = function(e){
    print(paste("Error in ", i))
    return(NA)
  })
}
