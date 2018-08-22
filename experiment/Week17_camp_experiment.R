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
idx <- grep("A", camp$cell.info$cell.type)
dat_subset <- dat[idx,]
dim(dat_subset)

#search for differentiation using scImpute's code
source("../experiment/em_gamma_normal.R")
diff_vec <- sapply(1:ncol(dat_subset), function(i){
  x <- dat_subset[,i]
  tryCatch({
    res <- get_mix(x, prop_init = 0.3)
    res[4] - res[2]/res[3]
  }, error = function(e) {
    NA
  }
  )
})

quantile(diff_vec, na.rm = T)
range(diff_vec, na.rm = T)
range(dat_subset)
