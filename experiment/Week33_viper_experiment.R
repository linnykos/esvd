rm(list=ls())
library(VIPER)
set.seed(10)

data(grun)
gene.expression <- t(as.matrix(gene.expression))
dat <- singlecell::downsample(gene.expression, downsample_rate = 0.5, dropoff_rate = 0.2)
tmp <- as.data.frame(t(dat$dat))
VIPER::VIPER(tmp, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 1,
             report = FALSE, outdir = NULL, prefix = NULL)
