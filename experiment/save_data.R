dat <- read.csv("Zeisel_expr.txt", sep = "\t")
rownames(dat) <- dat[,1]
dat <- dat[,-1]

save(dat, file = "Zeisel_expr.RData")


