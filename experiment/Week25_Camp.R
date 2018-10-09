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

png("../figure/experiment/25_camp_data_properties.png", height = 1000, width = 2400, res = 300, units = "px")
par(mfrow = c(1,3))
zz <- apply(dat, 1, function(x){length(which(x != 0))/length(x)})
plot(sort(zz), ylab = "Percentage non-zero", main = "Sparsity per cell (row)\n(Camp data)",
     xlab = "Row rank", pch = 16, col = rgb(0,0,0,0.25))
zz <- apply(dat, 2, function(x){length(which(x != 0))/length(x)})
plot(sort(zz), ylab = "Percentage non-zero", main = "Sparsity per gene (colum)\n(Camp data)",
     xlab = "Column rank", pch = 16, col = rgb(0,0,0,0.25))
plot(sort(dat[dat != 0]), xlab = "Sorted rank (non-zeros)", ylab = "Value",
     main = "Non-zero quantiles\n(Camp data)", pch = 16, col = rgb(0,0,0,0.25))
graphics.off()
