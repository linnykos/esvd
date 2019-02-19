rm(list=ls())

dat1 <- read.table("../../raw_data/Buettner/G1_singlecells_counts.txt",
                   stringsAsFactors = F)
dat2 <- read.table("../../raw_data/Buettner/G2M_singlecells_counts.txt",
                   stringsAsFactors = F)
dat3 <- read.table("../../raw_data/Buettner/S_singlecells_counts.txt",
                   stringsAsFactors = F)
