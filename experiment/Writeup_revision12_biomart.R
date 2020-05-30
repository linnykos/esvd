# https://www.bioconductor.org/help/course-materials/2015/BioC2015/Annotation_Resources.html

library("biomaRt")
head(listMarts())
ensembl <- useMart("ensembl")
ensembl
head(listDatasets(ensembl))

ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
ensembl
head(listAttributes(ensembl))

res <- getBM(attributes=c("hgnc_symbol"), mart = ensembl)
head(res)
