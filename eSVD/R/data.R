#' Zhang dataset of marker genes
#'
#' @name zhang_genes
#' @docType data
#' @format A data frame with two columns
#' \describe{
#'   \item{cell_type}{One of eight cell types: astrocytes, endothelial, microglia, MO (mature
#'   oligodendrocytes), neurons, NFO (newly-formed oligodendrocytes), OPCs (oligodendrocyte precursors),
#'   and pericytes}
#'   \item{enriched_genes}{Genes in mice that are enriched for the particular cell type. These gene
#'   names are as-given in the Zhang paper. One might want to use functions like \code{eSVD::convert_synonyms}
#'   to ensure these gene names are converted properly to their currently most-common synonym}
#' }
#' @author Kevin Lin \email{kevinl1@andrew.cmu.edu}
#' @source The dataset is Figure 4 of
#' "Zhang, Ye, et al. 'An RNA-sequencing transcriptome and splicing database of glia,
#' neurons, and vascular cells of the cerebral cortex.' Journal of Neuroscience 34.36 (2014):
#' 11929-11947," located at \url{https://www.jneurosci.org/content/jneuro/34/36/11929.full.pdf}.
#' @keywords data
NULL
