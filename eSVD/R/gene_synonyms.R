#' Convert genes (in mouse) to a standardized set of synonyms
#'
#' This uses the database in \code{org.Mm.eg.db} package. At the time of running this
#' designing this code, we are using the version 3.10.0 released on	October 30, 2019, but
#' this code should work for even newer versions of \code{org.Mm.eg.db}.
#'
#' If genes in \code{vec} are found to be in the database in \code{org.Mm.eg.db}, they are converted
#' to the their standardize symbol. If genes are not found in the database, they are left
#' untouched. Gene names in \code{vec} are replaced in-place.
#'
#' @param vec a vector of gene names
#'
#' @return a vector of length \code{length(vec)}
#' @export
convert_synonyms <- function(vec){
  stopifnot(is.character(vec))
  len <- length(vec)

  dbCon <- org.Mm.eg.db::org.Mm.eg_dbconn()
  sqlQuery <- 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
  aliasSymbol <- DBI::dbGetQuery(dbCon, sqlQuery)

  syn_vec <- sapply(1:len, function(i){
    bool <- any(c(vec[i] %in% aliasSymbol$alias_symbol, vec[i] %in% aliasSymbol$symbol))
    if(!bool | vec[i] %in% aliasSymbol$symbol) return(vec[i])

    idx <- which(aliasSymbol$alias_symbol %in% vec[i])[1]
    aliasSymbol$symbol[idx]
  })

  names(syn_vec) <- NULL
  syn_vec
}
