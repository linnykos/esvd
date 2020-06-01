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
