# rm(list=ls())
# source("parse_taxonomy_function.r")

parse_taxonomy <- function(data, column="classification", delim=";", 
                           levels=c("Domain","Phylum","Class","Order","Family","Genus","Species"),
                           write_to_file=F, file=NULL){
  #delim=";"
  #column="classification"
  #data=df
  #levels=c("Domain","Phylum","Class","Order","Family","Genus","Species")
  
  tax <- as.character(data[,column])
  parsed_tax <- strsplit(tax, split=delim)
  max_lev <- length(parsed_tax[[1]])
  
  tax_df <- do.call(rbind, lapply(parsed_tax, function(x) {
    length(x) <- max_lev  # pad with NA if needed
    names(x) <- levels
    return(x)
  }))
  
  # remove prefixes
  tax_df <- apply(tax_df, 2, function(x) {gsub(pattern = "[a-z]__",replacement = "", x)})
  
  # remove spaces
  tax_df <- apply(tax_df, 2, function(x) {gsub(pattern = "\\s+",replacement = "_", x)})
  
  # fill empties
  tax_df[tax_df == ""] <- NA  # First, treat empty strings as NA
  
  # Apply the logic row-wise
  tax_df <- as.data.frame(t(apply(tax_df, 1, function(r) {
    for (i in 2:length(r)) {
      if (is.na(r[i]) || r[i] == "") {
        r[i] <- paste(r[i - 1], tolower(levels[i]),sep = "_")
      }
    }
    return(r)
  })), stringsAsFactors = FALSE)
  
  if(write_to_file==T){
    if(is.null(file)){file=paste0(getwd(),"/parsed_taxonomy.csv",collapse = "")}
    write.csv(tax_df, file = file,quote = F,row.names = F)
  }
  
  return(tax_df)
}

deduplicate_species <- function(species_vector, type = c("numeric", "symbols"), symbol = "*") {
  # This sets the default to "numeric" (the first option)
  type <- match.arg(type)
  
  species_vector[species_vector == ""] <- NA
  
  dups <- unique(species_vector[duplicated(species_vector) & !is.na(species_vector)])
  
  dedup_data <- species_vector
  
  for (var in dups) {
    sub <- which(species_vector == var)
    l <- length(sub)
    
    if (type == "numeric") {
      dedup_sub <- paste0(var, "_", seq_len(l))
    } else if (type == "symbols") {
      dedup_sub <- paste0(var, sapply(seq_len(l), function(x) paste(rep(symbol, x), collapse = "")))
    }
    
    dedup_data[sub] <- dedup_sub
  }
  
  return(dedup_data)
}
