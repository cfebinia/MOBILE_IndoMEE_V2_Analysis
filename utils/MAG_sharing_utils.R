convert_to_prevalence <- function(species_matrix){
  # species matrix with samples as colums and species as rows
  x <- species_matrix
  x[x>0] <- 1
  return(x)
}

filter_empties <- function(species_matrix) {
  x <- species_matrix
  x <- x[which(rowSums(x, na.rm = T)>0),] 
  return(x)
}

get_MAGs_matrix <- function(raw_table, sample.type=c("DHF","DHP","DHV","DDF","DDH","DDV"), genomes="Genome"){
  if(length(sample.type)>1){
    ss <- paste0(sample.type, collapse = "|")
    ss <- grep(pattern = ss, substr(colnames(raw_table),1,6))
  } else if(length(sample.type)==1){
    ss <- grep(pattern = sample.type, substr(colnames(raw_table),1,6))
  }
  x <- raw_table[,ss]
  rownames(x) <- raw_table[,genomes]
  isna <- names(which(apply(x,2, function(matrix) any(is.na(matrix)))))
  x <- x[,!colnames(x) %in% isna]
  x <- filter_empties(x)
  message(paste("Result: ",ncol(x),"samples selected."))
  message(paste(nrow(x),"species retained at prevalence >0."))
  return(x)
}

remove_repeats <- function(species_matrix){
  # species matrix with samples as colums and species as rows
  # if repeats are detected, we'll select only one sample with the least lexicographical value (usually the id without the "R" label)
  data <- species_matrix
  with_repeats <- unique(substr(colnames(data)[duplicated(substr(colnames(data),1,6))],1,6))
  with_repeats <- data.frame(sampleid=colnames(data)[grep(pattern = paste0(with_repeats,collapse = "|"), colnames(data))])
  with_repeats$SAMPLE <- substr(with_repeats$sampleid, 1,6)
  to_keep <- c()
  for(sample in unique(with_repeats$SAMPLE)){
    all_s <- with_repeats$sampleid[with_repeats$SAMPLE==sample]
    first <- min(all_s)
    to_keep <- c(to_keep,first)
  }
  excluded <- with_repeats$sampleid[!with_repeats$sampleid%in%to_keep]
  data <- filter_empties(cbind(data[,!colnames(data) %in% with_repeats$sampleid],
                               data[,to_keep]))
  data <- data[,order(colnames(data))]
  data <- filter_empties(data)
  message(paste("The following samples are kept: ",
                paste0(to_keep, collapse = ", ")))
  message(paste("The following samples were removed: ",
                paste0(excluded, collapse = ", ")))
  message(paste("Result: ",ncol(data),"samples retained."))
  message(paste("Result: ",nrow(data),"species retained at prevalence >0."))
  return(data)
}
