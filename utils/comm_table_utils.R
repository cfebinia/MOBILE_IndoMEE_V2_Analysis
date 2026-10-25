# custom functions
convert_to_CPM <- function(comm.df){
  comm.df <- comm.df*10^6
  round(comm.df,0)}

prevalence_filter <- function(comm.df, n, sample_as_columns=T) {
  # if n is integer, it refers to the number of minimum samples 
  # if n is a fraction, it refers to the nearest integer number by the total number of samples
  # comm.df is a matrix of abundance with sampleid as columns
  
  if(sample_as_columns==T){x <- as.matrix(comm.df)}
  if(sample_as_columns==F){x <- t(as.matrix(comm.df))}
  if(n<1){n <- round(as.numeric(n)*ncol(x),0)}
  y <- x
  y[y>0] <- 1
  y <- y[rowSums(y)>n,]
  x <- x[rownames(x) %in% rownames(y),]
  x
}

headx <- function(data.frame,n=NULL){
  if(is.null(n)){
    data.frame[1:5,1:5]
  }else{
    data.frame[1:n,1:n]
  }
}
