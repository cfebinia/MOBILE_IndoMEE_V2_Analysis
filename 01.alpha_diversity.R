rm(list=ls())

library(ape)
library(phangorn)
library(picante)
library(abdiv)

# source("99.community_data_utils.R") # load custom scripts

basedir="Documents/Analysis/MOBILE/MOBILE_microbiome"
refdir="Documents/Analysis/MOBILE/INDOMEE_refdata"

setwd(basedir)
#################
## Input Files ##
#################
input_dir=file.path(refdir,"INDOMEE_V2_Abundance_Unified_draft3/")
input_file="bwa_counts_total_filtered_wMetadata_allUnified_draft3_weakFilter_relAbund.tsv"
output_dir="outfile/adiv"

dir.create(file.path(basedir,output_dir),showWarnings = F,recursive = T)

# sampleid
sampleid <- readRDS("input_files/sampleid_n116.rds")
sporeids <- readRDS("input_files/sampleid_spore.rds")

loop="no"
# dir <- paste0(input_dir,"/","original")
dir=input_dir
input_file <- file.path(input_dir,input_file)


if(loop=="no"){
  if(length(grep("comboRef", dir))==0){
    input_tree <- "input_files/taxonomy/hybrid/gtdbtk.bac120.user_msa.fasta.treefile"
  } else if (length(grep("comboRef", dir))==1){
    input_tree <- "input_files/taxonomy/comboRef/gtdbtk.bac120_wTax.unrooted.tree"
  }
  
  if(length(grep("comboRef", dir))==0){
    output_dir <- "output_files/adiv/hybrid"
  } else if (length(grep("comboRef", dir))==1){
    output_dir <- "output_files/adiv/comboRef"
  }
  
  assembly=basename(output_dir)
  
  if(length(grep("original", dir))==1){
    output_dir <- paste0(output_dir,"/original")
    dir.create(output_dir,recursive = T)
    iter.prefix=paste0(assembly,"_",basename(output_dir))
  } else if(length(grep("subsampled_iter0", dir))==1){
    output_dir <- paste0(output_dir,"/subsampled_iter0")
    dir.create(output_dir,recursive = T)
    iter.prefix=paste0(assembly,"_",basename(output_dir))
  } else if(length(grep("subsampled_iter1", dir))==1){
    output_dir <- paste0(output_dir,"/subsampled_iter1")
    dir.create(output_dir,recursive = T)
    iter.prefix=paste0(assembly,"_",basename(output_dir))
  } else if(length(grep("subsampled_iter2", dir))==1){
    output_dir <- paste0(output_dir,"/subsampled_iter2")
    dir.create(output_dir,recursive = T)
    iter.prefix=paste0(assembly,"_",basename(output_dir))
  }
  
  mypaths <- list()
  ipaths <- list(loop=c(list(input_file=input_file),
                        list(input_tree=input_tree),
                        list(output_dir=output_dir)))
  names(ipaths) <- iter.prefix
  mypaths <- append(mypaths, ipaths)
  
}

if(loop=="yes"){
  dir <- NULL
  input_file <- NULL
  mypaths <- list()
  folder_list<- list.files(input_dir)
  
  for(dir in folder_list){
    dir <- paste0(input_dir,"/",dir)
    file_list <- list.files(dir)
    n1 <- grep("strongFilter_relAbund",file_list)
    n2 <- grep("total",file_list)
    n <- intersect(n1, n2)
    if(length(grep("original", dir))>0){
      n3 <- grep("revised", file_list)
      n <- intersect(n, n3)
    }
    input_file <- paste0(dir,"/",file_list[n])
    
    if(length(grep("comboRef", dir))==0){
      input_tree <- "input_files/taxonomy/hybrid/gtdbtk.bac120.user_msa.fasta.treefile"
    } else if (length(grep("comboRef", dir))==1){
      input_tree <- "input_files/taxonomy/comboRef/gtdbtk.bac120_wTax.unrooted.tree"
    }
    
    if(length(grep("comboRef", dir))==0){
      output_dir <- "output_files/adiv/hybrid"
    } else if (length(grep("comboRef", dir))==1){
      output_dir <- "output_files/adiv/comboRef"
    }
    
    assembly=basename(output_dir)
    
    if(length(grep("original", dir))==1){
      output_dir <- paste0(output_dir,"/original")
      dir.create(output_dir,recursive = T,showWarnings = F)
      iter.prefix=paste0(assembly,"_",basename(output_dir))
    } else if(length(grep("subsampled_iter0", dir))==1){
      output_dir <- paste0(output_dir,"/subsampled_iter0")
      dir.create(output_dir,recursive = T,showWarnings = F)
      iter.prefix=paste0(assembly,"_",basename(output_dir))
    } else if(length(grep("subsampled_iter1", dir))==1){
      output_dir <- paste0(output_dir,"/subsampled_iter1")
      dir.create(output_dir,recursive = T,showWarnings = F)
      iter.prefix=paste0(assembly,"_",basename(output_dir))
    } else if(length(grep("subsampled_iter2", dir))==1){
      output_dir <- paste0(output_dir,"/subsampled_iter2")
      dir.create(output_dir,recursive = T,showWarnings = F)
      iter.prefix=paste0(assembly,"_",basename(output_dir))
    }
    
    ipaths <- list(loop=c(list(input_file=input_file),
                          list(input_tree=input_tree),
                          list(output_dir=output_dir)))
    names(ipaths) <- iter.prefix
    mypaths <- append(mypaths, ipaths)
    
  }
}


# Calculate alpha diversity and collate results
ndir=length(mypaths)

collate.adiv <- data.frame()
for(i in seq(ndir)){ 
  iter.prefix <- names(mypaths)[i]
  assembly <- unlist(strsplit(iter.prefix, split = "_"))[1]
  iteration <- unlist(strsplit(iter.prefix, split = "_"))[2]
  
  input_file=mypaths[[i]]$input_file
  input_tree=mypaths[[i]]$input_tree
  output_dir=mypaths[[i]]$output_dir
  
  print(iter.prefix)
  print(input_file)
  print(input_tree)
  print(output_dir)
  
  # load and filter abundance table
  df <- read.delim(input_file)
  reads <- df[,colnames(df) %in% c(sampleid,sporeids)]
  rownames(reads) <- df$original_bin
  comm.df <- convert_to_CPM(reads)
  table(colSums(comm.df)==0)
  dim(comm.df)
  
  # load tree
  tree <- read.tree(input_tree)
  
  # check if tree tip labels matches with rownames
  if(assembly=="hybrid"){
    rownames(comm.df) <- df$original_bin
  } else if(assembly=="comboRef"){
    rownames(comm.df) <- df$Original_bin_ID
  }
  
  if(all(tree$tip.label %in% rownames(comm.df)) == T){
    print("Tree tips match")
  } else if(all(tree$tip.label %in% rownames(comm.df)) == T){
    print("Tree tips match")
  } else {
    stop("Tree tips does not match. Stopping execution.")
  }
  
  # identify MAGs not present in trees (usually archaeas)
  table(rownames(comm.df) %in% tree$tip.label)
  not.present.table <- rownames(comm.df)[!rownames(comm.df) %in% tree$tip.label]
  comm.df.bacteria <- comm.df[rownames(comm.df) %in% tree$tip.label,]
  comm.df.bacteria <- comm.df.bacteria[tree$tip.label,] # reorder table to follow tree
  #comm.df.archaea <- comm.df[!rownames(comm.df) %in% tree$tip.label,] #  ignored for now
  
  # convert tree to usable format
  tree <- as.phylo(tree)
  tree <- midpoint(tree)
  
  ####################
  ## Alpha Diversiy ##
  ####################
  shannon <- as.numeric(apply(comm.df,2, function(x) abdiv::shannon(x,base = 10)))
  simpson <- as.numeric(apply(comm.df,2, function(x) abdiv::simpson(x)))
  invsimpson <- as.numeric(apply(comm.df,2, function(x) abdiv::invsimpson(x)))
  simpsonE <- as.numeric(apply(comm.df,2, function(x) abdiv::simpson_e(x)))
  dominance <- as.numeric(apply(comm.df,2, function(x) abdiv::strong(x)))
  observed <- vegan::specnumber(comm.df, MARGIN = 2)
  faith <- as.numeric(apply(comm.df.bacteria,2, function(x) abdiv::faith_pd(abdiv::match_to_tree(x, tree=tree), tree=tree)))
  
  adiv.df <- data.frame(sampleid = colnames(comm.df),
                        observed, shannon, simpson, invsimpson, simpsonE, dominance, faith)
  
  adiv.df$subsample <- iter.prefix
  collate.adiv <- rbind(collate.adiv,adiv.df)
}

# WRITE OUTPUT
write.table(collate.adiv, file = "input_files/alpha_diversity_subsampled_wFD.tsv", quote = F, row.names = F,sep = "\t")
saveRDS(mypaths, file="input_files/alpha_diversity_subsampled_paths.rds")
