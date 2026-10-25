#####################
# utility functions #
#####################
calc_weighted_sharing <- function(g1, g2, prev_df, binary_data) {
  # This function calculates the MAG sharing matrix with Inverse-Prevalence Weighting (rewards rare species)
  # g1, g2, are pairs of sample type (e.g. humanFecal, humanOral) 
  # prev_df data contains prevalence of specific MAG in each sample types (see site_preval above)
  # prev_df must have columns named after the 'groups' above
  # binary_data, the presence/absence matrix of species by sampleid
  # The count data is based in a subset MAG list found in the two groups,  to insure symmetry
  
  samples_g1 <- dat %>% filter(host_site == g1) %>% pull(seqID) %>% unique() %>% sort()
  samples_g2 <- dat %>% filter(host_site == g2) %>% pull(seqID) %>% unique() %>% sort()
  
  combo_binary <- binary_data[,unique(c(samples_g1,samples_g2))]
  combo_binary <- filter_empties(combo_binary)
  relevant_species <- unique(rownames(combo_binary))
  
  weight_df <- prev_df %>%
    filter(Genome %in% relevant_species)%>%
    select(all_of(c("Genome",g1,g2))) %>%
    mutate(W = 1 / ((.data[[g1]] + epsilon) * (.data[[g2]] + epsilon)))
  
  weight_vec <- weight_df$W[match(relevant_species, weight_df$Genome)]
  
  mat_g1 <- as.matrix(combo_binary[, samples_g1, drop=FALSE])
  mat_g2 <- as.matrix(combo_binary[, samples_g2, drop=FALSE])
  
  # Weighted pairwise sharing score
  sharing_score <- t(mat_g1) %*% (mat_g2 * weight_vec)
  
  return(as.data.frame(as.table(as.matrix(sharing_score)),stringsAsFactors = F))
}

