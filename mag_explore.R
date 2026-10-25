rm(list=ls())

refdir="Documents/Analysis/MOBILE/INDOMEE_refdata"
WD="Documents/Analysis/MOBILE/MOBILE_microbiome/"

setwd(refdir)

library(dplyr)
library(tidyr)
library(ggplot2)

source(file.path(WD,"utils/comm_table_utils.R"))
source(file.path(WD,"utils/parse_taxonomy_function.r"))

# input files
strongFilter <- file.path(refdir,"INDOMEE_V2_Abundance_Unified_draft3","bwa_counts_total_filtered_wMetadata_allUnified_draft3_strongFilter_relAbund.tsv")
weakFilter <- file.path(refdir,"INDOMEE_V2_Abundance_Unified_draft3","bwa_counts_total_filtered_wMetadata_allUnified_draft3_weakFilter_relAbund.tsv")
indv_sample_table <- file.path(WD,"richness/Indv_composite_allUnified_draft3_weakFilter_relAbund.tsv")

mapping_rates           <- read.delim(file.path(refdir, "MOBILE_Read_Counts","allUnified_draft3_strongFilter_mapStats.tsv"))
mapping_rates_strong    <- mapping_rates %>%
  rename(seqID = sample)

ids_wide0 <- read.delim(file.path(WD,"richness/ids_wide.tsv"), stringsAsFactors = F)
metadata <- read.delim("Consolidated_SampleID_combined_final_SeqIDmatched.tsv")

all_seqID <- strsplit(readLines(weakFilter, n = 1), "\t")[[1]][-c(1:25)]
controls <- all_seqID[grep("DHC|DCT|DEN|DOT|QHV|ZB|NTC|POS|MOCK|NTC|HC", all_seqID)]
outliers <- data.frame(read.delim("Outliers_list.txt") %>% filter(!Sample %in% controls))$Sample
replicates <- read.delim("replicate_list.tsv")$seqID
pilot_reseq <- all_seqID[grepl("MAL|BRU", substr(all_seqID, 1, 3)) & nchar(all_seqID) > 15]
pilot_seqID <- all_seqID[grepl("MAL|BRU|ASM|BA0|BAP|PBSHF", all_seqID)]

longitudinal_pairs <- read.delim(file.path("MOBILE_Social_Networks/processed/replicate_and_longitudinal_seqID_pairs.tsv"))  %>% filter(relation=="longitudinal")

bio_mobile <- read.delim("2023_MOBILE_AnthropometryData - CORRECTED.tsv") %>% 
  rename(Sample=sampleid, bmi=BMI) %>% 
  select(Sample, sex, age, w_kg, h_cm, bmi, population) %>%
  mutate(individual=Sample,
         population=substr(Sample,1,3)) %>%
  group_by(population, sex) %>%
  mutate(across(c(age, w_kg, h_cm), \(x) if_else(is.na(x), median(x, na.rm = TRUE), x)))%>%
  ungroup()

# load data
comm.df_clean <- list()
for(i in 1:2){
  if(i == 1){
    df <- read.delim(strongFilter)
    n = "strongFilter"
  } else {
    df <- read.delim(weakFilter)
    n = "weakFilter"
  }

  print(n)
  rownames(df) <- df$Genome
  df <- t(df[,-c(1:25)])

  message(paste("Raw n Sample =", nrow(df)))
  message(paste("Raw n MAGs =", ncol(df)))
  
  df[is.na(df)] <- 0
  df <- df[,colSums(df)>0]
  
  df <- df[!rownames(df) %in% controls,]
  
  message(paste("Raw n Sample (no ctrl) =", nrow(df)))
  message(paste("Raw n MAGs (no ctrl) =", ncol(df)))
  
  df <- df[!rownames(df) %in% replicates,]
  df <- df[!rownames(df) %in% pilot_reseq,]
  df <- df[rowSums(df)>0,]
  
  print(dim(df))
  message(paste("Filtered n Sample =", nrow(df)))
  message(paste("Filtered n MAGs =", ncol(df)))
  
  comm.df_clean[[i]] <- df
  names(comm.df_clean) <- n
}

df_ind <- read.delim(indv_sample_table)

ids_wide <- metadata %>%
  select(individual,population, host) %>%
  distinct() %>%
  rename(Sample=individual) %>%
  left_join(ids_wide0, by="Sample") %>%
  rowwise() %>%
  mutate(study="MOBILE",
         missing_sites = sum(is.na(c_across(c(Fecal, Skin, Oral))))) %>%
  ungroup() 

tmp <- data.frame(Sample=substr(pilot_seqID,1,3),
                  population=substr(pilot_seqID,1,3),
                  host="human",
                  Fecal=pilot_seqID, Skin=NA, Oral=NA, 
                  study="IndoMEE_V1", missing_sites=2) %>%
  filter(!Fecal %in% replicates) %>%
  distinct() %>%
  mutate(Sample = case_when(
    substr(Fecal,1,7) == "MALSUL3" ~ substr(gsub("MALSUL3","SUL003",Fecal),1,6),
    substr(Fecal,1,6) == "MALSUL" ~ substr(gsub("MAL","",Fecal),1,6),
    substr(Fecal,1,6) == "MALHUL" ~ substr(gsub("MAL","",Fecal),1,6),
    substr(Fecal,1,6) == "MALRPN" ~ substr(gsub("MAL","",Fecal),1,6),
    substr(Fecal,1,6) == "ASMDAI" ~ substr(gsub("DAI","",Fecal),1,6),
    substr(Fecal,1,6) == "BRUSBG" ~ substr(gsub("BRU","",Fecal),1,6),
    substr(Fecal,1,3) == "BAP" ~ substr(gsub("BAP","P",Fecal),1,6),
    substr(Fecal,1,3) == "BA0" ~ substr(gsub("BA0","DPS0",Fecal),1,6),
    substr(Fecal,1,3) == "PBS" ~ substr(gsub("PBSHF","PBS",Fecal),1,6),
    TRUE ~ Sample)) %>%
  mutate(population = substr(Sample,1,3)) %>%
  mutate(population = gsub("PBS","BTU", population),
         population = gsub("SUL","APT", population),
         population = gsub("RPN","TBU", population),
         population = gsub("SBG","BSP", population))

ids_wide <- rbind(ids_wide,tmp) %>%
  distinct() %>%
  select(-population) %>%
  mutate(host=factor(host, c("human","dog"))) %>%
  arrange(host, missing_sites, Sample) %>%
  filter(!Fecal %in% replicates) %>%
  filter(!Oral %in% replicates) %>%
  filter(!Skin %in% replicates)

colSums(is.na(ids_wide))
table(ids_wide$study, ids_wide$host)
colSums(!is.na(ids_wide[,c("Fecal","Oral","Skin")]))

write.table(ids_wide,"Consolidated_SampleID_combined_final_SeqIDmatched_wPilot_wide.tsv")


# Antrhop data
bio_pilot <- read.delim("../MOBILE_microbiome/richness/fecal_diversity_metadata_wPilot.tsv") %>% 
  select(individual, sex, age, w_kg, h_cm, bmi) %>% 
  distinct() %>% 
  filter(individual %in% pilot_seqID) %>%
  mutate(Sample = case_when(
    substr(individual,1,7) == "MALSUL3" ~ substr(gsub("MALSUL3","SUL003",individual),1,6),
    substr(individual,1,6) == "MALSUL" ~ substr(gsub("MAL","",individual),1,6),
    substr(individual,1,6) == "MALHUL" ~ substr(gsub("MAL","",individual),1,6),
    substr(individual,1,6) == "MALRPN" ~ substr(gsub("MAL","",individual),1,6),
    substr(individual,1,6) == "ASMDAI" ~ substr(gsub("DAI","",individual),1,6),
    substr(individual,1,6) == "BRUSBG" ~ substr(gsub("BRU","",individual),1,6),
    substr(individual,1,3) == "BAP" ~ substr(gsub("BAP","P",individual),1,6),
    substr(individual,1,3) == "BA0" ~ substr(gsub("BA0","DPS0",individual),1,6),
    substr(individual,1,3) == "PBS" ~ substr(gsub("PBSHF","PBS",individual),1,6),
    TRUE ~ individual))  %>%
  mutate(population = substr(Sample,1,3)) %>%
  mutate(population = gsub("PBS","BTU", population),
         population = gsub("SUL","APT", population),
         population = gsub("RPN","TBU", population),
         population = gsub("SBG","BSP", population)) %>%
  select(all_of(colnames(bio_mobile)))

bio <- rbind(bio_mobile, bio_pilot) %>%
  mutate(transition=case_when(
    population %in% c("ASM","BTU","ORS","ORT","HUL") ~ "earlyTransition",
    population %in% c("APT","TBU","BSP") ~ "lateTransition",
    population %in% c("LDY","PDW","DPS") ~ "Agri_Urban",
    TRUE ~ NA_character_
  )) %>%
  mutate(bmi=round((w_kg/((h_cm/100)^2)),1),
         is_obese=ifelse(bmi>=27.5,1,0))

colSums(is.na(bio))

write.table(bio,"anthrop_data_median_imputed_processed_wPilot.tsv")

table(metadata$type)
table(metadata$site)


# taxonomy 
full.taxonomy <- read.delim(weakFilter) %>% select("Genome","classification")
MAGs <- full.taxonomy$Genome
full.taxonomy <- cbind(Genome=full.taxonomy[,"Genome"], parse_taxonomy(full.taxonomy))
full.taxonomy$Species2 <- deduplicate_species(full.taxonomy$Species)

select_MAGs <- unique(sort(c(colnames(comm.df_clean[["weakFilter"]]), 
                         colnames(comm.df_clean[["weakFilter"]]))))

full.taxonomy <- full.taxonomy %>% filter(Genome %in%  select_MAGs)

apply(full.taxonomy, 2, function(x) length(unique(x)))

full.taxonomy %>%
  group_by(Domain, Phylum) %>%
  summarise(n=n_distinct(Genome),
            pct=round(100*n/length(MAGs),2)) %>%
  arrange(Domain, desc(n))

full.taxonomy %>%
  group_by(Family) %>%
  summarise(n=n_distinct(Genome),
            pct=round(100*n/length(MAGs),2)) %>%
  arrange(desc(n))

full.taxonomy %>%
  group_by(Genus) %>%
  summarise(n=n_distinct(Genome),
            pct=round(100*n/length(MAGs),2)) %>%
  arrange(desc(n))

full.taxonomy %>%
  group_by(Species) %>%
  summarise(n=n_distinct(Genome),
            pct=round(100*n/length(MAGs),2)) %>%
  arrange(desc(n))

# body site
setwd(WD)
df <- comm.df_clean[["weakFilter"]]
upset_matrix <- read.delim("input_files/upset_matrix_nofilter.tsv",row.names = 1) 
upset_matrix_filtered <- read.delim("input_files/upset_matrix_minAbund0.001.tsv", row.names = 1)

mags_count_noFilter <- colnames(upset_matrix)
mags_count_Filtered <- colnames(upset_matrix_filtered)

colnames(upset_matrix) <- c("humanFecal","humanOral","humanSkin", "dogFecal", "dogOral", "dogSkin")
colnames(upset_matrix_filtered) <- c("humanFecal","humanOral","humanSkin", "dogFecal", "dogOral", "dogSkin")

MAGs_byType <- list(
  DHF = rownames(upset_matrix)[upset_matrix$humanFecal==1],
  DHV = rownames(upset_matrix)[upset_matrix$humanOral==1],
  DHP = rownames(upset_matrix)[upset_matrix$humanSkin==1],
  DDF = rownames(upset_matrix)[upset_matrix$dogFecal==1],
  DDV = rownames(upset_matrix)[upset_matrix$dogOral==1],
  DDH = rownames(upset_matrix)[upset_matrix$dogSkin==1],
  Fecal = unique(rownames(upset_matrix)[upset_matrix$humanFecal==1|upset_matrix$dogFecal==1]),
  Oral = unique(rownames(upset_matrix)[upset_matrix$humanOral==1|upset_matrix$dogOral==1]),
  Skin = unique(rownames(upset_matrix)[upset_matrix$humanSkin==1|upset_matrix$dogSkin==1]),
  All = unique(rownames(upset_matrix))
)

lapply(MAGs_byType, length)

MAGs_byHost <- list(
  Human = unique(rownames(upset_matrix)[upset_matrix$humanFecal==1|upset_matrix$humanOral==1|upset_matrix$humanSkin==1]),
  Dog = unique(rownames(upset_matrix)[upset_matrix$dogFecal==1|upset_matrix$dogOral==1|upset_matrix$dogSkin==1])
)
lapply(MAGs_byHost, length)

seqID_byType <- lapply(names(MAGs_byType)[1:7], function(x) {
  isample <- rownames(df)[grepl(x,substr(rownames(df),1,3))]
  if(x=="DHF"){isample <- c(isample, pilot_seqID)}
  return(isample)
})
names(seqID_byType) <- names(MAGs_byType)[1:7]
seqID_byType[["Fecal"]] <- unique(c(seqID_byType[["DHF"]], seqID_byType[["DDF"]]))
seqID_byType[["Oral"]] <- unique(c(seqID_byType[["DHV"]], seqID_byType[["DDV"]]))
seqID_byType[["Skin"]] <- unique(c(seqID_byType[["DHP"]], seqID_byType[["DDH"]]))
seqID_byType[["All"]] <- unique(unlist(seqID_byType))

archeal_list <- full.taxonomy$Genome[full.taxonomy$Domain=="Archaea"]

archeal_MAGs <- lapply(names(seqID_byType), function(x) {
  isample <- seqID_byType[[x]]
  idf <- df[rownames(df) %in% isample, ]
  idf <- idf[,colnames(idf) %in% archeal_list]
  idf <- idf[rowSums(idf)>0,]
  idf <- idf[,colSums(idf)>0]
  imags <- colnames(idf)
  return(imags)
})
names(archeal_MAGs) <- names(seqID_byType)
archeal_MAGs[["Fecal"]] <- unique(c(archeal_MAGs[["DHF"]], archeal_MAGs[["DDF"]]))
archeal_MAGs[["Oral"]] <- unique(c(archeal_MAGs[["DHV"]], archeal_MAGs[["DDV"]]))
archeal_MAGs[["Skin"]] <- unique(c(archeal_MAGs[["DHP"]], archeal_MAGs[["DDH"]]))

top5 <- lapply(names(MAGs_byType), function(x) {
  imag <- MAGs_byType[[x]]
  amag <- archeal_MAGs[[x]]
  
  out <- full.taxonomy %>%
    filter(Genome %in% c(imag, amag)) %>%
    group_by(Domain, Phylum) %>%
    summarise(
      n = n_distinct(Genome),
      .groups = "drop_last"
      ) %>%
    mutate(pct = round(100 * n / length(imag), 2)) %>%
    arrange(Domain, desc(n)) %>%
    slice_max(pct, n = 5, with_ties = FALSE) %>%
    ungroup() %>%
    select(Phylum, pct)
  
  return(out)
})
names(top5) <- names(MAGs_byType)

library(purrr)
combined_top5 <- reduce(top5, full_join, by = "Phylum")
colnames(combined_top5) <- c("Phylum", names(MAGs_byType))
combined_top5 <- combined_top5 %>%
  left_join(distinct(full.taxonomy[,c("Domain","Phylum")]), by="Phylum") %>%
  select(Domain, Phylum, DHF, DHV, DHP, DDF, DDV, DDH, Fecal, Oral, Skin, All) %>%
  arrange(Domain)  

write.table(combined_top5, file.path(refdir, "MOBILE_Summary_Stats", "phylum_spread.tsv"),
            quote = F, sep = "\t",row.names = F)


mag_mean <- colMeans(comm.df_clean[["weakFilter"]]*100)

data.frame(Genome = names(mag_mean), meanRA = mag_mean) %>%
  left_join(full.taxonomy[, c("Genome", "Domain", "Phylum")], by = "Genome") %>%
  group_by(Domain, Phylum) %>%
  summarise(
    sumRA = round(sum(meanRA, na.rm = TRUE), 2),
    .groups = "drop_last"
  ) %>%
  slice_max(sumRA, n = 5, with_ties = FALSE) %>%
  arrange(Domain, desc(sumRA), Phylum) %>%
  ungroup()

# by average mag relab:
top5_RA <- lapply(names(MAGs_byType), function(x) {
  imag <- unique(c(MAGs_byType[[x]], archeal_MAGs[[x]]))
  isample <- seqID_byType[[x]]
  
  idf <- comm.df_clean[["weakFilter"]] 
  idf <- idf[rownames(idf) %in% isample, colnames(idf) %in% imag, drop = FALSE]
  idf <- idf*100
  mag_mean <- colMeans(idf)
  
  out <- data.frame(Genome = names(mag_mean), meanRA = mag_mean) %>%
    left_join(full.taxonomy[, c("Genome", "Domain", "Phylum")], by = "Genome") %>%
    group_by(Domain, Phylum) %>%
    summarise(sumRA = round(sum(meanRA, na.rm = TRUE), 5), .groups = "drop_last") %>%
    slice_max(sumRA, n = 5, with_ties = TRUE) %>%
    arrange(Domain, desc(sumRA), Phylum) %>%
    ungroup() %>%
    select(-Domain)
    
  return(out)
})
names(top5_RA) <- names(MAGs_byType)

library(purrr)
combined_top5_RA <- reduce(top5_RA, full_join, by = "Phylum")
colnames(combined_top5_RA) <- c("Phylum", names(MAGs_byType))
combined_top5_RA <- combined_top5_RA %>%
  left_join(distinct(full.taxonomy[,c("Domain","Phylum")]), by="Phylum") %>%
  select(Domain, Phylum, DHF, DHV, DHP, DDF, DDV, DDH, Fecal, Oral, Skin, All) %>%
  arrange(Domain)  

write.table(combined_top5_RA, file.path(refdir, "MOBILE_Summary_Stats", "phylum_abund.tsv"),
            quote = F, sep = "\t",row.names = F)

# genus
# by average mag relab:
top5_RA <- lapply(names(MAGs_byType), function(x) {
  imag <- unique(c(MAGs_byType[[x]], archeal_MAGs[[x]]))
  isample <- seqID_byType[[x]]
  
  idf <- comm.df_clean[["weakFilter"]] 
  idf <- idf[rownames(idf) %in% isample, colnames(idf) %in% imag, drop = FALSE]
  idf <- idf*100
  mag_mean <- colMeans(idf)
  
  out <- data.frame(Genome = names(mag_mean), meanRA = mag_mean) %>%
    left_join(full.taxonomy[, c("Genome", "Domain", "Genus")], by = "Genome") %>%
    group_by(Domain, Genus) %>%
    summarise(sumRA = round(sum(meanRA, na.rm = TRUE), 5), .groups = "drop_last") %>%
    #slice_max(sumRA, n = 5, with_ties = TRUE) %>%
    #filter(sumRA>0.01)%>%
    arrange(Domain, desc(sumRA), Genus) %>%
    ungroup() %>%
    select(-Domain)
  
  return(out)
})
names(top5_RA) <- names(MAGs_byType)
combined_top5_RA <- reduce(top5_RA, full_join, by = "Genus")
colnames(combined_top5_RA) <- c("Genus", names(MAGs_byType))
combined_top5_RA <- combined_top5_RA %>%
  left_join(distinct(full.taxonomy[,c("Domain","Genus")]), by="Genus") %>%
  select(Domain, Genus, DHF, DHV, DHP, DDF, DDV, DDH, Fecal, Oral, Skin, All) %>%
  arrange(Domain)  

write.table(combined_top5_RA, file.path(refdir, "MOBILE_Summary_Stats", "genus_abund.tsv"),
            quote = F, sep = "\t",row.names = F)


idf <- comm.df_clean[["weakFilter"]] 
idf <- idf*100
mag_mean <- colMeans(idf)

out <- data.frame(Genome = names(mag_mean), meanRA = mag_mean) %>%
  left_join(full.taxonomy[, c("Genome", "Domain", "Genus")], by = "Genome") %>%
  group_by(Domain, Genus) %>%
  summarise(sumRA = round(sum(meanRA, na.rm = TRUE), 5), .groups = "drop_last") %>%
  #slice_max(sumRA, n = 5, with_ties = TRUE) %>%
  arrange(Domain, desc(sumRA), Genus) %>%
  ungroup() %>%
  select(-Domain) %>%
  arrange(desc(sumRA))
out %>% filter(Genus=="Bacteroides")
out %>% filter(Genus=="Phocaeicola")
out %>% filter(Genus=="Bifidobacterium")


idf <- comm.df_clean[["weakFilter"]] 
idf[idf > 0] <- 1
prev <- colSums(idf) / nrow(idf)

out <- data.frame(Genome = names(prev), prev = prev) %>%
  left_join(full.taxonomy[, c("Genome", "Domain", "Phylum", "Genus")], by = "Genome") %>%
  group_by(Domain, Phylum, Genus) %>%
  summarise(maxPrev = max(prev), .groups = "drop") %>%
  arrange(desc(maxPrev))

out
out %>% filter(Genus=="Bifidobacterium")
out %>% filter(Phylum=="Actinomycetota")

top5 <- lapply(names(MAGs_byType), function(x) {
  imag <- unique(c(MAGs_byType[[x]], archeal_MAGs[[x]]))
  isample <- seqID_byType[[x]]
  
  idf <- comm.df_clean[["weakFilter"]] 
  idf <- idf[rownames(idf) %in% isample, colnames(idf) %in% imag, drop = FALSE]
  idf[idf>0] <- 1
  idf <- t(idf)
  
  out <- data.frame(Genome = rownames(idf), idf, check.names = FALSE) %>%
    left_join(full.taxonomy[, c("Genome", "Domain", "Phylum")], by = "Genome") %>%
    group_by(Domain, Phylum) %>%
    summarise(across(where(is.numeric), ~ max(.x)), .groups = "drop") %>%
    mutate(Prevalence = 100*rowSums(across(where(is.numeric)))/ncol(idf)) %>%
    select(Domain, Phylum, Prevalence) %>% 
    arrange(Domain, desc(Prevalence))  %>%
    group_by(Domain) %>%
    slice_max(order_by = Prevalence, n=5) %>%
    ungroup() %>%
    select(-Domain)
  
  return(out)
})
names(top5) <- names(MAGs_byType)
combined_top5_RA <- reduce(top5, full_join, by = "Phylum")
colnames(combined_top5_RA) <- c("Phylum", names(MAGs_byType))
combined_top5_RA <- combined_top5_RA %>%
  left_join(distinct(full.taxonomy[,c("Domain","Phylum")]), by="Phylum") %>%
  select(Domain, Phylum, DHF, DHV, DHP, DDF, DDV, DDH, Fecal, Oral, Skin, All) %>%
  arrange(Domain)  
write.table(combined_top5_RA, file.path(refdir, "MOBILE_Summary_Stats", "phyla_preval.tsv"),
            quote = F, sep = "\t",row.names = F)

top5 <- lapply(names(MAGs_byType), function(x) {
  imag <- unique(c(MAGs_byType[[x]], archeal_MAGs[[x]]))
  isample <- seqID_byType[[x]]
  
  idf <- comm.df_clean[["weakFilter"]] 
  idf <- idf[rownames(idf) %in% isample, colnames(idf) %in% imag, drop = FALSE]
  idf[idf>0] <- 1
  idf <- t(idf)
  
  out <- data.frame(Genome = rownames(idf), idf, check.names = FALSE) %>%
    left_join(full.taxonomy[, c("Genome", "Domain", "Genus")], by = "Genome") %>%
    group_by(Domain, Genus) %>%
    summarise(across(where(is.numeric), ~ max(.x)), .groups = "drop") %>%
    mutate(Prevalence = 100*rowSums(across(where(is.numeric)))/ncol(idf)) %>%
    select(Domain, Genus, Prevalence) %>% 
    arrange(Domain, desc(Prevalence))  %>%
    group_by(Domain) %>%
    slice_max(order_by = Prevalence, n=5) %>% 
    ungroup() %>%
    select(-Domain)
  
  return(out)
})
names(top5) <- names(MAGs_byType)
combined_top5_RA <- reduce(top5, full_join, by = "Genus")
colnames(combined_top5_RA) <- c("Genus", names(MAGs_byType))
combined_top5_RA <- combined_top5_RA %>%
  left_join(distinct(full.taxonomy[,c("Domain","Genus")]), by="Genus") %>%
  select(Domain, Genus, DHF, DHV, DHP, DDF, DDV, DDH, Fecal, Oral, Skin) %>%
  arrange(Domain)  
write.table(combined_top5_RA, file.path(refdir, "MOBILE_Summary_Stats", "genus_preval.tsv"),
            quote = F, sep = "\t",row.names = F)


## Prevotella vs bacteroides
lapply(names(MAGs_byType), function(x) {
  imag <- unique(c(MAGs_byType[[x]], archeal_MAGs[[x]]))
  isample <- seqID_byType[[x]]
  
  idf <- comm.df_clean[["weakFilter"]] 
  idf <- idf[rownames(idf) %in% isample, colnames(idf) %in% imag, drop = FALSE]
  idf <- idf*100
  
  prev_mags <- full.taxonomy$Genome[full.taxonomy$Genome=="Prevotella"]
  
