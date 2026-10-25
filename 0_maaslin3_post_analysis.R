rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)

# Set directories
basedir="Documents/Analysis/MOBILE/MOBILE_microbiome"
refdir="Documents/Analysis/MOBILE/INDOMEE_refdata"

setwd(basedir)

# Select analysis
with_pilot <- TRUE
is_genus <- FALSE

hostTypes <- c("humanFecal", "humanSkin", "humanOral", "dogFecal", "dogSkin", "dogOral")
# for non-iterative analysis
# humanFecal, humanSkin, humanOral, dogFecal, dogSkin, dogOral
hostType <- "humanFecal"

# select populations
pop.ord <- c("BTU","ORT","ORS","APT","TBU","BSP")
if(with_pilot == TRUE){
  pop.ord <- c("ASM","BTU","ORT","ORS","HUL","APT","TBU","BSP","LDY","PDW","DPS")
}

# load full taxon:
full.taxonomy <- read.delim(file.path(basedir,"input_files/unified_draft3_mobile_combined_taxonomy.tsv"))
full.maglist <- unique(full.taxonomy$Genome)
full.genuslist <- unique(full.taxonomy$Genus)
phylum.stats <- full.taxonomy %>% 
  group_by(Phylum) %>%
  summarise(nFamily=n_distinct(Family),
            nGenus=n_distinct(Genus),
            nSpecies=n_distinct(Species),
            nMAG=n_distinct(Genome)) %>%
  arrange(desc(nMAG))
phylum.ord <- as.character(unique(phylum.stats$Phylum))

family.stats <- full.taxonomy %>% 
  mutate(Phylum=factor(Phylum, phylum.ord)) %>%
  group_by(Phylum,Family) %>%
  summarise(nGenus=n_distinct(Genus),
            nSpecies=n_distinct(Species),
            nMAG=n_distinct(Genome)) %>%
  arrange(Phylum, desc(nMAG))
family.ord <- as.character(unique(family.stats$Family))

phylum_family_list <- distinct(family.stats[,c("Phylum","Family")])

dir.create(file.path(basedir,"outfile","MAASLIN3_cline"),recursive = T,showWarnings = F)
write.table(phylum.stats, "outfile/MAASLIN3_cline/phylum_n_species.tsv",quote = F,sep = "\t",row.names = F)
write.table(family.stats, "outfile/MAASLIN3_cline/family_n_species.tsv",quote = F,sep = "\t",row.names = F)
write.table(phylum_family_list, "outfile/MAASLIN3_cline/family_list_unique.tsv",quote = F,sep = "\t",row.names = F)


#################
# 1. Diet Assoc #
#################
inputdir=file.path(basedir, "input_files/MAASLIN3_onlyhumanFecal/collect/")
target_dir <- list.dirs(inputdir)
all_output <- list.files(path = target_dir, 
                         pattern = ".*all_results.tsv$", 
                         full.names = FALSE)
all_output <-  unlist(lapply(strsplit(all_output, split="/"), function(x) x))
if(is_genus==T){
  all_output <- all_output[grep("Genus", all_output)]
} else {
  all_output <- all_output[grep("MAGs", all_output)]
}

hostType <- "humanFecal"
ranef="Ranef" # "noranef" or "Ranef"
test="dietZ"
confounders <- c("clean_nReads","age","sex")

food.ord <- c("fish","white_meat","red_meat","honey","rice",
              "legumes","tuber","vegetables","fruits","eggs","noodles",
              "cooking_oil","coffee_tea","sugar_g","added_sugar",
              "lard","sago")

dropped.food <- c("added_sugar")
food.ord <- food.ord[!food.ord %in% dropped.food]

select_files <- all_output[grep(ranef, all_output)]
diet_files <- select_files[grep(test, select_files)]

infile <- diet_files[grep(hostType, diet_files)]

df <- read.delim(file.path(inputdir,infile)) %>%
  filter(is.na(error)) %>%
  filter(!is.na(as.numeric(null_hypothesis))) %>%
  filter(!is.na(as.numeric(coef))) %>%
  mutate(coef=as.numeric(coef)) %>%
  filter(value %in% food.ord)

# for(imodel in c("abundance","prevalence","joint")){
  imodel = "joint"
  
  if(imodel != "joint"){tmp <- df[grep(imodel, df$model),]}
  
  if(imodel == "joint"){
    tmp <- df %>%
      group_by(feature, metadata, value) %>%
      slice_min(order_by = pval_individual, n = 1, with_ties = FALSE) %>%
      mutate(pval_individual=pval_joint) %>%
      ungroup()
  }
  
  print(paste0(c(imodel,"valid MAGs:",length(unique(tmp$feature))),collapse = " "))
  
  signif <- tmp[which(tmp$pval_individual<0.05),]
  table(signif$value)
  
  heatdata <- data.frame(hostType, 
                         level=ifelse(is_genus==T, "Genus","MAGs"),
                         tmp[!is.na(tmp$pval_individual),c("feature","Species2","value","coef","stderr","pval_individual","qval_individual")])
  heatdata <- subset(heatdata, value %in% food.ord)
  heatdata$value <- factor(heatdata$value, food.ord)
  heatdata$t_statistic <- heatdata$coef/heatdata$stderr
  heatdata$inv_variance <- 1/(heatdata$stderr^2)
  heatdata$weighted_coef <- heatdata$coef*heatdata$inv_variance
  heatdata$sign <- sign(heatdata$coef)
  heatdata$z_score <- qnorm(1 - (heatdata$pval_individual / 2)) * sign(heatdata$coef)
  heatdata$direction <- ifelse(heatdata$sign >= 0, "Positive (+)", "Negative (-)")
  heatdata$direction <- factor(heatdata$direction, levels = c("Negative (-)", "Positive (+)"))
  heatdata$absZ <- abs(heatdata$z_score)
  
  ###########
  # Density #
  ###########
  ratio_data <- heatdata %>%
    count(value, direction) %>%
    pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
    mutate(
      Neg = `Negative (-)` / (`Negative (-)` + `Positive (+)`),
      Pos = `Positive (+)` / (`Negative (-)` + `Positive (+)`)
    )
  
  tmp <- heatdata %>%
    left_join(select(ratio_data, value, Neg, Pos), by = "value")
  
  density_data <- tmp %>%
    group_by(value, direction, Neg, Pos) %>%
    reframe(
      x = density(z_score)$x,
      y = density(z_score)$y
    ) %>%
    mutate(
      weighted_y = case_when(
        direction == "Positive (+)" ~ y * Pos,
        direction == "Negative (-)" ~ y * Neg,
        TRUE ~ NA_real_
      )
    ) %>%
    select(-Neg, -Pos)
  
  library(dplyr)
  
  peak_data <- density_data %>%
    group_by(value, direction) %>%
    summarise(
      peak_x = x[which.max(weighted_y)],
      peak_y = max(weighted_y, na.rm = TRUE),
      y_signif = if(direction[1] == "Negative (-)") max(weighted_y[which(x < -1.64)], na.rm = TRUE) else min(weighted_y[which(x > 1.64)], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      x_signif = ifelse(direction == "Positive (+)", 1.64, -1.64),
      nudge_x_val = ifelse(direction == "Negative (-)", -0.2, 0.2),
      hjust_val = ifelse(direction == "Negative (-)", 1, 0)
    ) %>%
    left_join(count(heatdata, value, direction), by = c("value", "direction"))
  # critical value for x_signif and y_signif at Zscore=1.64 for one-tailed test. 1.96 for two-tailed test
  
  p <- ggplot(density_data, aes(x = x, y = weighted_y, fill = direction)) +
    facet_wrap(~value, ncol = 4) +
    geom_area(alpha = 0.6, colour = "grey30", linewidth = 0.5, position = "identity") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.6) +
    #geom_segment(data = peak_data, 
    #             aes(x=x_signif, xend= x_signif, yend = y_signif, y=0), colour = "grey40",
    #             linetype = "solid", linewidth = 0.6, show.legend = FALSE) +
    geom_segment(data = peak_data, 
                 aes(x = peak_x, xend = peak_x, y = 0, yend = peak_y, colour = direction), 
                 linetype = "dashed", linewidth = 0.6, show.legend = FALSE) +
    geom_text(data = peak_data, 
              aes(x = peak_x + nudge_x_val, 
                  y = ifelse(peak_y > 2.5, peak_y - 1, peak_y - 0.5), 
                  label = round(peak_x, 2), 
                  colour = direction, 
                  hjust = hjust_val), 
              vjust = -0.5, size = 3, show.legend = FALSE) +
    scale_fill_manual(values = c("Negative (-)" = "deepskyblue4", "Positive (+)" = "firebrick3")) +
    scale_colour_manual(values = c("Negative (-)" = "deepskyblue4", "Positive (+)" = "firebrick3")) +
    labs(
      title = "Distribution of MaAsLin3 Z-Scores by Direction",
      subtitle = paste("Test:", imodel),
      x = "Z-score (from p-value)",
      y = "Weighted Density",
      fill = "Effect Direction"
    ) +
    #scale_x_continuous(limits = c(-1, 1) * max(abs(heatdata$z_score), na.rm = TRUE)) +
    theme_minimal(base_size = 14) +
    theme(
      plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"),
      panel.spacing = unit(0.2, "cm"),
      legend.position = "bottom",
      legend.title = element_text(size = 10, vjust = 1, face = "bold"),
      legend.text = element_text(size = 10),
      legend.key.size = unit(0.4, "cm"),
      panel.grid.minor = element_blank()
    )
  
  print(p)
  
  outpath=file.path(basedir,"figout/MAASLIN3")
  outfile=paste0(c(hostType,test,imodel,"ZscoreDensity.pdf"),collapse = "_")
  ggsave(plot=p, filename=file.path(outpath, outfile), device = "pdf",
         width = 15,height = 15,units = "cm",dpi = 1200, scale=1.5)
  print(file.path(outpath, outfile))

  ###########
  # Heatmap #
  ###########
  ggdata <- full.taxonomy %>%
    rename(feature = Genome) %>%
    filter(feature %in% unique(heatdata$feature)) %>%
    select(feature, Phylum, Family, Genus, Species) %>%
    right_join(heatdata, by="feature",keep = FALSE)%>%
    arrange(z_score,value,Species)
  
  # selecting for
  to_display=50
  if(imodel=="abundance"){
    p_tresh=0.001
    fam_cutoff=5
    spec_cutoff=1
    gen_cutoff=3
    value_cutoff=5}
  if(imodel=="prevalence"){
    p_tresh=0.001
    fam_cutoff=5
    spec_cutoff=1
    gen_cutoff=3
    value_cutoff=5}
  if(imodel=="joint"){
    p_tresh=0.001
    fam_cutoff=5
    spec_cutoff=1
    gen_cutoff=3
    value_cutoff=5}

  top_selected <- ggdata %>% 
    filter(pval_individual < p_tresh) %>%
    group_by(value,  direction, Species) %>% 
    slice_max(order_by = absZ, n = spec_cutoff, with_ties = FALSE) %>% ungroup() %>% 
    group_by(value, direction, Family) %>% 
    slice_max(order_by = absZ, n = fam_cutoff, with_ties = FALSE) %>% ungroup() %>% 
    group_by(value, direction, Genus) %>% 
    slice_max(order_by = absZ, n = gen_cutoff, with_ties = FALSE) %>% ungroup() %>% 
    group_by(value, direction) %>% 
    slice_max(order_by = absZ, n = value_cutoff, with_ties = FALSE) %>% ungroup() %>% 
    slice_max(order_by = absZ, n = to_display, with_ties = FALSE) %>% 
    ungroup() %>% 
    arrange(hostType,value, direction, Phylum, Family, desc(absZ))
  top_selected <- unique(top_selected$feature)
  
  selected_features <- top_selected
    
  ggdata <- ggdata %>%
    filter(feature %in% selected_features) %>%
    select(feature, Phylum, Family, Genus, Species, Species2, value, pval_individual, qval_individual, z_score, direction) 
  
  unique(ggdata$feature) %>% length()
  
  ggdata$lab_signif <- ifelse(ggdata$pval_individual < p_tresh, "*","" )
  
  ggdata_complete <- ggdata %>%
    tidyr::complete(value, nesting(Phylum, Family, Genus, Species2)) %>%
    mutate(
      value = droplevels(as.factor(value)),
      Phylum = droplevels(factor(Phylum, levels = phylum.ord)),
      Family = droplevels(factor(Family, levels = family.ord)),
      Genus = factor(Genus, levels = unique(Genus)),
      Species2 = factor(Species2, levels = unique(Species2)),
      lab_signif = ifelse(is.na(lab_signif), "", as.character(lab_signif))
    )
  
  p <- ggplot(ggdata_complete) +
    geom_tile(aes(x = value, y = Species2, fill = z_score)) +
    facet_grid(Family ~ ., scales = "free_y", space = "free_y") +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_discrete(limits = rev, expand = c(0, 0)) +
    scale_fill_gradient2(
      low = "#1E90FF",      
      mid = "#FFFFFF",      
      high = "#FF6347",     
      midpoint = 0,
      na.value = "grey",
      name = "Z-Score"
    ) +
    geom_text(aes(x = value, y = Species2, label = lab_signif), size = 4) +
    labs(
      title = "Significant Dietary Associations",
      subtitle = paste("Test:", imodel)
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 10, vjust = 1, face = "bold"),
      legend.text = element_text(size = 10),
      legend.key.size = unit(0.4, "cm"),
      panel.border = element_rect(colour = "grey40", linewidth = 0.2),
      panel.spacing = unit(0, "cm"),
      strip.text.y = element_text(angle = 0, hjust = 0, face = "bold", size = 10),
      axis.title = element_blank(),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(
        angle = 90, hjust = 0, vjust = 0.5, size = 10,
        colour = ifelse(levels(ggdata_complete$value) %in% c("sago", "lard"), "grey50", "black"),
        face = ifelse(levels(ggdata_complete$value) %in% c("sago", "lard"), "italic", "bold")
      )
    )
  print(p)
  
  outpath=file.path(basedir,"figout/MAASLIN3")
  outfile=paste0(c(hostType,test,imodel,"Heatmap.pdf"),collapse = "_")
  ggsave(plot=p, filename=file.path(outpath, outfile), device = "pdf",
         width = 15,height = 15,units = "cm",dpi = 1200, scale=1.5)
#}

# diet vs transition/community
inputdir="Documents/Analysis/MOBILE/MOBILE_microbiome/input_files/MAASLIN3_onlyhumanFecal"
read_maaslin_result <- function(infile, filter=NULL){
  x <- read.delim(file.path(inputdir,infile)) %>%
    filter(is.na(error)) %>%
    filter(!is.na(as.numeric(null_hypothesis))) %>%
    filter(!is.na(as.numeric(coef))) %>%
    mutate(coef=as.numeric(coef))
  if(!is.null(filter)){
    x <- x %>% filter(value %in% filter)
  }
  
  # select top by p_joint
  x <- x %>%
    filter(pval_joint < 0.05) %>%
    group_by(feature, value) %>%
    slice_min(order_by = pval_individual, n=1) %>%
    ungroup()
      
  return(x)
  }

diet_infile <- "diet_age_sex_noranef/all_results.tsv"
transition_infile <- "cline_age_sex_noranef_p0.05/all_results.tsv"
population_infile <- "population_age_sex_noranef_p0.05/all_results.tsv"
obese_infile <- "obese_age_sex_noranef_p0.05/all_results.tsv"

df_diet <- read_maaslin_result(diet_infile, filter=food.ord)
df_cline <- read_maaslin_result(transition_infile, filter="cline")
df_pop <- read_maaslin_result(population_infile, filter=pop.ord)
df_obese <- read_maaslin_result(obese_infile, filter="obese")

assoc_MAGs=lapply(
  list(Diet=df_diet, Cline=df_cline, Comm=df_pop, Obesity=df_obese),
  function(x) {return(unique(x$feature))})

library(ggvenn)

p <- ggvenn(
  assoc_MAGs,
  fill_color = c("#4DAF4A", "maroon", "royalblue","tan2"),
  stroke_size = 0.5,
  set_name_size = 4,
  text_size = 3
  )

length(unique(unlist(assoc_MAGs)))
unlist(lapply(assoc_MAGs, length)) * 100 /3260
1-(c(165,27,531,17)/unlist(lapply(assoc_MAGs, length)))
unlist(lapply(assoc_MAGs, length)) - c(165,27,531,17)
outpath=file.path(basedir,"figout/MAASLIN3")
outfile="overlap_maaslin3.pdf"
ggsave(plot=p, filename = file.path(outpath, outfile),
       width = 15, height = 10, device = "pdf", dpi = 150, units = "cm")

#####################
# 2. Obesity Assoc #
#####################
inputdir=file.path(basedir, "input_files/MAASLIN3_cline/")
target_dir <- list.dirs(inputdir)
all_output <- list.files(path = target_dir, 
                         pattern = ".*all_results.tsv$", 
                         full.names = FALSE)
all_output <-  unlist(lapply(strsplit(all_output, split="/"), function(x) x))
if(is_genus==T){
  all_output <- all_output[grep("Genus", all_output)]
} else {
  all_output <- all_output[grep("MAGs", all_output)]
}

hostType <- "humanFecal"
ranef="noranef" # "noranef" or "Ranef"
# test="Obesity_cline_wPilot"
test="obese_wPilot"
# test="obese_age"
# default is noranef + obese_wPilot

select_files <- all_output[grep(test, all_output)]
select_files <- select_files[grep(ranef, select_files)]

infile <- select_files[grep(hostType, select_files)]

confounders <- c("clean_nReads","study")
value.ord <- c("obese","cline","age","sex")

df <- read.delim(file.path(inputdir,infile)) %>%
  filter(is.na(error)) %>%
  filter(!is.na(as.numeric(null_hypothesis))) %>%
  filter(!is.na(as.numeric(coef))) %>%
  filter(metadata %in% value.ord) %>%
  mutate(coef=as.numeric(coef)) %>%
  mutate(coef=ifelse(metadata=="cline", -1*coef, coef)) # adjusted because BSP is the baselin in Pop

for(imodel in c("abundance","prevalence", "joint")){
    #imodel="joint"
  
    if(imodel != "joint"){tmp <- df[grep(imodel, df$model),]}
    
    if(imodel == "joint"){
      tmp <- df %>%
        group_by(feature, metadata, value) %>%
        slice_min(order_by = pval_individual, n = 1, with_ties = FALSE) %>%
        mutate(pval_individual=pval_joint) %>%
        ungroup()
    }
    
    print(paste0(c(imodel,"valid MAGs:",length(unique(tmp$feature))),collapse = " "))
    
    heatdata <- data.frame(level=ifelse(is_genus==T, "Genus","MAGs"),
                           tmp[!is.na(tmp$pval_individual),c("feature","Species2","metadata","value","coef","stderr","pval_individual","qval_individual")])
    
    heatdata <- heatdata %>%
      mutate(value = factor(metadata, value.ord),
             p_log = -log10(pval_individual),
             q_log = -log10(qval_individual),
             sign = sign(as.numeric(coef)))
    
    heatdata$t_statistic <- heatdata$coef/heatdata$stderr
    heatdata$inv_variance <- 1/(heatdata$stderr^2)
    heatdata$weighted_coef <- heatdata$coef*heatdata$inv_variance
    
    heatdata$pval_individual[heatdata$pval_individual<1e-6] <- 1e-6
    heatdata$qval_individual[heatdata$qval_individual<1e-6] <- 1e-6
    
    heatdata$z_score <- qnorm(1 - (heatdata$pval_individual / 2)) * sign(heatdata$coef)
    heatdata$absZ <- abs(heatdata$z_score)
    heatdata$direction <- ifelse(heatdata$sign >= 0, "Positive (+)", "Negative (-)")
    heatdata$direction <- factor(heatdata$direction, levels = c("Negative (-)", "Positive (+)"))
    rownames(heatdata) <- NULL
    
    p_tresh = 0.05
    species_cutoff = 1
    family_cutoff = 10
    top_selected <- full.taxonomy %>%
        rename(feature=Genome) %>%
        select(feature, Phylum, Family, Genus, Species) %>%
        distinct() %>%
        full_join(heatdata, by="feature",keep = F) %>%
        filter(pval_individual<p_tresh) %>%
        mutate(absZ = abs(z_score)) %>%
        group_by(Species, direction, metadata) %>% 
        slice_max(order_by = absZ, n = species_cutoff, with_ties = FALSE) %>%
        ungroup() %>% 
        group_by(Phylum, Family, metadata) %>% 
        slice_max(order_by = absZ, n = family_cutoff, with_ties = FALSE) %>%
        ungroup() %>% 
        arrange(value,desc(absZ))
      
    selected_features <- unique(top_selected$feature)
    heatdata <- subset(heatdata, feature %in% selected_features)
    
    bmi_signif <- unique(heatdata$feature[heatdata$metadata==value.ord[1]&heatdata$pval_individual<p_tresh])
    cline_signif <- unique(heatdata$feature[heatdata$metadata==value.ord[2]&heatdata$pval_individual<p_tresh])
    age_signif <- unique(heatdata$feature[heatdata$metadata==value.ord[3]&heatdata$pval_individual<p_tresh])
    sex_signif <- unique(heatdata$feature[heatdata$metadata==value.ord[4]&heatdata$pval_individual<p_tresh])
    
    bmi_only <- bmi_signif[!bmi_signif %in% c(cline_signif,age_signif,sex_signif)]
    
    bmi_cline <- intersect(bmi_signif,cline_signif)
    bmi_age <- intersect(bmi_signif,age_signif)
    bmi_sex <- intersect(bmi_signif,sex_signif)
    
    cline_age <- intersect(cline_signif,age_signif)
    cline_sex <- intersect(cline_signif,sex_signif)
    
    age_sex <- intersect(age_signif,sex_signif)
    
    heatdata <- heatdata %>%
        mutate(assoc=case_when(
        metadata %in% value.ord[1] & pval_individual<p_tresh ~ "Obese",
        metadata %in% value.ord[2] & pval_individual<p_tresh ~ "Lifestyle",
        metadata %in% value.ord[3] & pval_individual<p_tresh ~ "Age",
        metadata %in% value.ord[4] & pval_individual<p_tresh ~ "Male",
        TRUE ~ "Not significant"
        ))
    
    heatdata$assoc <- factor(heatdata$assoc, c("Obese","Lifestyle","Age", "Male", "Not significant"))
    
    count_assoc <- heatdata %>%
      group_by(assoc) %>%
      summarise(n_signif=n_distinct(feature))
    
    assoc_labs <- paste0(count_assoc$assoc, "\n(", count_assoc$n_signif, ")")
    names(assoc_labs) <-count_assoc$assoc
    assoc_labs
    
    # plot co-factors separately
    ggdata <- heatdata %>%
      arrange(assoc, desc(round(p_log,3)), desc(abs(coef)))
    
    y_tresh <- c("BMI"=0.85,"Lifestyle"=0.95,"Age"=0.95,"Sex"=0.95)
    if(imodel=="joint"){y_tresh <- c("BMI"=0.80,"Lifestyle"=0.90,"Age"=0.90,"Sex"=0.95)}
    
    quantile_bmi <- quantile(heatdata$coef[heatdata$assoc == "Obese"&heatdata$pval_individual<p_tresh], c(0.05,0.95), na.rm = TRUE)
    quantile_ls <- quantile(heatdata$coef[heatdata$assoc == "Lifestyle" & heatdata$pval_individual < p_tresh], c(0.025,0.975), na.rm = TRUE)
    quantile_age <- quantile(heatdata$coef[heatdata$assoc == "Age"&heatdata$pval_individual<p_tresh], c(0.025,0.975), na.rm = TRUE)
    quantile_sex <- quantile(heatdata$coef[heatdata$assoc == "Male"&heatdata$pval_individual<p_tresh], c(0.025,0.975), na.rm = TRUE)
    quantile_y_bmi <- quantile(heatdata$p_log[heatdata$assoc == "Obese"&heatdata$pval_individual<p_tresh], y_tresh["BMI"], na.rm = TRUE)
    quantile_y_ls <- quantile(heatdata$p_log[heatdata$assoc == "Lifestyle"&heatdata$pval_individual<p_tresh],y_tresh["Lifestyle"], na.rm = TRUE)
    quantile_y_age <- quantile(heatdata$p_log[heatdata$assoc == "Age"&heatdata$pval_individual<p_tresh], y_tresh["Age"], na.rm = TRUE)
    quantile_y_sex <- quantile(heatdata$p_log[heatdata$assoc == "Male"&heatdata$pval_individual<p_tresh], y_tresh["Sex"], na.rm = TRUE)
    
    ggdata <- ggdata %>%
      mutate(
        repel_labs = case_when(
          assoc == "Not significant" ~ "",
          assoc == "Obese" & p_log > quantile_y_bmi ~ as.character(Species2),
          assoc == "Lifestyle" & p_log > quantile_y_ls  ~ as.character(Species2),
          assoc == "Age" & p_log > quantile_y_age ~ as.character(Species2),
          assoc == "Male" & p_log > quantile_y_sex  ~ as.character(Species2),
          assoc == "Obese" & (coef < quantile_bmi[1] | coef > quantile_bmi[2]) ~ as.character(Species2),
          assoc == "Lifestyle" & (coef < quantile_ls[1] | coef > quantile_ls[2]) ~ as.character(Species2),
          assoc == "Age" & (coef < quantile_age[1] | coef > quantile_age[2]) ~ as.character(Species2),
          assoc == "Male" & (coef < quantile_sex[1] | coef > quantile_sex[2]) ~ as.character(Species2),
          TRUE ~ ""
          ))
    
    ggdata %>% group_by(assoc) %>% summarise(n_labs=n_distinct(repel_labs))
    
    # Annotate
    ggdata$repel_labs[ggdata$feature%in%bmi_cline & ggdata$repel_labs!=""] <- paste0(ggdata$repel_labs[ggdata$feature%in%bmi_cline & ggdata$repel_labs!=""], "\u1D38")
    ggdata$repel_labs[ggdata$feature%in%bmi_age & ggdata$repel_labs!=""] <- paste0(ggdata$repel_labs[ggdata$feature%in%bmi_age & ggdata$repel_labs!=""], "\u1D2C")
    ggdata$repel_labs[ggdata$feature%in%bmi_sex & ggdata$repel_labs!=""] <- paste0(ggdata$repel_labs[ggdata$feature%in%bmi_sex & ggdata$repel_labs!=""], "\u1D39")
    ggdata$repel_labs[ggdata$feature%in%cline_age & ggdata$repel_labs!=""] <- paste0(ggdata$repel_labs[ggdata$feature%in%cline_age & ggdata$repel_labs!=""], "\u1D2C")
    ggdata$repel_labs[ggdata$feature%in%cline_sex & ggdata$repel_labs!=""] <- paste0(ggdata$repel_labs[ggdata$feature%in%cline_sex & ggdata$repel_labs!=""], "\u1D39")
    ggdata$repel_labs[ggdata$feature%in%age_sex & ggdata$repel_labs!=""] <- paste0(ggdata$repel_labs[ggdata$feature%in%age_sex & ggdata$repel_labs!=""], "\u1D39")
    ggdata$repel_labs <- gsub("\u1D39\u1D39","\u1D39",ggdata$repel_labs)
    ggdata$repel_labs <- gsub("\u1D2C\u1D2C","\u1D2C",ggdata$repel_labs)
    
    # Shorten
    ggdata$repel_labs <- gsub(pattern = "species","sp", ggdata$repel_labs,ignore.case = F)
    
    lab_ord <- ggdata %>%
      filter(repel_labs != "") %>%
      arrange(desc(abs(z_score)))
    lab_ord <- c(unique(lab_ord$repel_labs),"")
    ggdata$repel_labs <- factor(ggdata$repel_labs, lab_ord)
    
    ggdata$metadata <- factor(ggdata$metadata, value.ord)
    
    library(ggrepel)
    library(Cairo)
    
    xlim <- round(max(abs(ggdata$coef)),0) +1
    ylim <- round(max(abs(ggdata$p_log)),0) +1
    
    p_volcano <- ggplot(ggdata, aes(x = coef, y = p_log)) +
      facet_wrap(~ metadata, ncol = 2, axes = "all", scales="free_y",
                 labeller = labeller(metadata = c("obese"="Obesity", "cline"="Lifestyle","age"="Age","sex"="Male")))+
      geom_point(aes(color = assoc), alpha = 0.7, size = 2.5) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.5) +
      geom_vline(xintercept = 0, linetype = "dotted", color = "grey40", linewidth = 0.5) +
      scale_color_manual(
        values = c(
          "Obese" = "#B22222",
          "Lifestyle" = "#458B74", 
          "Age" = "#E066FF",
          "Male" = "#777777",
          "Not significant" = "grey85"),
        name = "Association", 
        labels=assoc_labs) +
      labs(
        x = "Coefficient",
        y = expression(-log[10] ~ (P-Value))
      ) +
      theme_bw(base_size = 14) +
      theme(
        strip.background = element_blank(),
        strip.text = element_text(face="bold"),
        legend.position = "bottom",
        legend.text = element_text(size=12),
        panel.grid.minor = element_blank()
      )
    
    q <- p_volcano + 
      geom_text_repel(seed = 1,direction = "both",
                      aes(label = as.character(repel_labs), col=assoc), 
                      min.segment.length = 0,
                      size = 3.5, fontface = "italic", 
                      max.overlaps = ifelse(imodel=="abundance",100, 20),
                      box.padding = unit(2,"mm"),
                      point.padding = 0, 
                      segment.size = 0.2, 
                      segment.color = "grey50",
                      show.legend = FALSE) +
      scale_y_continuous(expand = c(0.15,0.15), breaks=seq(0,ylim, by=1))+
      scale_x_continuous(breaks = seq(-xlim,xlim,by=5),expand = c(0.1,0.1),
                         limits = c(-xlim,xlim))
    
    print(q)
    
    
    outpath=file.path(basedir,"figout/MAASLIN3")
    dir.create(outpath, showWarnings = F,recursive = T)
    outfile=paste0(c(hostType,ifelse(is_genus==F,"MAGs","Genus"),
                     paste(test,ranef,sep="_"),
                     imodel,
                     "BMI_Volcano.pdf"),
                   collapse = "_")
    ggsave(plot=q, filename=file.path(outpath, outfile), device = cairo_pdf,
           width = 15,height = 15,units = "cm",dpi = 1200, scale=2)
}

# count overlap
# Define significance threshold (e.g., p < 0.05 or q < 0.05)
p_thresh <- 0.05

obese_features <- unique(ggdata$feature[ggdata$assoc == "Obese" & ggdata$pval_individual < p_tresh])
lifestyle_features <- unique(ggdata$feature[ggdata$assoc == "Lifestyle" & ggdata$pval_individual < p_tresh])
age_features <- unique(ggdata$feature[ggdata$assoc == "Age" & ggdata$pval_individual < p_tresh])
sex_features <- unique(ggdata$feature[ggdata$assoc == "Male" & ggdata$pval_individual < p_tresh])

total_obese <- length(obese_features)

overlap_lifestyle <- length(intersect(obese_features, lifestyle_features)) / total_obese * 100
overlap_age <- length(intersect(obese_features, age_features)) / total_obese * 100
overlap_sex <- length(intersect(obese_features, sex_features)) / total_obese * 100

obese_only <- setdiff(obese_features, c(lifestyle_features, age_features, sex_features))
pct_obese_only <- length(obese_only) / total_obese * 100

arr <- c("Obesity-only", "Obesity+Lifestyle", "Obesity+Age", "Obesity+Sex")
overlap_data <- data.frame(
  Category = factor(arr, arr),
  Percentage = c(pct_obese_only, overlap_lifestyle, overlap_age,  overlap_sex)
)

p <- ggplot(overlap_data, aes(x = Category, y = Percentage, fill = Category)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("Obesity-only" = "#B22222",
                               "Obesity+Lifestyle" = "#458B74",
                               "Obesity+Age"= "#E066FF",
                               "Obesity+Sex:M"="#777777")) +
  labs(
    x = "Association Category",
    y = "Overlapping Features (%)",
    title = "Microbial Associations and Overlap with Obesity"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  expand_limits(y = c(0, 75))

print(p)

outpath=file.path(basedir,"figout/MAASLIN3")
outfile="obesity_overlap_barplot_from_data.pdf"
ggsave(file.path(outpath,outfile), plot = p, width = 15, height = 6, units = "cm", dpi = 300, scale = 1.6)


##############################
# 3. Lifestyle vs Population #
##############################
setwd(basedir)
inputdir <- file.path(basedir, "input_files/MAASLIN3_cline/")
hostType <- NULL

all_output <- list.files(path = inputdir, 
                         pattern = ".*all_results.tsv$", 
                         full.names = FALSE)
if(is_genus==T){
  all_output <- all_output[grep("Genus", all_output)]
} else {
  all_output <- all_output[grep("MAGs", all_output)]
}

ranef="noranef" # or: "Ranef"
test_1="cline"
test_2="population"

select_files <- all_output[grep(ranef, all_output)]
ls_files <- select_files[grep(test_1, select_files)]
pop_files <- select_files[grep(test_2, select_files)]

# load input
output_cline <- lapply(ls_files, function(x) {
  tmp <- read.delim(file.path(inputdir,x)) %>%
    filter(is.na(error)) %>%
    filter(!is.na(as.numeric(null_hypothesis))) %>%
    filter(!is.na(as.numeric(coef))) %>%
    filter(metadata==test_1) %>%
    mutate(coef=-1*as.numeric(coef)) %>% # adjusted because BSP is the baseline for Pop
    group_by(feature, metadata, value) %>%
    slice_min(order_by = pval_individual, n = 1, with_ties = FALSE) %>%
    ungroup()
  return(tmp)}
  )

output_Pop <- lapply(pop_files, function(x) {
  tmp <- read.delim(file.path(inputdir,x)) %>%
    filter(is.na(error)) %>%
    filter(!is.na(as.numeric(null_hypothesis))) %>%
    filter(!is.na(as.numeric(coef))) %>%
    filter(metadata==test_2) %>%
    mutate(coef=as.numeric(coef)) %>%
    group_by(feature, metadata, value) %>%
    slice_min(order_by = pval_individual, n = 1, with_ties = FALSE) %>%
    ungroup()
  return(tmp)}
)

hostTypes <- unlist(lapply(ls_files, function(x) unlist(strsplit(x,"_"))[1]))
hostTypes[grep("wPilot",ls_files)] <- paste0("humanFecal_wPilot")

names(output_cline) <- hostTypes
names(output_Pop) <- hostTypes

# load taxonomy
MAG_list <- readRDS("input_files/MAASLIN3_cline/MAASLIN3_cline_MAG_list.rds")
MAG_list_cline <- MAG_list[grep(test_1, names(MAG_list))]
MAG_list_Pop <- MAG_list[grep(test_2, names(MAG_list))]

names(MAG_list_cline) <- unlist(lapply(strsplit(gsub("_age",";",gsub("_cline|_population", "", names(MAG_list_cline))),split = ";"),function(x) x[[1]]))
names(MAG_list_Pop) <- unlist(lapply(strsplit(gsub("_age",";",gsub("_cline|_population", "", names(MAG_list_Pop))),split = ";"),function(x) x[[1]]))


taxonomy <- lapply(hostTypes, function(x) {
  i <- MAG_list_cline[[grep(x, names(MAG_list_cline),fixed = T)[1]]]
  j <- MAG_list_Pop[[grep(x, names(MAG_list_Pop),fixed = T)[1]]]
  k <- unique(c(i,j))
  
  tmp <- full.taxonomy %>% 
    rename(feature=Genome) %>%
    filter(feature %in% k)
  
  return(tmp)
  })
names(taxonomy) <- hostTypes
lapply(taxonomy, nrow)

# count samples
sampleid_list <- read.delim("input_files/MAASLIN3_cline/MAASLIN3_cline_input_stats.tsv")
sampleid_list <- sampleid_list[grep(paste(c(test_1,test_2), collapse = "|"),sampleid_list$test_name),]
sampleid_list$type[grep("wPilot",sampleid_list$test_name)] <- paste("humanFecal","wPilot",sep = "_")
sample_count <- sampleid_list %>%
  select(type, n_sample) %>%
  distinct() %>%
  mutate(labs = paste0(type, "\n(", n_sample, ")"))

sample_labs <- sample_count$labs
names(sample_labs) <- sample_count$type


## Filter significance
signif_cline <- lapply(hostTypes, function(x) {
  tmp <- subset(output_cline[[x]], pval_joint<0.05)
  tmp <- tmp %>% filter(metadata == test_1)
  return(tmp)
  })

signif_Pop <- lapply(hostTypes, function(x) {
  tmp <- subset(output_Pop[[x]], pval_joint<0.05)
  tmp <- tmp %>% filter(metadata == test_2)
  return(tmp)
})

names(signif_cline) <- hostTypes
names(signif_Pop) <- hostTypes

# count significance
signifMAG_cline <- lapply(hostTypes, function(x) unique(signif_cline[[x]]$feature))
signifMAG_Pop <- lapply(hostTypes, function(x) unique(signif_Pop[[x]]$feature))

names(signifMAG_cline) <- hostTypes
names(signifMAG_Pop) <- hostTypes

# Intersect Analysis and plot venn
library(ggvenn)

magList <- data.frame()
magCount <- data.frame()
vennPlots <- list()
for(type in hostTypes){
  #type <- "humanFecal"
  venndf <- list(signifMAG_cline[[type]],signifMAG_Pop[[type]])  
  features <- taxonomy[[type]]$feature
  nfeature <- length(features)
  outside <- features[!features %in% unique(unlist(venndf))]
  noutside <- length(outside)
  names(venndf) <- c("Transition","Population")
  annot <- paste(type, "(not_signif=",noutside, "; nTotal=",nfeature,")")
      pven <- ggvenn(venndf, # show_percentage = F,
                     fill_color = c("maroon", "royalblue"),
                     stroke_size = 0.5,text_size = 5, set_name_size = 5) + 
        theme(plot.background = element_rect(fill = "white", colour = "NA"),
              plot.margin = margin(10,10,10,10,unit = "mm"),
              plot.title = element_text(margin = margin(b = 5))) + 
        labs(title = annot)
      print(pven)
      
      outdir=file.path(basedir,"figout/MAASLIN3/venn")
      dir.create(outdir,recursive = T,showWarnings = F)
      
      figout=paste(type,"Venn.pdf",sep = "_")
      ggsave(plot = pven, filename = file.path(outdir, figout), device = "pdf",
             width = 15,height = 10,units = "cm",dpi = 600, scale = 1)
      vennPlots[[type]] <- pven
      
      # list intersect
      tr <- venndf[["Transition"]]
      pop <- venndf[["Population"]]
      inter <- intersect(tr,pop)
      
      taxout <- data.frame(hostSite=type, taxonomy[[type]], association="none")
      taxout[taxout$feature %in% tr,"association"] <- "transition_only"
      taxout[taxout$feature %in% pop,"association"] <- "population_only"
      taxout[taxout$feature %in% inter,"association"] <- "intersect"
      taxout$association <- factor(taxout$association, c("transition_only", "population_only", "intersect","none"))
      magList <- rbind(magList, taxout)
      
      # output counts
      out <- data.frame(table(taxout$association), hostType=type)
      colnames(out)[1] <- "association"
      magCount <- rbind(magCount, out)
      }

# write outputs
outdir=file.path(basedir,"figout/MAASLIN3/venn")
dir.create(outdir,recursive = T,showWarnings = F)
write.table(magList, file.path(outdir,"collated_signif_maglist.tsv"), quote = F,sep = "\t", row.names = F)
write.table(magCount, file.path(outdir, "signif_count.tsv"), quote = F,sep = "\t", row.names = F)
saveRDS(vennPlots,file.path(outdir, "venn.rds"))


## plot barplot
library(ggplot2)
ggdata <- magCount
ggdata$association <- factor(ggdata$association, c("transition_only", "intersect","population_only", "none"))

nTotal <- unlist(lapply(1:length(taxonomy), function(x) length(unique(taxonomy[[x]][,"feature"]))))
nTotal <- data.frame(hostType=names(taxonomy), nTotal)

ggdata <- merge.data.frame(ggdata, nTotal, by="hostType")

ggdata$Prop <- round(100 * ggdata$Freq / ggdata$nTotal,1)

ggdata$hostType <- factor(ggdata$hostType ,
                          c("humanFecal", "humanFecal_wPilot", "humanOral", "humanSkin", 
                            "dogFecal", "dogOral", "dogSkin"))

p <- ggdata %>%
  filter(hostType != "humanFecal") %>%
  ggplot(aes(x = hostType, y = Prop, fill = association)) +
  geom_col(position = "stack") +
  geom_text(aes(label = ifelse(Prop > 0, paste0(Freq, "\n(", Prop, "%)"), "")),
            position = position_stack(vjust = 0.5),size = 4) +
  scale_fill_manual(values = c("maroon","mediumpurple","royalblue","grey")) +
  scale_x_discrete(labels=sample_labs)+
  labs(caption="MAG ~ Lifestyle_cline + Age + Sex + Seq_Depth VS.  MAG ~ population + Age + Sex + Seq_Depth") +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(),
        legend.position = "bottom")

p
outdir=file.path(basedir,"figout/MAASLIN3/venn")
dir.create(outdir,recursive = T,showWarnings = F)
ggsave(plot = p, filename = file.path(outdir, "barplot.pdf"), device = "pdf", 
         width = 15, height = 9,units = "cm", dpi = 600, scale = 1.5)

###################
## 4. UPSET PLOT ##
###################
setwd(basedir)
library(dplyr)
library(tidyr)
library(ComplexUpset)

# prepare data
all_output <- list.files(path = inputdir, 
                         pattern = ".*all_results.tsv$", 
                         full.names = FALSE)

if(is_genus==T){
  all_output <- all_output[grep("Genus", all_output)]
} else {
  all_output <- all_output[grep("MAGs", all_output)]
}

ranef="Ranef" # or: "Ranef"
test_1="cline"

select_files <- all_output[grep(ranef, all_output)]
ls_files <- select_files[grep(test_1, select_files)]

hostTypes <- unlist(lapply(ls_files, function(x) unlist(strsplit(x,"_"))[1]))
hostTypes[grep("wPilot",ls_files)] <- paste0("humanFecal_wPilot")

output_cline <- lapply(ls_files, function(x) {
  tmp <- read.delim(file.path(inputdir,x)) %>%
    filter(is.na(error)) %>%
    filter(!is.na(as.numeric(null_hypothesis))) %>%
    filter(!is.na(as.numeric(coef))) %>%
    filter(metadata==test_1) %>%
        group_by(feature, metadata, value) %>%
    slice_min(order_by = pval_individual, n = 1, with_ties = FALSE) %>%
    ungroup()
  return(tmp)}
)
names(output_cline) <- hostTypes

signif_cline <- lapply(hostTypes, function(x) {
  tmp <- subset(output_cline[[x]], pval_joint<0.05)
  tmp <- tmp %>% filter(metadata == test_1)
  return(tmp)
})
names(signif_cline) <- hostTypes

signif_Tranef <- signif_cline[names(signif_cline) != "humanFecal"] # drop the MOBILE-only dataset
tmp <- do.call(rbind, lapply(names(signif_Tranef), function(x) data.frame(hostSite = x, signif_Tranef[[x]])))
tmp <- subset(tmp, metadata == "cline")

upsetdf <- as.data.frame(table(tmp$feature, tmp$hostSite), stringsAsFactors = FALSE)
colnames(upsetdf) <- c("feature", "hostSite", "Freq")

upsetdf_wide <- upsetdf %>%
  mutate(
    feature  = as.character(feature),
    hostSite = trimws(as.character(hostSite)),
    Freq     = as.integer(as.numeric(Freq) > 0)
  ) %>%
  group_by(feature, hostSite) %>%
  summarise(Freq = max(Freq, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from  = hostSite,
    values_from = Freq,
    values_fill = list(Freq = 0)
  )

upsetdf$hostSite <- factor(upsetdf$hostSite,hostTypes)
set_cols <- setdiff(names(upsetdf_wide), "feature")

# plot
p <- ComplexUpset::upset(
  upsetdf_wide,
  intersect = set_cols,
  base_annotations = list(
    'Intersection size' = intersection_size(
      counts = TRUE, text = list(colour = "black"),
      fill = "#FF6347"
    )
  )
) +
  labs(caption = "MAG ~ Lifestyle + Age + Sex + clean_nReads + (1|community)")

print(p)

figdir="figout/MAASLIN3/venn/"
ggsave(plot = p, filename = paste0(figdir,"upset.pdf", collapse = ""), device = "pdf",
         width = 15, height = 7, units = "cm", dpi = 600, scale = 1.5)


top5_signif <- data.frame()
ivar <- c("humanFecal_wPilot", "humanOral", "humanSkin", "dogFecal", "dogOral", "dogSkin")
for(i in seq_along(ivar)){
  nm <- ivar[[i]]
  
  tmp <- signif_cline[[nm]] %>%
    left_join(
      full.taxonomy %>%
        rename(feature = Genome) %>%
        select(feature, Species),
      by = "feature"
    ) %>%
    mutate(direction = ifelse(coef > 0, "Positive", "Negative")) %>%
    group_by(feature) %>%
    slice_min(order_by = pval_joint, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    group_by(Species) %>%
    slice_min(order_by = pval_joint, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    group_by(direction, model) %>%
    filter(qval_individual < 0.05) %>%
    slice_min(order_by = qval_individual, n = 10) %>%
    ungroup() %>%
    mutate(source = nm) %>%
    select(source, feature, Species2, coef, pval_individual, qval_individual, model, direction)
    
  
  top5_signif <- rbind(top5_signif, tmp)
}

top5_signif
write.table(top5_signif,file = "outfile/MAASLIN3_cline/top_features_across_sites.tsv", quote = F, row.names = F, sep="\t")

################################
# 5. Heatmap Top 5 Cpop Zscore #
################################
setwd(basedir)
inputdir <- file.path(basedir, "input_files/MAASLIN3_cline/")
hostType <- NULL

all_output <- list.files(path = inputdir, 
                         pattern = ".*all_results.tsv$", 
                         full.names = FALSE)
if(is_genus==T){
  all_output <- all_output[grep("Genus", all_output)]
} else {
  all_output <- all_output[grep("MAGs", all_output)]
}

test_1="cline"
test_2="population"

ls_files <- all_output[intersect(grep(test_1, all_output), grep("Ranef",all_output))]
pop_files <- all_output[intersect(grep(test_2, all_output), grep("noranef",all_output))]

hostTypes <- unlist(lapply(ls_files, function(x) unlist(strsplit(x,"_"))[1]))
hostTypes[grep("wPilot",ls_files)] <- paste0("humanFecal_wPilot")

# load input
output_cline <- lapply(ls_files, function(x) {
  tmp <- read.delim(file.path(inputdir,x)) %>%
    filter(is.na(error)) %>%
    filter(!is.na(as.numeric(null_hypothesis))) %>%
    filter(!is.na(as.numeric(coef))) %>%
    filter(metadata==test_1) %>%
    mutate(coef=-1*as.numeric(coef)) # adjusted because BSP is the baseline for pop
  return(tmp)}
)

output_Pop <- lapply(pop_files, function(x) {
  tmp <- read.delim(file.path(inputdir,x)) %>%
    filter(is.na(error)) %>%
    filter(!is.na(as.numeric(null_hypothesis))) %>%
    filter(!is.na(as.numeric(coef))) %>%
    filter(metadata==test_2) %>%
    mutate(coef=1*as.numeric(coef)) 
  return(tmp)}
)

names(output_cline) <- hostTypes
names(output_Pop) <- hostTypes

all_Tranef <- lapply(hostTypes, function(x) {
  tmp <- as.data.frame(output_cline[[x]]) %>%
    mutate(coef=as.numeric(coef))
  tmp$t_statistic <- tmp$coef/tmp$stderr
  tmp$inv_variance <- 1/(tmp$stderr^2)
  tmp$weighted_coef <- tmp$coef*tmp$inv_variance
  tmp$sign <- sign(tmp$coef)
  tmp$z_score <- qnorm(1 - (tmp$pval_individual / 2)) * sign(tmp$coef)
  cPop <- subset(tmp, metadata == test_1)
  
  tmp <- as.data.frame(output_Pop[[x]]) %>%
    mutate(coef=as.numeric(coef))
  tmp$t_statistic <- tmp$coef/tmp$stderr
  tmp$inv_variance <- 1/(tmp$stderr^2)
  tmp$weighted_coef <- tmp$coef*tmp$inv_variance
  tmp$sign <- sign(tmp$coef)
  tmp$z_score <- qnorm(1 - (tmp$pval_individual / 2)) * sign(tmp$coef)
  Pop <- subset(tmp, metadata == test_2)
  
  tmp <- rbind(cPop, Pop)
  tmp$direction <- ifelse(tmp$sign >= 0, "Positive (+)", "Negative (-)")
  tmp$direction <- factor(tmp$direction, levels = c("Negative (-)", "Positive (+)"))
  
  if(is_genus == FALSE){tmp <- tmp %>% select(-Species2, -classification)}
  
  tmp <- tmp %>%
    group_by(feature, metadata, value, model) %>%
    mutate(
      dominant_coef = coef[which.min(pval_individual)],
      z_score_joint = qnorm(1 - (pval_joint / 2)) * sign(dominant_coef)
    ) %>%
    ungroup() %>%
    select(-dominant_coef) %>%
    arrange(feature,metadata,value)
  
  # take the min pval of the two models (abund|prev)
  tmp <- tmp %>%
    group_by(feature, metadata, value) %>%
    slice_min(order_by = pval_individual, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  tmp %>%
    count(feature, metadata, value, model) %>%
    filter(n > 1)
  
  tmp <- full.taxonomy %>%
    rename(feature=Genome) %>%
    select(feature,Phylum,Family,Genus,Species,Species2) %>%
    right_join(tmp, by="feature",keep = FALSE) %>%
    arrange(feature,metadata,value)
  
  with_pilot <- length(x[grep("Pilot",x)]) >0
  
  value.ord <- c("cline","ASM","BTU","ORT","ORS","HUL","APT","TBU","BSP","LDY","PDW","DPS")
  
  tmp$value <- factor(tmp$value, value.ord)
  
  tmp <- data.frame(hostType=x, tmp)
})

all_Tranef <- do.call(rbind, all_Tranef)

# PLOT HEATMAP
# selecting for
p_tresh=0.05
spec_cutoff=1
fam_cutoff=15
gen_cutoff=5
value_cutoff=5
to_display=25 # select the top n associated feature

pop_heatmap <- list()
for(i in 1:length(hostTypes)){
  type = hostTypes[i]
  
  ggdata <- all_Tranef %>% 
    filter(hostType == type) %>%
    mutate(absZ=abs(z_score_joint))
  
  if(type=="humanFecal"){fam_cutoff=10; to_display=10}
  if(type=="humanFecal_wPilot"){fam_cutoff=10;value_cutoff=3; to_display=5}
  
  top_selected <- ggdata %>% 
    filter(pval_joint < p_tresh) %>%
    group_by(value,  direction, Species) %>% 
    slice_max(order_by = absZ, n = spec_cutoff, with_ties = FALSE) %>% ungroup() %>% 
    group_by(value, direction, Family) %>% 
    slice_max(order_by = absZ, n = fam_cutoff, with_ties = FALSE) %>% ungroup() %>% 
    group_by(value, direction, Genus) %>% 
    slice_max(order_by = absZ, n = gen_cutoff, with_ties = FALSE) %>% ungroup() %>% 
    group_by(value,direction) %>% 
    slice_max(order_by = absZ, n = value_cutoff, with_ties = FALSE) %>% ungroup() %>% 
    group_by(value) %>%
    slice_max(order_by = absZ, n = to_display, with_ties = FALSE) %>% 
    ungroup() %>% 
    arrange(hostType,value, direction, Phylum, Family, desc(absZ))
  top_selected <- unique(top_selected$feature)
  
  length(top_selected)
  
  ggdata <- ggdata %>% 
    filter(feature %in% top_selected) %>%
    arrange(desc(absZ)) %>%
    mutate(lab_signif=ifelse(z_score_joint>2,"*",""))%>%
    mutate(Phylum=factor(Phylum, phylum.ord),
           Family=factor(Family, family.ord),
           Genus=factor(Genus, unique(Genus)),
           Species2=factor(Species2, unique(Species2)))
  
  # cap zscore to 95% percentile value
  p95 <- quantile(ggdata$z_score_joint[!is.infinite(ggdata$z_score_joint)], c(0.05,0.95))
  ggdata[is.infinite(ggdata$z_score_joint)&ggdata$direction=="Negative (-)","z_score_joint"] <- p95[1]
  ggdata[is.infinite(ggdata$z_score_joint)&ggdata$direction=="Positive (+)","z_score_joint"] <- p95[2]
  ggdata[ggdata$z_score_joint<p95[1],"z_score_joint"] <- p95[1]
  ggdata[ggdata$z_score_joint>p95[2],"z_score_joint"] <- p95[2]
  
  # flatten non-significant (this replace low z scores with 1 and preserve the direction)
  ggdata[ggdata$pval_joint>0.1,"z_score_joint"] <- (ggdata[ggdata$pval_joint>0.1,"z_score_joint"]/ggdata[ggdata$pval_joint>0.1,"z_score_joint"]) * sign(ggdata[ggdata$pval_joint>0.1,"z_score_joint"])
  ggdata_complete <- ggdata %>%
    mutate(value = droplevels(value)) %>%
    tidyr::complete(value, nesting(Species2, Family))
  
  z_limit <- ceiling(max(abs(ggdata_complete$z_score_joint)))
  z_signif <- round(min(ggdata$absZ[ggdata$pval_joint<0.05]),0)

  p <- ggplot(ggdata_complete) +
    geom_tile(aes(x=value, y=Species2, fill=z_score_joint)) +
    facet_grid(Family ~ ., scales="free_y", space="free_y") +
    scale_x_discrete(position = "top", expand = c(0,0)) +
    scale_y_discrete(limits = rev, expand = c(0,0))+
    scale_fill_gradient2(
      low = "#1E90FF",      
      mid = "#FFFFFF",      
      high = "#FF6347",     
      midpoint = 0,
      na.value = "grey",
      name = "Z-Score (joint)",
      limits=c(-z_limit,z_limit),
      breaks=c(0,-z_limit,-z_signif,z_signif,z_limit))+
    #geom_text(aes(x = value, y = Species2, label = lab_signif), size=4)+
    labs(title = type) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size=10, vjust=1, face="bold"),
      legend.text = element_text(size=10),
      legend.key.size = unit(0.4, "cm"),
      panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.2),
      panel.spacing = unit(0, "cm"),
      strip.text.y = element_text(angle=0, hjust=0, face="bold", size=10),
      axis.title = element_blank(),
      axis.text.y= element_text(size=10),
      axis.text.x= element_text(angle=90, hjust=0, vjust=0.5, size=10, face="bold"))
  print(p)
  
  pop_heatmap[[i]] <- p
  names(pop_heatmap)[[i]] <- type
}

names(pop_heatmap)

outpath=file.path(basedir,"figout/MAASLIN3/Heatmap")
dir.create(outpath, showWarnings = F,recursive = T)
for(type in names(pop_heatmap)){
  outfile=paste0(c(type,"Heatmap.pdf"),collapse = "_")
  ggsave(plot=pop_heatmap[[type]], filename=file.path(outpath, outfile), device = "pdf",
         width = 10,height = 12,units = "cm",dpi = 600, scale=2)
}

library(patchwork)
sites <- c("Fecal","Oral", "Skin")
for(site in sites){
  n <- grep(site, names(pop_heatmap))
  n <- names(pop_heatmap)[n]
  n <- c(n[grep("human",n)],n[grep("dog",n)])

  if(site=="Fecal"){
    layout <- pop_heatmap[[n[2]]] + pop_heatmap[[n[3]]] +
      plot_layout(widths = c(2, 1.2))
    print(layout)
    
    outfile=paste0(c(site,"Heatmap.pdf"),collapse = "_")
    ggsave(plot=layout, filename=file.path(outpath, outfile), device = "pdf",
           width = 15,height = 12,units = "cm",dpi = 600, scale=2.3)
    
    layout <- pop_heatmap[[n[1]]] + pop_heatmap[[n[2]]] +
      plot_layout(widths = c(1, 2))
    print(layout)
    
    outfile=paste0(c("human",site,"vs_Pilot_Heatmap.pdf"),collapse = "_")
    ggsave(plot=layout, filename=file.path(outpath, outfile), device = "pdf",
           width = 15,height = 12,units = "cm",dpi = 600, scale=2.3)
    
    } else {
      
      layout <- pop_heatmap[[n[1]]] + pop_heatmap[[n[2]]]
      print(layout)
      
      outfile=paste0(c(site,"Heatmap.pdf"),collapse = "_")
      ggsave(plot=layout, filename=file.path(outpath, outfile), device = "pdf",
             width = 15,height = 12,units = "cm",dpi = 600, scale=2)
  }
}

######
print("all completed")



