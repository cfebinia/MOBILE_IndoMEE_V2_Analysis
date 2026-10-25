rm(list=ls())

library(dplyr)

refdir="Documents/Analysis/MOBILE/INDOMEE_refdata"

setwd(refdir)

sampleid <- read.delim("Consolidated_SampleID_combined_final_SeqIDmatched.tsv")

sampleid %>%
  group_by(host, type) %>%
  summarise(n_indivdual = n_distinct(individual),
            n_biological= n_distinct(labid),
            n_sequenced = n_distinct(substr(seqID,1,6)))

n_distinct(sampleid$seqID)

sampleid %>%
  group_by(host) %>%
  summarise(n_indivdual = n_distinct(individual),
            n_biological= n_distinct(labid),
            n_sequenced = n_distinct(substr(seqID,1,6)))

lib_qc <- read.delim("Sequencing_Library_QC_Stats_All_samples.tsv",row.names = 1)
lib_qc <- lib_qc %>%
  rename(labid=Sample.Name,
         concentration=Conc...ng.ul.,
         yield_mcg=Total.Amount..ug.,
         pass=Pass.Fail,
         warning=Result.,
         prep_kit=Library.Prep) %>%
  select(-X) %>%
  mutate(sampleid=substr(labid,1,6)) %>%
  mutate(type=substr(sampleid,1,3)) %>%
  mutate(site=case_when(
    type %in% c("DHF","DDF") ~ "Fecal",
    type %in% c("DHV","DDV") ~ "Oral",
    type %in% c("DHP","DDH") ~ "Skin",
    TRUE ~ "Others"
  )) %>%
  mutate(host=case_when(
    type %in% c("DHF","DHV","DHP") ~ "human",
    type %in% c("DDF","DDV","DDH") ~ "dog",
    TRUE ~ "Other"
  ))  

library(stringr)

repeat_seq <- lib_qc %>%
  mutate(across(c(type, prep_kit), str_trim),
         prep_kit = tolower(prep_kit)) %>%
  group_by(type, prep_kit, sampleid) %>%
  summarise(n = n_distinct(labid), .groups = "drop") 
repeat_samples <- unique(repeat_seq$sampleid)

repeat_seq <- repeat_seq %>%
  filter(n > 1) %>%
  group_by(type, prep_kit) %>%
  summarise(n_repeat = n(), .groups = "drop")


controls_cat <- data.frame(type=c("DCT","NTC","ZBD",
                                  "ZBS","ZBV","ZBL"),
                           s_type=c("Environmental NTC","Technical NTC","DNA Control",
                                    "Mock Control - Faecal", "Mock Control - Saliva", "Mock Control - Skin"))
sample_cat <- data.frame(type=c("DHF","DHV","DHP",
                                "DDF","DDV","DDH"),
                         s_type=c("humanFecal","humanOral","humanSkin",
                                  "dogFecal","dogOral","dogSkin"))
  
lib_qc %>%
  mutate(across(c(type, prep_kit), str_trim)) %>%
  mutate(prep_kit=tolower(prep_kit))%>%
  group_by(type, host, site, prep_kit) %>%
  summarise(n = n_distinct(sampleid),
            n_labid = n_distinct(labid),
            min_yield = min(yield_mcg),
            max_yield = max(yield_mcg),
            mean_yield = mean(yield_mcg),
            sem_yield = sd(yield_mcg) / sqrt(n),
            .groups = "drop") %>%
  filter(type %in% c(sample_cat$type,controls_cat$type)) %>%
  merge.data.frame(rbind(sample_cat, controls_cat), by="type") %>%
  mutate(type=factor(type, c(sample_cat$type,controls_cat$type))) %>%
  arrange(type) -> prep_table

prep_table %>% filter(n_labid>1) %>% View()

lib_qc %>% 
  mutate(across(c(type, prep_kit), str_trim)) %>%
  mutate(prep_kit=tolower(prep_kit)) %>%
  group_by(sampleid,type,host,site, prep_kit) %>%
  summarise(n=n_distinct(sampleid),.groups = "drop") %>%
  filter(n>1) %>%
  filter(type %in% c("DHF","DHV","DHP","DDF","DDV","DDH")) -> repeat_seq

ntc <- lib_qc %>%
  filter(!type %in% c("DHF","DHV","DHP","DDF","DDV","DDH","DEN","ZBD"))

ntc %>%  summarise(n = n_distinct(labid),
                   total_yield=mean(yield_mcg), min=min(yield_mcg), max=max(yield_mcg),
                   min_conc=min(concentration), max_conc=max(concentration),
                   min_DIN=min(DIN), max_DIN=max(DIN)) 

# host reads
nReads = "MOBILE_Read_Counts/allUnified_draft3_weakFilter_mapStats.tsv"
nReads <- read.delim(nReads) %>%
  rename(seqID=sample)

ids <- "Consolidated_SampleID_combined_final_SeqIDmatched.tsv"
ids <- read.delim(ids)

df <- ids %>%
  select(seqID, individual, sampleid, population, site, host, type) %>%
  right_join(nReads[,c("seqID","bwa_counts_unique_pass_pct","bwa_counts_total_pass_pct",
                       "fastq_raw_reads","fastq_cleaned_reads",
                       "host_pct","unmapped_pct","duplication_pct","est_unique_reads")],
             by="seqID", keep = FALSE) %>%
  filter(!is.na(type))

intersect(lib_qc$sampleid, df$sampleid) -> matched_sid

lib_qc[!lib_qc$sampleid %in% matched_sid,]
df[!df$sampleid %in% matched_sid,]


df_controls <- nReads[,c("seqID","bwa_counts_unique_pass_pct","bwa_counts_total_pass_pct",
                         "fastq_raw_reads","fastq_cleaned_reads",
                         "host_pct","unmapped_pct","duplication_pct","est_unique_reads")] %>%
  mutate(type=substr(seqID,1,3)) %>%
  filter(type %in% controls_cat$type) %>%
  merge.data.frame(controls_cat, by="type") %>%
  mutate(individual=substr(seqID,1,6)) %>%
  select(-type) %>%
  rename(type=s_type) %>%
  filter(seqID!="NTC052")

host_reads <- df %>%
  filter(seqID!="NTC052") %>%
  select(all_of(colnames(df_controls))) %>%
  rbind(df_controls) %>%
  mutate(type=factor(type, c(sample_cat$s_type, controls_cat$s_type))) %>%
  mutate(raw_Gb=fastq_raw_reads*150/1e9,
         clean_Gb=fastq_cleaned_reads*150/1e9) %>%
  mutate(fastq_raw_reads=fastq_raw_reads/1e6,
         fastq_cleaned_reads=fastq_cleaned_reads/1e6) %>%
  group_by(type) %>%
  summarise(n=n_distinct(individual),
            n_samples=n_distinct(seqID),
            mean_raw_Gb=mean(raw_Gb),
            mean_clean_Gb=mean(clean_Gb),
            mean_raw_mio=mean(fastq_raw_reads),
            mean_clean_mio=mean(fastq_cleaned_reads),
            sem_clean=sd(fastq_cleaned_reads)/sqrt(length(fastq_cleaned_reads)),
            med_clean=median(fastq_cleaned_reads),
            min_clean=min(fastq_cleaned_reads),
            max_clean=max(fastq_cleaned_reads),
            mean_host=mean(host_pct),
            sem_host=sd(host_pct)/sqrt(length(host_pct)),
            med_host=median(host_pct),
            min_host=min(host_pct),
            max_host=max(host_pct),
            .groups = "drop") %>%
  mutate(across(c(mean_raw_Gb, mean_clean_Gb, mean_raw_mio, mean_clean_mio, 
                  sem_clean, med_clean,min_clean, max_clean,
                  mean_host, sem_host, med_host, min_host, max_host), ~ round(., digits = 1))) %>%
  arrange(type)

df %>%
  filter(seqID!="NTC052") %>%
  select(all_of(colnames(df_controls))) %>%
  rbind(df_controls) %>%
  mutate(type=factor(type, c(sample_cat$s_type, controls_cat$s_type))) %>%
  filter(type=="Technical NTC") %>%
  mutate(raw_Gb=round(fastq_raw_reads*150/1e9, 3),
         clean_Gb=round(fastq_cleaned_reads*150/1e9, 3)) %>%
  mutate(fastq_raw_reads = round(fastq_cleaned_reads /1e6,3),
         fastq_cleaned_reads = round(fastq_cleaned_reads /1e6, 3)) %>%
  mutate(outlier=ifelse(seqID=="NTC052", "Outlier","Normal")) %>%
  group_by(outlier) %>%
  summarise(n=n_distinct(seqID), 
            mean_raw_Gb=mean(raw_Gb),
            mean_clean_Gb=mean(clean_Gb),
            mean_raw_reads=mean(fastq_raw_reads),
            mean_clean_reads=mean(fastq_cleaned_reads))


# mapping rates
weakF <- read.delim("MOBILE_Read_Counts/allUnified_draft3_weakFilter_mapStats.tsv") %>% rename(seqID=sample)
strgF <- read.delim("MOBILE_Read_Counts/allUnified_draft3_strongFilter_mapStats.tsv") %>% rename(seqID=sample)

mapping_df <- rbind(data.frame(method="weak filter", weakF),
                    data.frame(method="strong filter", strgF)) %>%
  mutate(sampleid=substr(seqID,1,6))

mapping_df <- ids %>%
  filter(seqID!="NTC052") %>%
  select(sampleid, individual, population, site, host, type) %>%
  distinct() %>%
  right_join(mapping_df, by="sampleid", keep = FALSE) %>%
  mutate(sampleid=substr(seqID, 1,6)) %>%
  filter(!is.na(type)) 

#dropped <- "Outliers_list.txt"
#dropped <- read.delim(dropped) 
#table(mapping_df$seqID %in% dropped$Sample)
#mapping_df$seqID [mapping_df$seqID %in% dropped$Sample]
#mapping_df <- mapping_df %>%
#  filter(!seqID %in% dropped$Sample)

mapping_rates <- mapping_df %>%
  filter(seqID!="NTC052") %>%
  group_by(method, type) %>%
  mutate(fastq_cleaned_reads=fastq_cleaned_reads/1e6,
         est_unique_reads=est_unique_reads/1e6) %>%
  summarise(n=n_distinct(sampleid),
            mean_clean=mean(fastq_cleaned_reads),
            mean_unique=mean(100*est_unique_reads/fastq_cleaned_reads),
            mean_bwaTotal=mean(bwa_counts_total_pass_pct),
            med_bwaTotal=median(bwa_counts_total_pass_pct),
            min_bwaTotal=min(bwa_counts_total_pass_pct),
            max_bwaTotal=max(bwa_counts_total_pass_pct),
            mean_mapped=mean(100-unmapped_pct),
            med_mapped=median(100-unmapped_pct),
            min_mapped=min(100-unmapped_pct),
            max_mapped=max(100-unmapped_pct),
            .groups = "drop") %>%
  mutate(across(c(mean_clean, mean_unique, 
                  mean_bwaTotal, med_bwaTotal, min_bwaTotal,max_bwaTotal, 
                  mean_mapped, med_mapped, min_mapped, max_mapped), ~ round(., digits = 1)))


# NTC052
ntc052_df <- as.matrix(read.delim("MOBILE_Summary_Stats/NTC052.txt",row.names = 1))
ntc052_df[is.na(ntc052_df)] <-  0
ntc052_df <- ntc052_df[,colSums(ntc052_df)>0]
n <- which(colnames(ntc052_df) == "NTC052")
ntc052_df <- cbind(NTC052=ntc052_df[,n],ntc052_df[,-n])
ntc052_corr <- cor(ntc052_df[,1],ntc052_df[,-1],method = "pearson")
ntc052_corr <- setNames(ntc052_corr[1,], colnames(ntc052_corr))
ggdata <- ids %>%
  select(seqID, type) %>%
  filter(seqID %in% names(ntc052_corr))

all(names(ntc052_corr) %in% ggdata$seqID)
setdi <- names(ntc052_corr)[! names(ntc052_corr)%in% ggdata$seqID]

tmp <- data.frame(seqID=setdiff(names(ntc052_corr), ggdata$seqID),
                  type=substr(setdiff(names(ntc052_corr), ggdata$seqID),1,3)) %>%
  mutate(type=case_when(
    type %in% c("ASM","BA0","BAP","MAL","BRU") ~ "humanFecal",
    type %in% c("DOT") ~ "Researcher",
    type %in% c("DEN") ~ "Environment",
    type %in% c("DDH") ~ "dogSkin",
    TRUE ~ "NTC"
  ))

ggdata <- rbind(ggdata, tmp) 
ggdata$Rho <- ntc052_corr[ggdata$seqID]
ggdata <- ggdata %>%
  arrange(desc(Rho)) %>%
  mutate(type=factor(type, unique(type)))

library(ggplot2)
ntc052_plot <- ggplot(ggdata, aes(x = type, y = Rho, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.7, show.legend = FALSE) +
  labs(
    title = "Pearson Correlation with NTC052",
    caption = "X-axis arranged in order of the highest Rho per category.",
    x = "Sample Type",
    y = "Rho Value"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ntc052_plot

# write out
outdir <- file.path(refdir,"MOBILE_Summary_Stats")
dir.create(outdir,showWarnings = F,recursive = T)

write.table(prep_table, quote = F, sep = "\t", row.names = F,
            file = file.path(outdir,"lib_prep_n_count.tsv"))
write.table(host_reads, quote = F, sep = "\t", row.names = F,
            file = file.path(outdir,"host_reads_pct.tsv"))
write.table(mapping_rates, quote = F, sep = "\t", row.names = F,
            file = file.path(outdir,"mapping_rates.tsv"))


write.table(ggdata, quote = F, sep = "\t", row.names = F,
            file = file.path(outdir,"ntc052_pearsonCor.tsv"))
ggsave(plot=ntc052_plot, filename = file.path(outdir,"ntc052_pearsonCor.pdf"),
       device = "pdf",width = 15,height = 10,units = "cm",dpi = 150, scale=1.5)
