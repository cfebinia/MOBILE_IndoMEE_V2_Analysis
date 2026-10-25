basedir="Documents/Analysis/MOBILE/Wellcome_poster"
refdir="Documents/Analysis/MOBILE/INDOMEE_refdata"
outdir="Documents/Analysis/MOBILE/Wellcome_poster/output_files/MAG_sharing"

# load metadata
popdf <- read.delim("input_files/Consolidated_SampleID_combined_final_SeqIDmatched.tsv")
popdf <- popdf %>%
  filter(grepl("DHF|DDF|DHV|DDV|DHP|DDH", sampleid)) %>%
  mutate(SAMPLE = substr(sampleid,1,6),
         individual=gsub("Dx|Hx","",individual)) %>%
  distinct()

# obtain sampleid
sampleid <- unique(popdf$seqID)
type <- unique(substr(sampleid,1,3)) # type <- c("DHF","DHP","DHV","DDF","DDH","DDV")

# Load MAG table
strgMAG <- read.delim(file.path(refdir,"INDOMEE_V2_Abundance_Unified_draft3/bwa_counts_total_filtered_wMetadata_allUnified_draft3_strongFilter_relAbund.tsv"))
weakMAG <- read.delim(file.path(refdir,"INDOMEE_V2_Abundance_Unified_draft3/bwa_counts_total_filtered_wMetadata_allUnified_draft3_weakFilter_relAbund.tsv"))

table(paste(strgMAG$Genome,strgMAG$original_bin) %in% paste(weakMAG$Genome,weakMAG$original_bin))
table(paste(weakMAG$Genome,weakMAG$original_bin) %in% paste(strgMAG$Genome,strgMAG$original_bin))

strgMAG <- get_MAGs_matrix(strgMAG, sample.type = type)
weakMAG <- get_MAGs_matrix(weakMAG, sample.type = type)

to_exclude <- c("DHV999","DDHC01","DDHC01","DDHC02","DDHC03","DDHC04","DDHC05","DDHC06")
strgMAG <- strgMAG[,!colnames(strgMAG)%in%to_exclude]
weakMAG <- weakMAG[,!colnames(weakMAG)%in%to_exclude]

strgMAG <- remove_repeats(strgMAG)
weakMAG <- remove_repeats(weakMAG)

summary(factor(substr(colnames(strgMAG),1,3)))
summary(factor(substr(colnames(weakMAG),1,3)))

missing_inStrong <- rownames(strgMAG)[!rownames(strgMAG) %in% rownames(weakMAG)]
missing_inWeak <- rownames(weakMAG)[!rownames(weakMAG) %in% rownames(strgMAG)]

length(missing_inStrong)
length(missing_inWeak)

matched_set <- unique(c(rownames(strgMAG)[rownames(strgMAG) %in% rownames(weakMAG)],
                        rownames(strgMAG)[rownames(strgMAG) %in% rownames(weakMAG)]))
length(matched_set)
