rm(list=ls())
setwd("Documents/Analysis/MOBILE/MOBILE_microbiome")
input_folder <- "input_files/humanFecal_draft1/summary"
abundance <- paste0(input_folder,"/","bwa_counts_total_filtered_wMetadata_humanFecal_strongFilter_relAbund.tsv")
metadata.raw <- "input_files/mobile_metadata.tsv"

sampleid <- read.table("input_files/sid.mobile.txt", header = T)[,1]

raw <- read.delim(abundance)
df <- data.frame(raw[,sampleid], row.names = raw$Genome)
colnames(df) <- gsub(x = colnames(df),pattern = "_POOLED",replacement = "")
df <- data.frame(t(df))

df.env <- read.delim(metadata.raw)
rownames(df.env) <- gsub(x = df.env$fastq_file,pattern = "_POOLED",replacement = "")
df.env$sampleid <- rownames(df.env)

# PCoA
library(vegan)

# add a small number to allow clr conversion
norm <- 10^-5
mydist <- vegdist(df+norm,method = "aitchison")

# run PCoA
cap <- capscale(mydist ~ 1,data = df.env)

# set plot variables
col.ord <- c("Basap","Lundayeh","PunanTubu","PunanAput","PunanBatu","OrangRimbaTerap","OrangRimbaSuban")
mycols <- c("#FFB90F","#2E8B57","#DF96B4","#B494D5","#5B3794","#FF6347","#A20000")
names(mycols) <- col.ord
print(mycols)

shape.ord <- unique(df.env$island)
myshape <- c(21,22)
names(myshape) <- shape.ord
print(myshape)


# plot
source("99.pretty_pcoa_plot_function.R")
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(svglite)
variance_explained <- summary(cap)$cont$importance
print(variance_explained)

myplot <- pretty_PCA(cap, 
           axes=c(1,2),
           metadata=df.env,
           colour.by="population",
           colour.scheme=mycols,
           colour.order=names(mycols),
           shape.by="island",
           shape.scheme=myshape,
           shape.order=names(myshape),
           output.path="figout/prelim_")

# taxonomy
source("99.parse_taxnomy.R")
tax.df <- parse_gtdbtk(raw,mag.colname="Genome",gtdbtk.colname="classification")                     
count_uniq_tax(tax.df)

count_unique_bylevel(tax.df,level="domain")
count_unique_bylevel(tax.df,level="phylum",order = T)
count_unique_bylevel(tax.df,level="class",order = T)
count_unique_bylevel(tax.df,level="order",order = T)
count_unique_bylevel(tax.df,level="family",order = T)
count_unique_bylevel(tax.df,level="genus", breakdown=F)
count_unique_bylevel(tax.df,level="species",breakdown=F)


################
# FAMILY ABUND #
################
library(reshape2)
library(dplyr)

taxfit.df <- cbind(tax.df, raw[,sampleid])
colnames(taxfit.df) <- gsub(pattern = "_POOLED",replacement = "", colnames(taxfit.df))
sid <- colnames(taxfit.df)[grep(pattern = "DHF", colnames(taxfit.df))]
taxfit.df %>% select("family",sid) %>% group_by(family) %>% summarise(across(everything(), sum)) -> tmp
tmp <- tmp[tmp$family!="",]
preval <- tmp[,-1]
preval[preval>0] <- 1
rownames(preval) <- tmp$family
select.genily <- which(rowSums(preval)>10)

family.abund <- t(tmp[select.genily,-1])
colnames(family.abund) <- tmp$family[select.genily]

source("99.clr_transform.R")
family.abund.clr <- clr_transform(family.abund, pseudocount = norm)

# fit.fam <- envfit(cap, family.abund.clr, permutations=9999)
# saveRDS(fit.fam,file = "rds/family.vectorfit.rds")
fit.fam <- readRDS("rds/family.vectorfit.rds")
signif.adjust <- data.frame(pval=fit.fam$vectors$pvals,
                            bonferroni=p.adjust(fit.fam$vectors$pvals,method = "bonferroni"),
                            fdr=p.adjust(fit.fam$vectors$pvals, method = "fdr"))
pmax.adj <- max(signif.adjust$pval[which(signif.adjust$bonferroni<0.05)])

plot(cap, type="none", main="Family") ; plot(fit.fam, p.max = pmax.adj, cex=0.8, col="rosybrown")

cbind(r=fit.fam$vectors$r, pval=fit.fam$vectors$pvals)[signif.adjust$bonferroni<0.05,]

signif.gen <- rownames(signif.adjust)[signif.adjust$bonferroni<0.05]

ggdata <- data.frame(sampleid=rownames(family.abund), family.abund[,signif.gen]+norm)
ggdata <- merge.data.frame(df.env[,c("sampleid","population")], ggdata, by="sampleid")
ggdata <- melt(ggdata, id.vars = c("sampleid","population"),variable.name = "family",value.name = "abund")

p <- ggplot(ggdata, aes(x=population, y=abund, colour=population))+
  geom_boxplot() + theme_pubr() +
  scale_colour_manual(values = mycols)+
  scale_y_continuous(transform = "log10")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.5, linetype=2, colour="grey"))+
  facet_wrap(~family, scale="free_x", ncol=4) +
  labs(caption=paste("Note: ", norm, "was added to allow log10 transformation"))

ggsave(plot=p, device = "pdf", filename="figout/family_boxplot_log10.pdf",width = 10,height = 40,dpi = 150,units = "in")

p <- ggplot(ggdata, aes(x=population, y=abund-norm, colour=population))+
  geom_boxplot() + theme_pubr() +
  scale_colour_manual(values = mycols)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.5, linetype=2, colour="grey"))+
  facet_wrap(~family, scale="free", ncol=4)
ggsave(plot=p, device = "pdf", filename="figout/family_boxplot.pdf",width = 15,height = 40,dpi = 200,units = "in")


###############
# GENUS ABUND #
###############
taxfit.df <- cbind(tax.df, raw[,sampleid])
colnames(taxfit.df) <- gsub(pattern = "_POOLED",replacement = "", colnames(taxfit.df))
sid <- colnames(taxfit.df)[grep(pattern = "DHF", colnames(taxfit.df))]
taxfit.df %>% select("genus",sid) %>% group_by(genus) %>% summarise(across(everything(), sum)) -> tmp
tmp <- tmp[tmp$genus!="",]
preval <- tmp[,-1]
preval[preval>0] <- 1
rownames(preval) <- tmp$genus
select.genus <- which(rowSums(preval)>10)

genus.abund <- t(tmp[select.genus,-1])
colnames(genus.abund) <- tmp$genus[select.genus]

genus.abund.clr <- clr_transform(genus.abund, pseudocount = norm)

# fit.gen <- envfit(cap, genus.abund.clr, permutations=9999)
# saveRDS(fit.gen,file = "rds/genus.vectorfit.rds")
fit.gen <- readRDS("rds/genus.vectorfit.rds")
signif.adjust <- data.frame(pval=fit.gen$vectors$pvals,
                            bonferroni=p.adjust(fit.gen$vectors$pvals,method = "bonferroni"),
                            fdr=p.adjust(fit.gen$vectors$pvals, method = "fdr"))
pmax.adj <- max(signif.adjust$pval[which(signif.adjust$bonferroni<0.05)])

plot(cap, type="none", main="Genus") ; plot(fit.gen, p.max = pmax.adj, cex=0.8, col="#836FFF")

cbind(r=fit.gen$vectors$r, pval=fit.gen$vectors$pvals)[signif.adjust$bonferroni<0.05,]

signif.gen <- rownames(signif.adjust)[signif.adjust$bonferroni<0.05]

ggdata <- data.frame(sampleid=rownames(genus.abund), genus.abund[,signif.gen]+norm)
ggdata <- merge.data.frame(df.env[,c("sampleid","population")], ggdata, by="sampleid")
ggdata <- melt(ggdata, id.vars = c("sampleid","population"),variable.name = "genus",value.name = "abund")

p <- ggplot(ggdata, aes(x=population, y=abund, colour=population))+
  geom_boxplot() + theme_pubr() +
  scale_colour_manual(values = mycols)+
  scale_y_continuous(transform = "log10")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.5, linetype=2, colour="grey"))+
  facet_wrap(~genus, scale="free_x", ncol=4) +
  labs(caption=paste("Note: ", norm, "was added to allow log10 transformation"))

ggsave(plot=p, device = "pdf", filename="figout/genus_boxplot_log10.pdf",width = 10,height = 40,dpi = 150,units = "in")

p <- ggplot(ggdata, aes(x=population, y=abund-norm, colour=population))+
  geom_boxplot() + theme_pubr() +
  scale_colour_manual(values = mycols)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.5, linetype=2, colour="grey"))+
  facet_wrap(~genus, scale="free", ncol=4)
ggsave(plot=p, device = "pdf", filename="figout/genus_boxplot.pdf",width = 15,height = 40,dpi = 200,units = "in")
