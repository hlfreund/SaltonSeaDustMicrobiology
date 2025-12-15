
#### Set WD & Load Libraries ####
getwd() # use setwd("path/to/files") if you are not in the right directory
suppressPackageStartupMessages({ # load packages quietly
  library(phyloseq)
  library(ggplot2)
  library(vegan)
  library(ggpubr)
  #library(scales)
  library(grid)
  library(ape)
  library(plyr)
  library(dplyr)
  library(viridis)
  library(readxl)
  library(metagenomeSeq)
  library(DESeq2)
  library(dplyr)
  library(magrittr)
  library(MASS)
  library(dendextend)
  library(tidyr)
  library(reshape)
  library(reshape2)
  library(wesanderson)
  library(nationalparkcolors)
  library(shades)
  library(ALDEx2)
  library(rstatix)
  library(devtools)
  library(decontam)
})

#### Import and Prepare Data for Analyses ####

## Import ALL env plate bacterial ASV count data
bac.ASV_counts<-data.frame(readRDS("data/Amplicon/EnvMiSeq_W23_16S.V3V4_ASVs_Counts_dada2_Robject.rds", refhook = NULL))
dim(bac.ASV_counts)
bac.ASV_counts[,1:4]

# Fix Sample IDs so that they match the sample IDs in the metadata
#colnames(bac.ASV_counts)<-gsub("_S[0-9]+", "", colnames(bac.ASV_counts)) # shorten sample names to match sample names in metadata file
#head(bac.ASV_counts)
colnames(bac.ASV_counts)<-gsub("_",".", colnames(bac.ASV_counts))
bac.ASV_counts$ASV_ID<-rownames(bac.ASV_counts)
head(bac.ASV_counts)
colnames(bac.ASV_counts)
#bac.ASV_counts<-bac.ASV_counts[, !duplicated(colnames(bac.ASV_counts))] # remove col duplicates
dim(bac.ASV_counts)

# save total pre-contamination counts
TotalSampleCounts_withContams<-data.frame(SampleID=colnames(bac.ASV_counts[,!(colnames(bac.ASV_counts) %in% "ASV_ID")]),
                                          TotalCounts=colSums(bac.ASV_counts[,!(colnames(bac.ASV_counts) %in% "ASV_ID")]))
write.table(TotalSampleCounts_withContams, "data/Amplicon/EnvMiSeq_W23_16S.V3V4_ASVs_TotalCounts_withContams.tsv",
            sep="\t", quote=F, col.names=T,row.names=F)

# Remove unwanted samples aka old Salton Seawater samples sequenced by Zymo
#remove_samples<-c("SS.OV.10m.seawater.0621", "SS.OV.2m.seawater.0621", "SS.OV.5m.seawater.0621")
#bac.ASV_counts<-bac.ASV_counts[,!(colnames(bac.ASV_counts) %in% remove_samples)]
#colnames(bac.ASV_counts)
#dim(bac.ASV_counts)

## Import ASV taxonomic data
bac.ASV_tax<-data.frame(readRDS("data/Amplicon/EnvMiSeq_W23_16S.V3V4_ASVs_Taxonomy_dada2_Robject.rds", refhook = NULL))
head(bac.ASV_tax)

# Turn NAs into Unknowns
bac.ASV_tax[is.na(bac.ASV_tax)]<- "Unknown" # turn all NAs into "Unkowns"
bac.ASV_tax$Species<-gsub("Unknown", "unknown", bac.ASV_tax$Species) # change uppercase Unkonwn to lowercase unknown for unknown species classification
head(bac.ASV_tax)
bac.ASV_tax$ASV_ID<-rownames(bac.ASV_tax) # create ASV ID column to use for merging data frames
head(bac.ASV_tax)

#save.image("data/Env_Seqs_All/env.seq_analysis.Rdata")
#### Import metadata ####

# Import all metadata
metadata<-as.data.frame(read_excel("data/Metadata_EnvMiSeqPlate_Winter23.xlsx", sheet="Metadata"), header=TRUE)
# ^ for csv: as.data.frame(read.csv("file.csv",header = TRUE, sep = ",", quote = "",))

head(metadata)
metadata$SampleID<-gsub("_",".", metadata$SampleID)

#metadata$SampleID<-gsub("(.*\\..*)\\..*","\\1", metadata$SampleID)
#metadata<-na.omit(metadata) # drop NAs from metadata
rownames(metadata)<-metadata$SampleID
#metadata<-subset(metadata, Project=="SaltonSea")
head(metadata)
#metadata<-subset(metadata, select=-c(Project))
#head(metadata)
dim(metadata)

unique(metadata$SampleID)

#### Identify & Remove Contaminants ####
ControlDF<-metadata[metadata$SampleType=="Control",] # pull out samples that are controls

vector_for_decontam<-metadata$Sample_or_Control # use for decontam package
# ^ tells us which are controls aka TRUE vs which are not aka FALSE

bac.ASV_counts[,-length(bac.ASV_counts)] <- as.data.frame(sapply(bac.ASV_counts[,-length(bac.ASV_counts)], as.numeric)) #convert data frame to numeric
bac.ASV_c2<-t(bac.ASV_counts[,-length(bac.ASV_counts)]) # transpose so that rows are Samples and columns are ASVs
contam_df <- isContaminant(bac.ASV_c2, neg=vector_for_decontam) # use isContaminant to identify ASVs that are contaminants
## ^ ASVs that are more prominant in negative controls are identified as contaminants because there is less actual DNA to compete with, allowing these contaminants to be amplifeid in neg. controls

table(contam_df$contaminant) # identify contaminants aka TRUE: 2156

contam_asvs <- (contam_df[contam_df$contaminant == TRUE, ]) # pull out ASV IDs for contaminating ASVs
dim(contam_asvs)

bac.ASV_tax[row.names(bac.ASV_tax) %in% row.names(contam_asvs),] # see which taxa are contaminants

Contam.Taxa<-bac.ASV_tax[row.names(bac.ASV_tax) %in% row.names(contam_asvs),] # see which taxa are contaminants
Contam.Taxa[1:4,]

write.table(Contam.Taxa, "data/Amplicon/EnvMiSeq_W23_16S.V3V4_ContaminationASVs_List.tsv",
            sep="\t", quote=F, col.names=T,row.names=F)

## Create new files that EXCLUDE contaminants!!!

# making new fasta file
#contam_indices <- which(asv_fasta %in% paste0(">", contam_asvs))
#dont_want <- sort(c(contam_indices, contam_indices + 1))
#asv_fasta_no_contam <- asv_fasta[- dont_want]

# making new count table
bac.ASV_counts[row.names(bac.ASV_counts) %in% row.names(contam_asvs), ] # number of reads in contam ASVs
ncol(bac.ASV_counts)
#colSums(bac.ASV_counts[row.names(bac.ASV_counts) %in% row.names(contam_asvs),-126])

bac.ASV_counts_no.contam <- bac.ASV_counts[!row.names(bac.ASV_counts) %in% row.names(contam_asvs), ] # drop ASVs found in contam_asvs
head(bac.ASV_counts_no.contam)
dim(bac.ASV_counts_no.contam)

# making new taxonomy table
bac.ASV_tax.no.contam <- bac.ASV_tax[!row.names(bac.ASV_tax) %in% row.names(contam_asvs), ] # drop ASVs found in contam_asvs
head(bac.ASV_tax.no.contam)

# Remove ASVs found in Controls from samples (in addition to contaminants previously ID'd by decontam)

Control_counts<-bac.ASV_counts_no.contam[,colnames(bac.ASV_counts_no.contam) %in% ControlDF$SampleID] # see which taxa are contaminants
Control_counts
Control_counts<-Control_counts[which(rowSums(Control_counts) > 0),] # drop ASVs that don't appear in Controls
dim(Control_counts)
head(Control_counts)

bac.ASV_counts_CLEAN<-bac.ASV_counts_no.contam[!bac.ASV_counts_no.contam$ASV_ID %in% row.names(Control_counts),!colnames(bac.ASV_counts_no.contam) %in% colnames(Control_counts)]
bac.ASV_taxa_CLEAN<-bac.ASV_tax.no.contam[!bac.ASV_tax.no.contam$ASV_ID %in% row.names(Control_counts),]

# sanity check
colnames(bac.ASV_counts_CLEAN) # check for control sample IDs

## and now writing them out to files
#write(asv_fasta_no_contam, "ASVs-no-contam.fa")
write.table(bac.ASV_counts_CLEAN, "data/Amplicon/EnvMiSeq_W23_16S.V3V4_ASVs_Counts_NoContam.tsv",
            sep="\t", quote=F, col.names=NA)
saveRDS(bac.ASV_counts_CLEAN, file = "data/Amplicon/EnvMiSeq_W23_16S.V3V4_ASVs_Counts_NoContam_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

write.table(bac.ASV_taxa_CLEAN, "data/Amplicon/EnvMiSeq_W23_16S.V3V4_ASVs_Taxa_NoContam.tsv",
            sep="\t", quote=F, col.names=NA)
saveRDS(bac.ASV_taxa_CLEAN, file = "data/Amplicon/EnvMiSeq_W23_16S.V3V4_ASVs_Taxa_NoContam_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

#### Update Metadata ####
# create color variable(s) to identify variables by colors
## color for sample type
unique(metadata$SampleType)
metadata$SampleType<-factor(metadata$SampleType, levels=c("Soil","Dust","Lung","Control"))

#colorset1 = reshape2::melt(c(Seawater="#1f547b",Soil="#c44536",Dust="#432818",Playa="#d00000",Fecal="#66615f",Lung="#47126b",Control="#b13d1e"))
colorset1 = reshape2::melt(c(Soil="#c44536",Dust="#432818",Fecal="#66615f",Lung="#47126b",Control="#b13d1e"))

colorset1$SampleType<-rownames(colorset1)
colnames(colorset1)[which(names(colorset1) == "value")] <- "Sample_Color"
colorset1

metadata<-merge(metadata, colorset1, by="SampleType")
head(metadata)
metadata$Sample_Color <- as.character(metadata$Sample_Color)
rownames(metadata)<-metadata$SampleID
head(metadata)
dim(metadata)
#metadata$Depth_m<-factor(metadata$Depth_m, levels=c("0","2","3","4","5","7","8","9","10","11"))

# unique(metadata$SampleMonth)
# metadata$SampleMonth<-factor(metadata$SampleMonth, levels=c("Januarhy","February","March","April","May","June","August","September", "October","November","December","NA"))
#
# cold2warm1<-get_palette(paste0("#",c("252A52", "66ADE5", "FFC465","BF1B0B")),k=10)
# names(cold2warm1) <- levels(metadata$Depth_m)
#
# fair_cols <- paste0("#",c("252A52", "66ADE5", "FFC465","BF1B0B"))
# names(fair_cols) <- letters[1:4]
# fair_ramp <- scales::colour_ramp(fair_cols)
# fair_sat <- saturation(fair_ramp, 1)

### Transform Data ####
# first we merge the ASV count object and the ASV taxonomy object together by column called "ASV_ID"
## then we need to melt the separate ASV & taxonomy so that we can rbind multiple data sets

# Drop singletons & zero count ASVs
dim(bac.ASV_counts_CLEAN)
bac.ASV_counts_CLEAN<-bac.ASV_counts_CLEAN[which(rowSums(bac.ASV_counts_CLEAN[,-length(bac.ASV_counts_CLEAN)]) > 0),]
dim(bac.ASV_counts_CLEAN)

# merge CLEAN aka contaminants/controls removed count & taxa tables
bac.ASV_all<-merge(bac.ASV_counts_CLEAN,bac.ASV_taxa_CLEAN, by="ASV_ID")
head(bac.ASV_all)
dim(bac.ASV_all)
bac.ASV_all<-bac.ASV_all[, !duplicated(colnames(bac.ASV_all))] # remove col duplicates
dim(bac.ASV_all)
#bac.ASV_dat<-left_join(bac.ASV_tax.no.contam,bac.ASV_dat2, by=c("ASV_ID","Kingdom","Phylum","Class","Order","Family","Genus","Species"),all=T) # all=T to keep all rows from the dataframes, not drop rows that are missing from one DF or the other

bac.dat<-reshape2::melt(bac.ASV_all)
head(bac.dat)
colnames(bac.dat)[which(names(bac.dat) == "variable")] <- "SampleID"
colnames(bac.dat)[which(names(bac.dat) == "value")] <- "Count"

# Drop all Zero counts & singletons ASVs
dim(bac.dat)
bac.dat<-bac.dat[which(bac.dat$Count > 1),]
dim(bac.dat)

#bac.ASV_dat<-bac.ASV_dat[, !duplicated(colnames(bac.ASV_dat))] # remove col duplicates
#dim(bac.ASV_dat)

#bac.asv.melt<-reshape2::melt(bac.ASV_dat, id.vars=c("ASV_ID","Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species"))
#colnames(bac.asv.melt)[which(names(bac.asv.melt) == "variable")] <- "SampleID"
#colnames(bac.asv.melt)[which(names(bac.asv.melt) == "value")] <- "Count"
#head(bac.asv.melt)

# Drop unknowns and eukaryotic hits
bac.dat<-subset(bac.dat, Kingdom!="Unknown") ## drop Unknowns from Kingdom
bac.dat<-subset(bac.dat, Phylum!="Unknown") ## drop Unknowns from Phylum
head(bac.dat)
dim(bac.dat)

# Create ASV count file that is filtered of eukaryotic taxa - for later use
bac.dat.with.euks<-bac.dat # contains eukaryotic sequence hits
colnames(bac.dat.with.euks)

# Drop chloroplast & mitochondria seqs
bac.dat<-subset(bac.dat, Class!="Chloroplast") ## exclude Chloroplast sequences
bac.dat<-subset(bac.dat, Order!="Chloroplast") ## exclude Chloroplast sequences
bac.dat<-subset(bac.dat, Family!="Mitochondria") ## exclude Mitochondrial sequences just in case

'Chloroplast' %in% bac.dat # check if Chloroplast counts are still in df, should be false because they've been removed
#'Mitochondria' %in% bac.dat # check if Mitochondria counts are still in df, should be false because they've been removed
#'Undetermined' %in% bac.dat # check if undetermined taxa in data frame
#NA %in% bac.ASV_dat

head(bac.dat)

# Create filtered ASV table (no contaminants & no ASVs found in controls!)
bac.ASV_table<-base::as.data.frame(dcast(bac.dat, SampleID~ASV_ID, value.var="Count", fun.aggregate=sum)) ###
bac.ASV_table[1:5,1:5]
bac.ASV_table[duplicated(rownames(bac.ASV_table))]
rownames(bac.ASV_table)<-bac.ASV_table$SampleID
#bac.ASV_table<-subset(bac.ASV_table, select=-c(SampleID))
rownames(bac.ASV_table)

# double check dimensions of metadata and ASV table
dim(metadata)
dim(bac.ASV_table)
# double check that the rownames exist + match
rownames(metadata)
rownames(bac.ASV_table)

# Find rows in metadata that are not in combined bacterial asv tables
setdiff(rownames(metadata), rownames(bac.ASV_table)) # check rows in metadata not in bac.ASV_table
setdiff(rownames(bac.ASV_table), rownames(metadata)) # check rows in bac.ASV_table not in metadata

bac.ASV_meta<-merge(bac.ASV_table,metadata, by.x="SampleID", by.y="SampleID") # drop samples from Sierra project
dim(bac.ASV_meta)
dim(metadata) # the difference between dimensions is that the original metadata contains info for the controls
dim(bac.ASV_table)
rownames(bac.ASV_meta)<-bac.ASV_meta$SampleID

# recreate ASV taxonomy table so that ASVs that were removed are not included in new taxa table
bac.tax.clean<-as.data.frame(unique(subset(bac.dat, select=-c(SampleID,Count))))
head(bac.tax.clean)

# recreate ASV table (excluding samples we don't have metadata or counts for)
colnames(bac.ASV_meta); rownames(bac.ASV_meta)
head(bac.ASV_meta)
bac.ASV_table<-subset(bac.ASV_meta, select=-c(SampleType,TubeID,OriginalSampleID,SorC1,Sample_or_Control,Sample_Color))
rownames(bac.ASV_table)<-bac.ASV_table$SampleID
head(bac.ASV_table)

rowSums(bac.ASV_table[,-1])

# save total post-contamination counts
TotalSampleCounts_NoContams<-data.frame(SampleID=rownames(bac.ASV_table),
                                          TotalCounts=rowSums(bac.ASV_table[,-1]))
write.table(TotalSampleCounts_NoContams, "data/Amplicon/EnvMiSeq_W23_16S.V3V4_ASVs_TotalCounts_NOContams.tsv",
            sep="\t", quote=F, col.names=T,row.names=F)

# triple check dimensions of metadata and ASV table
dim(metadata)
dim(bac.ASV_table)

# Find rows in metadata that are not in combined bacterial asv tables
setdiff(rownames(metadata), rownames(bac.ASV_table)) # check rows from metadata not in bac.ASV_table
# the difference between dimensions is that the original metadata contains info for the controls
setdiff(rownames(bac.ASV_table), rownames(metadata)) # check bac.ASV_table from metadata not in rows

# reorder metadata based off of ASV table
metadata=metadata[rownames(bac.ASV_table),] ## will drop rows that are not shared by both dataframes!
# here we are reordering our metadata by rows, using the rownames from our ASV table as a guide
# this indexing method will only work if the two dfs have the same # of rows AND the same row names!

# sanity check to see if this indexing step worked
head(metadata)
dim(metadata)

head(rownames(bac.ASV_table))
dim(bac.ASV_table)

## subset data by sample type
duplicated(rownames(bac.ASV_meta))

rownames(bac.ASV_meta)
head(bac.ASV_meta)
#rownames(bac.meta.melt)<-bac.meta.melt$SampleID

# Remove ASVs from taxa table
head(bac.ASV_taxa_CLEAN)
bac.tax<-bac.ASV_taxa_CLEAN[,-length(bac.ASV_taxa_CLEAN)]
head(bac.tax)

# create super metadata + taxa + counts data frame
bac.dat.meta<-merge(bac.dat,metadata,by="SampleID")

#write.table(sw.asv.meta, "data/Env_Seqs_All/EnvSeqsAll_9.14.2022/SaltonSeawater_16S_AllData.tsv",
#            sep="\t", quote=F, col.names=NA)
#saveRDS(sw.asv.meta, file = "data/Env_Seqs_All/EnvSeqsAll_9.14.2022/SaltonSeawater_16S_AllData_Robject.rds", ascii = FALSE, version = NULL,
#        compress = TRUE, refhook = NULL)

#### Separate Data by Project ####
unique(bac.dat.meta$SampleType)

# first separate combined metadata + ASV table
d.asv.meta<-subset(bac.dat.meta, SampleType=="Dust")
s.asv.meta<-subset(bac.dat.meta, SampleType=="Soil")
l.asv.meta<-subset(bac.dat.meta, SampleType=="Lung")

# separate metadata by project
Dust_meta<-metadata[metadata$SampleType=="Dust",] # pull out samples that are dust
dim(Dust_meta)
Soil_meta<-metadata[metadata$SampleType=="Soil",] # pull out samples that are soil
Lung_meta<-metadata[metadata$SampleType=="Lung",] # pull out samples that are lung

# find DUST samples in ASV table that are in metadata (example as safety check)
dim(bac.ASV_table)
bac.ASV_table[1:5,1:5]

# double check this indexing works
bac.ASV_table[rownames(bac.ASV_table) %in% row.names(Dust_meta),]
dim(bac.ASV_table[rownames(bac.ASV_table) %in% row.names(Dust_meta),])
dim(Dust_meta)

# separate data by projects
dust_b.ASV<-bac.ASV_table[rownames(bac.ASV_table) %in% row.names(Dust_meta),] # find DUST samples in ASV table that are in metadata
rownames(dust_b.ASV) %in% row.names(Dust_meta) # sanity check that this worked

soil_b.ASV<-bac.ASV_table[rownames(bac.ASV_table) %in% row.names(Soil_meta),] # find SOIL samples in ASV table that are in metadata
rownames(soil_b.ASV) %in% row.names(Soil_meta) # sanity check that this worked

lung_b.ASV<-bac.ASV_table[rownames(bac.ASV_table) %in% row.names(Lung_meta),] # find LUNG samples in ASV table that are in metadata
rownames(lung_b.ASV) %in% row.names(Lung_meta) # sanity check that this worked

#### Prep Dataframe for Relative Abundance ####
head(t(bac.ASV_table[,-1]))

bac.clean.counts<-as.data.frame(t(bac.ASV_table[,-1])) # all singletons, zeros, controls, contaminants dropped
head(bac.clean.counts)
bac.clean.counts$ASV_ID<-rownames(bac.clean.counts)

#bac.tax<-subset(bac.ASV_tax.no.contam, select=c("ASV_ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
bac.all<-merge(bac.clean.counts, bac.ASV_taxa_CLEAN, by="ASV_ID")
head(bac.all)
bac_melt<-reshape2::melt(bac.all, id.vars = c("ASV_ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
head(bac_melt)
names(bac_melt)[which(names(bac_melt) == "variable")] <- "SampleID"
names(bac_melt)[which(names(bac_melt) == "value")] <- "Counts"
head(bac_melt)

head(metadata)

all_bac<-merge(bac_melt, metadata, by = "SampleID")
head(all_bac) # contains metadata, ASV counts, and taxonomic IDs for ASVs

#### Save Files for Each Project ####

# save ASV table for dust; includes no contaminants, zeros, singletons, or eukaryotic hits
saveRDS(dust_b.ASV, file = "data/Amplicon/SaltonSeaDust_16S.V3V4_ASVTable_Robject.rds", ascii = FALSE, version = NULL,
       compress = TRUE, refhook = NULL)
#saveRDS(soil_b.ASV, file = "data/Amplicon/SaltonSeaDust_16S.V3V4_ASVTable_Robject.rds", ascii = FALSE, version = NULL,
#        compress = TRUE, refhook = NULL)
saveRDS(lung_b.ASV, file = "data/Amplicon/SaltonSea_Lungs_16S.V3V4_ASVTable_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

# save updated taxa file
saveRDS(bac.tax.clean, file = "data/Amplicon/EnvMiSeq_W23_16S.V3V4_ASVs_Taxonomy_dada2_Clean_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

#### Pull out ASV Sequences of Interest for MSA & Phylogenetic Tree ####
## pulling sequences from DADA2 results only for ASVs that we have in our cleaned up bac.ASV_table
## use this new set of ASVs of interest to construct phylogenetic tree, import the tree with ape package, and calculate Unifrac distance with rbiom
# giving our seq headers more manageable names (ASV_1, ASV_2...)
load(file = "data/Amplicon/mydada_16S.V3V4.Rdata") # load your DADA2 output for this project

asv_seqs <- colnames(seqtab.nochim) # pull out ASV sequences
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character") # pull out ASV headers

# add ASV in front of each # of each sequence
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste("ASV", i, sep="_")
}

# create table of ASVs and their respective sequences so that we can filter out ASVs from samples we aren't looking at
asv_all <- as.data.frame(cbind(asv_headers, asv_seqs))
length(unique(asv_all$asv_headers))

# only pull out ASVs and their sequences for ASVs in our taxa table for only samples of interest
new.asv_all<-asv_all[asv_all$asv_headers %in% colnames(bac.ASV_table[,-1]),]

# insert ">" before ASV IDs before we make the new FASTA file
new.asv_all$asv_headers<-gsub("ASV_",">ASV_",new.asv_all$asv_headers)
new.asv_all$asv_headers[1:5]

# create and save the new fasta file - will use this to generate Multiple Sequence Alignment & build phylogenetic tree
new.asv_fasta<-c(rbind(new.asv_all$asv_headers,new.asv_all$asv_seqs))
write(new.asv_fasta, "SSD_16SV3V4_ASVs_Updated.fa")

#### Save Global Env for Import into Other Scripts ####

save.image("data/Amplicon/EnvMiSeq_W23_16S.V3V4_DataReady.Rdata") # save global env to Rdata file

## ^ this will be loaded into other scripts for downstream analyses
