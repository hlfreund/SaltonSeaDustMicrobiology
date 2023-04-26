
#### Set WD & Load Libraries ####
getwd() # use setwd("path/to/files") if you are not in the right directory
setwd("/Volumes/HLF_SSD/Aronson_Lab_Data/Salton_Sea/SaltonSeaDust")
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

#load("data/EnvMiSeq_W23_Data.Rdata") # load Rdata to global env
#save.image("data/Env_Seqs_All/env.seq_analysis.Rdata") # save global env to Rdata file
load("data/MiSeq_16S.V3V4_W23_Data_Ready.Rdata")

#### Import and Prepare Data for Analyses ####

## Import ALL env plate bacterial ASV count data
bac.ASV_counts<-data.frame(readRDS("data/EnvMiSeq_W23_16S.V3V4_ASVs_Counts_dada2_Robject.rds", refhook = NULL))
dim(bac.ASV_counts)
head(bac.ASV_counts)

#colnames(bac.ASV_counts)<-gsub("_S[0-9]+", "", colnames(bac.ASV_counts)) # shorten sample names to match sample names in metadata file
#head(bac.ASV_counts)
colnames(bac.ASV_counts)<-gsub("_",".", colnames(bac.ASV_counts))
bac.ASV_counts$ASV_ID<-rownames(bac.ASV_counts)
head(bac.ASV_counts)
colnames(bac.ASV_counts)
#bac.ASV_counts<-bac.ASV_counts[, !duplicated(colnames(bac.ASV_counts))] # remove col duplicates
dim(bac.ASV_counts)

# Remove unwanted samples aka old Salton Seawater samples sequenced by Zymo
#remove_samples<-c("SS.OV.10m.seawater.0621", "SS.OV.2m.seawater.0621", "SS.OV.5m.seawater.0621")
#bac.ASV_counts<-bac.ASV_counts[,!(colnames(bac.ASV_counts) %in% remove_samples)]
#colnames(bac.ASV_counts)
#dim(bac.ASV_counts)

## Import ASV taxonomic data
bac.ASV_tax<-data.frame(readRDS("data/EnvMiSeq_W23_16S.V3V4_ASVs_Taxonomy_dada2_Robject.rds", refhook = NULL))
head(bac.ASV_tax)

bac.ASV_tax[is.na(bac.ASV_tax)]<- "Unknown" # turn all NAs into "Unkowns"
bac.ASV_tax$Species<-gsub("Unknown", "unknown", bac.ASV_tax$Species) # change uppercase Unkonwn to lowercase unknown for unknown species classification
head(bac.ASV_tax)
bac.ASV_tax$ASV_ID<-rownames(bac.ASV_tax) # create ASV ID column to use for merging data frames
head(bac.ASV_tax)

#### Import metadata ####

# Import all metadata
metadata<-as.data.frame(read_excel("data/Metadata_EnvMiSeqPlate_Winter23.xlsx", sheet="Metadata"), header=TRUE)
head(metadata)
metadata$SampleID<-gsub("_",".", metadata$SampleID)

#metadata$SampleID<-gsub("(.*\\..*)\\..*","\\1", metadata$SampleID)
#metadata<-na.omit(metadata) # drop NAs from metadata
rownames(metadata)<-metadata$SampleID
#metadata<-subset(metadata, Project=="SaltonSea")
head(metadata)
#metadata<-subset(metadata, select=-c(Project))
#head(metadata)


#### Identify & Remove Contaminants ####
ControlDF<-metadata[metadata$SampleType=="Control",] # pull out samples that are controls

vector_for_decontam<-metadata$Sample_or_Control # use for decontam package
# ^ tells us which are controls aka TRUE vs which are not aka FALSE

bac.ASV_counts[,-length(bac.ASV_counts)] <- as.data.frame(sapply(bac.ASV_counts[,-length(bac.ASV_counts)], as.numeric)) #convert data frame to numeric
bac.ASV_c2<-t(bac.ASV_counts[,-length(bac.ASV_counts)]) # transpose so that rows are Samples and columns are ASVs
contam_df <- isContaminant(bac.ASV_c2, neg=vector_for_decontam)

table(contam_df$contaminant) # identify contaminants aka TRUE: 2156

contam_asvs <- (contam_df[contam_df$contaminant == TRUE, ]) # pull out ASV IDs for contaminating ASVs

bac.ASV_tax[rownames(bac.ASV_tax) %in% rownames(contam_asvs),] # see which taxa are contaminants

## Create new files that EXCLUDE contaminants!!!

# making new fasta file
#contam_indices <- which(asv_fasta %in% paste0(">", contam_asvs))
#dont_want <- sort(c(contam_indices, contam_indices + 1))
#asv_fasta_no_contam <- asv_fasta[- dont_want]

# making new count table
bac.ASV_counts_no.contam <- bac.ASV_counts[!rownames(bac.ASV_counts) %in% rownames(contam_asvs), ] # drop ASVs found in contam_asvs
head(bac.ASV_counts_no.contam)

# making new taxonomy table
bac.ASV_tax.no.contam <- bac.ASV_tax[!rownames(bac.ASV_tax) %in% rownames(contam_asvs), ] # drop ASVs found in contam_asvs
head(bac.ASV_tax.no.contam)

# Remove ASVs found in Controls from samples (in addition to contaminants previously ID'd)

Control_counts<-bac.ASV_counts_no.contam[,colnames(bac.ASV_counts_no.contam) %in% ControlDF$SampleID] # see which taxa are contaminants
Control_counts
Control_counts<-Control_counts[which(rowSums(Control_counts) > 0),] # drop ASVs that don't appear in Controls
dim(Control_counts)
head(Control_counts)

# Create new DFs after removing contaminants and control samples
bac.ASV_counts_CLEAN<-bac.ASV_counts_no.contam[!bac.ASV_counts_no.contam$ASV_ID %in% rownames(Control_counts),!colnames(bac.ASV_counts_no.contam) %in% colnames(Control_counts)]
bac.ASV_taxa_CLEAN<-bac.ASV_tax.no.contam[!bac.ASV_tax.no.contam$ASV_ID %in% rownames(Control_counts),]
metadata_clean<-metadata[!metadata$SampleID %in% colnames(Control_counts),]

# sanity check
colnames(bac.ASV_counts_CLEAN) # check for control sample IDs

## and now writing them out to files
#write(asv_fasta_no_contam, "ASVs-no-contam.fa")
write.table(bac.ASV_counts_CLEAN, "data/EnvMiSeq_W23_16S.V3V4_ASVs_Counts_NoContam.tsv",
            sep="\t", quote=F, col.names=NA)
saveRDS(bac.ASV_counts_CLEAN, file = "data/EnvMiSeq_W23_16S.V3V4_ASVs_Counts_NoContam_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

write.table(bac.ASV_taxa_CLEAN, "data/EnvMiSeq_W23_16S.V3V4_ASVs_Taxa_NoContam.tsv",
            sep="\t", quote=F, col.names=NA)
saveRDS(bac.ASV_taxa_CLEAN, file = "data/EnvMiSeq_W23_16S.V3V4_ASVs_Taxa_NoContam_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

save.image("data/EnvMiSeq_W23_Data.Rdata")

#### Update Metadata ####
# create color variable(s) to identify variables by colors
## color for sample type
unique(metadata$SampleType)
metadata$SampleType<-factor(metadata$SampleType, levels=c("Soil","Dust","Lung","Control"))

metadata_clean$SampleType<-factor(metadata_clean$SampleType, levels=c("Soil","Dust","Lung"))

#colorset1 = melt(c(Seawater="#1f547b",Soil="#c44536",Dust="#432818",Playa="#d00000",Fecal="#66615f",Lung="#47126b",Control="#b13d1e"))
colorset1 = melt(c(Soil="#c44536",Dust="#432818",Lung="#47126b",Control="red"))

colorset1$SampleType<-rownames(colorset1)
colnames(colorset1)[which(names(colorset1) == "value")] <- "SampleType_Color"
colorset1

metadata_clean<-merge(metadata_clean, colorset1, by="SampleType")
head(metadata_clean)
metadata_clean$SampleType_Color <- as.character(metadata_clean$SampleType_Color)
rownames(metadata_clean)<-metadata_clean$SampleID
head(metadata_clean)

#cold2warm1<-get_palette(paste0("#",c("252A52", "66ADE5", "FFC465","BF1B0B")),k=10)
#names(cold2warm1) <- levels(metadata_clean$Depth_m)

fair_cols <- paste0("#",c("252A52", "66ADE5", "FFC465","BF1B0B"))
names(fair_cols) <- letters[1:4]
fair_ramp <- scales::colour_ramp(fair_cols)
fair_sat <- saturation(fair_ramp, 1)

#### Drop Singletons, Remove Eukaryotic Hits ####
# first we merge the ASV count object and the ASV taxonomy object together by column called "ASV_ID"
## then we need to melt the separate ASV & taxonomy so that we can rbind multiple data sets

# Drop singletons & zero count ASVs
dim(bac.ASV_counts_CLEAN)
colnames(bac.ASV_counts_CLEAN)
bac.ASV_counts_CLEAN<-bac.ASV_counts_CLEAN[which(rowSums(bac.ASV_counts_CLEAN[,-length(bac.ASV_counts_CLEAN)]) > 1),]
dim(bac.ASV_counts_CLEAN)

# merge CLEAN aka contaminants/controls removed count & taxa tables
bac.ASV_all<-merge(bac.ASV_counts_CLEAN,bac.ASV_taxa_CLEAN, by="ASV_ID")
head(bac.ASV_all)
dim(bac.ASV_all)
#bac.ASV_all<-bac.ASV_all[, !duplicated(colnames(bac.ASV_all))] # remove col duplicates

bac.dat<-melt(bac.ASV_all)
head(bac.dat)
colnames(bac.dat)[which(names(bac.dat) == "variable")] <- "SampleID"
colnames(bac.dat)[which(names(bac.dat) == "value")] <- "Count"

#bac.ASV_dat<-bac.ASV_dat[, !duplicated(colnames(bac.ASV_dat))] # remove col duplicates
#dim(bac.ASV_dat)

# Drop unknowns and eukaryotic hits
bac.dat<-subset(bac.dat, Kingdom!="Unknown") ## drop Unknowns from Kingdom
bac.dat<-subset(bac.dat, Phylum!="Unknown") ## drop Unknowns from Phylum
head(bac.dat)
dim(bac.dat)
bac.dat<-bac.dat[!bac.dat$SampleID=="Undetermined", ]
"Undetermined" %in% bac.dat$SampleID # sanity check

# Create ASV count file that is filtered of eukaryotic taxa - for later use
bac.dat.with.euks<-bac.dat

# Drop chloroplast & mitochondria seqs
bac.dat<-subset(bac.dat, Class!="Chloroplast") ## exclude Chloroplast sequences
bac.dat<-subset(bac.dat, Order!="Chloroplast") ## exclude Chloroplast sequences
bac.dat<-subset(bac.dat, Family!="Mitochondria") ## exclude Mitochondrial sequences just in case

# checking for eukaryotic hits takes a long time, caution!
#'Chloroplast' %in% bac.dat # check if Chloroplast counts are still in df, should be false because they've been removed
#'Mitochondria' %in% bac.dat # check if Mitochondria counts are still in df, should be false because they've been removed
#'Undetermined' %in% bac.dat # check if undetermined taxa in data frame
#NA %in% bac.ASV_dat

head(bac.dat) # No more eukaryotic hits here on out

#### Create Filtered ASV & Taxa Tables ####
# This ASV table has no contaminants & no ASVs found in controls, & no eukaryotic ASVs!
head(bac.dat) # No more eukaryotic hits here on out

bac.ASV_table<-base::as.data.frame(dcast(bac.dat, SampleID~ASV_ID, value.var="Count", fun.aggregate=sum)) ###
bac.ASV_table[1:5,1:5]
#bac.ASV_table[duplicated(rownames(bac.ASV_table))]
rownames(bac.ASV_table)<-bac.ASV_table$SampleID
bac.ASV_table[1:5,1:5]

# Create updated Taxa table (excluding eukaryotes)
b.taxa<-unique(subset(bac.dat,select=-c(SampleID, Count)))
b.taxa[duplicated(b.taxa$ASV_ID),] # check if there are duplicate ASV entries

#### Check if Metadata Exists for All Samples ####
# double check dimensions of metadata_clean and ASV table
dim(metadata_clean)
dim(bac.ASV_table)
# double check that the rownames exist + match
rownames(metadata_clean)
rownames(bac.ASV_table)

# Find rows in metadata_clean that are not in combined bacterial asv tables
setdiff(rownames(metadata_clean), rownames(bac.ASV_table)) # check rows in metadata_clean not in bac.ASV_table
setdiff(rownames(bac.ASV_table), rownames(metadata_clean)) # check rows in bac.ASV_table not in metadata_clean

meta_final<-metadata_clean[rownames(metadata_clean) %in% rownames(bac.ASV_table),]
"SN.SJER.M.42.50cm.03.01.22.C" %in% meta_final$SampleID # confirm the sample we don't have data for was dropped (FALSE means it's been dropped)
"SN.SJER.M.42.50cm.03.01.22.C" %in% rownames(meta_final)

#### Create Super DF w/ Taxa, Meta, & Count Data ####
bac.dat.all<-merge(bac.dat,meta_final,by="SampleID")

#### Separate Data by Project ####
bac.ASV_table[1:4,1:4]
head(meta_final)

DustDF<-meta_final[meta_final$SampleType=="Dust",] # pull out samples that are dust
dim(DustDF)
SoilDF<-meta_final[meta_final$SampleType=="Soil",] # pull out samples that are soil
LungDF<-meta_final[meta_final$SampleType=="Lung",] # pull out samples that are lung

# find DUST samples in ASV table that are in metadata
bac.ASV_table[rownames(bac.ASV_table) %in% rownames(DustDF),]
dim(bac.ASV_table[rownames(bac.ASV_table) %in% rownames(DustDF),])

# pull out dust sample ASV table
dust_b.ASV<-bac.ASV_table[rownames(bac.ASV_table) %in% rownames(DustDF),] # find DUST samples in ASV table that are in metadata
rownames(dust_b.ASV) %in% rownames(DustDF) # sanity check that this worked
dust_b.ASV[1:4,1:4]

write.table(dust_b.ASV, "data/SaltonSea_16S.V3V4_Dust_ASVs_Clean.tsv",
            sep="\t", quote=F, col.names=NA)
saveRDS(dust_b.ASV, file = "data/SaltonSea_16S.V3V4_Dust_ASVs_Clean.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

# pull out soil sample ASV table
soil_b.ASV<-bac.ASV_table[rownames(bac.ASV_table) %in% rownames(SoilDF),] # find DUST samples in ASV table that are in metadata
rownames(soil_b.ASV) %in% rownames(SoilDF) # sanity check that this worked
soil_b.ASV[1:4,1:4]

write.table(soil_b.ASV, "data/CZN_SaltonSea_SoilSamples_16S.V3V4_Clean_W23.tsv",
           sep="\t", quote=F, col.names=NA)
saveRDS(soil_b.ASV, file = "data/CZN_SaltonSea_SoilSamples_16S.V3V4_Clean_W23.tsv.rds", ascii = FALSE, version = NULL,
       compress = TRUE, refhook = NULL)

# pull out lung sample ASV table
lung_b.ASV<-bac.ASV_table[rownames(bac.ASV_table) %in% rownames(LungDF),] # find DUST samples in ASV table that are in metadata
rownames(lung_b.ASV) %in% rownames(LungDF) # sanity check that this worked
lung_b.ASV[1:4,1:4]

write.table(lung_b.ASV, "data/SaltonSea_LungSamples_16S.V3V4_Clean_W23.tsv",
           sep="\t", quote=F, col.names=NA)
saveRDS(lung_b.ASV, file = "data/SaltonSea_LungSamples_16S.V3V4_Clean_W23.tsv.rds", ascii = FALSE, version = NULL,
       compress = TRUE, refhook = NULL)

#### Upload/Update Dust Metadata ####
dust_meta<-as.data.frame(read_excel("data/Metadata_EnvMiSeqPlate_Winter23.xlsx", sheet="SSea_Dust_Metadata"), header=TRUE)
head(dust_meta)
unique(dust_meta$SampleMonth)
dust_meta$SampleMonth<-factor(dust_meta$SampleMonth, levels=c("July","August","September","October","November","December"))
colorset2 = melt(c(July="#2b9348",August="#ffd60a",September="#CA6702",October="#d00000",November="#6930c3",December="#03045e"))

colorset2$SampleMonth<-rownames(colorset2)
colnames(colorset2)[which(names(colorset2) == "value")] <- "SampMonth_Color"
colorset2

dust_meta<-merge(dust_meta, colorset2, by="SampleMonth")
head(dust_meta)
dust_meta$SampMonth_Color <- as.character(dust_meta$SampMonth_Color)
rownames(dust_meta)<-dust_meta$SampleID
head(dust_meta)

#### Scale Environmental Metadata ####
head(metadata)
meta_scaled<-metadata
meta_scaled[,8:15]<-scale(meta_scaled[,8:15],center=TRUE,scale=TRUE) # only scale chem env data
head(meta_scaled)


#### Save Global Env for Import into Other Scripts ####

save.image("data/MiSeq_16S.V3V4_W23_Data_Ready.Rdata") # save global env to Rdata file

## ^ this will be loaded into other scripts for downstream analyses
