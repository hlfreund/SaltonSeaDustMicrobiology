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
#load("data/MiSeq_16S.V3V4_W23_Data_Ready.Rdata")

#### Import and Prepare Data for Analyses ####

## Import ALL env plate bacterial ASV count data
bac.ASV_table<-as.data.frame(readRDS("data/SaltonSea_16S.V3V4_Dust_ASVs_Clean.rds", refhook = NULL))
dim(bac.ASV_table)
bac.ASV_table[1:5,1:5]
colnames(bac.ASV_table)
#bac.ASV_table<-bac.ASV_table[, !duplicated(colnames(bac.ASV_table))] # remove col duplicates

## Import ASV taxonomic data
bac.ASV_tax<-data.frame(readRDS("data/EnvMiSeq_W23_16S.V3V4_ASVs_Taxonomy_dada2_Robject.rds", refhook = NULL))
head(bac.ASV_tax)

bac.ASV_tax[is.na(bac.ASV_tax)]<- "Unknown" # turn all NAs into "Unkowns"
bac.ASV_tax$Species<-gsub("Unknown", "unknown", bac.ASV_tax$Species) # change uppercase Unkonwn to lowercase unknown for unknown species classification
head(bac.ASV_tax)
bac.ASV_tax$ASV_ID<-rownames(bac.ASV_tax) # create ASV ID column to use for merging data frames
head(bac.ASV_tax)

### Import & Update Metadata ####
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

#### Create Super DF with Count, Taxa, & Metadata ####
# merge CLEAN aka contaminants/controls removed count & taxa tables
bac.ASV_table[1:5,1:5]
bac.ASV_tax[1:5,1:5]

bac.ASV_t.table<-as.data.frame(t(bac.ASV_table[,-1]))
bac.ASV_t.table[1:5,1:5]
bac.ASV_t.table$ASV_ID<-rownames(bac.ASV_t.table)
bac.ASV_t.table[1:5,5-ncol(bac.ASV_t.table):ncol(bac.ASV_t.table)]

bac.ASV_all<-merge(bac.ASV_t.table,bac.ASV_tax, by="ASV_ID")
head(bac.ASV_all)
dim(bac.ASV_all)
#bac.ASV_all<-bac.ASV_all[, !duplicated(colnames(bac.ASV_all))] # remove col duplicates

bac.dat.dust<-melt(bac.ASV_all)
head(bac.dat.dust)
colnames(bac.dat.dust)[which(names(bac.dat.dust) == "variable")] <- "SampleID"
colnames(bac.dat.dust)[which(names(bac.dat.dust) == "value")] <- "Count"

b.dust.all<-merge(bac.dat.dust,dust_meta,by="SampleID")

#### Save Global Env for Import into Other Scripts ####

save.image("data/SSDust_16S.V3V4_W23_Data_Ready.Rdata") # save global env to Rdata file

## ^ this will be loaded into other scripts for downstream analyses
