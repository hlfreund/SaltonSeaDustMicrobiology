#### Set WD & Load Libraries ####
getwd() # use setwd("path/to/files") if you are not in the right directory
#setwd("/Volumes/HLF_SSD/Aronson_Lab_Data/Salton_Sea/SaltonSeaDust")
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
  library(fst)
})

## Info about these data
# dust samples were sequenced with some soil & lungs via UCR Core, MiSeq (targeting 16S V3V4)
# contaminants removed with decontam, also control ASVs removed as well as singletons, zero ASVs, and eukaryotic hits
# taxonomy file also had those contaminating ASVs removed

#### Import and Prepare Data for Analyses ####

## Import ALL env plate bacterial ASV count data
bac.ASV_table<-as.data.frame(readRDS("data/SaltonSeaDust_16S.V3V4_ASVTable_Robject.rds", refhook = NULL))
dim(bac.ASV_table)
bac.ASV_table[1:5,1:5]
#bac.ASV_table<-bac.ASV_table[, !duplicated(colnames(bac.ASV_table))] # remove col duplicates

## Import ASV taxonomic data
bac.ASV_tax<-data.frame(readRDS("data/EnvMiSeq_W23_16S.V3V4_ASVs_Taxonomy_dada2_Clean_Robject.rds", refhook = NULL))
head(bac.ASV_tax)

### Import & Update Metadata ####
dust_meta<-as.data.frame(read_excel("data/Metadata_EnvMiSeqPlate_Winter23.xlsx", sheet="SSea_Dust_Metadata_Updated"), header=TRUE)
head(dust_meta)

# confirm that categorical variables of interest are factors
dust_meta$CollectionYear<-factor(dust_meta$CollectionYear,levels=c("2020","2021"))
unique(dust_meta$SampleMonth)
dust_meta$SampleMonth<-factor(dust_meta$SampleMonth, levels=c("July","August","September","October","November","December"))

# create sample month palette
colorset2 = melt(c(July="#2b9348",August="#ffd60a",September="#CA6702",October="#d00000",November="#6930c3",December="#03045e"))

colorset2$SampleMonth<-rownames(colorset2)
colnames(colorset2)[which(names(colorset2) == "value")] <- "SampMonth_Color"
colorset2

dust_meta<-merge(dust_meta, colorset2, by="SampleMonth")
head(dust_meta)
dust_meta$SampMonth_Color <- as.character(dust_meta$SampMonth_Color)
rownames(dust_meta)<-dust_meta$SampleID
head(dust_meta)

# create Summer vs Fall palette
unique(dust_meta$Season_General)
colorset3 = melt(c(Summer="#4361ee",Fall="#9a031e"))

colorset3$Season_General<-rownames(colorset3)
colnames(colorset3)[which(names(colorset3) == "value")] <- "SeasonGen_Color"
colorset3

dust_meta<-merge(dust_meta, colorset3, by="Season_General")
head(dust_meta)
dust_meta$SeasonGen_Color <- as.character(dust_meta$SeasonGen_Color)
rownames(dust_meta)<-dust_meta$SampleID
dust_meta$Season_General<-factor(dust_meta$Season_General, levels=c("Summer","Fall"))
head(dust_meta)

# create more specific seasons palette
unique(dust_meta$Season_Specific)
colorset4 = melt(c(Early.Summer="#4cc9f0",Late.Summer="#5e60ce",Early.Fall="#e85d04",Late.Fall="#9a031e",Fall.Winter="#63003a"))

colorset4$Season_Specific<-rownames(colorset4)
colnames(colorset4)[which(names(colorset4) == "value")] <- "SeasonSpec_Color"
colorset4

dust_meta<-merge(dust_meta, colorset4, by="Season_Specific")
head(dust_meta)
dust_meta$SeasonSpec_Color <- as.character(dust_meta$SeasonSpec_Color)
rownames(dust_meta)<-dust_meta$SampleID
dust_meta$Season_Specific<-factor(dust_meta$Season_Specific, levels=c("Early.Summer","Late.Summer","Early.Fall","Late.Fall","Fall.Winter"))

head(dust_meta)

# create year palette
unique(dust_meta$CollectionYear)
colorset5 = melt(c("2020"="#751966","2021"="#135767"))

colorset5$CollectionYear<-rownames(colorset5)
colnames(colorset5)[which(names(colorset5) == "value")] <- "Year_Color"
colorset5

dust_meta<-merge(dust_meta, colorset5, by="CollectionYear")
head(dust_meta)
dust_meta$Year_Color <- as.character(dust_meta$Year_Color)
rownames(dust_meta)<-dust_meta$SampleID
dust_meta$CollectionYear<-factor(dust_meta$CollectionYear, levels=c("2020","2021"))

head(dust_meta)

# create site palette
unique(dust_meta$Site)
colorset6 = melt(c("BDC"="#390099","DP"="#ffbd00","PD"="#eb5e28","RHB"="#a11d33","SB"="#0077b6","WI"="#008000"))

colorset6$Site<-rownames(colorset6)
colnames(colorset6)[which(names(colorset6) == "value")] <- "Site_Color"
colorset6

dust_meta<-merge(dust_meta, colorset6, by="Site")
head(dust_meta)
dust_meta$Site_Color <- as.character(dust_meta$Site_Color)
rownames(dust_meta)<-dust_meta$SampleID
#dust_meta$Site<-factor(dust_meta$Site, levels=c("WI","RHB","SB","DP","PD","BDC"))

head(dust_meta)

# create SampDate variable
dust_meta$SampDate<-interaction(dust_meta$SampleMonth,dust_meta$CollectionYear,sep=".")
head(dust_meta)
unique(dust_meta$SampDate)
dust_meta$SampDate<-factor(dust_meta$SampDate, levels=c("July.2020","August.2020","October.2020","November.2020",
                                                        "July.2021","August.2021","September.2021","December.2021"))

#### Import Wind Back Trajectory Data ####
# data from Will Porter's group

# back.traj<-read.fst("data/SaltonSeaDust_PorterLab_TrajectoryData.fst")
# colnames(back.traj)[which(names(back.traj) == "site")] <- "Site"
# head(back.traj)
# ^^ individual data points for each time steps that it was running backwards
# _i --> starting point of wind trajectory (aka origin aka the site)
# _add --> additional back trajectory points, moving backwards in time (5 min increments)
# so back trajectory info was used to determine which surface the wind traveled over
# this raw model output data


perc.cov.bck.trj<-as.data.frame(read.csv("data/Porter_CollectorBackTrajectoriesData.csv",header=TRUE))
head(perc.cov.bck.trj)
# ^ predicted source material by site, by time --> % from shrubs, urban, etc
# describes what land was covered by wind trajectory, indicating how much dust can be attributed from that area based on wind trajectory
#look at all points passed over by each trajectory on its way to the dust collector, and weight them
## weighted by wind speed (higher wind = more contribution), estimated erodibility of the surface (based on published sediment availability maps), & proximity (closer surfaces are more likely to contribute than distant ones due to gravitational settling)

#### Scale Environmental Metadata ####
#head(metadata)
#meta_scaled<-metadata
#meta_scaled[,8:15]<-scale(meta_scaled[,8:15],center=TRUE,scale=TRUE) # only scale chem env data
#head(meta_scaled)

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
head(b.dust.all)

### Drop Duplicates from ASV Table ####
# sanity check that this indexing gets what we want
dim(bac.ASV_table[rownames(bac.ASV_table) %in% rownames(dust_meta),])

# update ASV table
bac.ASV_table<-bac.ASV_table[rownames(bac.ASV_table) %in% rownames(dust_meta),]

#### Drop Rep Letter from Sample Names ####

# v gsub() - (.*\\..*) means that we are keeping all .*.*., no matter how many there are in string
## (.*\\..*) & (.*\\..*\\..*) & (.*\\..*\\..*\\..*) ALL do the same thing
## whatever is outside of the parentheses will be removed; here it's last occurence of .*
dust_meta$SampleID<-gsub("(.*\\..*)\\..*","\\1", dust_meta$SampleID)
dust_meta$SampleID # sanity check
rownames(dust_meta)<-dust_meta$SampleID

bac.ASV_table$SampleID<-gsub("(.*\\..*)\\..*","\\1", bac.ASV_table$SampleID)
bac.ASV_table$SampleID # sanity check
rownames(bac.ASV_table)<-bac.ASV_table$SampleID

# check if sampleIDs match between metadata & ASV table
rownames(dust_meta) %in% rownames(bac.ASV_table)

# reorder metadata based off of ASV table
dust_meta=dust_meta[rownames(bac.ASV_table),] ## will drop rows that are not shared by both dataframes!

#### Save Global Env for Import into Other Scripts ####

save.image("data/SSDust_16S.V3V4_W23_Data_Ready.Rdata") # save global env to Rdata file

## ^ this will be loaded into other scripts for downstream analyses
