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
# ^ you can remove all of this steps in the script titled 1a_Prep_Amplicon_Data_All.R

#### Import and Prepare Data for Analyses ####

## Import ALL env plate bacterial ASV count data
bac.ASV_table_og<-as.data.frame(readRDS("data/Amplicon/SaltonSeaDust_16S.V3V4_ASVTable_Robject.rds", refhook = NULL))
dim(bac.ASV_table_og)
bac.ASV_table_og[1:5,1:5]
#bac.ASV_table<-bac.ASV_table[, !duplicated(colnames(bac.ASV_table))] # remove col duplicates

# drop SB & RHB samples due to not having sample sample sizes as other sites
#define patterns to search for
dropped_samps <- c('SB', 'RHB')

# use patterns separated by | to remove from SampleID column
bac.ASV_table=bac.ASV_table_og[!grepl(paste(dropped_samps, collapse='|'),bac.ASV_table_og$SampleID),]

# create separate ASV table to see if we are dropping all the samples we need to drop
dropped_samps_ASVs=bac.ASV_table_og[grepl(paste(dropped_samps, collapse='|'),bac.ASV_table_og$SampleID),]

# compare dimensions of original ASV table, new ASV table with unwanted samples dropped, and separate ASV table containing only dropped samples
dim(bac.ASV_table_og) # original ASV table with all samples
dim(bac.ASV_table) # new ASV table with unwanted samples excluded
dim(dropped_samps_ASVs) # table of only unwanted samples

# one final sanity check to see if unwanted samples are in new ASV table
dropped_samps_ASVs$SampleID %in% bac.ASV_table$SampleID

## Import ASV taxonomic data
bac.ASV_tax<-data.frame(readRDS("data/Amplicon/EnvMiSeq_W23_16S.V3V4_ASVs_Taxonomy_dada2_Clean_Robject.rds", refhook = NULL))
head(bac.ASV_tax)

### Import & Update Metadata ####
dust_meta<-as.data.frame(read_excel("data/Metadata_EnvMiSeqPlate_Winter23.xlsx", sheet="SSea_Dust_Metadata_Updated"), header=TRUE)
head(dust_meta)

# confirm that categorical variables of interest are factors
dust_meta$CollectionYear<-factor(dust_meta$CollectionYear,levels=c("2020","2021"))
unique(dust_meta$SampleMonth)
dust_meta$SampleMonth<-factor(dust_meta$SampleMonth, levels=c("July","August","September","October","November","December"))

# create sample month palette
colorset2 = as.data.frame(t(data.frame("July"="#2b9348","August"="#ffd60a","September"="#CA6702","October"="#d00000","November"="#6930c3","December"="#03045e")))

colorset2$SampleMonth<-rownames(colorset2)
colnames(colorset2)[which(names(colorset2) == "V1")] <- "SampMonth_Color"
colorset2

dust_meta<-merge(dust_meta, colorset2, by="SampleMonth")
head(dust_meta)
dust_meta$SampMonth_Color <- as.character(dust_meta$SampMonth_Color)
rownames(dust_meta)<-dust_meta$SampleID
head(dust_meta)

# create Summer vs Fall palette
unique(dust_meta$Season_General)
colorset3 = as.data.frame(t(data.frame("Summer"="#4361ee","Fall"="#9a031e")))

colorset3$Season_General<-rownames(colorset3)
colnames(colorset3)[which(names(colorset3) == "V1")] <- "SeasonGen_Color"
colorset3

dust_meta<-merge(dust_meta, colorset3, by="Season_General")
head(dust_meta)
dust_meta$SeasonGen_Color <- as.character(dust_meta$SeasonGen_Color)
rownames(dust_meta)<-dust_meta$SampleID
dust_meta$Season_General<-factor(dust_meta$Season_General, levels=c("Summer","Fall"))
head(dust_meta)

# create more specific seasons palette
unique(dust_meta$Season_Specific)
colorset4 = as.data.frame(t(data.frame("Early.Summer"="#4cc9f0","Late.Summer"="#5e60ce","Late.Fall"="#e85d04","Fall.Winter"="#63003a")))

colorset4$Season_Specific<-rownames(colorset4)
colnames(colorset4)[which(names(colorset4) == "V1")] <- "SeasonSpec_Color"
colorset4

dust_meta<-merge(dust_meta, colorset4, by="Season_Specific")
head(dust_meta)
dust_meta$SeasonSpec_Color <- as.character(dust_meta$SeasonSpec_Color)
rownames(dust_meta)<-dust_meta$SampleID
dust_meta$Season_Specific<-factor(dust_meta$Season_Specific, levels=c("Early.Summer","Late.Summer","Early.Fall","Late.Fall","Fall.Winter"))

head(dust_meta)

# create year palette
unique(dust_meta$CollectionYear)
colorset5 = as.data.frame(t(data.frame("2020"="#751966","2021"="#135767")))

colorset5$CollectionYear<-gsub("X","",rownames(colorset5))
colnames(colorset5)[which(names(colorset5) == "V1")] <- "Year_Color"
colorset5

dust_meta<-merge(dust_meta, colorset5, by="CollectionYear")
head(dust_meta)
dust_meta$Year_Color <- as.character(dust_meta$Year_Color)
rownames(dust_meta)<-dust_meta$SampleID
dust_meta$CollectionYear<-factor(dust_meta$CollectionYear, levels=c("2020","2021"))

head(dust_meta)

# create site palette
unique(dust_meta$Site) # dropped SB, RHB
# "RHB"="#a11d33","SB"="#0077b6"
colorset6 = as.data.frame(t(data.frame("BDC"="#390099","DP"="#ffbd00","PD"="#eb5e28","WI"="#008000")))

colorset6$Site<-rownames(colorset6)
colnames(colorset6)[which(names(colorset6) == "V1")] <- "Site_Color"
colorset6

dust_meta<-merge(dust_meta, colorset6, by="Site")
head(dust_meta)
dust_meta$Site_Color <- as.character(dust_meta$Site_Color)
rownames(dust_meta)<-dust_meta$SampleID
dust_meta$Site<-factor(dust_meta$Site, levels=c("PD","BDC","DP","WI"))

head(dust_meta)

# create collection palette
unique(dust_meta$Seas_Coll_Year)
colorset7 = as.data.frame(t(data.frame("S.1.2020"="#14c9cb","S.2.2020"="#2962ff","S.3.2020"="#9500ff","F.1.2020"="#ff0059",
                   "S.1.2021"="#ff8c00","S.2.2021"="#0B6623","F.1.2021"="#ffd500")))
colorset7$Seas_Coll_Year<-rownames(colorset7)
colnames(colorset7)[which(names(colorset7) == "V1")] <- "SCY_Color"
colorset7

dust_meta<-merge(dust_meta, colorset7, by="Seas_Coll_Year")
head(dust_meta)
dust_meta$SCY_Color <- as.character(dust_meta$SCY_Color)
rownames(dust_meta)<-dust_meta$SampleID
dust_meta$Seas_Coll_Year<-factor(dust_meta$Seas_Coll_Year, levels=c("S.1.2020","S.2.2020","S.3.2020","F.1.2020","S.1.2021","S.2.2021","F.1.2021"))

head(dust_meta)

# create SampDate variable
dust_meta$SampDate<-interaction(dust_meta$SampleMonth,dust_meta$CollectionYear,sep=".")
head(dust_meta)
unique(dust_meta$SampDate)
length(unique(dust_meta$SampDate)) # 8 different sampdates
dust_meta$SampDate<-factor(dust_meta$SampDate, levels=c("July.2020","August.2020","October.2020","November.2020",
                                                        "July.2021","August.2021","September.2021","December.2021"))

# check brewer palette for 8 different colors for SampDate variable, then create colorset8 for SampDate
# display.brewer.pal(n = 8, name = 'Dark2')
# colorset8=data.frame(SampDate_Color=(brewer.pal(n = 8, name = "Dark2")))
# rownames(colorset8)<-unique(dust_meta$SampDate[order(dust_meta$SampDate)])
# colorset8
# colorset8$SampDate<-rownames(colorset8)

# create collection palette
unique(dust_meta$SampDate)
colorset8 = as.data.frame(t(data.frame("July.2020"="green2","August.2020"="orange","October.2020"="red","November.2020"="purple",
                   "July.2021"="darkgreen","August.2021"="darkorange3","September.2021"="red4","December.2021"="purple4")))
colorset8$SampDate<-rownames(colorset8)
colnames(colorset8)[which(names(colorset8) == "V1")] <- "SampDate_Color"
colorset8

dust_meta<-merge(dust_meta, colorset8, by="SampDate")
head(dust_meta)
dust_meta$SampDate_Color <- as.character(dust_meta$SampDate_Color)
rownames(dust_meta)<-dust_meta$SampleID

# check if SampDate still has its factor levels
unique(dust_meta$SampDate)

# dust_meta$SampDate<-factor(dust_meta$SampDate, levels=c("July.2020","August.2020","October.2020","November.2020",
#                                                         "July.2021","August.2021","September.2021","December.2021"))



#### Import Surface Type Frequencies Data ####
# data from Will Porter's group
# predicted source material by site, by time --> % from shrubs, urban, etc
# describes what land was covered by wind trajectory, indicating how much dust can be attributed from that area based on wind trajectory
#look at all points passed over by each trajectory on its way to the dust collector, and weight them
## weighted by wind speed (higher wind = more contribution), estimated erodibility of the surface (based on published sediment availability maps), & proximity (closer surfaces are more likely to contribute than distant ones due to gravitational settling)


SurfTypFreq<-as.data.frame(read.csv("data/Climate/PorterLab_SSD_SurfaceTypeFrequencies_Amplicon_SOI_Only.csv",header = TRUE, sep = ",", quote = "",))
head(SurfTypFreq)
rownames(SurfTypFreq)<-SurfTypFreq$SampleID
head(SurfTypFreq) #sanity check

unique(SurfTypFreq$STF_Date)
SurfTypFreq$STF_Date<-factor(SurfTypFreq$STF_Date, levels=c("July.2020","August.2020","October.2020","November.2020",
                                                        "August.2021","October.2021","December.2021"))


# back.traj<-read.fst("data/SaltonSeaDust_PorterLab_TrajectoryData.fst")
# colnames(back.traj)[which(names(back.traj) == "site")] <- "Site"
# head(back.traj)
# ^^ individual data points for each time steps that it was running backwards
# _i --> starting point of wind trajectory (aka origin aka the site)
# _add --> additional back trajectory points, moving backwards in time (5 min increments)
# so back trajectory info was used to determine which surface the wind traveled over
# this raw model output data

#### Create Super DF with Count, Taxa, & Metadata ####
# merge CLEAN aka contaminants/controls removed count & taxa tables
bac.ASV_table[1:5,1:5]
bac.ASV_tax[1:5,1:5]

bac.ASV_t.table<-as.data.frame(t(bac.ASV_table[,-1])) # transpose ASV table
bac.ASV_t.table[1:5,1:5]
bac.ASV_t.table$ASV_ID<-rownames(bac.ASV_t.table)
bac.ASV_t.table[1:5,5-ncol(bac.ASV_t.table):ncol(bac.ASV_t.table)]

bac.ASV_all<-merge(bac.ASV_t.table,bac.ASV_tax, by="ASV_ID") # merge ASV table & taxa table by ASV_IDs
head(bac.ASV_all)
dim(bac.ASV_all)
#bac.ASV_all<-bac.ASV_all[, !duplicated(colnames(bac.ASV_all))] # remove col duplicates

# melt combined ASV + taxa table to then merge with metadata by SampleID
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
bac.ASV_table<-bac.ASV_table[rownames(bac.ASV_table) %in% rownames(dust_meta),] # drop technical reps of samples
dim(bac.ASV_table)
bac.ASV_table[1:6,1:6]

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

b.dust.all$SampleID<-gsub("(.*\\..*)\\..*","\\1", b.dust.all$SampleID)
b.dust.all$SampleID # sanity check
head(b.dust.all)

# check if sampleIDs match between metadata & ASV table
rownames(dust_meta) %in% rownames(bac.ASV_table)

# reorder metadata based off of ASV table
dust_meta=dust_meta[rownames(bac.ASV_table),] ## will drop rows that are not shared by both dataframes!

#### Merge Surface Type Freqs with Metadata ####

dust.meta.surf<-merge(dust_meta,SurfTypFreq,by=c("Site","SampleID"))
head(dust.meta.surf)
dim(dust.meta.surf)

rownames(dust.meta.surf)<-dust.meta.surf$SampleID

#### Check ASV Count Distribution Across Samples ####
rowSums(bac.ASV_table[,-1])

#### Import Synoptic Climate Data & Merge with Metadata ####

clim.data<-data.frame(readRDS("data/Climate/SaltonSea_SynopticClimateData_Robject.rds", refhook = NULL))
head(clim.data)

# which column names match in the dfs
colnames(clim.data)[which(colnames(clim.data) %in% colnames(dust.meta.surf))]

dust.meta.all<-merge(clim.data,dust.meta.surf,by=c("CollectionNum","STID","Deploy_dth","Collect_dth"))
rownames(dust.meta.all)<-dust.meta.all$SampleID
rownames(dust.meta.all) # sanity check

# reorder metadata based off of ASV table
dust.meta.all=dust.meta.all[rownames(bac.ASV_table),] ## will drop rows that are not shared by both dataframes!
rownames(dust.meta.all) # sanity check

#### Scale Environmental Metadata ####
#head(metadata)
meta.all.scaled<-dust.meta.all
head(meta.all.scaled)
meta.all.scaled[,c(5:10,33:42)]<-scale(meta.all.scaled[,c(5:10,33:42)],center=TRUE,scale=TRUE) # only scale chem env data
head(meta.all.scaled)
# meta.all.scaled should still have SampleIDs as rownames!

#### Save Global Env for Import into Other Scripts ####

save.image("data/Amplicon/SSDust_16S.V3V4_W23_Data_Ready.Rdata") # save global env to Rdata file

## ^ this will be loaded into other scripts for downstream analyses
