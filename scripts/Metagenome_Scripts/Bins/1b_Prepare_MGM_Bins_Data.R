#### Set WD & Load Libraries ####
getwd() # use setwd("path/to/files") if you are not in the right directory
#setwd("/Volumes/HLF_SSD/Aronson_Lab_data/Metagenomes/Salton_Sea/SaltonSeaDust")
suppressPackageStartupMessages({ # load packages quietly
  library(phyloseq)
  library(ggplot2)
  library(vegan)
  library(ggpubr)
  #library(scales)
  library(grid)
  library(data.table)
  library(ape)
  #library(apeglm)
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
  #library(wesanderson)
  #library(nationalparkcolors)
  library(fitdistrplus)
  library(logspline)
  #library(shades)
  #library(ALDEx2)
  library(rstatix)
  library(devtools)
  #library(decontam)
  #library(batchtools)
})

#load("SSD_mgm.bin_contigs_analysis_data.Rdata") # load Rdata to global env
#setwd("/bigdata/aronsonlab/shared/SaltonSea/Metagenomes/SeqCenter_3.30.2023/MGM_Analyses/Bins")

#### Import Custom Functions ####

counts_to_binary <- function(dataFrame){
  new_m <- matrix(nrow=dim(dataFrame)[1],ncol = dim(dataFrame)[2]) # create new matrix w/ same rows and cols as input dataframe
  ## dim(df)[1] gives you first dimensions (x aka rows), dim(df)[2] gives you second dimensions (y aka columns)

  for( currentRow in 1:nrow(dataFrame)){ # for every row
    for( currentCol in 1:ncol(dataFrame)){ # for every column

      if ( is.na(dataFrame[currentRow, currentCol]) & is.numeric(dataFrame[currentRow, currentCol])){ # if both row and col (specifies each cell) are NA, change val to 0
        new_m[currentRow, currentCol] = 0
        # is.numeric(df[currentRow,currentCol]) is to confirm each cell contains a numeric element
      } else if( is.numeric(dataFrame[currentRow, currentCol]) & dataFrame[currentRow, currentCol] > 0){ # if both row and col (specifies each cell) are > 0, change val to 1
        new_m[currentRow, currentCol] = 1
      } else if ( is.numeric(dataFrame[currentRow, currentCol]) & dataFrame[currentRow, currentCol] == 0){ # if both row and col (specifies each cell) == 0 , change val to 0
        new_m[currentRow, currentCol] = 0
      } else if ( is.character(dataFrame[currentRow, currentCol])){ # if both row and col (specifies each cell) == 0 , change val to 0
        new_m[currentRow, currentCol] = dataFrame[currentRow, currentCol]
      }
    }
  }
  new_df <- as.data.frame(new_m) #turns matrix into dataframe
  names(new_df) <- names(dataFrame) #names rows & cols of new dataframe to be same as row names and col names from input dataframe
  rownames(new_df) <- rownames(dataFrame)
  #  new_df2=new_df[,order(ncol(new_df):1)]
  new_df2=new_df[rownames(dataFrame),colnames(dataFrame)]
  return(new_df2) # ensures only output is the new dataframe
}

# function below is to convert large data.table NA to 0; https://stackoverflow.com/questions/7235657/fastest-way-to-replace-nas-in-a-large-data-table

remove_na <- function(x){
  dm <- as.matrix(x)
  dm[is.na(dm)] <- 0
  #dm<-as.matrix(sapply(dm[,-1], as.numeric))
  new_df <- as.data.frame(dm) #turns matrix into dataframe
  if ("SampleID" %in% colnames(new_df)){
    new_df[,!names(new_df) %in% c("SampleID")]<-as.data.frame(lapply(new_df[,!names(new_df) %in% c("SampleID")],as.numeric))
  } else{
    new_df<-as.data.frame(lapply(new_df,as.numeric))
  }
  #new_df[,-1]<-lapply(new_df[,-1],as.numeric)
  names(new_df) <- names(x) #names rows & cols of new dataframe to be same as row names and col names from input dataframe
  rownames(new_df) <- rownames(x)
  colnames(new_df) <- colnames(x)
  #  new_df2=new_df[,order(ncol(new_df):1)]
  return(new_df) # ensures only output is the new dataframe
  #data.table(dm)
}

remove_na_bins <- function(x){
  dm <- as.matrix(x)
  dm[is.na(dm)] <- 0
  #dm<-as.matrix(sapply(dm[,-1], as.numeric))
  new_df <- as.data.frame(dm) #turns matrix into dataframe
  if ("SampleID" %in% colnames(new_df) && "Bin_ID" %in% colnames(new_df)){
    new_df[,!names(new_df) %in% c("SampleID","Bin_ID")]<-as.data.frame(lapply(new_df[,!names(new_df) %in% c("SampleID","Bin_ID")],as.numeric))
  } else{
    new_df<-as.data.frame(lapply(new_df,as.numeric))
  }
  #new_df[,-1]<-lapply(new_df[,-1],as.numeric)
  names(new_df) <- names(x) #names rows & cols of new dataframe to be same as row names and col names from input dataframe
  rownames(new_df) <- rownames(x)
  colnames(new_df) <- colnames(x)
  #  new_df2=new_df[,order(ncol(new_df):1)]
  return(new_df) # ensures only output is the new dataframe
  #data.table(dm)
}

##save.image("SSD_mgm.bin_contigs_analysis_data.Rdata") # save global env to Rdata file

## Notes:
# code & info came from :
## https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#quick-start
## https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/07_practical.pdf
## https://www.reneshbedre.com/blog/deseq2.html
## https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html

#### Import Bin Coverage Data ####
mgm.bin_genes.cov<-fread(file = 'data/Metagenomes/SSD_Bins_Gene_Coverages_1.2.24.txt', sep='\t',header = TRUE)
dim(mgm.bin_genes.cov)
head(mgm.bin_genes.cov)

# Clean up SampleID names
unique(mgm.bin_genes.cov$SampleID) # make sure we do not have duplicate sample ID names
mgm.bin_genes.cov$SampleID<-gsub("_",".", mgm.bin_genes.cov$SampleID)
mgm.bin_genes.cov$SampleID<-gsub("(.*\\..*)\\..*","\\1", mgm.bin_genes.cov$SampleID) # drop rep letter at the end of SampleiD
mgm.bin_genes.cov$SampleID # sanity check
head(mgm.bin_genes.cov)

# Clean up BinID names
unique(mgm.bin_genes.cov$Bin_ID) # make sure we do not have duplicate sample ID names
mgm.bin_genes.cov$Bin_ID<-gsub("_",".", mgm.bin_genes.cov$Bin_ID)
mgm.bin_genes.cov$Bin_ID<-gsub("(.*\\..*)\\..*.(.bin.*)","\\1\\2", mgm.bin_genes.cov$Bin_ID) # drop the rep # for the Sample ID that is before the Bin ID, then save the SampleID and BinID
unique(mgm.bin_genes.cov$Bin_ID) # sanity check
head(mgm.bin_genes.cov)

# create unique list of bins & simpleIDs to set aside just in case
bin_ID_list<-unique(data.frame(Bin_ID=mgm.bin_genes.cov$Bin_ID,SampleID=mgm.bin_genes.cov$SampleID))

# count occurrences (frequency) of all traits across all mgms
n_kofxn_occur <- data.frame(table(mgm.bin_genes.cov$KO_Function)) # frequency of KO functions
n_kofxn_occur[n_kofxn_occur$Freq > 1,]

n_koid_occur1 <- data.frame(table(mgm.bin_genes.cov$KO_ID)) # frequency of KO IDs
n_koid_occur1[n_koid_occur1$Freq > 1,]

# Divide gene counts by gene length to account for sample differences in assembly
## we do this because we did not co-assemble contigs, so each KO assignments across samples may not come from genes with the same length
## need to account for gene length here first (since KO assignments can come from genes of different lengths)
mgm.bin_genes.cov$CovPerGene<-mgm.bin_genes.cov$ReadsPerGene/mgm.bin_genes.cov$Gene_Length
head(mgm.bin_genes.cov)

# Drop rows where we do not have a KO ID for the gene
mgm.bin_genes.cov.clean<-na.omit(mgm.bin_genes.cov, cols=c("KO_ID"))
head(mgm.bin_genes.cov.clean) # clean = no NAs
tail(mgm.bin_genes.cov.clean)

which(is.na(mgm.bin_genes.cov.clean$KO_ID)) # sanity check - are there NAs?

# create unique list of GeneIDs & KO IDs
gene_KOs<-unique(data.frame(Gene_ID=mgm.bin_genes.cov.clean$Gene_ID, KO_ID=mgm.bin_genes.cov.clean$KO_ID,KO_Function=mgm.bin_genes.cov.clean$KO_Function))
dim(gene_KOs)
head(gene_KOs)
which(is.na(gene_KOs)) # 0 NAs found

# how many unique KOs are there total, and how many genes are assigned to each KO?
KO_uniq_num_occur <- data.frame(table(gene_KOs$KO_ID))
KO_uniq_num_occur[KO_uniq_num_occur$Freq > 1,]

KOs_unique<-unique(data.frame(KO_ID=mgm.bin_genes.cov.clean$KO_ID,KO_Function=mgm.bin_genes.cov.clean$KO_Function))
dim(KOs_unique)
# ^ # of unique KO IDs is the same # as if you also include functions with KO IDs, meaning no KO IDs have been assigned to multiple functions

# create Bin ID x Gene ID count table, using coverage aka reads per gene that were divided by gene length
mgm.bin_gene.cov_table<-as.data.frame(dcast(mgm.bin_genes.cov.clean, Bin_ID~Gene_ID, value.var="CovPerGene"))
mgm.bin_gene.cov_table[,1:4] # sanity check
rownames(mgm.bin_gene.cov_table)<-mgm.bin_gene.cov_table$Bin_ID
mgm.bin_gene.cov_table[,1:4] # sanity check

# check for NAs in data.table
which(is.na(mgm.bin_gene.cov_table))

# convert all NAs to 0; will take a while
mgm.bin_gene.cov_table[,-1][is.na(mgm.bin_gene.cov_table[,-1])]<-0
mgm.bin_gene.cov_table[,1:6] # sanity check

# compare gene mean cov vs variance of the gene cov across samples
mean_gen.cov <- apply(mgm.bin_gene.cov_table[,-1], 2, mean)        #The second argument '2' of 'apply' function indicates the function being applied to columns Use '1' if applied to rows
variance_gen.cov <- apply(mgm.bin_gene.cov_table[,-1], 2, var)
df1 <- data.frame(mean_gen.cov, variance_gen.cov)

png(filename="figures/SSD_MAG_Bins_GeneMeanCov_vs_GeneCovVariance.png",width = 7, height = 7, units = "in",res = 800)
ggplot(df1) +
  geom_point(aes(x=mean_gen.cov, y=variance_gen.cov)) +
  scale_y_log10(limits = c(0.000000001,1e4)) +
  scale_x_log10(limits = c(0.000000001,1e4)) +
  geom_abline(intercept = 0, slope = 1, color="red")
dev.off()
# mean gen cov > variance gene cov until a certain point, then switches....

#### Create Bin x KO Matrix of Coverage Data ####
mgm.bin_gene.cov_table[1:4,1:4] # sanity check

mgm.gene.cov.melt<-reshape2::melt(mgm.bin_gene.cov_table, by="Bin_ID")
colnames(mgm.gene.cov.melt)[which(names(mgm.gene.cov.melt) == "variable")] <- "Gene_ID"
colnames(mgm.gene.cov.melt)[which(names(mgm.gene.cov.melt) == "value")] <- "CovPerGene"
head(mgm.gene.cov.melt)

# merge melted coverage data with KO IDs
head(gene_KOs)
mgm.gene.cov.fxn<-merge(mgm.gene.cov.melt, gene_KOs,by="Gene_ID")
tail(mgm.gene.cov.fxn) # confirmed that multiple genes are assigned to the same KOs

KO_num_occur_noNA <- data.frame(table(mgm.gene.cov.fxn$KO_ID))
KO_num_occur_noNA[KO_num_occur_noNA$Freq > 1,]
min(KO_num_occur_noNA$Freq)
max(KO_num_occur_noNA$Freq)

# recreate table so that it is Bin ID x KO ID instead of Gene ID
# create Bin ID x KO ID count table, using coverage aka reads per gene that were divided by gene length
mgm.bin_fxn.cov_table<-as.data.frame(dcast(mgm.gene.cov.fxn, Bin_ID~KO_ID, value.var="CovPerGene",fun.aggregate=sum))
mgm.bin_fxn.cov_table[,1:4] # sanity check
# TIP to make sure there are no NAs in KO IDs aka column names - type this into console: mgm.bin_fxn.cov_table$ and scroll through KO IDs
NA %in% colnames(mgm.bin_fxn.cov_table) # another sanity check that there are no NAs in KO IDs aka column names
"NA" %in% colnames(mgm.bin_fxn.cov_table) # yet another sanity check that there are no NAs in KO IDs aka column names

rownames(mgm.bin_fxn.cov_table)<-mgm.bin_gene.cov_table$Bin_ID
dim(mgm.bin_fxn.cov_table)
mgm.bin_fxn.cov_table[1:6,1:6]

# drop KOs with low coverage
mgm.bin.fxns.all<-mgm.bin_fxn.cov_table # save original Bin x KO table before you drop low KOs just in case
mgm.bin_fxn.cov_table<-mgm.bin_fxn.cov_table[,which(colSums(mgm.bin_fxn.cov_table[,-1])>=1)] # remove functions with less than 5 total counts across mgms
#mgm.bin_fxn.binary<-counts_to_binary(mgm.bin_fxn.cov_table[,-1]) # custom function to convert all counts to binary (presence/absence)
mgm.bin_fxn.cov_table[1:4,1:4] # sanity check
dim(mgm.bin_fxn.cov_table)

#save.image("SSD_MGM_Bins_analysis_data.Rdata")

# NOTE 1/9/24: mgm.bin_fxn.cov_table contains gene coverage summed by KOs

#### Import MAG Taxonomic Annotations ####

mgm.bin.tax<-as.data.frame(read.delim("data/Metagenomes/SSD_MAGs_TaxoAnnotation.tsv",header=TRUE,sep="\t", quote = "",na.strings=""))
mgm.bin.tax[is.na(mgm.bin.tax)]<-"Unknown"
head(mgm.bin.tax)

# Clean up BinID names
unique(mgm.bin.tax$Bin_ID) # make sure we do not have duplicate sample ID names
mgm.bin.tax$Bin_ID<-gsub("_",".", mgm.bin.tax$Bin_ID)
mgm.bin.tax$Bin_ID<-gsub("(.*\\..*)\\..*.(.bin.*)","\\1\\2", mgm.bin.tax$Bin_ID) # drop the rep # for the Sample ID that is before the Bin ID, then save the SampleID and BinID
unique(mgm.bin.tax$Bin_ID) # sanity check
head(mgm.bin.tax)

# Create SampleID column in taxo annotation df
mgm.bin.tax$SampleID<-gsub("\\.bin.*","", mgm.bin.tax$Bin_ID) # drop the rep # for the Sample ID that is before the Bin ID, then save the SampleID and BinID
head(mgm.bin.tax)

#### Drop MAGs > 90% Complete ####

# import CheckM results of only MAGs with >90% complete
highqual.MAGs<-as.data.frame(read_excel("data/Metagenomes/SSD_Samples_GoodBins.xlsx", sheet="SSD_Samples_HighQualityBins"), header=TRUE)
highqual.MAGs$Bin_ID<-gsub("(.*\\..*)\\..*.(.bin.*)","\\1\\2", highqual.MAGs$Bin_ID) # drop the rep # for the Sample ID that is before the Bin ID, then save the SampleID and BinID

# drop MAGs < 90% complete from taxa data
head(mgm.bin.tax)
mgm.bin.tax<-mgm.bin.tax[mgm.bin.tax$Bin_ID %in% highqual.MAGs$Bin_ID,]

# drop MAGs < 90% complete from functional coverage-per-bin data
mgm.bin_fxn.cov_table[1:4,1:4] # sanity check
mgm.bin_fxn.cov_table<-mgm.bin_fxn.cov_table[mgm.bin_fxn.cov_table$Bin_ID %in% highqual.MAGs$Bin_ID,]

#### Import Surface Type Frequencies ####
# upload surface type frequencies from Will Porter and Yanning Miao

SurfTypFreq_mgm<-as.data.frame(read.csv("data/Climate/PorterLab_SSD_SurfaceTypeFrequencies_MGM_SOI_Only.csv",header = TRUE, sep = ",", quote = "",))
head(SurfTypFreq_mgm)

#### Import & Update Metadata ####
dust_mgm_meta<-as.data.frame(read_excel("data/Metagenomes/SSD_MGM_Metadata.xlsx", sheet="Sheet1"), header=TRUE)
head(dust_mgm_meta)

# Clean up SampleID names
dust_mgm_meta$SampleID<-gsub("(.*\\..*)\\..*","\\1", dust_mgm_meta$SampleID)
dust_mgm_meta$SampleID # sanity check
head(dust_mgm_meta)

# confirm that categorical variables of interest are factors
dust_mgm_meta$CollectionYear<-factor(dust_mgm_meta$CollectionYear,levels=c("2020","2021"))
unique(dust_mgm_meta$SampleMonth)
dust_mgm_meta$SampleMonth<-factor(dust_mgm_meta$SampleMonth, levels=c("July","August","September","November","December"))

# create sample month palette
colorset2 = melt(c(July="#2b9348",August="#ffd60a",September="#CA6702",November="#6930c3",December="#03045e"))

colorset2$SampleMonth<-rownames(colorset2)
colnames(colorset2)[which(names(colorset2) == "value")] <- "SampMonth_Color"
colorset2

dust_mgm_meta<-merge(dust_mgm_meta, colorset2, by="SampleMonth")
head(dust_mgm_meta)
dust_mgm_meta$SampMonth_Color <- as.character(dust_mgm_meta$SampMonth_Color)
#rownames(dust_mgm_meta)<-dust_mgm_meta$SampleID
head(dust_mgm_meta)

# create Summer vs Fall palette
unique(dust_mgm_meta$Season_General)
colorset3 = melt(c(Summer="#4361ee",Fall="#9a031e"))

colorset3$Season_General<-rownames(colorset3)
colnames(colorset3)[which(names(colorset3) == "value")] <- "SeasonGen_Color"
colorset3

dust_mgm_meta<-merge(dust_mgm_meta, colorset3, by="Season_General")
head(dust_mgm_meta)
dust_mgm_meta$SeasonGen_Color <- as.character(dust_mgm_meta$SeasonGen_Color)
#rownames(dust_mgm_meta)<-dust_mgm_meta$SampleID
dust_mgm_meta$Season_General<-factor(dust_mgm_meta$Season_General, levels=c("Summer","Fall"))
head(dust_mgm_meta)

# create more specific seasons palette
unique(dust_mgm_meta$Season_Specific)
colorset4 = melt(c(Early.Summer="#4cc9f0",Late.Summer="#5e60ce",Late.Fall="#9a031e",Fall.Winter="#63003a"))

colorset4$Season_Specific<-rownames(colorset4)
colnames(colorset4)[which(names(colorset4) == "value")] <- "SeasonSpec_Color"
colorset4

dust_mgm_meta<-merge(dust_mgm_meta, colorset4, by="Season_Specific")
head(dust_mgm_meta)
dust_mgm_meta$SeasonSpec_Color <- as.character(dust_mgm_meta$SeasonSpec_Color)
#rownames(dust_mgm_meta)<-dust_mgm_meta$SampleID
dust_mgm_meta$Season_Specific<-factor(dust_mgm_meta$Season_Specific, levels=c("Early.Summer","Late.Summer","Late.Fall","Fall.Winter"))

head(dust_mgm_meta)

# create year palette
unique(dust_mgm_meta$CollectionYear)
colorset5 = melt(c("2020"="#751966","2021"="#135767"))

colorset5$CollectionYear<-rownames(colorset5)
colnames(colorset5)[which(names(colorset5) == "value")] <- "Year_Color"
colorset5

dust_mgm_meta<-merge(dust_mgm_meta, colorset5, by="CollectionYear")
head(dust_mgm_meta)
dust_mgm_meta$Year_Color <- as.character(dust_mgm_meta$Year_Color)
#rownames(dust_mgm_meta)<-dust_mgm_meta$SampleID
dust_mgm_meta$CollectionYear<-factor(dust_mgm_meta$CollectionYear, levels=c("2020","2021"))

head(dust_mgm_meta)

# create site palette
unique(dust_mgm_meta$Site)
colorset6 = melt(c("BDC"="#390099","DP"="#ffbd00","PD"="#eb5e28","WI"="#008000"))

colorset6$Site<-rownames(colorset6)
colnames(colorset6)[which(names(colorset6) == "value")] <- "Site_Color"
colorset6

dust_mgm_meta<-merge(dust_mgm_meta, colorset6, by="Site")
head(dust_mgm_meta)
dust_mgm_meta$Site_Color <- as.character(dust_mgm_meta$Site_Color)
#rownames(dust_mgm_meta)<-dust_mgm_meta$SampleID
dust_mgm_meta$Site<-factor(dust_mgm_meta$Site, levels=c("PD","BDC","DP","WI"))

head(dust_mgm_meta)

# create collection palette
unique(dust_mgm_meta$Seas_Coll_Year)
colorset7 = melt(c("S.1.2020"="#14c9cb","S.2.2020"="#2962ff","F.1.2020"="#ff0059","S.1.2021"="#ff8c00","S.2.2021"="#0B6623","F.1.2021"="#ffd500"))
colorset7$Seas_Coll_Year<-rownames(colorset7)
colnames(colorset7)[which(names(colorset7) == "value")] <- "SCY_Color"
colorset7

dust_mgm_meta<-merge(dust_mgm_meta, colorset7, by="Seas_Coll_Year")
head(dust_mgm_meta)
dust_mgm_meta$SCY_Color <- as.character(dust_mgm_meta$SCY_Color)
#rownames(dust_mgm_meta)<-dust_mgm_meta$SampleID
dust_mgm_meta$Seas_Coll_Year<-factor(dust_mgm_meta$Seas_Coll_Year, levels=c("S.1.2020","S.2.2020","F.1.2020","S.1.2021","S.2.2021","F.1.2021"))

head(dust_mgm_meta)

# create SampDate variable
dust_mgm_meta$SampDate<-interaction(dust_mgm_meta$SampleMonth,dust_mgm_meta$CollectionYear,sep=".")
head(dust_mgm_meta)
unique(dust_mgm_meta$SampDate)
dust_mgm_meta$SampDate<-factor(dust_mgm_meta$SampDate, levels=c("July.2020","August.2020","November.2020",
                                                        "July.2021","August.2021","September.2021","December.2021"))
head(dust_mgm_meta)

# are there any NAs in the metadata now? sanity check
which(is.na(dust_mgm_meta))

save.image("data/Metagenomes/SSD_MGM_Bins_analysis_data.Rdata")

#### Merge All Metadata ####
# merge the dust collection data + surface type frequency data
all_mgm_meta<-merge(dust_mgm_meta, SurfTypFreq_mgm, by=c("SampleID","Site"))
head(all_mgm_meta)
rownames(all_mgm_meta)<-all_mgm_meta$SampleID

# are there any NAs in the merged metadata now? sanity check
which(is.na(all_mgm_meta))

# pull out all Bin IDs, generate the SampleIDs, and merge with metadata
mgm.bin_fxn.cov_table[1:6,1:6]
mgm.bin.list<-data.frame(Bin_ID=mgm.bin_fxn.cov_table$Bin_ID)
# Create SampleID column in bin list df
mgm.bin.list$SampleID<-gsub("\\.bin.*","", mgm.bin.list$Bin_ID) # drop the rep # for the Sample ID that is before the Bin ID, then save the SampleID and BinID
head(mgm.bin.list)

# merge bin list with all the metadata
all_bin_meta<-merge(all_mgm_meta,mgm.bin.list,by="SampleID")
head(all_bin_meta)
rownames(all_bin_meta)<-all_bin_meta$Bin_ID
head(all_bin_meta)

# merge metadata with MAG taxonomic annotation
all_bin_meta.taxa<-merge(all_bin_meta,mgm.bin.tax,by=c("Bin_ID","SampleID"))
head(all_bin_meta.taxa)

#### Scale Chem Data in Metadata ####
# head(SurfTypFreq_mgm)
# #all_mgm_meta<-subset(SurfTypFreq_mgm, select=-c(Salinity_ppt)) # drop salinity, for now
# all_mgm_meta<-SurfTypFreq_mgm
#
# head(all_mgm_meta)
# all_mgm_meta[,8:17]<-as.data.frame(scale(all_mgm_meta[,8:17], center=TRUE, scale=TRUE)) #centering before scaling
# all_mgm_meta
#
# # create numeric version of depth column for later analyses
# all_mgm_meta$Depth.num<-as.numeric(as.character(all_mgm_meta$Depth_m))

#save.image("SSD_mgm.bin_contigs_analysis_data.Rdata")

#### Import Gene Info from KEGG ####

# all_goi.kegg<-read.table("/bigdata/Metagenomes/aronsonlab/shared/SaltonSea/Metagenomes/SeqCenter_3.30.2023/MGM_Analyses/KO_KEGG_Files/Genes_of_Interest_All_KOs_Pathway_Cycle_KEGG.txt", header = TRUE, quote = "", sep = "\t", dec = ".")

all_goi.kegg<-read.table("data/Metagenomes/KO_KEGG_Files/Genes_of_Interest_All_KOs_Pathway_Cycle_KEGG.txt", header = TRUE, quote = "", sep = "\t", dec = ".")
carb.kegg<-read.table("data/Metagenomes/KO_KEGG_Files/CarbonFixation_KOs_KEGG.txt", header = TRUE, quote = "", sep = "\t", dec = ".")
sulf.kegg<-read.table("data/Metagenomes/KO_KEGG_Files/Sulfur_KOs_KEGG.txt", header = TRUE, quote = "", sep = "\t", dec = ".")
nitro.kegg<-read.table("data/Metagenomes/KO_KEGG_Files/N_KOs_Pathway_Cycle_KEGG.txt", header = TRUE, quote = "", sep = "\t", dec = ".")
arsen.kegg<-read.table("data/Metagenomes/KO_KEGG_Files/Arsenic_KOs_KEGG.txt", header = TRUE, quote = "", sep = "\t", dec = ".")
tempshock.kegg<-read.table("data/Metagenomes/KO_KEGG_Files/TempShock_KOs_KEGG.txt", header = TRUE, quote = "", sep = "\t", dec = ".")
osmo.kegg<-read.table("data/Metagenomes/KO_KEGG_Files/Osmoprot_KOs_KEGG.txt", header = TRUE, quote = "", sep = "\t", dec = ".")
selen.kegg<-read.table("data/Metagenomes/KO_KEGG_Files/Selenium_KOs_KEGG.txt", header = TRUE, quote = "", sep = "\t", dec = ".")
metal.re.kegg<-read.table("data/Metagenomes/KO_KEGG_Files/MetalResistance_KOs_KEGG.txt", header = TRUE, quote = "", sep = "\t", dec = ".")
photo.kegg<-read.table("data/Metagenomes/KO_KEGG_Files/Photo_KO_KEGG.txt", header = TRUE, quote = "", sep = "\t", dec = ".")
aero.rep.kegg<-read.table("data/Metagenomes/KO_KEGG_Files/OxidativePhosphorylation_KO_KEGG.txt", header = TRUE, quote = "", sep = "\t", dec = ".")
lps.kegg<-read.table("data/Metagenomes/KO_KEGG_Files/LPS_Biosynthesis_KOs_KEGG.txt", header = TRUE, quote = "", sep = "\t", dec = ".")
UV.DNA.rep.kegg<-read.table("data/Metagenomes/KO_KEGG_Files/UVDamageRepair_KOs_KEGG.txt", header = TRUE, quote = "", sep = "\t", dec = ".")
spor.kegg<-read.table("data/Metagenomes/KO_KEGG_Files/Sporulation_KOs_KEGG.txt", header = TRUE, quote = "", sep = "\t", dec = ".")
QuorSens.kegg<-read.table("data/Metagenomes/KO_KEGG_Files/QuorumSensing_KOs_KEGG.txt", header = TRUE, quote = "", sep = "\t", dec = ".")

#save.image("SSD_MGM_Bins_analysis_data.Rdata")

#### Check Summed Gene Coverage by KO Distribution in MGMs ####
mgm.bin_fxn.cov_table[,1:4] # sanity check
mgm.bin_cov_t<-as.data.frame(t(mgm.bin_fxn.cov_table[,-1]))
class(mgm.bin_cov_t$PD.D.11.6.20.bin.7)
# mgm.bin_cov_t <- as.data.frame(sapply(mgm.bin_cov_t, as.numeric)) # convert integer to numeric across df
# class(mgm.bin_cov_t$PD.D.11.6.20.bin.7) # sanity check

head(mgm.bin_cov_t)
#descdist(mgm.bin_cov_t$SSD.4.13.22.5m, discrete = TRUE)
#descdist(mgm.bin_cov_t$SSD.12.22.21.5m, discrete = TRUE)
#descdist(mgm.bin_cov_t$SSD.8.24.21.0m, discrete = TRUE)

# check distribution of gene # x gene cov with two samples
png(filename="figures/SSD_PD.D.11.6.20.bin.7_KOCov_histogram.png",width = 9, height = 9, units = "in",res = 800)

ggplot(mgm.bin_cov_t, aes(PD.D.11.6.20.bin.7)) +
  geom_histogram(color = "#000000", fill = "#0099F8")
## most KOs have very low coverage, few have higher coverage
dev.off()

ggplot(mgm.bin_cov_t, aes(WI.D.9.18.21.bin.9)) +
  geom_histogram(color = "#000000", fill = "#0099F8")

ggplot(mgm.bin_cov_t, aes(DP.D.12.8.21.bin.24)) +
  geom_histogram(color = "#000000", fill = "#0099F8")

ggplot(mgm.bin_cov_t) +
  geom_histogram(aes(x = PD.D.11.6.20.bin.7)) +
  xlab("KO summed gene coverage") +
  ylab("Number of KOs")
dev.off()

# # histogram of gene cov
# hist(mgm.bin_cov_t$PD.D.11.6.20.bin.7, col="blue")
# # visualize Q-Q plot for
# qqnorm(mgm.bin_cov_t$PD.D.11.6.20.bin.7, pch = 1, frame = FALSE)
# qqline(mgm.bin_cov_t$PD.D.11.6.20.bin.7, col = "red", lwd = 2)
#
# # histogram of gene cov
# hist(mgm.bin_cov_t$WI.D.9.18.21.bin.9, col="blue")
# # visualize Q-Q plot
# qqnorm(mgm.bin_cov_t$WI.D.9.18.21.bin.9, pch = 1, frame = FALSE)
# qqline(mgm.bin_cov_t$WI.D.9.18.21.bin.9, col = "red", lwd = 2)

# compare mean KO summed gene cov vs variance of the KO summed gene cov across samples
mean_KO_cov <- apply(mgm.bin_cov_t, 1, mean)        #The second argument '1' of 'apply' function indicates the function being applied to rows. Use '2' if applied to columns
variance_KO_cov <- apply(mgm.bin_cov_t, 1, var)
df_ko <- data.frame(mean_KO_cov, variance_KO_cov)

png(filename="figures/SSD_Contigs_Mean_KO_SummedCov_vs_Variance_KO_SummedCov.png",width = 7, height = 7, units = "in",res = 800)
ggplot(df_ko) +
  geom_point(aes(x=mean_KO_cov, y=variance_KO_cov)) +
  scale_y_log10(limits = c(0.00001,1e9)) +
  scale_x_log10(limits = c(0.00001,1e9)) +
  geom_abline(intercept = 0, slope = 1, color="red")
dev.off()
#save.image("SSD_MGM_Bins_analysis_data.Rdata")

#### Drop KOs with Lowest Summed Gene Coverages Before Transformation ####
# drop KO IDs with lowest summed coverages before we transform
keepKOs <- colSums(round(mgm.bin_fxn.cov_table[,-1])) >= 3
mgm.bin_fxn.cov_table.no_lows <- mgm.bin_fxn.cov_table[,-1][,keepKOs] # exclude Bin_ID column, then use keepKOs to keep the KO columns you want -- otherwise, this indexing will screw up your output if you keep Bin_ID column in...
mgm.bin_fxn.cov_table.no_lows[1:4,1:4]

length(mgm.bin_fxn.cov_table.no_lows) # we've dropped almost half the KOs becuase they had less than 3 summed gene coverage (so pretty low representation)
length(mgm.bin_fxn.cov_table[,-1])

### NOTE ^ the table mgm.bin_fxn.cov_table.no_lows will be used for transformations & normalizations that are NOT done using DESeq2
### this is because the DESeq2 object prep step allows us to drop the KOs with the lowest summed gene coverages

#### Prepare Contig Feature Count Data for Normalization w/ DESeq2 ####

head(all_bin_meta)

#convert data frame to numeric (must be numeric before normalization step)
class(mgm.bin_fxn.cov_table$K18890) #is data numeric?
#mgm.bin_fxn.cov_table[,-1] <- as.data.frame(sapply(mgm.bin_fxn.cov_table[,-1], as.numeric))

# convert coverage table into matrix for DESeq2 normalization
sum.cov_matrix<-as.matrix(mgm.bin_fxn.cov_table[,!names(mgm.bin_fxn.cov_table) %in% c("Bin_ID")]) # convert count table into matrix & exclude column called Bin_ID (not sure where it is in df)
sum.cov_matrix_transp<-t(sum.cov_matrix) # transpose matrix so KO_IDs are rows, samples are columns
# ^ will be used in DESeq2 functions
sum.cov_matrix_transp[1:3,1:3] #sanity check

# check if rownames in metadata (Bin_ID) match column names in count data
rownames(all_bin_meta) %in% colnames(sum.cov_matrix_transp)
dim(all_bin_meta)
dim(sum.cov_matrix_transp)

# create the DESeq DataSet object for DGE analysis
# DESeq2 needs whole # data, so need raw read counts, NOT coverage for these tables...questionable
# mgm.bin_dds has the scaled coverage that was calculated by dividing reads from featureCounts by gene lengths
mgm.bin_dds<-DESeqDataSetFromMatrix(countData=(round(sum.cov_matrix_transp)+1), colData = all_bin_meta, design = ~ 1)
# REMEMBER: we are using gene coverages summed up per KO, and are normalizing those summed coverages
# this function is rounding the summed gene coverages per KO into whole #s
# added pseudocount of 1 so that we could get the geometric means without an error later on

# design = ~ 1 means no design
head(counts(mgm.bin_dds)) # check if summed coverages below match the rounded counts in mgm.bin_dds
head(sum.cov_matrix_transp)

colSums(counts(mgm.bin_dds)) %>% barplot
rowSums(counts(mgm.bin_dds)) %>% barplot

# add rownames from KO IDs to DESeq Data set object
length(mgm.bin_dds) # # of elements in object mgm.bin_dds
rownames(mgm.bin_dds) # rownames of DESeq2 object are KO IDs
featureData <- data.frame(KO_ID=rownames(mgm.bin_dds)) # create object of KO IDs
mcols(mgm.bin_dds) <- DataFrame(mcols(mgm.bin_dds), featureData)
mcols(mgm.bin_dds) #
# mcols notes: Get or set the metadata columns. If use.names=TRUE and the metadata columns are not NULL, then the names of x are propagated as the row names of the returned DataFrame object.
# basically, mcols(mgm.bin_dds) should be the same as the rownames(mgm.bin_dds)

# drop KO IDs with lowest summed coverages (that have been rounded)
keepKOs <- rowSums(counts(mgm.bin_dds)) >= 3
mgm.bin_dds <- mgm.bin_dds[keepKOs,]
length(mgm.bin_dds) # we've dropped almost half the KOs becuase they had less than 3 summed gene coverage (so pretty low representation)

# Estimate size factor - aka normalization factor, divide all read counts by each size factor per sample
mgm.bin_dds <- estimateSizeFactors(mgm.bin_dds,type="ratio")
## To calculate size factor in DESeq2:
# calculates geometric mean of each gene in each sample and uses this as a pseudoreference
# calculates ratio of each sample by dividing each gene count by it's pseudoreference in each sample
# The median value of all ratios for a given sample is taken as the normalization factor (size factor)
# The differentially expressed genes should not affect the median value
# median of ratios method makes the assumption that not ALL genes are differentially expressed; therefore, the normalization factors should account for sequencing depth and RNA composition of the sample
## (large outlier genes will not represent the median ratio values)
# more here: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html

sizeFactors(mgm.bin_dds) # check the size factors

#plot(sizeFactors(mgm.bin_dds), colSums(counts(mgm.bin_dds)))

# Does sequencing depth influence normalization?
# par(mfrow=c(1,2)) # to plot the two box plots next to each other
# boxplot(log2(counts(mgm.bin_dds[1:100,1:10])), notch=TRUE,
#         main = "Non-normalized read counts\n(log-transformed)",
#         ylab="read counts")
# boxplot(log2(counts(mgm.bin_dds[1:100,1:10], normalize= TRUE)), notch=TRUE,
#         main = "Size-factor-normalized read counts\n(log-transformed)",
#         ylab="read counts")
# dev.off()

#save.image("SSD_MGM_Bins_analysis_data.Rdata")

#### Median-Ratio Normalized - Gene in Contigs ####
# NOTES on median ratio normalization
# Briefly, the size factor is calculated by first dividing the observed counts for each sample by its geometric mean.
## Remember: Geometric mean: multiply all the numbers together and take the nth root, where n = number of values multiplied together
## if we had 100 genes, we'd take the 100th root of all gene counts multiplied together to get geometric mean for each gene
# The size factor is then calculated as the median of this ratio for each sample.
# This size factor then used for normalizing raw count data for each sample.
# more at these links: https://scienceparkstudygroup.github.io/research-data-management-lesson/median_of_ratios_manual_normalization/index.html, https://www.reneshbedre.com/blog/expression_units.html#deseq-or-deseq2-normalization-median-of-ratios-method,
## https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html)

counts(mgm.bin_dds) #sanity check, what data are we starting with here

sizeFactors(mgm.bin_dds) # will be used to normalize counts
## to normalize counts, we divide each raw count value in a given sample by that sample’s normalization factor (aka size factor) to generate normalized count values.
### This is performed for all count values (every gene in every sample)
mgm.bin_fxn_mr <- counts(mgm.bin_dds, normalized=TRUE) ## median-ratio transformation is how DESeq2 Normalizes counts!
mgm.bin_fxn_mr[1:6,1:10]
#write.table(mgm.bin_fxn_counts.norm, file="data/Metagenomes/Metagenomes/mgm.bin_NoBins_MedianRatio_GeneCounts_2.27.23.txt", sep="\t", quote=F, col.names=T)

# turn this normalized table of summed gene coverages by KOs into a Sample x KO ID table (for vegan functions)
mgm.bin.mr<-as.data.frame(t(mgm.bin_fxn_mr))
mgm.bin.mr[1:4,1:4]

#mgm.bin.mr$Bin_ID<-rownames(mgm.bin.mr)
mgm.bin.mr[1:4,1:4]

#### Variance Stabilizing Transformation - Gene in Contigs ####
# using DESea2 for VST -- so if you have not done this yet, go to "Prepare Contig Feature Count Data for Normalization w/ DESeq2" section and create mgm.bin_dds
# notes on VST done by DESeq2 here: https://rdrr.io/bioc/DESeq2/man/varianceStabilizingTransformation.html
# more info on VST here: https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-data-transformations
# we are running blind = FALSE because size factors for mgm.bin_dds were already calculated previously

# you should be able to use matrix or DESeq2 object for this next function, but matrix was not working?
# variance stabilizing transformation
mgm.bin_fxn_vst <- varianceStabilizingTransformation(mgm.bin_dds, blind = FALSE)
# note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.

assay(mgm.bin_fxn_vst)[1:4,1:4] #see output of VST

# turn this VST output of summed gene coverages by KOs into a Sample x KO ID table (for vegan functions)
mgm.bin.vst<-as.data.frame(t(assay(mgm.bin_fxn_vst)))
mgm.bin.vst[1:4,1:4]

#### Centered Log Ratio Transformation - Gene in Contigs ####
mgm.bin_fxn.cov_table[1:4,1:4]
# ^ table contains gene coverages (reads mapped to each gene divided by gene length) summed by KO across samples
mgm.bin_fxn.cov_table.no_lows[1:4,1:4] # << table with lowest summed gene coverage KOs dropped

# use table that has dropped KOs with summed gene coverage < 3
# df must have rownames are Bin_IDs, columns are ASV IDs for vegan functions below\
mgm.bin.clr<-decostand(mgm.bin_fxn.cov_table.no_lows,method = "clr",pseudocount=1) #CLR transformation
mgm.bin.clr[1:4,1:4]
# NOTE: CLR transformation does not treat all 0s equally, it has to do with the pseudocount that's added before transformation
# The method can operate only with positive data; a common way to deal with zeroes is to add pseudocount, either by adding it manually to the input data, or by using the argument pseudocount as in decostand(x, method = "clr", pseudocount = 1).
# Adding pseudocount will inevitably introduce some bias; see the rclr method for one available solution

# check rownames of CLR transformed KO-summed coverage data & metadata
rownames(mgm.bin.clr) %in% rownames(all_bin_meta)

#### Robust Centered Log Ratio Transformation - Gene in Contigs ####
mgm.bin_fxn.cov_table[1:4,1:4]
# ^ table contains gene coverage, Sample IDs as rows & genes as columns
mgm.bin_fxn.cov_table.no_lows[1:4,1:4] # << table with lowest summed gene coverage KOs dropped

# use table that has dropped KOs with summed gene coverage < 3
# df must have rownames are Bin_IDs, columns are ASV IDs for vegan functions below\
mgm.Rclr<-decostand(mgm.bin_fxn.cov_table.no_lows,method = "rclr") #CLR transformation
mgm.Rclr[1:4,1:4]
# NOTE: Robust CLR just excludes 0s and performs CLR transformation without pseudocount
# robust clr ("rclr") is similar to regular clr (see above) but allows data that contains zeroes.
# This method does not use pseudocounts, unlike the standard clr. Robust clr divides the values by geometric mean of the observed features; zero values are kept as zeroes, and not taken into account.
#In high dimensional data, the geometric mean of rclr is a good approximation of the true geometric mean

# check rownames of CLR transformed ASV data & metadata
rownames(mgm.Rclr) %in% rownames(all_bin_meta)

#### Copies Per Million Transformation - Gene in Contigs ####
mgm.bin_fxn.cov_table[1:4,1:4] # sanity check
mgm.bin_fxn.cov_table.no_lows[1:4,1:4] # << table with lowest summed gene coverage KOs dropped

# use table that has dropped KOs with summed gene coverage < 3
mgm.bin_fxn.cov_t.table<-as.data.frame(t(mgm.bin_fxn.cov_table.no_lows))
mgm.bin_fxn.cov_t.table[1:4,1:4]

mgm.bin_fxn_cpm<-(mgm.bin_fxn.cov_t.table/colSums(mgm.bin_fxn.cov_t.table))*10^6
mgm.bin_fxn_cpm[1:4,1:4]
colSums(mgm.bin_fxn_cpm)

mgm.bin.cpm<-as.data.frame(t(mgm.bin_fxn_cpm))
mgm.bin.cpm[1:4,1:4]

#write.table(mgm.bin_fxn_cpm, file="data/Metagenomes/Metagenomes/mgm.bin_NoBins_CopiesPerMillion_GeneCounts_2.27.23.txt", sep="\t", quote=F, col.names=T)

#### Create Presence/Absence Table of Functions in MAGs ####
mgm.bin_fxn.cov_table[,-c(1)][1:4,1:4]
mgm.bin_fxn.cov_table.no_lows[1:4,1:4] # << table with lowest summed gene coverage KOs dropped

mgm.bin_fxn.binary<-counts_to_binary(mgm.bin_fxn.cov_table.no_lows) # custom function to convert all counts to binary (presence/absence)
# sanity check that function worked below
mgm.bin_fxn.binary[1:5,1:5]
mgm.bin_fxn.cov_table.no_lows[1:5,1:5]

# #### Create MR-transformed Coverage Table w/ NAs for Absent Functions ####
#
# # are NA table and mgm.bin.mr in the same order?
# rownames(mgm.bin_fxn.binary)
# rownames(mgm.bin.mr)
#
# colnames(mgm.bin_fxn.binary)
# colnames(mgm.bin.mr)
#
# colnames(mgm.bin_fxn.binary)[which(colnames(mgm.bin_fxn.binary) %in% colnames(mgm.bin.mr) == FALSE)] #
# colnames(mgm.bin.mr)[which(colnames(mgm.bin.mr) %in% colnames(mgm.bin_fxn.binary) == FALSE)]
#
# which(colnames(mgm.bin.mr) %in% colnames(mgm.bin_fxn.binary)=="FALSE")
#
# # then create logical table saying which functions are NA in this ko.cov.na table AND have a negative coverage value in mgm.bin.mr
# # TRUE means they are absent, FALSE means they are present
# NA.fxns <- (mgm.bin_fxn.binary==0)
# NA.fxns[1:5,1:5]
# mgm.bin_fxn.binary[1:5,1:5] #compare to NA.fxns -- are all absent functions (0s) given TRUE, and all present functions (1s) given FALSE
# mgm.bin.mr[1:5,1:5] # remember, MR transforms all elements of df including those with 0s, so this is why we are going through this process of labeling absent KOs with NAs
#
# # create data frame that will be MR transformed sum coverages + NA values
# mgm.bin.mr.na<-mgm.bin.mr
# mgm.bin.mr.na[NA.fxns] <- NA # could be because we dropped functions that are absent?
# #mgm.bin.mr.na[NA.fxns == TRUE]<- NA
#
# # did this work?
# mgm.bin.mr[1:4,1:4] # original MR transformed summed KO coverages
# mgm.bin.mr.na[1:4,1:4] # new MR transformed df where absent KOs given NA
# mgm.bin_fxn.binary[1:4,1:4] # presence/absence matrix of KOs
# NA.fxns[1:4,1:4] # present KOs == false, absent KOs == true
#
# # create sample ID column
# mgm.bin.mr.na$Bin_ID<-rownames(mgm.bin.mr.na)
# # * use mgm.bin.mr.na for heatmaps of functions and coverage!
# # to get rid of Bin_ID column later for mgm.bin.mr.na, use the following code
# ## mgm.bin.mr.na[,!names(mgm.bin.mr.na) %in% c("Bin_ID")]
#
# #### Create CLR-transformed Coverage Table w/ NAs for Absent Functions ####
#
# # are NA table and mgm.bin.clr in the same order?
# rownames(mgm.bin_fxn.binary)
# rownames(mgm.bin.clr)
#
# colnames(mgm.bin_fxn.binary)
# colnames(mgm.bin.clr)
#
# colnames(mgm.bin_fxn.binary)[which(colnames(mgm.bin_fxn.binary) %in% colnames(mgm.bin.clr) == FALSE)]
# colnames(mgm.bin.clr)[which(colnames(mgm.bin.clr) %in% colnames(mgm.bin_fxn.binary) == FALSE)]
#
# which(colnames(mgm.bin.clr) %in% colnames(mgm.bin_fxn.binary)=="FALSE")
#
# # then create logical table saying which functions are NA in this ko.cov.na table AND have a negative coverage value in mgm.bin.clr
# # TRUE means they are absent, FALSE means they are present
# NA.fxns <- (mgm.bin_fxn.binary==0)
# NA.fxns[1:5,1:5]
# mgm.bin_fxn.binary[1:5,1:5] #compare to NA.fxns -- are all absent functions (0s) given TRUE, and all present functions (1s) given FALSE
# mgm.bin.clr[1:5,1:5] # remember, CLR transforms all elements of df including those with 0s, so this is why we are going through this process of labeling absent KOs with NAs
#
# # create data frame that will be CLR transformed sum coverages + NA values
# mgm.bin.clr.na<-mgm.bin.clr
# mgm.bin.clr.na[NA.fxns] <- NA # could be because we dropped functions that are absent?
# #mgm.bin.clr.na[NA.fxns == TRUE]<- NA
#
# # did this work?
# mgm.bin.clr[1:4,1:4] # original CLR transformed summed KO coverages
# mgm.bin.clr.na[1:4,1:4] # new CLR transformed df where absent KOs given NA
# mgm.bin_fxn.binary[1:4,1:4] # presence/absence matrix of KOs
# NA.fxns[1:4,1:4] # present KOs == false, absent KOs == true
#
# # create sample ID column
# mgm.bin.clr.na$Bin_ID<-rownames(mgm.bin.clr.na)
# # * use mgm.bin.clr.na for heatmaps of functions and coverage!
# # to get rid of Bin_ID column later for mgm.bin.clr.na, use the following code
# ## mgm.bin.clr.na[,!names(mgm.bin.clr.na) %in% c("Bin_ID")]

#### Compare Sequencing Depth Across Samples ####

# median ratio normalizaztion
mgm.bin.mr[1:4,(ncol(mgm.bin.mr)-4):(ncol(mgm.bin.mr))]
total_mr_counts<-rowSums(mgm.bin.mr[,-(ncol(mgm.bin.mr))])

# copies per million
mgm.bin_fxn_cpm[1:4,1:4]
total_cpm_counts<-colSums(mgm.bin_fxn_cpm)

# centered log ratio transformation
mgm.bin.clr[1:4,1:4]
total_clr_counts<-rowSums(mgm.bin.clr)

# total raw counts
total_counts<-rowSums(mgm.bin_gene.cov_table[,-1])
#total_counts %>% barplot

# variance stabilizing transformation
total_vst_counts<-rowSums(mgm.bin.vst)

par(mfrow=c(1,4)) # to plot the three box plots next to each other (1 row, 3 columns)
total_counts %>% barplot(main = "Total Counts per Sample")
total_mr_counts  %>% barplot(main = "Total Median-Ratio Transformed Counts per Sample")
total_vst_counts %>% barplot(main = "Total Variance-Stabilized Transformed Counts per Sample")
total_cpm_counts  %>% barplot(main = "Total Copies per Million (CPM) per Sample")

dev.off()

### Pull out traits of interest ####
# create unique list of KO ID and functions (did this above, pasting the same code below)
# check for duplicates to make sure each KO_ID has a unique function assignment
# create unique list of GeneIDs & KO IDs
gene_KOs<-unique(data.frame(Gene_ID=mgm.bin_genes.cov.clean$Gene_ID, KO_ID=mgm.bin_genes.cov.clean$KO_ID,KO_Function=mgm.bin_genes.cov.clean$KO_Function))
dim(gene_KOs)
head(gene_KOs)
which(is.na(gene_KOs)) # 0 NAs found

# how many unique KOs are there total, and how many genes are assigned to each KO?
KO_uniq_num_occur <- data.frame(table(gene_KOs$KO_ID))
KO_uniq_num_occur[KO_uniq_num_occur$Freq > 1,]

KOs_unique<-unique(data.frame(KO_ID=mgm.bin_genes.cov.clean$KO_ID,KO_Function=mgm.bin_genes.cov.clean$KO_Function))
dim(KOs_unique)
# ^ # of unique KO IDs is the same # as if you also include functions with KO IDs, meaning no KO IDs have been assigned to multiple functions
KOs_unique[1:10,]

n_occur <- data.frame(table(KOs_unique$KO_ID)) # see how many duplicates there are of KO IDs, compare duplicates
n_occur[n_occur$Freq > 1,] # what traits appear more than once? should be none!

## pull out functions of interest
sulfur.bin.fxns<-KOs_unique[which(KOs_unique$KO_ID %in% sulf.kegg$KO_ID),]
nitro.bin.fxns<-KOs_unique[which(KOs_unique$KO_ID %in% nitro.kegg$KO_ID),]
carb.bin.fxns<-KOs_unique[which(KOs_unique$KO_ID %in% carb.kegg$KO_ID),]
All_GOI.bin.fxns<-KOs_unique[which(KOs_unique$KO_ID %in% all_goi.kegg$KO_ID),]
osmo.bin.fxns<-KOs_unique[which(KOs_unique$KO_ID %in% osmo.kegg$KO_ID),]
selen.bin.fxns<-KOs_unique[which(KOs_unique$KO_ID %in% selen.kegg$KO_ID),]
arsen.bin.fxns<-KOs_unique[which(KOs_unique$KO_ID %in% arsen.kegg$KO_ID),]
tempshock.bin.fxns<-KOs_unique[which(KOs_unique$KO_ID %in% tempshock.kegg$KO_ID),]
metal.bin.fxns<-KOs_unique[which(KOs_unique$KO_ID %in% metal.re.kegg$KO_ID),]
photo.bin.fxns<-KOs_unique[which(KOs_unique$KO_ID %in% photo.kegg$KO_ID),]
aero.rep.bin.fxns<-KOs_unique[which(KOs_unique$KO_ID %in% aero.rep.kegg$KO_ID),]
lps.bin.fxns<-KOs_unique[which(KOs_unique$KO_ID %in% lps.kegg$KO_ID),]
UVrep.bin.fxns<-KOs_unique[which(KOs_unique$KO_ID %in% UV.DNA.rep.kegg$KO_ID),]
spor.bin.fxns<-KOs_unique[which(KOs_unique$KO_ID %in% spor.kegg$KO_ID),]
QuorSens.bin.fxns<-KOs_unique[which(KOs_unique$KO_ID %in% QuorSens.kegg$KO_ID),]

### Merge Metadata & MR Normalized, Summed KO Coverages Together ####
mgm.bin.mr[1:4,1:4]
head(all_bin_meta)
head(KOs_unique)

# melt data with normalized feature counts to merge with all traits & traits of interest
#mgm.bin.mr$Bin_ID<-rownames(mgm.bin.mr)
mgm.bin_mr_melt<-melt(mgm.bin.mr, by="Bin_ID")
head(mgm.bin_mr_melt)
colnames(mgm.bin_mr_melt)[which(names(mgm.bin_mr_melt) == "variable")] <- "KO_ID"
colnames(mgm.bin_mr_melt)[which(names(mgm.bin_mr_melt) == "value")] <- "MRSumCovPerKO"
mgm.bin_mr_melt[1:4,]

# Merge mr-transformed, summed gene coverage per KO w/ KO function info
mgm.bin.mr.all<-as.data.frame(merge(KOs_unique, mgm.bin_mr_melt, by=c("KO_ID"),allow.cartesian = TRUE))
head(mgm.bin.mr.all)

NA %in% mgm.bin.mr.all$KO_ID
dim(mgm.bin.mr.all)

mgm.bin.mr.all<-mgm.bin.mr.all[!is.na(mgm.bin.mr.all$KO_ID),] # only looking at functions we have KO IDs for
head(mgm.bin.mr.all)
dim(mgm.bin.mr.all)

mgm.bin_mr_melt[1:4,] #sanity check

# Merge median-ratio normalized, summed gene coverage per KO w/ KO functions of interest
mgm.bin.mr.sulf<-merge(sulfur.bin.fxns, mgm.bin_mr_melt, by=c("KO_ID"),allow.cartesian = TRUE)
head(mgm.bin.mr.sulf)

mgm.bin.mr.ars<-merge(arsen.bin.fxns, mgm.bin_mr_melt, by=c("KO_ID"),allow.cartesian = TRUE)
head(mgm.bin.mr.ars)

mgm.bin.mr.osmo<-merge(osmo.bin.fxns, mgm.bin_mr_melt, by=c("KO_ID"),allow.cartesian = TRUE)
head(mgm.bin.mr.osmo)

#save.image("SSD_MGM_Bins_analysis_data.Rdata")


### Merge Metadata & CLR Transformed, Summed KO Coverages Together ####
mgm.bin.clr[1:4,1:4]
head(all_mgm_meta)
head(KOs_unique)

# melt data with normalized feature counts to merge with all traits & traits of interest
mgm.bin.clr$Bin_ID<-rownames(mgm.bin.clr)
mgm.bin_clr_melt<-melt(mgm.bin.clr, by="Bin_ID")
head(mgm.bin_clr_melt)
colnames(mgm.bin_clr_melt)[which(names(mgm.bin_clr_melt) == "variable")] <- "KO_ID"
colnames(mgm.bin_clr_melt)[which(names(mgm.bin_clr_melt) == "value")] <- "CLRSumCovPerKO"
mgm.bin_clr_melt[1:4,]

# Merge CLR-transformed, summed gene coverage per KO w/ KO function info
mgm.bin.clr.all<-as.data.frame(merge(KOs_unique, mgm.bin_clr_melt, by=c("KO_ID"),allow.cartesian = TRUE))
head(mgm.bin.clr.all)

NA %in% mgm.bin.clr.all$KO_ID

mgm.bin.clr.all<-mgm.bin.clr.all[!is.na(mgm.bin.clr.all$KO_ID),] # only looking at functions we have KO IDs for
head(mgm.bin.clr.all)
#NA %in% mgm.bin.clr.all

mgm.bin_clr_melt[1:4,] #sanity check

# Merge CLR-transformed, summed gene coverage per KO w/ KO functions of interest
mgm.bin.clr.sulf<-merge(sulfur.bin.fxns, mgm.bin_clr_melt, by=c("KO_ID"),allow.cartesian = TRUE)
head(mgm.bin.clr.sulf)

mgm.bin.clr.ars<-merge(arsen.bin.fxns, mgm.bin_clr_melt, by=c("KO_ID"),allow.cartesian = TRUE)
head(mgm.bin.clr.ars)

mgm.bin.clr.osmo<-merge(osmo.bin.fxns, mgm.bin_clr_melt, by=c("KO_ID"),allow.cartesian = TRUE)
head(mgm.bin.clr.osmo)


#save.image("SSD_MGM_Bins_analysis_data.Rdata")

### Export Global Env for Other Scripts ####
save.image("data/Metagenomes/SSD_MGM_Bins_analysis_data.Rdata")
# ^ includes all data combined in object bac.dat.all, ASV table (samples are rows, ASVs are columns), SurfTypFreq_mgm, and an ASV count table (where ASVs are rows, not columns)

# NOTES for downstream analysis...
# DESeq2 creators say that for visualization and clustering, they recommend using transformed counts like VST or regularized logarithm (rlog)
# previously I have used CLR for visualization + stats (SSW project)
# not sure why they do not suggestion median ratio normalization? maybe because it normalizes the data which is not what transformations necessarily do...
# could use normalized data for all things (models, ordinations, PERMANOVA)

# Version Information
sessionInfo()
