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
  library(data.table)
  library(ape)
  library(apeglm)
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
  library(decontam)
  library(batchtools)
})

#load("SSD_mgm_contigs_analysis_data.Rdata") # load Rdata to global env

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

##save.image("SSD_mgm_contigs_analysis_data.Rdata") # save global env to Rdata file

## Notes:
# code & info came from :
## https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#quick-start
## https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/07_practical.pdf
## https://www.reneshbedre.com/blog/deseq2.html
## https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html

#### Import Contig Coverage Data ####
mgm_fxns.cov<-fread(file = 'SSD_Samples_NoBins_Gene_Coverages_1.2.24.txt', sep='\t',header = TRUE)
dim(mgm_fxns.cov)
head(mgm_fxns.cov)

# update SampleID names so they match w/ metadata
# mgm_fxns.cov$SampleID<-gsub("_",".", mgm_fxns.cov$SampleID)
# head(mgm_fxns.cov)
# mgm_fxns.cov$SampleID<-gsub("m\\..*","m",mgm_fxns.cov$SampleID)
# mgm_fxns.cov$PlotID<-gsub("^SSD.","",mgm_fxns.cov$SampleID)
# head(mgm_fxns.cov)
#
# # Drop June 2021 samples
# unique(mgm_fxns.cov$SampleID)
# mgm_fxns.cov<-mgm_fxns.cov[!grepl("6.15.21", mgm_fxns.cov$SampleID),]
# unique(mgm_fxns.cov$SampleID) #sanity check

# count occurrences of all traits across all mgms
n_occur <- data.frame(table(mgm_fxns.cov$KO_Function))
n_occur[n_occur$Freq > 1,]

# Divide gene counts by gene length to account for sample differences in assembly
## we do this because we did not co-assemble contigs, so each KO assignments across samples may not come from genes with the same length
## need to account for gene length here first (since KO assignments can come from genes of different lengths)
mgm_fxns.cov$CovPerGene<-mgm_fxns.cov$ReadsPerGene/mgm_fxns.cov$Gene_Length
head(mgm_fxns.cov)

# Drop rows where we do not have a KO ID for the gene
mgm_fxns.cov.clean<-na.omit(mgm_fxns.cov, cols=c("KO_ID"))
head(mgm_fxns.cov.clean)

# create list of GeneIDs & KO IDs
gene_KOs<-unique(data.frame(Gene_ID=mgm_fxns.cov.clean$Gene_ID, KO_ID=mgm_fxns.cov.clean$KO_ID,KO_Function=mgm_fxns.cov.clean$KO_Function))
dim(gene_KOs)
head(gene_KOs)

# how many unique KOs are there total, and how many genes are assigned to each KO?
n_occur <- data.frame(table(gene_KOs$KO_ID))
n_occur[n_occur$Freq > 1,]

KOs_unique<-unique(data.frame(KO_ID=mgm_fxns.cov.clean$KO_ID,KO_Function=mgm_fxns.cov.clean$KO_Function))
dim(KOs_unique)
# ^ # of unique KO IDs is the same # as if you also include functions with KO IDs, meaning no KO IDs have been assigned to multiple functions

# create Sample ID x Gene ID count table, using coverage aka reads per gene that were divided by gene length
mgm_fxn.cov_table1<-as.data.frame(dcast(mgm_fxns.cov.clean, SampleID~Gene_ID, value.var="CovPerGene"))
mgm_fxn.cov_table1[,1:4] # sanity check
rownames(mgm_fxn.cov_table1)<-mgm_fxn.cov_table1$SampleID

# check for NAs in data.table
which(is.na(mgm_fxn.cov_table1))

# convert all NAs to 0; will take a while
# mgm_fxn.cov_table<-remove_na(mgm_fxn.cov_table1)
# mgm_fxn.cov_table[,1:6] # sanity check

#### Create Sample x KO Matrix of Coverage Data ####
mgm_fxn.cov_table1[,1:4] # sanity check

mgm.gene.cov.m<-melt(mgm_fxn.cov_table1, by="SampleID")
colnames(mgm.gene.cov.m)[which(names(mgm.gene.cov.m) == "variable")] <- "Gene_ID"
colnames(mgm.gene.cov.m)[which(names(mgm.gene.cov.m) == "value")] <- "CovPerGene"
head(mgm.gene.cov.m)

# merge melted coverage data with KO IDs
head(gene_KOs)
mgm.fxn.cov<-merge(mgm.gene.cov.m, gene_KOs,by="Gene_ID")
head(mgm.fxn.cov) # confirmed that multiple genes are assigned to the same KOs

# recreate table so that it is Sample ID x KO ID instead of Gene ID
# create Sample ID x Gene ID count table, using coverage aka reads per gene that were divided by gene length
mgm_fxn.cov_table<-as.data.frame(dcast(mgm_fxns.cov, SampleID~KO_ID, value.var="CovPerGene",fun.aggregate=sum))
mgm_fxn.cov_table[,1:4] # sanity check
rownames(mgm_fxn.cov_table)<-mgm_fxn.cov_table1$SampleID
## Note: multiple genes have same KO assignment, so their coverage values are all over the place -- should scale these down to make them more legible?
## ^ CPM transformation as Mike suggested?? Wondering what best transformation is for stats?

# drop KOs with low coverage
mgm.fxn.nolows<-mgm_fxn.cov_table[,which(colSums(mgm_fxn.cov_table[,-1])>=5)] # remove functions with less than 5 total counts across mgms
#mgm_fxn.binary<-counts_to_binary(mgm_fxn.cov_table[,-1]) # custom function to convert all counts to binary (presence/absence)
mgm.fxn.nolows[1:4,1:4] # sanity check

save.image("SSD_mgm_analysis_data.Rdata")

#### Update Metadata ####
# upload UNSCALED geochem data from Lyons lab
# ^^ scale env variable data in respective scripts

mgm_meta<-as.data.frame(read_xlsx("data/Metagenomes/SSD_Lyons_Aronson_Metadata_MGM_All.xlsx", sheet="Metagenomes_Metadata"))
head(mgm_meta)
mgm_meta$SampleID<-gsub("m\\..*","m",mgm_meta$SampleID) # remove rep # after m
mgm_meta$PlotID<-gsub("^SSD.","",mgm_meta$SampleID)
head(mgm_meta)

# drop data from June 2021
mgm_june<-mgm_meta # separate DF containing June 2021 data
mgm_meta<-mgm_meta[mgm_meta$SampleMonth!="June",]
"June" %in% mgm_meta$SampleMonth # confirm there are no Junes left
unique(mgm_meta$SampleMonth) # another sanity check to be safe

# create factor levels for certain groups
unique(mgm_meta$Depth_m)
mgm_meta$Depth_m<-factor(mgm_meta$Depth_m, levels=c("0","5","10"))

unique(mgm_meta$SampleMonth)
mgm_meta$SampleMonth<-factor(mgm_meta$SampleMonth, levels=c("August","December","April"))

mgm_meta$SampDate<-interaction(mgm_meta$SampleMonth, mgm_meta$SampleYear)
head(mgm_meta)
mgm_meta$SampDate<-factor(mgm_meta$SampDate, levels=c("August.2021", "December.2021", "April.2022"))
head(mgm_meta)

# Create color palettes to be used for gradient colors

cold2warm1<-get_palette(paste0("#",c("BF1B0B","66ADE5","252A52")),k=3)
names(cold2warm1) <- levels(mgm_meta$Depth_m)
cold2warm1

fair_cols <- paste0("#",c("252A52", "66ADE5", "FFC465","BF1B0B"))
names(fair_cols) <- letters[1:4]
fair_ramp <- scales::colour_ramp(fair_cols)
fair_sat <- saturation(fair_ramp, 1)

# Add colors for specific variables

colorset1 = melt(c(August.2021="#ef781c",December.2021="#03045e",April.2022="#059c3f"))

colorset1$SampDate<-rownames(colorset1)
colorset1
colnames(colorset1)[which(names(colorset1) == "value")] <- "SampDate_Color"
colorset1

mgm_meta<-merge(mgm_meta, colorset1, by="SampDate")
head(mgm_meta)
mgm_meta$SampDate_Color <- as.character(mgm_meta$SampDate_Color)
rownames(mgm_meta)<-mgm_meta$SampleID
# #save.image("data/Metagenomes/SSD_mgm_analysis.Rdata")

# create color palette for binary heat maps of function data
binary.cols<-c("1"="red","0"="aliceblue")

#### Scale Chem Data in Metadata ####
head(mgm_meta)
#meta_scaled<-subset(mgm_meta, select=-c(Salinity_ppt)) # drop salinity, for now
meta_scaled<-mgm_meta

head(meta_scaled)
meta_scaled[,8:17]<-as.data.frame(scale(meta_scaled[,8:17], center=TRUE, scale=TRUE)) #centering before scaling
meta_scaled

# create numeric version of depth column for later analyses
meta_scaled$Depth.num<-as.numeric(as.character(meta_scaled$Depth_m))

#save.image("SSD_mgm_contigs_analysis_data.Rdata")

#### Import Gene Info from KEGG ####

all_goi.kegg<-read.table("/bigdata/aronsonlab/shared/SaltonSea/Metagenomes/SeqCenter_3.30.2023/MGM_Analyses/KO_KEGG_Files/Genes_of_Interest_All_KOs_Pathway_Cycle_KEGG.txt", header = TRUE, sep = "\t", dec = ".")
carb.kegg<-read.table("/bigdata/aronsonlab/shared/SaltonSea/Metagenomes/SeqCenter_3.30.2023/MGM_Analyses/KO_KEGG_Files/CarbonFixation_KOs_KEGG.txt", header = TRUE, sep = "\t", dec = ".")
sulf.kegg<-read.table("/bigdata/aronsonlab/shared/SaltonSea/Metagenomes/SeqCenter_3.30.2023/MGM_Analyses/KO_KEGG_Files/Sulfur_KOs_KEGG.txt", header = TRUE, sep = "\t", dec = ".")
nitro.kegg<-read.table("/bigdata/aronsonlab/shared/SaltonSea/Metagenomes/SeqCenter_3.30.2023/MGM_Analyses/KO_KEGG_Files/N_KOs_Pathway_Cycle_KEGG.txt", header = TRUE, sep = "\t", dec = ".")
arsen.kegg<-read.table("/bigdata/aronsonlab/shared/SaltonSea/Metagenomes/SeqCenter_3.30.2023/MGM_Analyses/KO_KEGG_Files/Arsenic_KOs_KEGG.txt", header = TRUE, sep = "\t", dec = ".")
heatshock.kegg<-read.table("/bigdata/aronsonlab/shared/SaltonSea/Metagenomes/SeqCenter_3.30.2023/MGM_Analyses/KO_KEGG_Files/HeatShock_KOs_KEGG.txt", header = TRUE, sep = "\t", dec = ".")
osmo.kegg<-read.table("/bigdata/aronsonlab/shared/SaltonSea/Metagenomes/SeqCenter_3.30.2023/MGM_Analyses/KO_KEGG_Files/Osmoprot_KOs_KEGG.txt", header = TRUE, sep = "\t", dec = ".")
selen.kegg<-read.table("/bigdata/aronsonlab/shared/SaltonSea/Metagenomes/SeqCenter_3.30.2023/MGM_Analyses/KO_KEGG_Files/Selenium_KOs_KEGG.txt", header = TRUE, sep = "\t", dec = ".")
metal.re.kegg<-read.table("/bigdata/aronsonlab/shared/SaltonSea/Metagenomes/SeqCenter_3.30.2023/MGM_Analyses/KO_KEGG_Files/MetalResistance_KOs_KEGG.txt", header = TRUE, sep = "\t", dec = ".")
photo.kegg<-read.table("/bigdata/aronsonlab/shared/SaltonSea/Metagenomes/SeqCenter_3.30.2023/MGM_Analyses/KO_KEGG_Files/Photo_KO_KEGG.txt", header = TRUE, sep = "\t", dec = ".")
aero.kegg<-read.table("/bigdata/aronsonlab/shared/SaltonSea/Metagenomes/SeqCenter_3.30.2023/MGM_Analyses/KO_KEGG_Files/OxidativePhosphorylation_KO_KEGG.txt", header = TRUE, sep = "\t", dec = ".")
lps.kegg<-read.table("/bigdata/aronsonlab/shared/SaltonSea/Metagenomes/SeqCenter_3.30.2023/MGM_Analyses/KO_KEGG_Files/LPS_Biosynthesis_KOs_KEGG.txt.txt", header = TRUE, sep = "\t", dec = ".")

#### Check Gene Distribution in MGMs ####
mgm_fxn.cov_table[,1:4] # sanity check
mgm_cov_t<-as.data.frame(t(mgm_fxn.cov_table[,-1]))
class(mgm_cov_t$SSD.12.22.21.5m)
mgm_cov_t <- as.data.frame(sapply(mgm_cov_t, as.numeric)) # convert integer to numeric across df
class(mgm_cov_t$SSD.12.22.21.5m) # sanity check

#descdist(mgm_cov_t$SSD.4.13.22.5m, discrete = TRUE)
#descdist(mgm_cov_t$SSD.12.22.21.5m, discrete = TRUE)
#descdist(mgm_cov_t$SSD.8.24.21.0m, discrete = TRUE)

# check distribution of gene # x gene cov with two samples
ggplot(mgm_cov_t) +
  geom_histogram(aes(x = SSD.4.13.22.5m), stat = "bin", bins = 200) +
  xlab("Gene coverage") +
  ylab("Number of genes")

ggplot(mgm_cov_t) +
  geom_histogram(aes(x = SSD.12.22.21.0m), stat = "bin", bins = 200) +
  xlab("Gene coverage") +
  ylab("Number of genes")
## data is not normally distributed

# histogram of gene cov
hist(mgm_cov_t$SSD.12.22.21.5m, col="blue")
# visualize Q-Q plot for
qqnorm(mgm_cov_t$SSD.12.22.21.5m, pch = 1, frame = FALSE)
qqline(mgm_cov_t$SSD.12.22.21.5m, col = "red", lwd = 2)

# histogram of gene cov
hist(mgm_cov_t$SSD.4.13.22.5m, col="blue")
# visualize Q-Q plot
qqnorm(mgm_cov_t$SSD.4.13.22.5m, pch = 1, frame = FALSE)
qqline(mgm_cov_t$SSD.4.13.22.5m, col = "red", lwd = 2)

# compare mean cov vs variance of the cov across samples
mean_cov <- apply(mgm_cov_t, 1, mean)        #The second argument '1' of 'apply' function indicates the function being applied to rows. Use '2' if applied to columns
variance_cov <- apply(mgm_cov_t, 1, var)
df <- data.frame(mean_cov, variance_cov)

ggplot(df) +
  geom_point(aes(x=mean_cov, y=variance_cov)) +
  scale_y_log10(limits = c(1,1e9)) +
  scale_x_log10(limits = c(1,1e9)) +
  geom_abline(intercept = 0, slope = 1, color="red")

#save.image("data/Metagenomes/SSD_mgm_analysis.Rdata")

#### Sum Gene Coverage by KO ID in Contigs Before Transformations ####
# mgm_fxns.cov[1:4,]
# mgm_fxns.cov_noNA<-as.data.frame(mgm_fxns.cov[!is.na(mgm_fxns.cov$KO_ID),]) # drop genes with KOs given NA as assignment (aka no KO ID assigned at all)
#
# ko.cov.sum_table<-as.data.frame(dcast(mgm_fxns.cov_noNA, SampleID~KO_ID, value.var="CovPerGene", fun.aggregate=sum)) ###
# ko.cov.sum_table[1:4,1:4]
# rownames(ko.cov.sum_table)<-ko.cov.sum_table$SampleID
# ko.cov.sum_table[1:4,1:4]
#
# # check rownames of summed transformed feature coverage data & metadata
# rownames(ko.cov.sum_table) %in% rownames(meta_scaled)
#
# ko.cov.sum_table[1:4,1:4] # contains the sum of coverages per gene per KO -- featureCounts was normalized by gene length across samples first to get coverage, then summed up per KO ID

#### Prepare Contig Feature Count Data for Normalization w/ DESeq2 ####
# make sure count data & mgm_meta are in the same order
ko.cov.sum_table[1:4,1:5]

# scale function coverage by 100 so that you can round to create DESeq2 object
## DESeq2 requires whole integers
mgm_scale.sum.cov_table<-ko.cov.sum_table
mgm_fxn.cov_table[,2:6]*100 # sanity check that this method will work
mgm_scale.sum.cov_table[,-1]<-mgm_scale.sum.cov_table[,-1]*100

mgm_scale.sum.cov_table[,1:10]

rownames(mgm_scale.sum.cov_table) # sanity check
rownames(mgm_meta)

rownames(mgm_meta) %in% rownames(mgm_scale.sum.cov_table)

mgm_meta=mgm_meta[rownames(mgm_scale.sum.cov_table),] ## reorder mgm_meta to have same rows as original OTU table

# convert all NAs to 0; will take a while
mgm_scale.sum.cov_table[1:6,1:6]

#convert data frame to numeric
#mgm_scale.sum.cov_table[,-1] <- as.data.frame(sapply(mgm_scale.sum.cov_table[,-1], as.numeric))

scale.sum.cov_matrix<-as.matrix(mgm_scale.sum.cov_table[,!names(mgm_scale.sum.cov_table) %in% c("SampleID")]) # convert count table into matrix & exclude column called SampleID
scale.sum.cov_matrix2<-t(scale.sum.cov_matrix) # transpose matrix so KO_IDs are rows, samples are columns
# ^ will be used in DESeq2 functions
rownames(mgm_meta) %in% colnames(scale.sum.cov_matrix2) # check if rownames in metadata (SampleID) match column names in count data
dim(mgm_meta)
dim(scale.sum.cov_matrix2)

# create the DESeq DataSet object for DGE analysis
# DESeq2 needs whole # data, so need raw read counts, NOT coverage for these tables...questionable
# mgm_dds has the scaled coverage that was calculated by dividing reads from featureCounts by gene lengths
mgm_dds<-DESeqDataSetFromMatrix(countData=round(scale.sum.cov_matrix2), colData = mgm_meta, design = ~ 1)

# design = ~ 1 means no design
head(counts(mgm_dds))
colSums(counts(mgm_dds)) %>% barplot

# add rownames from KO IDs to DESeq Data set object
featureData <- data.frame(KO_ID=rownames(mgm_dds)) # create object of KO IDs
mcols(mgm_dds) <- DataFrame(mcols(mgm_dds), featureData) # add KO IDs to rownames of DeSEQ dataset object
mcols(mgm_dds)

# drop KO IDs with lowest counts
keep <- rowSums(counts(mgm_dds)) >= 5
mgm_dds <- mgm_dds[keep,]

# Estimate size factor - aka normalization factor, divide all read counts by each size factor per sample
mgm_dds <- estimateSizeFactors(mgm_dds,type="ratio")
## To calculate size factor in DESeq2:
# calculates geometric mean of each gene in each sample and uses this as a pseudoreference
# calculates ratio of each sample by dividing each gene count by it's pseudoreference in each sample
# The median value of all ratios for a given sample is taken as the normalization factor (size factor)
# The differentially expressed genes should not affect the median value
# median of ratios method makes the assumption that not ALL genes are differentially expressed; therefore, the normalization factors should account for sequencing depth and RNA composition of the sample
## (large outlier genes will not represent the median ratio values)
# more here: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html

sizeFactors(mgm_dds)

plot(sizeFactors(mgm_dds), colSums(counts(mgm_dds)))

# Does sequencing depth influence normalization?
par(mfrow=c(1,2)) # to plot the two box plots next to each other
boxplot(log2(counts(mgm_dds)), notch=TRUE,
        main = "Non-normalized read counts\n(log-transformed)",
        ylab="read counts")
boxplot(log2(counts(mgm_dds, normalize= TRUE)), notch=TRUE,
        main = "Size-factor-normalized read counts\n(log-transformed)",
        ylab="read counts")
dev.off()

#save.image("data/Metagenomes/SSD_mgm_analysis.Rdata")

#### Median-Ratio Normalized - Gene in Contigs ####
#
sizeFactors(mgm_dds) # will be used to normalize counts
## to normalize counts, we divide each raw count value in a given sample by that sample’s normalization factor (aka size factor) to generate normalized count values.
### This is performed for all count values (every gene in every sample)
mgm_fxn_mr <- counts(mgm_dds, normalized=TRUE) ## median-ratio transformation is how DESeq2 Normalizes counts!
#write.table(mgm_fxn_counts.norm, file="data/Metagenomes/MGM_NoBins_MedianRatio_GeneCounts_2.27.23.txt", sep="\t", quote=F, col.names=T)

mgm.mr<-as.data.frame(t(mgm_fxn_mr))
mgm.mr$SampleID<-rownames(mgm.mr)
mgm.mr[1:4,1:4]

#### Variance Stabilizing Transformation - Gene in Contigs ####

# code below is also in previous section - if you ran previous section, skip to VST step
ko.cov.sum_table[1:4,1:4]

scale.sum.cov_matrix<-as.matrix(ko.cov.sum_table[,!names(ko.cov.sum_table) %in% c("SampleID")]) # convert count table into matrix
scale.sum.cov_matrix2<-t(scale.sum.cov_matrix) # transpose matrix so KO_IDs are rows, samples are columns
# ^ will be used in DESeq2 functions
rownames(mgm_meta) %in% colnames(scale.sum.cov_matrix2) # check if rownames in metadata (SampleID) match column names in count data
# sanity check
dim(mgm_meta)
dim(scale.sum.cov_matrix2)

# you should be able to use matrix or DESeq2 object for this next function, but matrix was not working?
# variance stabilizing transformation
mgm_fxn_vst <- varianceStabilizingTransformation(mgm_dds, blind = TRUE, fitType = "parametric")
assay(mgm_fxn_vst) #see output of VST

mgm.vst<-assay(mgm_fxn_vst)
#total_fxn_vst<-colSums(mgm_fxn_vst)

#### Centered Log Ratio Transformation - Gene in Contigs ####
ko.cov.sum_table[1:4,1:4]
# ^ table contains gene coverage, Sample IDs as rows & genes as columns

# df must have rownames are SampleIDs, columns are ASV IDs for vegan functions below\
mgm.clr<-decostand(ko.cov.sum_table[,-1],method = "clr",pseudocount=1) #CLR transformation
mgm.clr[1:4,1:4]
# NOTE: CLR transformation does not treat all 0s equally, it has to do with the pseudocount that's added before transformation
# The method can operate only with positive data; a common way to deal with zeroes is to add pseudocount, either by adding it manually to the input data, or by using the argument pseudocount as in decostand(x, method = "clr", pseudocount = 1).
# Adding pseudocount will inevitably introduce some bias; see the rclr method for one available solution

# check rownames of CLR transformed ASV data & metadata
rownames(mgm.clr) %in% rownames(meta_scaled)

#### Robust Centered Log Ratio Transformation - Gene in Contigs ####
ko.cov.sum_table[1:4,1:4]
# ^ table contains gene coverage, Sample IDs as rows & genes as columns

# df must have rownames are SampleIDs, columns are ASV IDs for vegan functions below\
mgm.Rclr<-decostand(ko.cov.sum_table[,-1],method = "rclr") #CLR transformation
mgm.Rclr[1:4,1:4]
# NOTE: Robust CLR just excludes 0s and performs CLR transformation without pseudocount
# robust clr ("rclr") is similar to regular clr (see above) but allows data that contains zeroes.
# This method does not use pseudocounts, unlike the standard clr. Robust clr divides the values by geometric mean of the observed features; zero values are kept as zeroes, and not taken into account.
#In high dimensional data, the geometric mean of rclr is a good approximation of the true geometric mean

# check rownames of CLR transformed ASV data & metadata
rownames(mgm.Rclr) %in% rownames(meta_scaled)


#### Copies Per Million Transformation - Gene in Contigs ####
ko.cov.sum_table[1:4,1:4] # sanity check
mgm_fxn.sum.cov_t.table<-as.data.frame(t(as.data.frame(ko.cov.sum_table[,-1])))
mgm_fxn.sum.cov_t.table[1:4,1:4]

mgm_fxn_cpm<-(mgm_fxn.sum.cov_t.table/colSums(mgm_fxn.sum.cov_t.table))*10^6
mgm_fxn_cpm[1:4,1:4]
colSums(mgm_fxn_cpm)

#write.table(mgm_fxn_cpm, file="data/Metagenomes/MGM_NoBins_CopiesPerMillion_GeneCounts_2.27.23.txt", sep="\t", quote=F, col.names=T)

#### Create Presence/Absence Table of Functions in MAGs ####
ko.cov.sum_table[,-c(1)][1:4,1:4]

mgm_fxn.binary<-counts_to_binary(ko.cov.sum_table[,-1]) # custom function to convert all counts to binary (presence/absence)
# sanity check that function worked below
mgm_fxn.binary[1:5,1:5]
ko.cov.sum_table[,-c(1)][1:5,1:5]

#### Create CLR-transformed Coveraged Table w/ NAs for Absent Functions ####

# are NA table and mgm.clr in the same order?
rownames(mgm_fxn.binary)
rownames(mgm.clr)

colnames(mgm_fxn.binary)
colnames(mgm.clr)
colnames(mgm_fxn.binary)[which(colnames(mgm_fxn.binary) %in% colnames(mgm.clr) == FALSE)] # mgm_fxn.binary does not have SampleID column

identical(mgm.clr,mgm_fxn.binary) # SampleID columns are in different places in these two data frames...

# then create logical table saying which functions are NA in this ko.cov.na table AND have a negative coverage value in mgm.clr
# TRUE means they are absent, FALSE means they are present
NA.fxns <- (mgm_fxn.binary==0)

# create data frame that will be CLR transformed sum coverages + NA values
mgm.clr.na<-mgm.clr
mgm.clr.na[NA.fxns] <- NA
#mgm.clr.na[NA.fxns == TRUE]<- NA

# did this work?
mgm.clr[1:4,1:4]
mgm.clr.na[1:4,1:4]
mgm_fxn.binary[1:4,1:4]

# create sample ID column
mgm.clr.na$SampleID<-rownames(mgm.clr.na)
# * use mgm.clr.na for heatmaps of functions and coverage!
# to get rid of sampleID column later for mgm.clr.na, use the following code
## mgm.clr.na[,!names(mgm.clr.na) %in% c("SampleID")]


#### Compare Sequencing Depth Across Samples ####
mgm.mr[1:4,(ncol(mgm.mr)-4):(ncol(mgm.mr))]
total_mr_counts<-rowSums(mgm.mr[,-(ncol(mgm.mr))])

mgm_fxn_cpm[1:4,1:4]
total_cpm_counts<-colSums(mgm_fxn_cpm)

mgm.clr[1:4,1:4]
total_clr_counts<-rowSums(mgm.clr)

total_counts<-rowSums(mgm_fxn.cov_table1[,-1])
#total_counts %>% barplot

total_vst_counts<-colSums(mgm.vst)

par(mfrow=c(1,3)) # to plot the three box plots next to each other (1 row, 3 columns)
#total_counts %>% barplot(main = "Total Counts per Sample")
total_mr_counts  %>% barplot(main = "Total Median-Ratio Transformed Counts per Sample")
total_vst_counts %>% barplot(main = "Total Variance-Stabilized Transformed Counts per Sample")
#total_cpm_counts  %>% barplot(main = "Total Copies per Million (CPM) per Sample")

dev.off()

### Pull out traits of interest ####
# create unique list of KO ID and functions
# check for duplicates to make sure each KO_ID has a unique function assignment
head(mgm_fxns.cov)
NA %in% mgm_fxns.cov$CovPerGene # just to ensure there are no NAs for genes in this df

ko_fxns1<-unique(mgm_fxns.cov[,1:3]) # first subset out data based on unique KO functions

n_occur <- data.frame(table(ko_fxns1$KO_ID)) # see how many duplicates there are of KO IDs, compare duplicates
n_occur[n_occur$Freq > 1,] # what traits appear more than once?

ko_ID<-unique(data.frame(KO_ID=ko_fxns1$KO_ID)) # get a list of unique KO IDs in data

ko_fxns<-as.data.frame(ko_fxns1[!is.na(ko_fxns1$KO_ID),])  # use unique KO ID list to subset out KO function data
head(ko_fxns)

## pull out functions of interest
#sulfur.fxns<-as.data.frame(ko_fxns[grep("sulf|thio", ko_fxns$KO_Function), ]) # pull out sulfur functions
sulfur.fxns<-ko_fxns[which(ko_fxns$KO_ID %in% sulf.kegg$KO_ID),]
nitro.fxns<-ko_fxns[which(ko_fxns$KO_ID %in% nitro.kegg$KO_ID),]
carb.fxns<-ko_fxns[which(ko_fxns$KO_ID %in% carb.kegg$KO_ID),]
All_GOI.fxns<-ko_fxns[which(ko_fxns$KO_ID %in% all_goi.kegg$KO_ID),]
osmo.fxns<-ko_fxns[which(ko_fxns$KO_ID %in% osmo.kegg$KO_ID),]
selen.fxns<-ko_fxns[which(ko_fxns$KO_ID %in% selen.kegg$KO_ID),]
arsen.fxns<-ko_fxns[which(ko_fxns$KO_ID %in% arsen.kegg$KO_ID),]
HS.fxns<-ko_fxns[which(ko_fxns$KO_ID %in% heatshock.kegg$KO_ID),]
metal.fxns<-ko_fxns[which(ko_fxns$KO_ID %in% metal.re.kegg$KO_ID),]
photo.fxn<-ko_fxns[which(ko_fxns$KO_ID %in% photo.kegg$KO_ID),]
aero.fxn<-ko_fxns[which(ko_fxns$KO_ID %in% aero.kegg$KO_ID),]

### Merge Metadata & Contig Count Data & Taxa Data Together ####
mgm.clr[1:4,1:4]
head(mgm_meta)
head(ko_fxns)

# melt data with normalized feature counts to merge with all traits & traits of interest
mgm.clr$SampleID<-rownames(mgm.clr)
mgm_clr_melt<-melt(mgm.clr, by="SampleID")
head(mgm_clr_melt)
colnames(mgm_clr_melt)[which(names(mgm_clr_melt) == "variable")] <- "KO_ID"
colnames(mgm_clr_melt)[which(names(mgm_clr_melt) == "value")] <- "SumCovPerKO"
mgm_clr_melt[1:4,]

# Merge CLR-transformed, summed gene coverage per KO w/ KO function info
mgm.clr.all<-as.data.frame(merge(ko_fxns, mgm_clr_melt, by=c("KO_ID"),allow.cartesian = TRUE))
head(mgm.clr.all)

NA %in% mgm.clr.all$KO_ID

mgm.clr.all<-mgm.clr.all[!is.na(mgm.clr.all$KO_ID),] # only looking at functions we have KO IDs for
head(mgm.clr.all)
#NA %in% mgm.clr.all

mgm_clr_melt[1:4,] #sanity check

# Merge CLR-transformed, summed gene coverage per KO w/ KO functions of interest
mgm.clr.sulf<-merge(sulfur.fxns, mgm_clr_melt, by=c("KO_ID"),allow.cartesian = TRUE)
head(mgm.clr.sulf)

mgm.clr.ars<-merge(arsen.fxns, mgm_clr_melt, by=c("KO_ID"),allow.cartesian = TRUE)
head(mgm.clr.ars)

mgm.clr.osmo<-merge(osmo.fxns, mgm_clr_melt, by=c("KO_ID"),allow.cartesian = TRUE)
head(mgm.clr.osmo)

#save.image("data/Metagenomes/SSD_mgm_analysis.Rdata")

### Export Global Env for Other Scripts ####
save.image("data/Metagenomes/SSD_mgm_analysis.Rdata")
# ^ includes all data combined in object bac.dat.all, ASV table (samples are rows, ASVs are columns), mgm_meta, and an ASV count table (where ASVs are rows, not columns)
# Version Information
sessionInfo()
