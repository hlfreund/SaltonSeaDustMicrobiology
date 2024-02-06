#### Set WD & Load Libraries ####
getwd() # use setwd("path/to/files") if you are not in the right directory
#setwd("/Volumes/HLF_SSD/Aronson_Lab_Data/Salton_Sea/SaltonSeaDust")
suppressPackageStartupMessages({ # load packages quietly
  library(devtools)
  library(phyloseq)
  library(ggplot2)
  library(vegan)
  library(ggpubr)
  library(lme4)
  #library(scales)
  library(grid)
  library(ape)
  library(plyr)
  library(dplyr)
  library(viridis)
  library(ggbiplot)
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
  library(decontam)
  library(ggvegan)
  library(microbiome)
  library(pairwiseAdonis)
  library(corrplot)
  library(lmtest)
  #library(car) # car is imported w/ rstatix so you do not need to load library individually?
})

#### Load Global Env to Import Count/ASV Tables ####
load("data/Amplicon/SSDust_16S.V3V4_W23_Data_Ready.Rdata") # save global env to Rdata file

#save.image("data/Env_Seqs_All/env.seq_analysis.Rdata") # save global env to Rdata file
b.dust.all[1:6,1:6]
bac.ASV_table[1:4,1:4]
bac.ASV_table[(nrow(bac.ASV_table)-4):(nrow(bac.ASV_table)),(ncol(bac.ASV_table)-4):(ncol(bac.ASV_table))] # last 4 rows & cols
head(dust_meta)

## DO NOT RUN THIS LINE, THIS IS YOUR COLOR REFERENCE!!!!
#(WI="#ef781c",DP="#03045e",BDC="#059c3f")

#### Rarefaction Curves & Species Accumulation Curves ####
# bacteria/archaea

# Species Accumulation Curve
sc2<-specaccum(bac.ASV_table[,-1],"random")
plot(sc2, ci.type="poly", col="darkgreen", lwd=2, ci.lty=0, ci.col="lightgreen")
boxplot(sc2, col="yellow", add=TRUE, pch=20)

# Prep for Rarefaction Curve
rowSums(bac.ASV_table[,-1]) # total # ASVs per sample, excluding SampleID from calculation
sort(colSums(bac.ASV_table[,-1]))

# Create Rarefaction curve
png('figures/AlphaDiversity/SSD_16S_rarecurve.png')
rarecurve(as.matrix(bac.ASV_table[,-1]),col=dust_meta$SampDate_Color, step=1000, label=T,ylab="ASVs")
# to show sampel labels per curve, change label=T
dev.off()

#### Good's Coverage ####

# Good's coverage is defined as 1 - (F1/N) where F1 is the number of singleton OTUs and N is the sum of counts for all OTUs.
# If a sample has a Good's coverage == . 98, this means that 2% of reads in that sample are from OTUs that appear only once in that sample.
# NOTE: Good's Coverage is not the best estimate of coverage because DADA2 removes singletons during process - so we are basing this on singeltons that remain (??)
bac.ASV.melt<-melt(bac.ASV_table,by="SampleID")
colnames(bac.ASV.melt)[which(colnames(bac.ASV.melt) == "variable")] <- "ASV_ID"
colnames(bac.ASV.melt)[which(colnames(bac.ASV.melt) == "value")] <- "Count"

# calculate # of singletons per sample
sings<-data.frame(Singletons=rowSums(bac.ASV_table[,-1]==1),SampleID=bac.ASV_table$SampleID) # how many ASVs have a count of only 1
# no singletons?

# find total counts per sample
totseq<-data.frame(TotalSeqs=rowSums(bac.ASV_table[,-1]),SampleID=bac.ASV_table$SampleID) # total # of counts per sample

# Calculate Good's Coverage
good.cov<-data.frame(GoodsCov=(100*(1-(sings$Singletons/totseq$TotalSeqs))),SampleID=sings$SampleID)

# Merge Total Seqs and Goods DFs together for plotting
good.df<-merge(good.cov,totseq,by="SampleID")

ggplot(good.df, aes(x=TotalSeqs,y=GoodsCov))+geom_point()+xlab("Total Seqs per Sample")+ylab("Good's Coverage (%)")

#### Rarefaction of Raw Counts ####
# this section explores different ways of utilizing rarefaction and rarefaction functions in vegan

# in vegan ROWS need to be SITES/samples; COLUMNS are SPECIES (OTUs, ASVs)
min.rar<-min(rowSums(bac.ASV_table[,-1])) ## seeing min sum of OTUs so we can see what min is for rarefaction
min.rar

bac.ASV.rar1<-rarefy(bac.ASV_table[,-1],min.rar) ## be cognizant of min for rarefaction
bac.ASV.rar1 # rarefy gives the expected species richness in random subsamples of size sample from the community

bac.ASV.rar<-rrarefy(bac.ASV_table[,-1],min.rar)
# rrarefy generates one randomly rarefied community data frame or vector of given sample size
# The random rarefaction is made without replacement so that the variance of rarefied communities is rather related to rarefaction proportion than to the size of the sample
rowSums(bac.ASV.rar) # should all be the same after rarefaction

# Calculating Average Shannon Diversity w/ repeated rarefaction with rrarefy()
sdiv.rarefy<-function(df,min){
  # rownames must be samples, ASVs must by columns
  # min = min(rowSums(ASV table)) aka minimum sum of ASVs in each sample
  rare<-rrarefy(df,min)
  # rrarefy generates one randomly rarefied community data frame or vector of given sample size
  # The random rarefaction is made without replacement so that the variance of rarefied communities is rather related to rarefaction proportion than to the size of the sample

  s.ent<-vegan::diversity(rare, index="shannon")
  # calculate Shannon entropy from rarefied data

  sdiv<- exp(s.ent) # Shannon Diversity aka Hill number 1
  # calculate Shannon diversity from rarefied data

  return(sdiv)
}

ave.sdiv <- data.frame(AveShanDiv=rowMeans(data.frame(lapply(as.list(1:100), function(x) sdiv.rarefy(bac.ASV_table[,-1], min.rar)))))
# ^ contains average Shannon Diversity calculations after rarefaction and Shan div calculations 100 times
# sdiv.rarefy function --> Rrarefy ASV warnings() table with given min, calculate Shannon diversity
# lapply(as.list(1:100), function(x) sdiv.rarefy(bac.ASV_table[,-1], min.rar))) --> creates list by running sdiv.rarefy() function 100 times & storing results in list each time
# data.frame(lapply(...)) --> saves lapply() output as data frame, not list
# rowMeans(df(lapply(...))) --> get rowMeans aka sample means (aka the average) of Shannon diversity after 100 calculations of Shan div.

# compare our results with just one iteration of Shannon Div calculation
## calculate Shannon entropy from rarefied data
bac.S.Ent<-vegan::diversity(bac.ASV.rar, index="shannon")

## calculations Shannon diversity from entropy
bac.S.Div<- exp(bac.S.Ent) # Shannon Diversity aka Hill number 1

# Results from bac.S.Div are similar to ave.sdiv results which is great!

# do the same thing with Shannon entropy just so you have those data
# calculate the average Shannon entropy of rarefied data (taking average of 100 calculations of Shannon entropy)
ave.s.ent<-data.frame(AveShanEnt=rowMeans(data.frame(lapply(as.list(1:100), function(x) vegan::diversity(bac.ASV.rar, index="shannon")))))

bac.ASV.probs<-drarefy(bac.ASV_table[,-1],min.rar)
# drarefy returns probabilities that species occur in a rarefied community of size sample

#### Variance Stabilizing Transformation VST of Raw (?) counts ####
#Prepare Contig Feature Count Data for Normalization w/ DESeq2
# make sure count data & mgm_meta are in the same order
bac.ASV_table[1:5,1:5]

#bac.ASV_matrix<-as.matrix(bac.ASV_table[,!names(bac.ASV_table) %in% c("SampleID")]) # convert count table into matrix & exclude column called SampleID
bac.ASV_matrix2<-t(as.matrix(bac.ASV_table[,-1])) # transpose matrix so ASVs are rows, samples are columns
# ^ will be used in DESeq2 functions
rownames(dust_meta) %in% colnames(bac.ASV_matrix2) # check if rownames in dust_meta (SampleID) match column names in count data
dim(dust_meta)
dim(bac.ASV_matrix2)

# Reorder matrix to match order of dust_meta
bac.ASV_matrix2=bac.ASV_matrix2[,rownames(dust_meta)] ## reorder ASV matrix by column name to match order of rownames in dust_meta
colnames(bac.ASV_matrix2) # sanity check that this reordering worked
rownames(dust_meta)

# create the DESeq DataSet object for DGE analysis
# DESeq2 needs whole # data, so need raw read counts, NOT coverage for these tables...questionable
# b_dds has the scaled coverage that was calculated by dividing reads from featureCounts by gene lengths
b_dds<-DESeqDataSetFromMatrix(countData=round(bac.ASV_matrix2), colData = dust_meta, design = ~ 1)

# design = ~ 1 means no design
head(counts(b_dds))
colSums(counts(b_dds)) %>% barplot

# Estimate size factor - aka normalization factor, divide all read counts by each size factor per sample
b_dds <- estimateSizeFactors(b_dds,type="ratio")
## To calculate size factor in DESeq2:
# calculates geometric mean of each gene in each sample and uses this as a pseudoreference
# calculates ratio of each sample by dividing each gene count by it's pseudoreference in each sample
# The median value of all ratios for a given sample is taken as the normalization factor (size factor)
# The differentially expressed genes should not affect the median value
# median of ratios method makes the assumption that not ALL genes are differentially expressed; therefore, the normalization factors should account for sequencing depth and RNA composition of the sample
## (large outlier genes will not represent the median ratio values)
# more here: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html

sizeFactors(b_dds)

plot(sizeFactors(b_dds), colSums(counts(b_dds)))

### Variance Stabilizing Transformation
# you should be able to use matrix or DESeq2 object for this next function, but matrix was not working?
b_vst1 <- varianceStabilizingTransformation(b_dds) # add pseudocount
assay(b_vst1) #see output of VST

b.vst<-assay(b_vst1)

#### Compare Transformed ASVs vs Raw Counts ####

bac.ASV_table$SampleID = factor(bac.ASV_table$SampleID, levels=unique(bac.ASV_table$SampleID[order(rowSums(bac.ASV_table[,-1]))]), ordered=TRUE)

raw.tot.counts<-ggplot(data=bac.ASV_table, aes(x=bac.ASV_table$SampleID, y=rowSums(bac.ASV_table[,-1]))) +
  geom_bar(stat="identity",colour="black",fill="dodgerblue")+theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(title="Total ASVs per Sample",subtitle="Based on Raw ASV Counts")+ylab("Total ASVs")+xlab("SampleID")

ggsave(raw.tot.counts,filename = "figures/AlphaDiversity/SSD_16S_Total_ASVs_per_Sample_barplot.png", width=13, height=10, dpi=600)

# Calculate rarefied ASVs per Sample

ggplot(data=as.data.frame(bac.ASV.rar), aes(x=rownames(bac.ASV.rar), y=rowSums(bac.ASV.rar))) +
  geom_bar(stat="identity",colour="black",fill="dodgerblue")+theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(title="Total Rarefied ASVs per Sample",subtitle="Based on Rarefied ASV Counts")+ylab("Total Rarefied ASVs")+xlab("SampleID")

# Calculate VST ASVs per Sample
total_vst_asvs<-data.frame(ASV_Total=colSums(b.vst),dust_meta)
total_vst_asvs$SampleID = factor(total_vst_asvs$SampleID, levels=unique(total_vst_asvs$SampleID[order(total_vst_asvs$ASV_Total)]), ordered=TRUE)

ggplot(data=total_vst_asvs, aes(x=SampleID, y=ASV_Total,fill=Sample_Type)) +
  geom_bar(stat="identity",colour="black")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## average ASV per sample month & depth
aggregate(bac.ASV_all$Count, list(bac.ASV_all$SampleMonth), FUN=mean)
aggregate(bac.ASV_all$Count, list(bac.ASV_all$Depth_m), FUN=mean)

#### Alpha Diversity & Species Richness - Rarefied Data ####

## Calculate Shannon Diversity (abundance + richness considered in diversity calculation)
# if you have another package loaded that has a diversity function, you can specify that you want to use vegan's diversity function as shown below
# Shan_ent.16s.rar<-vegan::diversity(bac.ASV.rar, index="shannon") # Shannon entropy
# Shan_div.16s.rar<- exp(Shan_ent.16s.rar) # Shannon Diversity aka Hill number 1

# create data frame with Shannon entropy and Shannon diversity values
div_16s.rar<-data.frame(Bac_Shannon_Entropy=ave.s.ent,AveShanDiv=ave.sdiv)
class(div_16s.rar)
head(div_16s.rar)

div_16s.rar$SampleID<-rownames(div_16s.rar)
head(div_16s.rar)

# Calculate species richness (number of species per sample)
specnumber(bac.ASV.rar)

# Create a DF with Species Richness
S_16s.rar<-data.frame(Bac_Species_Richness=specnumber(bac.ASV.rar), SampleID=rownames(bac.ASV.rar)) # finds # of species per sample using RAW count data; if MARGIN = 2 it finds frequencies of species

# merge richness and diversity dataframes together
d.r_16s.rar<-merge(div_16s.rar, S_16s.rar, by.x="SampleID", by.y="SampleID")

# merge w/ dust_meta
bac.div.metadat.rar <- merge(d.r_16s.rar,dust_meta, by.x="SampleID", by.y="SampleID")
head(bac.div.metadat.rar)
class(bac.div.metadat.rar) # want data frame

unique(bac.div.metadat.rar$SampDate) # see how many elements there are in the Group variable
unique(bac.div.metadat.rar$Site) # see how many elements there are in the Group variable

# drop the outliers
#bac.div.metadat.rar<-bac.div.metadat.rar[bac.div.metadat.rar$AveShanDiv<300 & bac.div.metadat.rar$Bac_Species_Richness>100,]

# create numeric variable for depth to be used for models later

# Find highest/lowest values of Shannon div per sample date
max(bac.div.metadat.rar$AveShanDiv[bac.div.metadat.rar$Site=="WI"]) # max div WI
min(bac.div.metadat.rar$AveShanDiv[bac.div.metadat.rar$Site=="WI"]) # min div WI
bac.div.metadat.rar[bac.div.metadat.rar$Site=="WI",]

max(bac.div.metadat.rar$AveShanDiv[bac.div.metadat.rar$Site=="DP"]) # max div DP
min(bac.div.metadat.rar$AveShanDiv[bac.div.metadat.rar$Site=="DP"]) # min div DP
bac.div.metadat.rar[bac.div.metadat.rar$Site=="DP",]

max(bac.div.metadat.rar$AveShanDiv[bac.div.metadat.rar$Site=="BDC"]) # max div BDC
min(bac.div.metadat.rar$AveShanDiv[bac.div.metadat.rar$Site=="BDC"]) # min div BDC
bac.div.metadat.rar[bac.div.metadat.rar$Site=="BDC",]

max(bac.div.metadat.rar$AveShanDiv[bac.div.metadat.rar$Site=="PD"]) # max div BDC
min(bac.div.metadat.rar$AveShanDiv[bac.div.metadat.rar$Site=="PD"]) # min div BDC
bac.div.metadat.rar[bac.div.metadat.rar$Site=="PD",]

# Find highest/lowest values of Species richness per sample date
max(bac.div.metadat.rar$Bac_Species_Richness[bac.div.metadat.rar$Site=="WI"]) # max sr WI
min(bac.div.metadat.rar$Bac_Species_Richness[bac.div.metadat.rar$Site=="WI"]) # min sr WI
bac.div.metadat.rar[bac.div.metadat.rar$Site=="WI",]

max(bac.div.metadat.rar$Bac_Species_Richness[bac.div.metadat.rar$Site=="DP"]) # max sr DP
min(bac.div.metadat.rar$Bac_Species_Richness[bac.div.metadat.rar$Site=="DP"]) # min sr DP
bac.div.metadat.rar[bac.div.metadat.rar$Site=="DP",]

max(bac.div.metadat.rar$Bac_Species_Richness[bac.div.metadat.rar$Site=="BDC"]) # max sr BDC
min(bac.div.metadat.rar$Bac_Species_Richness[bac.div.metadat.rar$Site=="BDC"]) # min sr BDC
bac.div.metadat.rar[bac.div.metadat.rar$Site=="BDC",]

max(bac.div.metadat.rar$Bac_Species_Richness[bac.div.metadat.rar$Site=="PD"]) # max sr BDC
min(bac.div.metadat.rar$Bac_Species_Richness[bac.div.metadat.rar$Site=="PD"]) # min sr BDC
bac.div.metadat.rar[bac.div.metadat.rar$Site=="PD",]

# save diversity data
save.image("data/SSeaDust_AlphaDiv_Data_Rarefied.Rdata")

#### Using Shapiro-Wilk test for Normality - Rarefied Data ####
shapiro.test(bac.div.metadat.rar$AveShanDiv) # what is the p-value?
# p-value = 0.7112
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(bac.div.metadat.rar$AveShanDiv, col="blue") # with outliars

# visualize Q-Q plot for alpha div
# The Q-Q plot, or quantile-quantile plot, is a graphical tool to help us assess if a set of data plausibly came from some theoretical distribution such as a normal or exponential.
# For example, if we run a statistical analysis that assumes our residuals are normally distributed, we can use a normal Q-Q plot to check that assumption
# more on Q-Q plots here: https://data.library.virginia.edu/understanding-q-q-plots/
# more here too: https://grodri.github.io/glms/notes/c2s9#:~:text=8%20The%20Q%2DQ%20Plot,versus%20quantiles%20of%20a%20distribution.
qqnorm(bac.div.metadat.rar$AveShanDiv, pch = 1, frame = FALSE)
qqline(bac.div.metadat.rar$AveShanDiv, col = "red", lwd = 2)

shapiro.test(bac.div.metadat.rar$Bac_Species_Richness) # what is the p-value?
# p-value = 0.009454
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(bac.div.metadat.rar$Bac_Species_Richness, col="blue")

# visualize Q-Q plot for species richness
qqnorm(bac.div.metadat.rar$Bac_Species_Richness, pch = 1, frame = FALSE) # with outliars
qqline(bac.div.metadat.rar$Bac_Species_Richness, col = "red", lwd = 2)

qqnorm(bac.div.metadat.rar$Bac_Species_Richness, pch = 1, frame = FALSE) # without outliars
qqline(bac.div.metadat.rar$Bac_Species_Richness, col = "red", lwd = 2)

#### Compare Means of Shannon Diversity - Rarefied Data ####
head(bac.div.metadat.rar)
WI.divmeta<-bac.div.metadat.rar[bac.div.metadat.rar$Site=="WI",]
DP.divmeta<-bac.div.metadat.rar[bac.div.metadat.rar$Site=="DP",]
BDC.divmeta<-bac.div.metadat.rar[bac.div.metadat.rar$Site=="BDC",]
PD.divmeta<-bac.div.metadat.rar[bac.div.metadat.rar$Site=="PD",]

# First run individual t-tests by time point
# Aug 2021 vs DP
t.test.WI.DP<-t.test(WI.divmeta$AveShanDiv,DP.divmeta$AveShanDiv)
t.test.WI.DP$p.value

# Aug 2021 vs BDC
t.test.WI.BDC<-t.test(WI.divmeta$AveShanDiv,BDC.divmeta$AveShanDiv)
t.test.WI.BDC$p.value

# DP vs BDC
t.test.DP.BDC<-t.test(DP.divmeta$AveShanDiv,BDC.divmeta$AveShanDiv)
t.test.DP.BDC$p.value

# Combine the p-values
ttest.Div.pvals<-c(t.test.WI.DP$p.value, t.test.WI.BDC$p.value, t.test.DP.BDC$p.value)

# Adjust the p-values based on the # of comparisons you did
p.adjust(ttest.Div.pvals, method="bonferroni",n=3)
# ^ matches findings on the figure, which uses the t_test function from the rstatix package (see geom_pwc())

#### Compare Means of Species Richness - Rarefied Data ####
head(bac.div.metadat.rar)
#WI.divmeta<-bac.div.metadat.rar[bac.div.metadat.rar$SampDate=="WI",]
#DP.divmeta<-bac.div.metadat.rar[bac.div.metadat.rar$SampDate=="DP",]
#BDC.divmeta<-bac.div.metadat.rar[bac.div.metadat.rar$SampDate=="BDC",]

# First run individual wilcox-tests by time point
# Aug 2021 vs DP
w.test.WI.DP<-wilcox.test(WI.divmeta$Bac_Species_Richness,DP.divmeta$Bac_Species_Richness)
w.test.WI.DP$p.value

# Aug 2021 vs BDC
w.test.WI.BDC<-wilcox.test(WI.divmeta$Bac_Species_Richness,BDC.divmeta$Bac_Species_Richness)
w.test.WI.BDC$p.value

# DP vs BDC
w.test.DP.BDC<-wilcox.test(DP.divmeta$Bac_Species_Richness,BDC.divmeta$Bac_Species_Richness)
w.test.DP.BDC$p.value

# Combine the p-values
wilc.SR.pvals<-c(w.test.WI.DP$p.value, w.test.WI.BDC$p.value, w.test.DP.BDC$p.value)

# Adjust the p-values based on the # of comparisons you did
p.adjust(wilc.SR.pvals, method="bonferroni",n=3)
# ^ matches findings on the figure, which uses the wilcox_test function from the rstatix package (see geom_pwc())

#### Visualize Alpha Diversity & Species Richness - from Rarefied Data ####

## Shannon Diversity by Site & Sample Date
bac.a.div.rar<-ggplot(bac.div.metadat.rar, aes(x=Site, y=AveShanDiv)) +geom_jitter(aes(color=Site), size=3, width=0.15, height=0) +
  scale_color_manual(name ="Site",values=unique(bac.div.metadat.rar$Site_Color[order(bac.div.metadat.rar$Site)]),labels=c("PD","BDC","DP","WI")) +
  geom_boxplot(fill=NA, outlier.color=NA)+scale_x_discrete(labels=c("PD","BDC","DP","WI"))+theme_bw()+theme_classic()+
  labs(title = "Dust Bacterial Shannon Diversity by Site", subtitle="Using Rarefied Counts", x="Site", y="Shannon Diversity", color="Depth (m)")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))
  #geom_pwc(method = "t_test", label = "p.adj.format",p.adjust.method = "bonferroni")

ggsave(bac.a.div.rar,filename = "figures/AlphaDiversity/RarefiedCounts/SSD_16S_rarefied_alpha_diversity_site_boxplot.png", width=13, height=10, dpi=600)

bac.a.div.rarB<-ggplot(bac.div.metadat.rar, aes(x=Site, y=AveShanDiv)) +geom_jitter(aes(color=SampDate,shape=Site), size=3, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date",values=unique(bac.div.metadat.rar$SampDate_Color[order(bac.div.metadat.rar$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  scale_shape_manual(values = c(7,10, 15,16)) +  geom_boxplot(fill=NA, outlier.color=NA)+scale_x_discrete(labels=c("PD","BDC","DP","WI"))+theme_bw()+theme_classic()+
  labs(title = "Dust Bacterial Shannon Diversity by Collection Date & Site", subtitle="Using Rarefied Counts", x="Site", y="Shannon Diversity", color="Collection Date")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))
#geom_pwc(method = "t_test", label = "p.adj.format",p.adjust.method = "bonferroni")

ggsave(bac.a.div.rarB,filename = "figures/AlphaDiversity/RarefiedCounts/SSD_16S_rarefied_alpha_diversity_sampledate_site_boxplot.png", width=13, height=10, dpi=600)

## Species Richness by Site & Sample Date
bac.a.sr.rar<-ggplot(bac.div.metadat.rar, aes(x=Site, y=Bac_Species_Richness)) +geom_jitter(aes(color=Site), size=3, width=0.15, height=0) +
  scale_color_manual(name ="Site",values=unique(bac.div.metadat.rar$Site_Color[order(bac.div.metadat.rar$Site)]),labels=c("PD","BDC","DP","WI")) +
  geom_boxplot(fill=NA, outlier.color=NA)+scale_x_discrete(labels=c("PD","BDC","DP","WI"))+theme_bw()+theme_classic()+
  labs(title = "Dust Bacterial Species Richness by Site", subtitle="Using Rarefied Counts", x="Site", y="Species Richness", color="Depth (m)")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))
#geom_pwc(method = "t_test", label = "p.adj.format",p.adjust.method = "bonferroni")

ggsave(bac.a.sr.rar,filename = "figures/SpeciesRichness/RarefiedCounts/SSD_16S_rarefied_species_richness_site_boxplot.png", width=13, height=10, dpi=600)

bac.a.sr.rarB<-ggplot(bac.div.metadat.rar, aes(x=Site, y=Bac_Species_Richness)) +geom_jitter(aes(color=SampDate,shape=Site), size=3, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date",values=unique(bac.div.metadat.rar$SampDate_Color[order(bac.div.metadat.rar$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  scale_shape_manual(values = c(7,10, 15,16)) +  geom_boxplot(fill=NA, outlier.color=NA)+scale_x_discrete(labels=c("PD","BDC","DP","WI"))+theme_bw()+theme_classic()+
  labs(title = "Dust Bacterial Species Richness by Collection Date & Site", subtitle="Using Rarefied Counts", x="Site", y="Species Richness", color="Collection Date")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))
#geom_pwc(method = "t_test", label = "p.adj.format",p.adjust.method = "bonferroni")

ggsave(bac.a.sr.rarB,filename = "figures/SpeciesRichness/RarefiedCounts/SSD_16S_rarefied_species_richness_sampledate_site_boxplot.png", width=13, height=10, dpi=600)

#### Compare Variance w/ Shannon Diversity ####
# shannon diversity is normally distributed; used rarefied counts to calculate ShanDiv
# use the following statisitcal tests for variance comparisons
## ANOVA: are variances significantly different between groups
## Tukey test: which groups' variances are significant different from one another
## Levene's test: is variance homogenous aka equal across samples?

fit1<-aov(AveShanDiv ~ Site, data=bac.div.metadat.rar)
# ANOVA is basically a regression but w/ categorical variables more info here https://www.statology.org/anova-vs-regression/
#pairwise.adonis(bac.div.metadat.rar$AveShanDiv, bac.div.metadat.rar$SampDate, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(fit1)
#Df           Sum Sq Mean Sq    F value   Pr(>F)
#SampDate     2   3648  1824.0   6.677 0.0057 **
# Residuals   21   5737   273.2

p.adjust(summary(fit1)[[1]][["Pr(>F)"]][1],method="bonferroni")

# Tukey test - tells us which groups are significantly different from each other (more here: https://www.r-bloggers.com/2013/06/anova-and-tukeys-test-on-r/)
Tuk1<-TukeyHSD(fit1)
Tuk1$SampDate
#                             diff        lwr      upr      p adj
# DP-WI 25.076433   4.245680 45.90719 0.016639615 *
# BDC-WI    27.111615   6.280862 47.94237 0.009571804 **
# BDC-DP   2.035182 -18.795571 22.86593 0.967173383

# fit.0<-aov(DustComplexity ~ as.factor(Elevation), data=bac.div.metadat.rar)
# summary(fit.0)
# TukeyHSD(fit.0)

# Levene's test - test for homogeneity of variance
## Levene’s test is an inferential statistic used to check if the variances of a variable obtained for two or more groups are equal or not when data comes from a non-normal distribution
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## t test assumes that variances are the same, so Levene's test needs to be non significant
car::leveneTest(AveShanDiv ~ Site, data = bac.div.metadat.rar, center=mean)
# Levene's Test for Homogeneity of Variance (center = mean)
#       Df F value Pr(>F)
# group  2  1.3552 0.5035
#       21
# ^ p value is > 0.05 so we cannot reject the Null hypothesis -- variances are equal

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
# fligner.test(AveShanDiv ~ SampDate, data = bac.div.metadat.rar)
# Fligner-Killeen:med chi-squared = 1.9504, df = 2, p-value = 0.3771
# Which shows that the data do not deviate significantly from homogeneity.

compare_means(AveShanDiv ~ Site, data=bac.div.metadat.rar, method="anova",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input

#### Compare Variance w/ Species Richness ####
# species richness is NOT normally distributed; used rarefied counts to calculate SR
# use the following statisitcal tests for variance comparisons
## Kruskal: are variances significantly different between groups
## Dunn test: which groups' variances are significant different from one another
## Fligner test: is variance homogenous aka equal across samples?

# Kruskal-Wallis test is an ANOVA for non-normal data
fit2<-kruskal.test(Bac_Species_Richness ~ Site, data=bac.div.metadat.rar)

fit2

# Instead of using Tukey test, we can use Dunn's test to see which groups significantly vary if Kruskal-Wallis test is significant
# ANOVA + Tukey for normally distributed data, Kruskal-Wallis + Dunn's test for non-normal data
rstatix::dunn_test(bac.div.metadat.rar, Bac_Species_Richness ~ Site, p.adjust.method = "bonferroni", detailed = TRUE)

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(Bac_Species_Richness ~ Site, data = bac.div.metadat.rar)
# Fligner-Killeen:med chi-squared = 0.6397, df = 3, p-value = 0.8873
# Which shows that the data do NOT deviate significantly from homogeneity.

compare_means(Bac_Species_Richness ~ Site, data=bac.div.metadat.rar, method="anova",p.adjust.method = "bonferroni")

#### Linear Regression - Shannon Diversity ####
## here the focus is comparing env variables of interest to see if they can predict diversity and richness
head(bac.div.metadat.rar)

# just look at everything at once in step-wise fashion
step1<-step(glm(formula = AveShanDiv ~ ., data=bac.div.metadat.rar[,c(3,11,13:14,18:20)]))
summary(step1)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)      86.842      3.002  28.926  < 2e-16 ***
#   ORP_mV           -6.478      3.664  -1.768 0.092285 .
# Temp_DegC       -21.212      4.471  -4.744 0.000124 ***
#   Sulfate_milliM  -10.994      3.886  -2.829 0.010374 *

# Run GLM w/ all env variables of interest
div.glm.all<-glm(formula = AveShanDiv ~ ., data=bac.div.metadat.rar[,c(3,11,13:14,18:20)])
summary(div.glm.all)

# NOTE: GLMs use maximum likelihood and not ordinary least squares, so R^2 is not provided... (more here https://stats.stackexchange.com/questions/232471/why-does-the-glm-function-does-not-return-an-r2-value#:~:text=Since%20the%20coefficient%20of%20determination,reason%20to%20display%20this%20measure.)
# BUT we can use McFadden's R-squared, which uses the following formula..
## 1 - (log likelihood of current model / log likelihood null model)
## null model is a model with only the intercept, whereas the current model is the one you are fitting
## more here: https://www.statology.org/glm-r-squared/

# McFadden's R-squared
with(summary(div.glm.all), 1 - div.glm.all$deviance/div.glm.all$null.deviance)
# 0.5734642 --> model explains a lot of variance in data

div.glm.p<-coef(summary(div.glm.all))[,4] # p-values
Div.GLM.Pval<-data.frame(Div.GLM.AdjPval=p.adjust(div.glm.p, method="bonferroni",n=length(div.glm.p)),Div.GLM.Pval=div.glm.p)
Div.GLM.Pval
#                               Div.GLM.AdjPval Div.GLM.Pval
# (Intercept)                    8.353610e-15 1.193373e-15
# DO_Percent_Local               1.000000e+00 6.614790e-01
# ORP_mV                         1.000000e+00 3.249102e-01
# Temp_DegC                      4.793721e-02 6.848173e-03 *
# Dissolved_OrganicMatter_RFU    1.000000e+00 6.835361e-01
# Sulfate_milliM                 2.762433e-01 3.946333e-02 *
# Sulfide_microM                 1.000000e+00 6.239275e-01

# [Example Code for Different Models]

s.div.glm.fit1<-glm(formula = AveShanDiv ~ DO_Percent_Local, data=bac.div.metadat.rar)%>%
  adjust_pvalue(method="bonferroni")

# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.div.glm.fit1)

# sanity check that lm() vs glm(familiy=Gaussian) is the same thing - and it is!
summary(lm(formula = AveShanDiv ~ DO_Percent_Local, data=bac.div.metadat.rar)%>%
          adjust_pvalue(method="bonferroni"))

# code for mixed effects model
mixed1 = lmer(AveShanDiv ~ DO_Percent_Local+ (1 | Depth_m), data = bac.div.metadat.rar)
summary(mixed1)

# Shan Div ~ DO%
plot(AveShanDiv ~ DO_Percent_Local, data=bac.div.metadat.rar,col=SampDate_Color)

# Shan Div ~ ORP

plot(AveShanDiv ~ ORP_mV, data=bac.div.metadat.rar,col=SampDate_Color)

# Shan Div ~ Temp (C)

plot(AveShanDiv ~ Temp_DegC, data=bac.div.metadat.rar,col=SampDate_Color)

# Shan Div ~ DOM

plot(AveShanDiv ~ Dissolved_OrganicMatter_RFU, data=bac.div.metadat.rar,col=SampDate_Color)

# Shan Div ~ Sulfate

plot(AveShanDiv ~ Sulfate_milliM, data=bac.div.metadat.rar,col=SampDate_Color)

# Shan Div ~ Sulfide

plot(AveShanDiv ~ Sulfide_microM, data=bac.div.metadat.rar,col=SampDate_Color)

# Shan Div ~ Depth
plot(AveShanDiv ~ Depth.num, data=bac.div.metadat.rar,col=SampDate_Color)

#### Linear Regression - Species Richness ####
# do env variables predict species richness?
# NOTES: species richness (using rarefied counts) is NOT normally distributed [but shannon div is...]
# need to decide if we should use the following models (designed for non-normal discrete count data): Poisson, quasi-Poisson, or Negative Binomial
# below we run these three models, calculate their residuals, and compare the dispersion aka variance of their residuals after each model
# we can also run a likelihood ratio test that will compare the goodness to fit of two models
# reminder: A goodness-of-fit is a statistical test that tries to determine whether a set of observed values match those expected under the applicable model [aka do our observations reflect what we would expect to see?]
## Goodness of Fit of a linear regression model attempts to get at how well a model fits a given set of data, or how well it will predict a future set of observations.
head(bac.div.metadat.rar)

# Generate poisson, quasiposson, and negative binomial (NB) regressions to see which env variables can predict species richness
sr.pois.all<-glm(formula = Bac_Species_Richness ~ ., family=poisson,data=bac.div.metadat.rar[,c(4,11,13:14,18:20)])
sr.qpois.all<-glm(formula = Bac_Species_Richness ~ ., family=quasipoisson,data=bac.div.metadat.rar[,c(4,11,13:14,18:20)])
sr.NB.all<-glm.nb(formula = Bac_Species_Richness ~ ., data=bac.div.metadat.rar[,c(4,11,13:14,18:20)])

# Likeliehood ratio test between Poisson model and NB model
# Null H: both models fit data equally well; alt H: one model fits data better than other model
lrtest(sr.pois.all, sr.NB.all) # p < 2.2e-16 -> indicates that NB model has much better fit than Poisson model
lrtest(sr.qpois.all, sr.NB.all) # no p value, not quite sure why?
lrtest(sr.qpois.all, sr.pois.all) # no p value

# We can also plot the residuals for each model and see their dispersion/spread
#Residual plot for Poisson regression
p_res <- resid(sr.pois.all) # Poisson residuals
qp_res <- resid(sr.qpois.all) #Quasipoisson Residuals
nb_res <- resid(sr.NB.all) # Negative Binomial Residuals

par(mfrow = c(1, 3))

#Residual plot for poisson regression
plot(fitted(sr.pois.all), p_res, col='blue', pch=16,
     xlab='Predicted Offers', ylab='Standardized Residuals', main='Poisson')
abline(0,0)

#Residual plot for quasipoisson regression
plot(fitted(sr.qpois.all), qp_res, col='blue', pch=16,
     xlab='Predicted Offers', ylab='Standardized Residuals', main='Quasi-Poisson')
abline(0,0)

#Residual plot for negative binomial regression
plot(fitted(sr.NB.all), nb_res, col='blue', pch=16,
     xlab='Predicted Offers', ylab='Standardized Residuals', main='Negative Binomial')
abline(0,0)

# NB shows least amount of spread in residuals, suggests this is the better fitting model!

summary(sr.NB.all)

# NOTE: GLMs use maximum likelihood and not ordinary least squares, so R^2 is not provided... (more here https://stats.stackexchange.com/questions/232471/why-does-the-glm-function-does-not-return-an-r2-value#:~:text=Since%20the%20coefficient%20of%20determination,reason%20to%20display%20this%20measure.)
# BUT we can use McFadden's R-squared, which uses the following formula..
## 1 - (log likelihood of current model / log likelihood null model)
## null model is a model with only the intercept, whereas the current model is the one you are fitting
## more here: https://www.statology.org/glm-r-squared/

# McFadden's R-squared
with(summary(sr.NB.all), 1 - sr.NB.all$deviance/sr.NB.all$null.deviance)

sr.NB.p<-coef(summary(sr.NB.all))[,4] # p-values
SR.NB.Pval<-data.frame(SR.NB.AdjPval=p.adjust(sr.NB.p, method="bonferroni",n=length(sr.NB.p)),SR.NB.Pval=sr.NB.p)
SR.NB.Pval
#                            SR.NB.AdjPval SR.NB.Pval
# (Intercept)                     0.0000000 0.00000000
# DO_Percent_Local                1.0000000 0.54606633
# ORP_mV                          0.2310918 0.03301311 *
# Temp_DegC                       1.0000000 0.35482885
# Dissolved_OrganicMatter_RFU     1.0000000 0.32657918
# Sulfate_milliM                  0.3101118 0.04430168 *
# Sulfide_microM                  0.4688277 0.06697538

# Species Richness ~ DO%

plot(Bac_Species_Richness ~ DO_Percent_Local, data=bac.div.metadat.rar,col=SampDate_Color)

# Species Richness ~ ORP

plot(Bac_Species_Richness ~ ORP_mV, data=bac.div.metadat.rar,col=SampDate_Color)

# Species Richness ~ Temp

plot(Bac_Species_Richness ~ Temp_DegC, data=bac.div.metadat.rar,col=SampDate_Color)

# Species Richness ~ DOM

plot(Bac_Species_Richness ~ Dissolved_OrganicMatter_RFU, data=bac.div.metadat.rar,col=SampDate_Color)

# Species Richness ~ Sulfate

plot(Bac_Species_Richness ~ Sulfate_milliM, data=bac.div.metadat.rar,col=SampDate_Color)

# Species Richness ~ Sulfide

plot(Bac_Species_Richness ~ Sulfide_microM, data=bac.div.metadat.rar,col=SampDate_Color)

# Species Richness ~ Depth

plot(Bac_Species_Richness ~ Depth.num, data=bac.div.metadat.rar,col=SampDate_Color)

#### Prep Data for Linear Regressions within Timepoints ####
## here the focus is comparing dust complexity to alpha diversity, species richness, & elevation
head(bac.div.metadat.rar)

# create the dataframes
aug21.div<-subset(bac.div.metadat.rar, bac.div.metadat.rar$SampDate=="WI")
dec21.div<-subset(bac.div.metadat.rar, bac.div.metadat.rar$SampDate=="DP")
apr22.div<-subset(bac.div.metadat.rar, bac.div.metadat.rar$SampDate=="BDC")

#### August - Shannon Diversity ####
# WI
aug21.div.glm.fit1<-glm(formula = AveShanDiv ~ DO_Percent_Local, data=aug21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(aug21.div.glm.fit1)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)       0.0112706  0.0004664  24.166   <2e-16 ***
#DO_Percent_Local -0.0009547  0.0005092  -1.875   0.0687 .

aug21.div.glm.fit2<-glm(formula = AveShanDiv ~ ORP_mV, data=aug21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(aug21.div.glm.fit2)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)   0.0112343  0.0004731  23.745   <2e-16 ***
#ORP_mV      -0.0001512  0.0004968  -0.304    0.763

aug21.div.glm.fit3<-glm(formula = AveShanDiv ~ Temp_DegC, data=aug21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(aug21.div.glm.fit3)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept) 0.0113225  0.0004451  25.439   <2e-16 ***
#Temp_DegC   0.0011947  0.0004861   2.458   0.0188 *

aug21.div.glm.fit4<-glm(formula = AveShanDiv ~ Dissolved_OrganicMatter_RFU, data=aug21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(aug21.div.glm.fit4)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)                 0.0112493  0.0004708  23.896   <2e-16 ***
#Dissolved_OrganicMatter_RFU 0.0004269  0.0004659   0.916    0.365

aug21.div.glm.fit5<-glm(formula = AveShanDiv ~ Sulfate_milliM, data=aug21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(aug21.div.glm.fit5)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)     0.0112325  0.0004772  23.536   <2e-16 ***
#Sulfate_milliM -0.0002664  0.0004893  -0.545    0.589

aug21.div.glm.fit6<-glm(formula = AveShanDiv ~ Sulfide_microM, data=aug21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(aug21.div.glm.fit6)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)    0.0112379  0.0004723   23.80   <2e-16 ***
#Sulfide_microM 0.0002944  0.0005160    0.57    0.572

fit1<-aov(AveShanDiv ~ as.factor(Depth_m), data=aug21.div)
# ANOVA is basically a regression but w/ categorical variables more info here https://www.statology.org/anova-vs-regression/
#pairwise.adonis(aug21.div$AveShanDiv, aug21.div$Depth_m, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(fit1)
#Df           Sum Sq Mean Sq    F value   Pr(>F)
#Depth_m      7   4097   585.3   1.114   0.38
#Residuals   31  16294   525.6
Tuk1<-TukeyHSD(fit1)
Tuk1$Depth_m

#plot(AveShanDiv ~ Depth_m, data=aug21.div)
#abline(aov(DustComplexity ~ Elevation, data=aug21.div))

# fit.0<-aov(DustComplexity ~ as.factor(Elevation), data=aug21.div)
# summary(fit.0)
# TukeyHSD(fit.0)
# Levene's test with one independent variable
## Levene's tests whether variances of 2 samples are equal
## we want variances to be the same -- want NON SIGNIFICANCE!
## t test assumes that variances are the same, so Levene's test needs to be non significant
## Fligner's test is a Levene's test for data that are not normally distributed
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(AveShanDiv ~ Depth_m, data = aug21.div)
# Fligner-Killeen:med chi-squared = 4.091, df = 7, p-value = 0.7692
# Which shows that the data do not deviate significantly from homogeneity.
compare_means(AveShanDiv ~ Depth_m, data=aug21.div, method="anova",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input


#### August - Species Richness ####
# WI
aug21.sr.glm.fit1<-glm(formula = Bac_Species_Richness ~ DO_Percent_Local, data=aug21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(aug21.sr.glm.fit1)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)       0.0112706  0.0004664  24.166   <2e-16 ***
#DO_Percent_Local -0.0009547  0.0005092  -1.875   0.0687 .

aug21.sr.glm.fit2<-glm(formula = Bac_Species_Richness ~ ORP_mV, data=aug21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(aug21.sr.glm.fit2)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)   0.0112343  0.0004731  23.745   <2e-16 ***
#ORP_mV      -0.0001512  0.0004968  -0.304    0.763

aug21.sr.glm.fit3<-glm(formula = Bac_Species_Richness ~ Temp_DegC, data=aug21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(aug21.sr.glm.fit3)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept) 0.0113225  0.0004451  25.439   <2e-16 ***
#Temp_DegC   0.0011947  0.0004861   2.458   0.0188 *

aug21.sr.glm.fit4<-glm(formula = Bac_Species_Richness ~ Dissolved_OrganicMatter_RFU, data=aug21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(aug21.sr.glm.fit4)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)                 0.0112493  0.0004708  23.896   <2e-16 ***
#Dissolved_OrganicMatter_RFU 0.0004269  0.0004659   0.916    0.365

aug21.sr.glm.fit5<-glm(formula = Bac_Species_Richness ~ Sulfate_milliM, data=aug21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(aug21.sr.glm.fit5)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)     0.0112325  0.0004772  23.536   <2e-16 ***
#Sulfate_milliM -0.0002664  0.0004893  -0.545    0.589

aug21.sr.glm.fit6<-glm(formula = Bac_Species_Richness ~ Sulfide_microM, data=aug21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(aug21.sr.glm.fit6)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)    0.0112379  0.0004723   23.80   <2e-16 ***
#Sulfide_microM 0.0002944  0.0005160    0.57    0.572

fit1<-aov(Bac_Species_Richness ~ as.factor(Depth_m), data=aug21.div)
# ANOVA is basically a regression but w/ categorical variables more info here https://www.statology.org/anova-vs-regression/
#pairwise.adonis(aug21.div$Bac_Species_Richness, aug21.div$Depth_m, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(fit1)
#Df           Sum Sq Mean Sq    F value   Pr(>F)
#Depth_m      7   4097   585.3   1.114   0.38
#Residuals   31  16294   525.6
Tuk1<-TukeyHSD(fit1)
Tuk1$Depth_m

#plot(Bac_Species_Richness ~ Depth_m, data=aug21.div)
#abline(aov(DustComplexity ~ Elevation, data=aug21.div))

# fit.0<-aov(DustComplexity ~ as.factor(Elevation), data=aug21.div)
# summary(fit.0)
# TukeyHSD(fit.0)
# Levene's test with one independent variable
## Levene's tests whether variances of 2 samples are equal
## we want variances to be the same -- want NON SIGNIFICANCE!
## t test assumes that variances are the same, so Levene's test needs to be non significant
## Fligner's test is a Levene's test for data that are not normally distributed
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(Bac_Species_Richness ~ Depth_m, data = aug21.div)
# Fligner-Killeen:med chi-squared = 4.091, df = 7, p-value = 0.7692
# Which shows that the data do not deviate significantly from homogeneity.
compare_means(Bac_Species_Richness ~ Depth_m, data=aug21.div, method="anova",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input


#### December - Shannon Diversity ####
# DP
dec21.div.glm.fit1<-glm(formula = AveShanDiv ~ DO_Percent_Local, data=dec21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(dec21.div.glm.fit1)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)       0.0112706  0.0004664  24.166   <2e-16 ***
#DO_Percent_Local -0.0009547  0.0005092  -1.875   0.0687 .

dec21.div.glm.fit2<-glm(formula = AveShanDiv ~ ORP_mV, data=dec21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(dec21.div.glm.fit2)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)   0.0112343  0.0004731  23.745   <2e-16 ***
#ORP_mV      -0.0001512  0.0004968  -0.304    0.763

dec21.div.glm.fit3<-glm(formula = AveShanDiv ~ Temp_DegC, data=dec21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(dec21.div.glm.fit3)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept) 0.0113225  0.0004451  25.439   <2e-16 ***
#Temp_DegC   0.0011947  0.0004861   2.458   0.0188 *

dec21.div.glm.fit4<-glm(formula = AveShanDiv ~ Dissolved_OrganicMatter_RFU, data=dec21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(dec21.div.glm.fit4)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)                 0.0112493  0.0004708  23.896   <2e-16 ***
#Dissolved_OrganicMatter_RFU 0.0004269  0.0004659   0.916    0.365

dec21.div.glm.fit5<-glm(formula = AveShanDiv ~ Sulfate_milliM, data=dec21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(dec21.div.glm.fit5)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)     0.0112325  0.0004772  23.536   <2e-16 ***
#Sulfate_milliM -0.0002664  0.0004893  -0.545    0.589

dec21.div.glm.fit6<-glm(formula = AveShanDiv ~ Sulfide_microM, data=dec21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(dec21.div.glm.fit6)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)    0.0112379  0.0004723   23.80   <2e-16 ***
#Sulfide_microM 0.0002944  0.0005160    0.57    0.572

fit1<-aov(AveShanDiv ~ as.factor(Depth_m), data=dec21.div)
# ANOVA is basically a regression but w/ categorical variables more info here https://www.statology.org/anova-vs-regression/
#pairwise.adonis(dec21.div$AveShanDiv, dec21.div$Depth_m, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(fit1)
#Df           Sum Sq Mean Sq    F value   Pr(>F)
#Depth_m      7   4097   585.3   1.114   0.38
#Residuals   31  16294   525.6
Tuk1<-TukeyHSD(fit1)
Tuk1$Depth_m

#plot(AveShanDiv ~ Depth_m, data=dec21.div)
#abline(aov(DustComplexity ~ Elevation, data=dec21.div))

# fit.0<-aov(DustComplexity ~ as.factor(Elevation), data=dec21.div)
# summary(fit.0)
# TukeyHSD(fit.0)
# Levene's test with one independent variable
## Levene's tests whether variances of 2 samples are equal
## we want variances to be the same -- want NON SIGNIFICANCE!
## t test assumes that variances are the same, so Levene's test needs to be non significant
## Fligner's test is a Levene's test for data that are not normally distributed
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(AveShanDiv ~ Depth_m, data = dec21.div)
# Fligner-Killeen:med chi-squared = 4.091, df = 7, p-value = 0.7692
# Which shows that the data do not deviate significantly from homogeneity.
compare_means(AveShanDiv ~ Depth_m, data=dec21.div, method="anova",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input

#### December - Species Richness ####
# DP
dec21.sr.glm.fit1<-glm(formula = Bac_Species_Richness ~ DO_Percent_Local, data=dec21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(dec21.sr.glm.fit1)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)       0.0112706  0.0004664  24.166   <2e-16 ***
#DO_Percent_Local -0.0009547  0.0005092  -1.875   0.0687 .

dec21.sr.glm.fit2<-glm(formula = Bac_Species_Richness ~ ORP_mV, data=dec21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(dec21.sr.glm.fit2)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)   0.0112343  0.0004731  23.745   <2e-16 ***
#ORP_mV      -0.0001512  0.0004968  -0.304    0.763

dec21.sr.glm.fit3<-glm(formula = Bac_Species_Richness ~ Temp_DegC, data=dec21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(dec21.sr.glm.fit3)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept) 0.0113225  0.0004451  25.439   <2e-16 ***
#Temp_DegC   0.0011947  0.0004861   2.458   0.0188 *

dec21.sr.glm.fit4<-glm(formula = Bac_Species_Richness ~ Dissolved_OrganicMatter_RFU, data=dec21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(dec21.sr.glm.fit4)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)                 0.0112493  0.0004708  23.896   <2e-16 ***
#Dissolved_OrganicMatter_RFU 0.0004269  0.0004659   0.916    0.365

dec21.sr.glm.fit5<-glm(formula = Bac_Species_Richness ~ Sulfate_milliM, data=dec21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(dec21.sr.glm.fit5)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)     0.0112325  0.0004772  23.536   <2e-16 ***
#Sulfate_milliM -0.0002664  0.0004893  -0.545    0.589

dec21.sr.glm.fit6<-glm(formula = Bac_Species_Richness ~ Sulfide_microM, data=dec21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(dec21.sr.glm.fit6)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)    0.0112379  0.0004723   23.80   <2e-16 ***
#Sulfide_microM 0.0002944  0.0005160    0.57    0.572

fit1<-aov(Bac_Species_Richness ~ as.factor(Depth_m), data=dec21.div)
# ANOVA is basically a regression but w/ categorical variables more info here https://www.statology.org/anova-vs-regression/
#pairwise.adonis(dec21.div$Bac_Species_Richness, dec21.div$Depth_m, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(fit1)
#Df           Sum Sq Mean Sq    F value   Pr(>F)
#Depth_m      7   4097   585.3   1.114   0.38
#Residuals   31  16294   525.6
Tuk1<-TukeyHSD(fit1)
Tuk1$Depth_m

#plot(Bac_Species_Richness ~ Depth_m, data=dec21.div)
#abline(aov(DustComplexity ~ Elevation, data=dec21.div))

# fit.0<-aov(DustComplexity ~ as.factor(Elevation), data=dec21.div)
# summary(fit.0)
# TukeyHSD(fit.0)
# Levene's test with one independent variable
## Levene's tests whether variances of 2 samples are equal
## we want variances to be the same -- want NON SIGNIFICANCE!
## t test assumes that variances are the same, so Levene's test needs to be non significant
## Fligner's test is a Levene's test for data that are not normally distributed
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(Bac_Species_Richness ~ Depth_m, data = dec21.div)
# Fligner-Killeen:med chi-squared = 4.091, df = 7, p-value = 0.7692
# Which shows that the data do not deviate significantly from homogeneity.
compare_means(Bac_Species_Richness ~ Depth_m, data=dec21.div, method="anova",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input




#### April - Shannon Diversity ####
# BDC
apr22.div.glm.fit1<-glm(formula = AveShanDiv ~ DO_Percent_Local, data=apr22.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(apr22.div.glm.fit1)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)       0.0112706  0.0004664  24.166   <2e-16 ***
#DO_Percent_Local -0.0009547  0.0005092  -1.875   0.0687 .

apr22.div.glm.fit2<-glm(formula = AveShanDiv ~ ORP_mV, data=apr22.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(apr22.div.glm.fit2)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)   0.0112343  0.0004731  23.745   <2e-16 ***
#ORP_mV      -0.0001512  0.0004968  -0.304    0.763

apr22.div.glm.fit3<-glm(formula = AveShanDiv ~ Temp_DegC, data=apr22.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(apr22.div.glm.fit3)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept) 0.0113225  0.0004451  25.439   <2e-16 ***
#Temp_DegC   0.0011947  0.0004861   2.458   0.0188 *

apr22.div.glm.fit4<-glm(formula = AveShanDiv ~ Dissolved_OrganicMatter_RFU, data=apr22.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(apr22.div.glm.fit4)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)                 0.0112493  0.0004708  23.896   <2e-16 ***
#Dissolved_OrganicMatter_RFU 0.0004269  0.0004659   0.916    0.365

apr22.div.glm.fit5<-glm(formula = AveShanDiv ~ Sulfate_milliM, data=apr22.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(apr22.div.glm.fit5)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)     0.0112325  0.0004772  23.536   <2e-16 ***
#Sulfate_milliM -0.0002664  0.0004893  -0.545    0.589

apr22.div.glm.fit6<-glm(formula = AveShanDiv ~ Sulfide_microM, data=apr22.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(apr22.div.glm.fit6)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)    0.0112379  0.0004723   23.80   <2e-16 ***
#Sulfide_microM 0.0002944  0.0005160    0.57    0.572

fit1<-aov(AveShanDiv ~ as.factor(Depth_m), data=apr22.div)
# ANOVA is basically a regression but w/ categorical variables more info here https://www.statology.org/anova-vs-regression/
#pairwise.adonis(apr22.div$AveShanDiv, apr22.div$Depth_m, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(fit1)
#Df           Sum Sq Mean Sq    F value   Pr(>F)
#Depth_m      7   4097   585.3   1.114   0.38
#Residuals   31  16294   525.6
Tuk1<-TukeyHSD(fit1)
Tuk1$Depth_m

#plot(AveShanDiv ~ Depth_m, data=apr22.div)
#abline(aov(DustComplexity ~ Elevation, data=apr22.div))

# fit.0<-aov(DustComplexity ~ as.factor(Elevation), data=apr22.div)
# summary(fit.0)
# TukeyHSD(fit.0)
# Levene's test with one independent variable
## Levene's tests whether variances of 2 samples are equal
## we want variances to be the same -- want NON SIGNIFICANCE!
## t test assumes that variances are the same, so Levene's test needs to be non significant
## Fligner's test is a Levene's test for data that are not normally distributed
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(AveShanDiv ~ Depth_m, data = apr22.div)
# Fligner-Killeen:med chi-squared = 4.091, df = 7, p-value = 0.7692
# Which shows that the data do not deviate significantly from homogeneity.
compare_means(AveShanDiv ~ Depth_m, data=apr22.div, method="anova",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input



#### April - Species Richness ####
# BDC
apr22.sr.glm.fit1<-glm(formula = Bac_Species_Richness ~ DO_Percent_Local, data=apr22.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(apr22.sr.glm.fit1)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)       0.0112706  0.0004664  24.166   <2e-16 ***
#DO_Percent_Local -0.0009547  0.0005092  -1.875   0.0687 .

apr22.sr.glm.fit2<-glm(formula = Bac_Species_Richness ~ ORP_mV, data=apr22.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(apr22.sr.glm.fit2)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)   0.0112343  0.0004731  23.745   <2e-16 ***
#ORP_mV      -0.0001512  0.0004968  -0.304    0.763

apr22.sr.glm.fit3<-glm(formula = Bac_Species_Richness ~ Temp_DegC, data=apr22.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(apr22.sr.glm.fit3)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept) 0.0113225  0.0004451  25.439   <2e-16 ***
#Temp_DegC   0.0011947  0.0004861   2.458   0.0188 *

apr22.sr.glm.fit4<-glm(formula = Bac_Species_Richness ~ Dissolved_OrganicMatter_RFU, data=apr22.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(apr22.sr.glm.fit4)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)                 0.0112493  0.0004708  23.896   <2e-16 ***
#Dissolved_OrganicMatter_RFU 0.0004269  0.0004659   0.916    0.365

apr22.sr.glm.fit5<-glm(formula = Bac_Species_Richness ~ Sulfate_milliM, data=apr22.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(apr22.sr.glm.fit5)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)     0.0112325  0.0004772  23.536   <2e-16 ***
#Sulfate_milliM -0.0002664  0.0004893  -0.545    0.589

apr22.sr.glm.fit6<-glm(formula = Bac_Species_Richness ~ Sulfide_microM, data=apr22.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(apr22.sr.glm.fit6)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)    0.0112379  0.0004723   23.80   <2e-16 ***
#Sulfide_microM 0.0002944  0.0005160    0.57    0.572

fit1<-aov(Bac_Species_Richness ~ as.factor(Depth_m), data=apr22.div)
# ANOVA is basically a regression but w/ categorical variables more info here https://www.statology.org/anova-vs-regression/
#pairwise.adonis(apr22.div$Bac_Species_Richness, apr22.div$Depth_m, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(fit1)
#Df           Sum Sq Mean Sq    F value   Pr(>F)
#Depth_m      7   4097   585.3   1.114   0.38
#Residuals   31  16294   525.6
Tuk1<-TukeyHSD(fit1)
Tuk1$Depth_m

#plot(Bac_Species_Richness ~ Depth_m, data=apr22.div)
#abline(aov(DustComplexity ~ Elevation, data=apr22.div))

# fit.0<-aov(DustComplexity ~ as.factor(Elevation), data=apr22.div)
# summary(fit.0)
# TukeyHSD(fit.0)
# Levene's test with one independent variable
## Levene's tests whether variances of 2 samples are equal
## we want variances to be the same -- want NON SIGNIFICANCE!
## t test assumes that variances are the same, so Levene's test needs to be non significant
## Fligner's test is a Levene's test for data that are not normally distributed
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(Bac_Species_Richness ~ Depth_m, data = apr22.div)
# Fligner-Killeen:med chi-squared = 4.091, df = 7, p-value = 0.7692
# Which shows that the data do not deviate significantly from homogeneity.
compare_means(Bac_Species_Richness ~ Depth_m, data=apr22.div, method="anova",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input

#### Visualize Richness, Diversity vs Env Variables ####

## Shannon Diversity & Environmental Variables
# note: R (correlation coefficient) vs R^2 (coefficient of determination): https://towardsdatascience.com/r%C2%B2-or-r%C2%B2-when-to-use-what-4968eee68ed3

ggplot(bac.div.metadat.rar, aes(x = DO_Percent_Local, y = AveShanDiv)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampDate), size=3) + theme_classic() +
  stat_smooth(method = "glm", col = "black", se=FALSE, size=1)+ labs(title="Dissolved Oxygen x 16S Shannon Diversity", color="Depth (m)")+ylab("Shannon Diversity")+xlab("Dissolved Oxygen (%)")+
  scale_colour_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) +
  scale_shape_discrete(labels=c("June 2021","WI","DP","BDC"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  stat_cor(label.y = 150, label.x=3) +
  stat_regline_equation(aes(label=paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),label.y = 160,label.x=3)

ggplot(bac.div.metadat.rar, aes(x = DO_Percent_Local, y = AveShanDiv)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampDate), size=3) + theme_classic() +
  labs(title="Dissolved Oxygen x 16S Shannon Diversity", color="Depth (m)")+ylab("Shannon Diversity")+xlab("Dissolved Oxygen (%)")+
  scale_colour_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) +
  scale_shape_discrete(labels=c("June 2021","WI","DP","BDC"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))

ggplot(bac.div.metadat.rar, aes(x = ORP_mV, y = AveShanDiv)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampDate), size=3) + theme_classic() +
  stat_smooth(method = "glm", col = "black", se=FALSE, size=1)+ labs(title="Oxidation-Reduction Potential x 16S Shannon Diversity", color="Depth (m)")+ylab("Shannon Diversity")+xlab("Redox Potential (mV)")+
  scale_colour_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) +
  scale_shape_discrete(labels=c("June 2021","WI","DP","BDC"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  stat_cor(label.y = 3, label.x=1) +
  stat_regline_equation(aes(label=paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),label.y = 3.2,label.x=1)

ggplot(bac.div.metadat.rar, aes(x = Temp_DegC, y = AveShanDiv)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampDate), size=3) + theme_classic() + labs(title="Temperature x 16S Shannon Diversity", color="Depth (m)")+ylab("16S Shannon Diversity")+xlab("Temperature (C)")+
  scale_colour_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) +
  scale_shape_discrete(labels=c("June 2021","WI","DP","BDC"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))

## Species Richness & Environmental Variables
# note: R (correlation coefficient) vs R^2 (coefficient of determination): https://towardsdatascience.com/r%C2%B2-or-r%C2%B2-when-to-use-what-4968eee68ed3

ggplot(bac.div.metadat.rar, aes(x = DO_Percent_Local, y = Bac_Species_Richness)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampDate), size=3) + theme_classic() +
  stat_smooth(method = "glm", col = "black", se=FALSE, size=1)+ labs(title="Dissolved Oxygen x 16S Species Richness", color="Depth (m)")+ylab("Species Richness")+xlab("Dissolved Oxygen (%)")+
  scale_colour_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) +
  scale_shape_discrete(labels=c("June 2021","WI","DP","BDC"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  stat_cor(label.y = 3, label.x=1) +
  stat_regline_equation(aes(label=paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),label.y = 3.1,label.x=1)

ggplot(bac.div.metadat.rar, aes(x = ORP_mV, y = Bac_Species_Richness)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampDate), size=3) + theme_classic() +
  stat_smooth(method = "glm", col = "black", se=FALSE, size=1)+ labs(title="Oxidation-Reduction Potential x 16S Species Richness", color="Depth (m)")+ylab("Species Richness")+xlab("Redox Potential (mV)")+
  scale_colour_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) +
  scale_shape_discrete(labels=c("June 2021","WI","DP","BDC"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  stat_cor(label.y = 3, label.x=1) +
  stat_regline_equation(aes(label=paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),label.y = 3.2,label.x=1)

ggplot(bac.div.metadat.rar, aes(x = Temp_DegC, y = Bac_Species_Richness)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampDate), size=3) + theme_classic() + labs(title="Temperature x 16S Species Richness", color="Depth (m)")+ylab("16S Species Richness")+xlab("Temperature (C)")+
  scale_colour_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) +
  scale_shape_discrete(labels=c("June 2021","WI","DP","BDC"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))


#### Save Everything ####
save.image("data/SSeaDust_AlphaDiv_Data_Rarefied.Rdata")
