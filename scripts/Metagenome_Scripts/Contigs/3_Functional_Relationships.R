#### Set WD & Load Libraries ####
getwd() # use setwd("path/to/files") if you are not in the right directory
#setwd("/bigdata/aronsonlab/shared/SaltonSea/Metagenomes/SeqCenter_3.30.2023/MGM_Analyses/Contigs")
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
  library(htmlwidgets)
  #library(heatmaply)
  #library(rgl)
  library(plotly)
  library(DESeq2)
  library(dplyr)
  library(magrittr)
  library(MASS)
  library(dendextend)
  library(plotly)
  library(tidyr)
  library(reshape)
  library(reshape2)
  #library(wesanderson)
  #library(nationalparkcolors)
  library(fitdistrplus)
  library(logspline)
  library(shades)
  #library(ALDEx2)
  library(rstatix)
  library(devtools)
  library(pairwiseAdonis)
  library(ggcorrplot)
  library(ggpmisc)
  library(lmtest)
  library(car)
})

#### Load Data & See Info About Data ####
#load("data/SSD_MGM_Fxn_BetaDiv.Rdata") # load Rdata to global env
load("data/Metagenomes/SSD_mgm_Fxn_GLMs_Ready.Rdata") # load Rdata to global env

head(meta.all.scaled)
arsen.fxns[1:4,]

# ABOUT THE DATA:
# Before normalizations (i.e., VST, MR, etc) were done, the following was performed:
# contigs were co-assembled, so genes are found in same assembly (same set of contigs)
# non-normalized reads were mapped to genes in contigs (reads mapped by metagenome aka sample)
# featureCounts counted reads that mapped to genes on contigs
# Reads mapped to genes were divided by gene length for all genes across all samples because same KO can be assigned to multiple genes
# Gene coverage was then added together for each KO ID, since multiple genes were assigned the same KO ID
# Summed coverage per KO was then transformed and/or normalized via median-ratio normalization, vst, and MR

# Notes about which objects are which..
mgm_genes.cov.clean # has all read counts, coverages, scaled coverages, and scaled up coverages by gene
mgm_gene.cov_table[,1:4] # uses scaled up, relative coverage per genes (which were also scaled by deployment duration before being scaled up)
mgm_fxn.cov_table[,1:4] # summed scaled up relative coverages per KO
mgm.mr[1:4,1:4] # median of ratio normalized summed KO coverages

## For pathway analyses -- after gene coverage was calculated and added together per KO ID, they were added together for each pathway
## summed coverages per KO ID, then per pathway were transformed by MR

#### Load Functions for Generating Plots ####

# scatter plot function with preset df, y var, color var, shape var
scat.plot.fxn<-function(df,x_var,y_var){
  scat.plt<-ggplot(df,aes(x=x_var, y=y_var, col=SampDate,shape=Site))+
    geom_jitter(aes(color=factor(SampDate),shape=Site), size=5, width=0.15, height=0) +
    scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                       labels=c("July 2020","August 2020","November 2020",
                                "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
    scale_shape_manual(name="Site",values = c(0,1,16,15)) +
    theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
          axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
    labs(x = "",
         y = "",
         color = "R",
         title = paste(as.character(y_var),"~",as.character(x_var),sep=" "))

  return(scat.plt)
}

date.mgm.scat.plot.fxn<-function(df, y_var){
  date.scat.plt<-ggplot(df,aes(x=SampDate, y=y_var, col=SampDate,shape=Site))+
    geom_jitter(aes(color=factor(SampDate),shape=Site), size=5, width=0.15, height=0) +
    scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                       labels=c("July 2020","August 2020","November 2020",
                                "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
    scale_shape_manual(name="Site",values = c(0,1,16,15)) +
    theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
          axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
    labs(x = "",
         y = "",
         shape="",
         color = "R")

  return(date.scat.plt)
}

#### Functional Beta Diversity - MRN data ####
# MR = median-ratio normalization
mgm.mr[1:4,1:4] # sample IDs are rows, genes (KOs) are columns
mgm_fxn.cov_table[1:4,1:4] # sanity check --> KOs with low coverage are still in this df

# check rownames of MR & Mr transformed feature count data & metadata
rownames(mgm.mr) %in% rownames(meta.all.scaled)

## PCOA with VST transformed data first
# calculate our Euclidean distance matrix using VST data
mgm.euc_dist.mr <- dist(mgm.mr, method = "euclidean")

# creating our hierarcical clustering dendrogram
mgm.euc.mr_clust <- hclust(mgm.euc_dist.mr, method="ward.D2")

# let's make it a little nicer...
mgm.euc.mr_dend <- as.dendrogram(mgm.euc.mr_clust, hang=0.2)
mgm.dend_cols <- as.character(meta.all.scaled$SampDate_Color[order.dendrogram(mgm.euc.mr_dend)])
labels_colors(mgm.euc.mr_dend) <- mgm.dend_cols

par(mar=c(2,2,2,2))
plot(mgm.euc.mr_dend, ylab="Median-Ratio Normalized, Euclidean Distance",cex = 0.5,horiz=TRUE) + title(main = "Bacteria Functional Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 3)
#legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c( "#ef781c","#03045e","#059c3f"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# let's use our Euclidean distance matrix from before
mgm.pcoa.mr <- pcoa(mgm.euc_dist.mr) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
##save.image("data/ssd_mr.euc.dist1_3.7.23.Rdata")

# The proportion of variances explained is in its element values$Relative_eig
mgm.pcoa.mr$values
#     Eigenvalues Relative_eig Broken_stick Cumul_eig Cumul_br_stick
# 1  10012830.6   0.50406277   0.16236050 0.5040628      0.1623605
# 2   5190786.9   0.54271296   0.11888224 0.7653757      0.2812427

# extract principal coordinates
mgm.pcoa.mr.vectors<-data.frame(mgm.pcoa.mr$vectors)
mgm.pcoa.mr.vectors$SampleID<-rownames(mgm.pcoa.mr$vectors)

# merge pcoa coordinates w/ metadata
mgm.pcoa.mr.meta<-merge(mgm.pcoa.mr.vectors, meta.all.scaled, by.x="SampleID", by.y="SampleID")
mgm.pcoa.mr.meta$SampleMonth
mgm.pcoa.mr.meta$SampDate

head(mgm.pcoa.mr.meta)

head(mgm.pcoa.mr$values) # pull out Relative (Relative_eig) variation % to add to axes labels

#### PERMANOVAs to Env Variables Across Groups ####

## The currently preferred analysis for evaluating differences among groups is PERMANOVA.
## This analysis partitions sums of squares using dissimilarities,
##  evaluating differences in the centroids of groups in multivariate space.
##  The vegan functions “adonis” and “adonis2” are used to compute PERMANOVA in R.

help(adonis)

## can specify dataframes for analysis, or we can alternatively specify a dissimilarity matrix:

#Other advantages of using PERMANOVA are that we can test for interactions between predictor variables,
## and we can use both categorical and continuous predictor variables.
## An advantage of adonis2 is that we can test for overall model fit, setting by=NULL, or by individual terms (w/ by="terms")
## w/ distance matrices - The adonis2 tests are identical to anova.cca of dbrda. With Euclidean distances, the tests are also identical to anova.cca of rda.

# drop SampleID column from mgm.mr
mgm.mr<-mgm.mr[,!names(mgm.mr) %in% c("SampleID")]
mgm.mr[1:4,1:4]

# First make sure your data frames you're comparing are in the same exact order!!
rownames(mgm.mr) %in% rownames(meta.all.scaled)

meta.all.scaled=meta.all.scaled[rownames(mgm.mr),] ## reorder metadata to match order of MR data
perm <- with(meta.all.scaled, how(nperm = 1000))

pnova0<-adonis2(mgm.mr ~ Site,data=meta.all.scaled,method = "euclidean",by=NULL,permutations=perm)
pnova0
# Df SumOfSqs      R2     F Pr(>F)
# Site      3   747070 0.18647 1.528 0.1179
# Residual 20  3259390 0.81353
# Total    23  4006460 1.00000

### SELF REMINDER FOR R^2
### Coefficient of Determination, denoted R2 or r2
### is the proportion of the variance in the dependent variable that is predictable from the independent variable(s)

### Pseudo F stat for PERMANOVA
### pseudo F-ratio: It compares the total sum of squared dissimilarities (or ranked dissimilarities) among objects belonging to different groups to that of objects belonging to the same group.
### Larger F-ratios indicate more pronounced group separation, however, the significance of this ratio is usually of more interest than its magnitude.

#### Do LPS Biosynthesis Fxns of Interest Correlate? ####

mr.cov.sum.lps.ko[1:4,]
grep("wecP;", colnames(mr.cov.sum.lps.ko)) # find index for column we want
grep("wbgP;", colnames(mr.cov.sum.lps.ko)) # find index for column we want

lpsOI.cor<-cor(mr.cov.sum.lps.ko[,c(65,76)])
lps.cor.p <- cor_pmat(mr.cov.sum.lps.ko[,c(65,76)])

ggcorrplot(lpsOI.cor,method="square",lab = T,p.mat=lps.cor.p,sig.level = 0.05,
           hc.order=TRUE,outline.color="white",type = "lower")

par(mar=c(1,1,1,1))
plot(mr.cov.sum.lps.ko$`wecP; UDP-GalNAc:undecaprenyl-phosphate GalNAc-1-phosphate transferase [EC:2.7.8.40]`,mr.cov.sum.lps.ko$`wbgP; O55-antigen biosynthesis glycosyltransferase`)

#### Create LPS Coverage Table, & Merge with Metadata ####

head(mr.lps.ko)

## create table of Sample x KO_ID
mr.lps.ko_table<-as.data.frame(dcast(mr.lps.ko, SampleID~KO_Function.KEGG, value.var="MR_SumCovPerKO"))
mr.lps.ko_table[1:5,1:5]

## merge LPS functions with scaled metadata

mr.lps.ko.all<-merge(meta.all.scaled,mr.lps.ko_table,by="SampleID")
mr.lps.ko.all[1:5,1:20]


#### LPS Fxns & Env Vars - Scatterplots ####

# prep env data for downstream visualizations/analyses
env.vars<-mr.lps.ko.all[,c(6,9:10,12:14,36:45)] # subset env var data
env.vars.names<-names(mr.lps.ko.all[,c(6,9:10,12:14,36:45)]) # create vector list of variable names

## create function name variable
lps.fxn.names<-names(mr.lps.ko.all)[names(mr.lps.ko.all) %in% c("wecP; UDP-GalNAc:undecaprenyl-phosphate GalNAc-1-phosphate transferase [EC:2.7.8.40]",
                                                                "wbgP; O55-antigen biosynthesis glycosyltransferase")] # pull out names of columns in df that contain "Axis" in name

# first let's do some basic plotting and see what looks important...
# reminder:

ggplot(data=mr.lps.ko.all,aes(x=ave.air_temp, y=`wecP; UDP-GalNAc:undecaprenyl-phosphate GalNAc-1-phosphate transferase [EC:2.7.8.40]`, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=3, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Ave Air Temp", y="wecP KO Coverage", title="wecP Coverage x Average Air Temp",subtitle="Using Scaled Climate Data")

ggplot(data=mr.lps.ko.all,aes(x=ave.air_temp, y=`wbgP; O55-antigen biosynthesis glycosyltransferase`, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=3, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Ave Air Temp", y="wbgP KO Coverage", title="wbgP Coverage x Average Air Temp",subtitle="Using Scaled Climate Data")

## Loop to Generate Scatter Plots
### comparing y ~ x where y = lps.fxnilia relative coverage and x = whatever environmental variable
head(mr.lps.ko.all)
lps.scatter.plot.list<-list() # create empty list for each plot to be stored in
lps.fxn.names
env.vars.names # sanity check
plotnum <- 1 # counting our model numbers for indexes purposes in the loop

for (i in lps.fxn.names) {
  for (j in env.vars.names){
    lps.fxn.scat.plot=scat.plot.fxn(mr.lps.ko.all,mr.lps.ko.all[,j],mr.lps.ko.all[,i])
    plot.titled = lps.fxn.scat.plot + ggtitle(paste(as.character(i),"Relative Coverage ~",as.character(j),sep=" "),subtitle="Using Scaled Climate Data")
    lps.scatter.plot.list[[plotnum]]=plot(plot.titled)
    ggsave(plot.titled,filename = paste("figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/LPS_Biosynthesis/SSD_",as.character(i),"~",as.character(j),"_scatterplot.png",sep=""), width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).
    #print(names(clsOI.meta.all[i]))
    #print(clsOI.meta.all[j])
    plotnum<-1+plotnum
  }
}


# check if our loop + functions worked!
lps.scatter.plot.list[[2]]

ggplot(data=mr.lps.ko.all,aes(x=ave.air_temp, y=`wbgP; O55-antigen biosynthesis glycosyltransferase`, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Ave Air Temp", y="wbgP KO Coverage", title="wbgP Coverage x Average Air Temp",subtitle="Using Scaled Climate Data")

lps.scatter.plot.list[[27]]

ggplot(data=mr.lps.ko.all,aes(x=ave.air_temp, y=`wecP; UDP-GalNAc:undecaprenyl-phosphate GalNAc-1-phosphate transferase [EC:2.7.8.40]`, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Ave Air Temp", y="wecP KO Coverage", title="wecP Coverage x Average Air Temp",subtitle="Using Scaled Climate Data")


# now visualize genes by sampdate
lps.date.plot.list<-list() # create empty list for each plot to be stored in

lps.fxn.names
env.vars.names
plotnum <- 1 # counting our model numbers for indexes purposes in the loop

for (i in lps.fxn.names) {
  lps.date.scat.plot=date.mgm.scat.plot.fxn(mr.lps.ko.all,mr.lps.ko.all[,i])
  date.plot.titled = lps.date.scat.plot + ggtitle(paste(as.character(i),"Relative Coverage by Collection Date",sep=" "))
  lps.date.plot.list[[plotnum]]=plot(date.plot.titled)
  ggsave(date.plot.titled,filename = paste("figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/LPS_Biosynthesis/SSD_",as.character(i),"_by_SampDate_scatterplot.png",sep=""), width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

  #print(names(clsOI.meta.all[i]))
  #print(clsOI.meta.all[j])
  plotnum<-1+plotnum

}

lps.date.plot.list[[1]]
ggplot(mr.lps.ko.all,aes(x=SampDate, y=`wbgP; O55-antigen biosynthesis glycosyltransferase`, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x = "",
       y = "",
       shape="",
       color = "R")

lps.date.plot.list[[2]]
ggplot(mr.lps.ko.all,aes(x=SampDate, y=`wecP; UDP-GalNAc:undecaprenyl-phosphate GalNAc-1-phosphate transferase [EC:2.7.8.40]`, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x = "",
       y = "",
       shape="",
       color = "R")

#### What Env Vars Predict Distribution of LPS Biosynthesis Genes ####
head(mr.lps.ko.all)

# first with wecP
# is wecP normally distributed?
shapiro.test(mr.lps.ko.all$`wecP; UDP-GalNAc:undecaprenyl-phosphate GalNAc-1-phosphate transferase [EC:2.7.8.40]`) # what is the p-value?
# W = 0.97668, p-value = 0.8277, suggests that this function is normally distributed...
hist(mr.lps.ko.all$`wecP; UDP-GalNAc:undecaprenyl-phosphate GalNAc-1-phosphate transferase [EC:2.7.8.40]`)

# visualize Q-Q plot for alpha div
# The Q-Q plot, or quantile-quantile plot, is a graphical tool to help us assess if a set of data plausibly came from some theoretical distribution such as a normal or exponential.
# For example, if we run a statistical analysis that assumes our residuals are normally distributed, we can use a normal Q-Q plot to check that assumption
# more on Q-Q plots here: https://data.library.virginia.edu/understanding-q-q-plots/
# more here too: https://grodri.github.io/glms/notes/c2s9#:~:text=8%20The%20Q%2DQ%20Plot,versus%20quantiles%20of%20a%20distribution.
qqnorm(mr.lps.ko.all$`wecP; UDP-GalNAc:undecaprenyl-phosphate GalNAc-1-phosphate transferase [EC:2.7.8.40]`, pch = 1, frame = FALSE)
qqline(mr.lps.ko.all$`wecP; UDP-GalNAc:undecaprenyl-phosphate GalNAc-1-phosphate transferase [EC:2.7.8.40]`, col = "red", lwd = 2)

# what column is wecP?
grep("wecP;", colnames(mr.lps.ko.all)) # find index for column we want

# next let's run some models
lps.glm.step1<-step(glm(formula = round(`wecP; UDP-GalNAc:undecaprenyl-phosphate GalNAc-1-phosphate transferase [EC:2.7.8.40]`) ~ ., data=mr.lps.ko.all[,c(6,9:10,12:14,36:45,123)],family="gaussian"))
summary(lps.glm.step1)

summary(glm(formula = `wecP; UDP-GalNAc:undecaprenyl-phosphate GalNAc-1-phosphate transferase [EC:2.7.8.40]` ~ ave.air_temp+ave.u_E.W.wind, data=mr.lps.ko.all,family="gaussian"))
# ^ seems to be the best model with the lowest AIC

# final model:
lps.glm1<-lm(formula = `wecP; UDP-GalNAc:undecaprenyl-phosphate GalNAc-1-phosphate transferase [EC:2.7.8.40]` ~ ave.air_temp+ave.u_E.W.wind, data=mr.lps.ko.all)
lps.glm1.sum<-summary(lps.glm1, corr = TRUE)
lps.glm1.sum
lps.glm1.sum$correlation # true correlation between the predictor variables

par(mar=c(5, 4, 4, 2) + 0.1) # change plot margins
avPlots(lps.glm1, id=FALSE, pt.wts=FALSE)

# confirm that this is the best model by plotting the residuals and seeing their distribution
par(mar=c(1,1,1,1)) # change plot margins
plot(fitted(lps.glm1), lps.glm1$residuals) # produce residual vs. fitted plot
abline(0,0) # add a horizontal line at 0

# Q-Q norm plot of residuals
qqnorm(lps.glm1$residuals)
qqline(lps.glm1$residuals)

# density of residuals (how many points have certain residuals)
plot(density(lps.glm1$residuals))

# a better way to plot the residuals and their densities
ggplot(data=mr.lps.ko.all,aes(lps.glm1$residuals)) +
  geom_histogram(aes(y=..density..),bindwidth=5,color='black',fill='grey') +
  geom_density(alpha=.2,fill='lightpink')

#plot(lps.glm1) # another way to plot residuals

# now let's plot the model!
wecP.plot1<-ggplot(data=mr.lps.ko.all,aes(x=ave.air_temp, y=`wecP; UDP-GalNAc:undecaprenyl-phosphate GalNAc-1-phosphate transferase [EC:2.7.8.40]`))+
  stat_smooth(aes(x=ave.air_temp, y=`wecP; UDP-GalNAc:undecaprenyl-phosphate GalNAc-1-phosphate transferase [EC:2.7.8.40]`), method=stats::glm,method.args = list(family = "gaussian")) +
  geom_jitter(aes(shape=Site,color=SampDate),size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Ave Air Temp", y="wecP KO Coverage", title="wecP Coverage x Ave Air Temp",subtitle="Using Scaled Climate Data") +
  stat_poly_eq(use_label(c("adj.R2", "P")), formula=y~x,small.p = TRUE)

ggsave(wecP.plot1,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/LPS_Biosynthesis/SSD_wecp_AveAirTemp_scatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

wecP.plot2<-ggplot(data=mr.lps.ko.all,aes(x=ave.u_E.W.wind, y=`wecP; UDP-GalNAc:undecaprenyl-phosphate GalNAc-1-phosphate transferase [EC:2.7.8.40]`))+
  stat_smooth(aes(x=Others, y=`wecP; UDP-GalNAc:undecaprenyl-phosphate GalNAc-1-phosphate transferase [EC:2.7.8.40]`), method=stats::glm,method.args = list(family = "gaussian")) +
  geom_jitter(aes(shape=Site,color=SampDate),size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="West to East Wind (u)", y="wecP KO Coverage", title="wecP Coverage x West to East Winds",subtitle="Using Scaled Climate Data") +
  stat_poly_eq(use_label(c("adj.R2", "P")), formula=y~x,small.p = TRUE)

ggsave(wecP.plot2,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/LPS_Biosynthesis/SSD_wecp_WindVector_u_scatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

# next with wbgP
# is wbgP normally distributed?
shapiro.test(mr.lps.ko.all$`wbgP; O55-antigen biosynthesis glycosyltransferase`) # what is the p-value?
# W = 0.79087, p-value = 0.0002097, suggests that this function is NOT normally distributed...
hist(mr.lps.ko.all$`wbgP; O55-antigen biosynthesis glycosyltransferase`)

# visualize Q-Q plot for alpha div
# The Q-Q plot, or quantile-quantile plot, is a graphical tool to help us assess if a set of data plausibly came from some theoretical distribution such as a normal or exponential.
# For example, if we run a statistical analysis that assumes our residuals are normally distributed, we can use a normal Q-Q plot to check that assumption
# more on Q-Q plots here: https://data.library.virginia.edu/understanding-q-q-plots/
# more here too: https://grodri.github.io/glms/notes/c2s9#:~:text=8%20The%20Q%2DQ%20Plot,versus%20quantiles%20of%20a%20distribution.
qqnorm(mr.lps.ko.all$`wbgP; O55-antigen biosynthesis glycosyltransferase`, pch = 1, frame = FALSE)
qqline(mr.lps.ko.all$`wbgP; O55-antigen biosynthesis glycosyltransferase`, col = "red", lwd = 2)

grep("wbgP;", colnames(mr.lps.ko.all)) # find index for column we want

# not sure if we should use a Poisson distribution or a negative binomial distribution for our GLM...
# so we are going to run both models, compare the results, and then calculate the log likelihood to see which is the better model
# for more info on the log likelihood, watch this: https://www.youtube.com/watch?v=8nogLkirA3I
## basically we are trying to find out which distribution fits the data (poisson or negative binomial in this case)
# more on why we use log likelihood and not just the likelihood here: https://math.stackexchange.com/questions/892832/why-we-consider-log-likelihood-instead-of-likelihood-in-gaussian-distribution
# and here: https://www.reddit.com/r/learnmachinelearning/comments/vewczb/why_do_we_maximize_loglikelihood_instead_of/
## likelihoods are a product, & usually you're mulitplying small values together which makes them even smaller and harder to measure
## log likelihoods turn these products into sums, and help to scale our small values without losing any information

# first a glm w/ Poisson distribution
lps.glm.step1<-step(glm(formula = round(`wbgP; O55-antigen biosynthesis glycosyltransferase`) ~ ., data=mr.lps.ko.all[,c(6,9:10,12:14,36:45,112)],family="poisson"))
summary(lps.glm.step1)
lps.glm.step1$anova

# then a  NB glm...
lps.glm.step2<-step(glm.nb(formula = round(`wbgP; O55-antigen biosynthesis glycosyltransferase`) ~ ., data=mr.lps.ko.all[,c(6,9:10,12:14,36:45,112)]))
summary(lps.glm.step2)
lps.glm.step2$anova

# both show same sig variables, and it seems Salton Sea & Shrub are the two best vars to explain distribution of wbgP
# so now let's re run the models and do some comparisons
lps.nb<-glm.nb(formula = round(`wbgP; O55-antigen biosynthesis glycosyltransferase`) ~ SaltonSea*Shrub, data=mr.lps.ko.all)
lps.nb.sum<-summary(lps.nb)
lps.nb.sum

par(mar=c(1,1,1,1)) # change plot margins
avPlots(lps.nb, id=FALSE, pt.wts=FALSE)

lps.pois<-glm(formula = round(`wbgP; O55-antigen biosynthesis glycosyltransferase`) ~ SaltonSea*Shrub, data=mr.lps.ko.all, family="poisson")
lps.pois.sum<-summary(lps.pois)
lps.pois.sum

# now we compare the loglikehoods * LRT only works for nested models
# In summary, the LRT tells us if it is beneficial to add parameters to our model, or if we should stick with our simpler model.
# In their most basic form, the hypotheses for the LRT are:
# H0: You should use the nested model.
# Ha: You should use the complex model.
# Thus, if you reject the H0, you can conclude that the complex model is significantly more accurate than the nested model, and you would choose to use the complex model.
# If you fail to reject the H0, you can conclude that the complex model is NOT significantly more accurate than the nested model, so you would choose to use the nested model instead.

## the higher the value of loglikehood, the better the model fits the data
# the smaller the AIC the better it's predictive power is
logLik(lps.nb) # negative binomial
AIC(lps.nb)

logLik(lps.pois) # poisson
AIC(lps.nb,lps.pois)

# lastly let's try a likelihood ratio test to compare
lrtest(lps.nb,lps.pois) # lps.pois is significantly better model based on likelihood ratio test...

# final model for wbgP:
summary(glm(formula = round(`wbgP; O55-antigen biosynthesis glycosyltransferase`) ~ SaltonSea*Shrub, data=mr.lps.ko.all, family="poisson"), corr = TRUE)
lps.pois.sum<-summary(glm(formula = round(`wbgP; O55-antigen biosynthesis glycosyltransferase`) ~ SaltonSea*Shrub, data=mr.lps.ko.all, family="poisson"), corr = TRUE)
lps.pois.sum
lps.pois.sum$correlation # true correlation between the predictor variables


# calculate McFadden's Pseudo R-squared for model for GLM aka Analysis of Deviance
# more here: https://bookdown.org/egarpor/PM-UC3M/glm-deviance.html
with(summary(lps.pois), 1 - deviance/null.deviance)
# 0.4629555

# confirm that this is the best model by plotting the residuals and seeing their distribution
par(mar=c(1,1,1,1)) # change plot margins
plot(fitted(lps.pois), lps.pois$residuals) # produce residual vs. fitted plot
abline(0,0) # add a horizontal line at 0

# Q-Q norm plot of residuals
qqnorm(lps.pois$residuals)
qqline(lps.pois$residuals)

# density of residuals (how many points have certain residuals)
plot(density(lps.pois$residuals))

# a better way to plot the residuals and their densities
ggplot(data=mr.lps.ko.all,aes(lps.pois$residuals)) +
  geom_histogram(aes(y=..density..),bindwidth=30,color='black',fill='grey') +
  geom_density(alpha=.2,fill='lightpink')

#plot(lps.pois) # another way to plot residuals

# now let's plot the model!

ggplot(data=mr.lps.ko.all,aes(x=SaltonSea, y=`wbgP; O55-antigen biosynthesis glycosyltransferase`))+
  geom_smooth(aes(x=SaltonSea, y=`wbgP; O55-antigen biosynthesis glycosyltransferase`), method=stats::glm,method.args = list(family = "poisson")) +
  geom_jitter(aes(shape=Site,color=SampDate),size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Salton Sea STF", y="wbgP KO Coverage", title="wbgP Coverage x Salton Sea STF",subtitle="Using Scaled Climate Data") +
  stat_fit_tidy(method = "glm",method.args = list(family=poisson,formula=y~x),
                mapping = aes(label = sprintf("p = %.3g",after_stat(x_p.value))))

#### Compare Variance of LPS Biosynthesis/Modification Genes Across Sites ####

# wecP is normally distributed; used rarefied counts to calculate ShanDiv
# use the following statisitcal tests for variance comparisons
## ANOVA: are variances significantly different between groups
## Tukey test: which groups' variances are significant different from one another
## Levene's test: is variance homogenous aka equal across samples?

head(mr.lps.ko.all)

fit1<-aov(`wecP; UDP-GalNAc:undecaprenyl-phosphate GalNAc-1-phosphate transferase [EC:2.7.8.40]`~ Site, data=mr.lps.ko.all)
# ANOVA is basically a regression but w/ categorical variables more info here https://www.statology.org/anova-vs-regression/
#pairwise.adonis(mr.lps.ko.all$`spoIVCA; site-specific DNA recombinase`, mr.lps.ko.all$SampDate, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(fit1)
# p = 0.508

# Tukey test - tells us which groups are significantly different from each other (more here: https://www.r-bloggers.com/2013/06/anova-and-tukeys-test-on-r/)
Tuk1<-TukeyHSD(fit1)
Tuk1$Site
# no sig differences

# fit.0<-aov(DustComplexity ~ as.factor(Elevation), data=mr.lps.ko.all)
# summary(fit.0)
# TukeyHSD(fit.0)

# Levene's test - test for homogeneity of variance
## Levene’s test is an inferential statistic used to check if the variances of a variable obtained for two or more groups are equal or not when data comes from a non-normal distribution
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## t test assumes that variances are the same, so Levene's test needs to be non significant
car::leveneTest(`wecP; UDP-GalNAc:undecaprenyl-phosphate GalNAc-1-phosphate transferase [EC:2.7.8.40]`~ Site, data = mr.lps.ko.all, center=mean)
# Levene's Test for Homogeneity of Variance (center = mean)
#       Df F value Pr(>F)
# group  3  1.1249 0.3627
# 20
# ^ p value is > 0.05 so we cannot reject the Null hypothesis -- variances are equal

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
# fligner.test(`wecP; UDP-GalNAc:undecaprenyl-phosphate GalNAc-1-phosphate transferase [EC:2.7.8.40]`~ SampDate, data = mr.lps.ko.all)
# Fligner-Killeen:med chi-squared = 1.9504, df = 2, p-value = 0.3771
# Which shows that the data do not deviate significantly from homogeneity.

#compare_means(`wecP; UDP-GalNAc:undecaprenyl-phosphate GalNAc-1-phosphate transferase [EC:2.7.8.40]`~ Site, data=mr.lps.ko.all, method="anova",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input


#### Do Sporulation Fxns of Interest Correlate? ####

mr.cov.sum.spor.ko[1:4,]
grep("spoIVCA;", colnames(mr.cov.sum.spor.ko)) # find index for column we want
grep("spmA;", colnames(mr.cov.sum.spor.ko)) # find index for column we want

sporOI.cor<-cor(mr.cov.sum.spor.ko[,c(2,22)])
spor.cor.p <- cor_pmat(mr.cov.sum.spor.ko[,c(2,22)])

ggcorrplot(sporOI.cor,method="square",lab = T,p.mat=spor.cor.p,sig.level = 0.05,
           hc.order=TRUE,outline.color="white",type = "lower")

par(mar=c(1,1,1,1))
plot(mr.cov.sum.spor.ko$`spmA; spore maturation protein A`,mr.cov.sum.spor.ko$`spoIVCA; site-specific DNA recombinase`)

# NOTE: not all of the genes we looked at are required for sporulation, so it's not worrisome necessarily that these genes do not correlate in mgms
## it's also possible that we missed the microbes in their endospore state since they are hard to sequence from

#### Create Sporulation Coverage Table, then Merge with Metadata ####

head(mr.spor.ko)

## create table of Sample x KO_ID
mr.spor.ko_table<-as.data.frame(dcast(mr.spor.ko, SampleID~KO_Function.KEGG, value.var="MR_SumCovPerKO"))
mr.spor.ko_table[1:5,1:5]
rownames(mr.spor.ko_table)<-mr.spor.ko_table$SampleID

## merge sporulation functions with scaled metadata

mr.spor.ko.all<-merge(meta.all.scaled,mr.spor.ko_table,by="SampleID")
mr.spor.ko.all[1:5,1:20]


#### Sporulation Fxns & Env Vars - Scatterplots ####
## heatmaps of traits of interest

# prep env data for downstream visualizations/analyses
env.vars<-mr.spor.ko.all[,c(6,9:10,12:14,36:45)] # subset env var data
env.vars.names<-names(mr.spor.ko.all[,c(6,9:10,12:14,36:45)]) # create vector list of variable names

## create function name variable
spor.fxn.names<-names(mr.spor.ko.all)[names(mr.spor.ko.all) %in% c("spoIVCA; site-specific DNA recombinase",
                                                                "spmA; spore maturation protein A")] # pull out names of columns in df that contain "Axis" in name


## Loop to Generate Scatter Plots
### comparing y ~ x where y = spore relative coverage and x = whatever environmental variable
head(mr.spor.ko.all)
spor.scatter.plot.list<-list() # create empty list for each plot to be stored in
spor.fxn.names
env.vars.names # sanity check
plotnum <- 1 # counting our model numbers for indexes purposes in the loop

for (i in spor.fxn.names) {
  for (j in env.vars.names){
    spor.fxn.scat.plot=scat.plot.fxn(mr.spor.ko.all,mr.spor.ko.all[,j],mr.spor.ko.all[,i])
    plot.titled = spor.fxn.scat.plot + ggtitle(paste(as.character(i),"Relative Coverage ~",as.character(j),sep=" "),subtitle="Using Scaled Climate Data")
    spor.scatter.plot.list[[plotnum]]=plot(plot.titled)
    ggsave(plot.titled,filename = paste("figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/Sporulation/SSD_",as.character(i),"~",as.character(j),"_scatterplot.png",sep=""), width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).
    #print(names(clsOI.meta.all[i]))
    #print(clsOI.meta.all[j])
    plotnum<-1+plotnum
  }
}


# check if our loop + functions worked!
spor.scatter.plot.list[[1]]

ggplot(data=mr.spor.ko.all,aes(x=ave.air_temp, y=`spmA; spore maturation protein A`, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Ave Air Temp", y="spmA KO Coverage", title="spmA Coverage x Average Air Temp",subtitle="Using Scaled Climate Data")

spor.scatter.plot.list[[17]]

ggplot(data=mr.spor.ko.all,aes(x=ave.relative_humidity, y=`spoIVCA; site-specific DNA recombinase`, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Ave Air Temp", y="spoIVCA KO Coverage", title="spoIVCA Coverage x Average Air Temp",subtitle="Using Scaled Climate Data")


# now visualize genes by sampdate
spor.date.plot.list<-list() # create empty list for each plot to be stored in

spor.fxn.names
env.vars.names
plotnum <- 1 # counting our model numbers for indexes purposes in the loop

for (i in spor.fxn.names) {
  spor.date.scat.plot=date.mgm.scat.plot.fxn(mr.spor.ko.all,mr.spor.ko.all[,i])
  date.plot.titled = spor.date.scat.plot + ggtitle(paste(as.character(i),"Relative Coverage by Collection Date",sep=" "))
  spor.date.plot.list[[plotnum]]=plot(date.plot.titled)
  ggsave(date.plot.titled,filename = paste("figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/Sporulation/SSD_",as.character(i),"_by_SampDate_scatterplot.png",sep=""), width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

  #print(names(clsOI.meta.all[i]))
  #print(clsOI.meta.all[j])
  plotnum<-1+plotnum

}

spor.date.plot.list[[1]]
ggplot(mr.spor.ko.all,aes(x=SampDate, y=`spmA; spore maturation protein A`, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x = "",
       y = "",
       shape="",
       color = "R")

spor.date.plot.list[[2]]
ggplot(mr.spor.ko.all,aes(x=SampDate, y=`spoIVCA; site-specific DNA recombinase`, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x = "",
       y = "",
       shape="",
       color = "R")

#### What Env Vars Predict Distribution of Sporulation Genes ####
# first with spoIVCA
# is spoIVCA normally distributed?
shapiro.test(mr.spor.ko.all$`spoIVCA; site-specific DNA recombinase`) # what is the p-value?
# W = 0.95675, p-value = 0.3766, suggests that this function is normally distributed...
hist(mr.spor.ko.all$`spoIVCA; site-specific DNA recombinase`)

# visualize Q-Q plot for alpha div
# The Q-Q plot, or quantile-quantile plot, is a graphical tool to help us assess if a set of data plausibly came from some theoretical distribution such as a normal or exponential.
# For example, if we run a statistical analysis that assumes our residuals are normally distributed, we can use a normal Q-Q plot to check that assumption
# more on Q-Q plots here: https://data.library.virginia.edu/understanding-q-q-plots/
# more here too: https://grodri.github.io/glms/notes/c2s9#:~:text=8%20The%20Q%2DQ%20Plot,versus%20quantiles%20of%20a%20distribution.
qqnorm(mr.spor.ko.all$`spoIVCA; site-specific DNA recombinase`, pch = 1, frame = FALSE)
qqline(mr.spor.ko.all$`spoIVCA; site-specific DNA recombinase`, col = "red", lwd = 2)

grep("spoIVCA;", colnames(mr.spor.ko.all)) # find index for column we want

spor.glm.step1<-step(glm(formula = round(`spoIVCA; site-specific DNA recombinase`) ~ ., data=mr.spor.ko.all[,c(6,9:10,12:14,36:45,69)],family="gaussian"))
summary(spor.glm.step1)
spor.glm.step1$anova

# final model
spor.lm.1<-lm(formula = `spoIVCA; site-specific DNA recombinase` ~ precip_24hr_accum*ave.relative_humidity, data=mr.spor.ko.all)
summary(spor.lm.1, corr = TRUE)
summary(spor.lm.1, corr = TRUE)$correlation # true correlation between the predictor variables


# confirm that this is the best model by plotting the residuals and seeing their distribution
par(mar=c(1,1,1,1)) # change plot margins
plot(fitted(spor.lm.1), spor.lm.1$residuals) # produce residual vs. fitted plot
abline(0,0) # add a horizontal line at 0

# Q-Q norm plot of residuals
qqnorm(spor.lm.1$residuals)
qqline(spor.lm.1$residuals)

# density of residuals (how many points have certain residuals)
plot(density(spor.lm.1$residuals))

# a better way to plot the residuals and their densities
ggplot(data=mr.spor.ko.all,aes(spor.lm.1$residuals)) +
  geom_histogram(aes(y=..density..),bindwidth=30,color='black',fill='grey') +
  geom_density(alpha=.2,fill='lightpink')

#plot(spor.lm.1) # another way to plot residuals

# now let's plot the model!

spoivca.plot1<-ggplot(data=mr.spor.ko.all,aes(x=precip_24hr_accum, y=`spoIVCA; site-specific DNA recombinase`))+
  stat_smooth(aes(x=precip_24hr_accum, y=round(`spoIVCA; site-specific DNA recombinase`)), method=stats::lm) +
  geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Accum. Precipitation (24 hrs)", y="spoIVCA KO Coverage", title="spoIVCA Coverage x Accumulated Precipitation",subtitle="Using Scaled Climate Data") +
  stat_poly_eq(use_label(c("adj.R2", "P")), formula=y~x,small.p = TRUE)

ggsave(spoivca.plot1,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/Sporulation/SSD_spoIVCA_AccumPrecip_scatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

spoivca.plot2<-ggplot(data=mr.spor.ko.all,aes(x=ave.relative_humidity, y=`spoIVCA; site-specific DNA recombinase`))+
  stat_smooth(aes(x=ave.relative_humidity, y=round(`spoIVCA; site-specific DNA recombinase`)), method=stats::lm) +
  geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Ave Relative Humidity", y="spoIVCA KO Coverage", title="spoIVCA Coverage x Ave Relative Humidity",subtitle="Using Scaled Climate Data") +
  stat_poly_eq(use_label(c("adj.R2", "P")), formula=y~x,small.p = TRUE)

ggsave(spoivca.plot2,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/Sporulation/SSD_spoIVCA_AveRelativeHumidity_scatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).


# next with spmA
# is spmA normally distributed?
shapiro.test(mr.spor.ko.all$`spmA; spore maturation protein A`) # what is the p-value?
# W = 0.80544, p-value = 0.000362, suggests that this function is NOT normally distributed...
hist(mr.spor.ko.all$`spmA; spore maturation protein A`)

# visualize Q-Q plot for alpha div
# The Q-Q plot, or quantile-quantile plot, is a graphical tool to help us assess if a set of data plausibly came from some theoretical distribution such as a normal or exponential.
# For example, if we run a statistical analysis that assumes our residuals are normally distributed, we can use a normal Q-Q plot to check that assumption
# more on Q-Q plots here: https://data.library.virginia.edu/understanding-q-q-plots/
# more here too: https://grodri.github.io/glms/notes/c2s9#:~:text=8%20The%20Q%2DQ%20Plot,versus%20quantiles%20of%20a%20distribution.
qqnorm(mr.spor.ko.all$`spmA; spore maturation protein A`, pch = 1, frame = FALSE)
qqline(mr.spor.ko.all$`spmA; spore maturation protein A`, col = "red", lwd = 2)

grep("spmA;", colnames(mr.spor.ko.all)) # find index for column we want

# not sure if we should use a Poisson distribution or a negative binomial distribution for our GLM...
# so we are going to run both models, compare the results, and then calculate the log likelihood to see which is the better model
# for more info on the log likelihood, watch this: https://www.youtube.com/watch?v=8nogLkirA3I
## basically we are trying to find out which distribution fits the data (poisson or negative binomial in this case)
# more on why we use log likelihood and not just the likelihood here: https://math.stackexchange.com/questions/892832/why-we-consider-log-likelihood-instead-of-likelihood-in-gaussian-distribution
# and here: https://www.reddit.com/r/learnmachinelearning/comments/vewczb/why_do_we_maximize_loglikelihood_instead_of/
## likelihoods are a product, & usually you're mulitplying small values together which makes them even smaller and harder to measure
## log likelihoods turn these products into sums, and help to scale our small values without losing any information

# first a glm w/ Poisson distribution
spor.glm.step1<-step(glm(formula = round(`spmA; spore maturation protein A`) ~ ., data=mr.spor.ko.all[,c(6,9:10,12:14,36:45,49)],family="poisson"))
summary(spor.glm.step1)
spor.glm.step1$anova

# then a NB glm...
spor.glm.step2<-step(glm.nb(formula = round(`spmA; spore maturation protein A`) ~ ., data=mr.spor.ko.all[,c(6,9:10,12:14,36:45,49)]))
summary(spor.glm.step2)
spor.glm.step2$anova

# now we compare the loglikehoods * LRT only works for nested models
## the higher the value of loglikehood, the better the model fits the data
# In summary, the LRT tells us if it is beneficial to add parameters to our model, or if we should stick with our simpler model.
# In their most basic form, the hypotheses for the LRT are:
# H0: You should use the nested model.
# Ha: You should use the complex model.
# Thus, if you reject the H0, you can conclude that the complex model is significantly more accurate than the nested model, and you would choose to use the complex model.
# If you fail to reject the H0, you can conclude that the complex model is NOT significantly more accurate than the nested model, so you would choose to use the nested model instead.

# the smaller the AIC the better it's predictive power is
AIC(spor.glm.step1,spor.glm.step2) # compare their AICs
lrtest(spor.glm.step1,spor.glm.step2) #
# results between Poisson and NB GLM are the same, go with Poisson

# now, we have two models (both Poisson) that are significant -- we need to figure out which is better
spo.pois1<-glm(formula = round(`spmA; spore maturation protein A`) ~ precip_24hr_accum*SaltonSea, data=mr.spor.ko.all, family="poisson")
spo.pois.sum1<-summary(spo.pois1)
spo.pois.sum1

spo.pois2<-glm(formula = round(`spmA; spore maturation protein A`) ~ precip_24hr_accum*ave.v_N.S.wind, data=mr.spor.ko.all, family="poisson")
spo.pois.sum2<-summary(spo.pois2)
spo.pois.sum2

anova(spo.pois2,spo.pois1,test="Chisq") # second model seems to explain more variation?

# calculate McFadden's Pseudo R-squared for model for GLM aka aka Analysis of Deviance
# more here: https://bookdown.org/egarpor/PM-UC3M/glm-deviance.html
with(summary(spo.pois1), 1 - deviance/null.deviance)
with(summary(spo.pois2), 1 - deviance/null.deviance) # seems like second Poisson model is better

# final model:
summary(glm(formula = round(`spmA; spore maturation protein A`) ~ precip_24hr_accum*ave.v_N.S.wind, data=mr.spor.ko.all, family="poisson"), corr = TRUE)
spo.pois2<-glm(formula = round(`spmA; spore maturation protein A`) ~ precip_24hr_accum*ave.v_N.S.wind, data=mr.spor.ko.all, family="poisson")
summary(spo.pois2, corr = TRUE)$correlation # true correlation between the predictor variables

# now we can see the loglikehoods and see which one is greater
## the higher the value, the better the model fits the data
logLik(spo.pois2) # poisson
# the smaller the AIC the better it's predictive power is
AIC(spo.pois2)

# calculate McFadden's Pseudo R-squared for model for GLM aka aka Analysis of Deviance
# more here: https://bookdown.org/egarpor/PM-UC3M/glm-deviance.html
with(summary(spo.pois2), 1 - deviance/null.deviance)

# confirm that this is the best model by plotting the residuals and seeing their distribution
par(mar=c(1,1,1,1)) # change plot margins
plot(fitted(spo.pois2), spo.pois2$residuals) # produce residual vs. fitted plot
abline(0,0) # add a horizontal line at 0

# Q-Q norm plot of residuals
qqnorm(spo.pois2$residuals)
qqline(spo.pois2$residuals)

# density of residuals (how many points have certain residuals)
plot(density(spo.pois2$residuals))

# a better way to plot the residuals and their densities
ggplot(data=mr.spor.ko.all,aes(spo.pois2$residuals)) +
  geom_histogram(aes(y=..density..),bindwidth=30,color='black',fill='grey') +
  geom_density(alpha=.2,fill='lightpink')

#plot(spo.pois) # another way to plot residuals

# now let's plot the model!
spmA.plot1<-ggplot(data=mr.spor.ko.all,aes(x=precip_24hr_accum, y=`spmA; spore maturation protein A`))+
  stat_smooth(aes(x=precip_24hr_accum, y=round(`spmA; spore maturation protein A`)), method=stats::glm,method.args = list(family = "poisson")) +
  geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Accum Precip (24 hr)", y="spmA KO Coverage", title="spmA Coverage x Accum. Precipitation (24 hrs)",subtitle="Using Scaled Climate Data") +
  stat_fit_tidy(method = "glm",method.args = list(family=poisson,formula=y~x),
                mapping = aes(label = sprintf("p = %.3g",after_stat(x_p.value))),label.x=0.02,label.y=26)

ggsave(spmA.plot1,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/Sporulation/SSD_spmA_AccumPrecip_scatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

spmA.plot2<-ggplot(data=mr.spor.ko.all,aes(x=ave.v_N.S.wind, y=`spmA; spore maturation protein A`))+
  stat_smooth(aes(x=ave.v_N.S.wind, y=round(`spmA; spore maturation protein A`)), method=stats::glm,method.args = list(family = "poisson")) +
  geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="North to South Wind Vector Component (v)", y="spmA KO Coverage", title="spmA Coverage x North to South Wind",subtitle="Using Scaled Climate Data") +
  stat_fit_tidy(method = "glm",method.args = list(family=poisson,formula=y~x),
                mapping = aes(label = sprintf("p = %.3g",after_stat(x_p.value))),label.x=0.02,label.y=26)+
  stat_fit_glance(method = "glm",method.args = list(family=poisson,formula=y~x),
                  mapping = aes(label = sprintf("'r^2~=~%.3f~~italic(P)~=~%.2g'", after_stat(..r.squared..), after_stat(..p.value..))),label.x=0.02,label.y=28)

ggsave(spmA.plot2,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/Sporulation/SSD_spmA_SouthtoNorth_WindVector_v_scatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

#### Compare Variance of Sporulation Genes Across Sites ####

# spoIVCA is normally distributed; used rarefied counts to calculate ShanDiv
# use the following statisitcal tests for variance comparisons
## ANOVA: are variances significantly different between groups
## Tukey test: which groups' variances are significant different from one another
## Levene's test: is variance homogenous aka equal across samples?

head(mr.spor.ko.all)

fit1<-aov(`spoIVCA; site-specific DNA recombinase` ~ Site, data=mr.spor.ko.all)
# ANOVA is basically a regression but w/ categorical variables more info here https://www.statology.org/anova-vs-regression/
#pairwise.adonis(mr.spor.ko.all$`spoIVCA; site-specific DNA recombinase`, mr.spor.ko.all$SampDate, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(fit1)
# p = 0.334

# Tukey test - tells us which groups are significantly different from each other (more here: https://www.r-bloggers.com/2013/06/anova-and-tukeys-test-on-r/)
Tuk1<-TukeyHSD(fit1)
Tuk1$Site
# no sig differences

# fit.0<-aov(DustComplexity ~ as.factor(Elevation), data=mr.spor.ko.all)
# summary(fit.0)
# TukeyHSD(fit.0)

# Levene's test - test for homogeneity of variance
## Levene’s test is an inferential statistic used to check if the variances of a variable obtained for two or more groups are equal or not when data comes from a non-normal distribution
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## t test assumes that variances are the same, so Levene's test needs to be non significant
car::leveneTest(`spoIVCA; site-specific DNA recombinase` ~ Site, data = mr.spor.ko.all, center=mean)
# Levene's Test for Homogeneity of Variance (center = mean)
#       Df F value Pr(>F)
# group  3  1.5424 0.2345
# 20
# ^ p value is > 0.05 so we cannot reject the Null hypothesis -- variances are equal

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
# fligner.test(`spoIVCA; site-specific DNA recombinase` ~ SampDate, data = mr.spor.ko.all)
# Fligner-Killeen:med chi-squared = 1.9504, df = 2, p-value = 0.3771
# Which shows that the data do not deviate significantly from homogeneity.

#compare_means(`spoIVCA; site-specific DNA recombinase` ~ Site, data=mr.spor.ko.all, method="anova",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input

#### Do Temp Shock Fxns of Interest Correlate? ####

mr.cov.sum.tempshock.ko[1:4,]
grep("cspA;", colnames(mr.cov.sum.tempshock.ko)) # find index for column we want
grep("htpX;", colnames(mr.cov.sum.tempshock.ko)) # find index for column we want

tempshockOI.cor<-cor(mr.cov.sum.tempshock.ko[,c(2,10)])
tempshock.cor.p <- cor_pmat(mr.cov.sum.tempshock.ko[,c(2,10)])

ggcorrplot(tempshockOI.cor,method="square",lab = T,p.mat=tempshock.cor.p,sig.level = 0.05,
           hc.order=TRUE,outline.color="white",type = "lower")

par(mar=c(1,1,1,1))
plot(mr.cov.sum.tempshock.ko$`cspA; cold shock protein`,mr.cov.sum.tempshock.ko$`htpX; heat shock protein HtpX [EC:3.4.24.-]`)

#### Create Temp Shock Coverage Table, then Merge with Metadata ####

head(mr.tempshock.ko)

## create table of Sample x KO_ID
mr.tempshock.ko_table<-as.data.frame(dcast(mr.tempshock.ko, SampleID~KO_Function.KEGG, value.var="MR_SumCovPerKO"))
mr.tempshock.ko_table[1:5,1:5]

## merge tempshock functions with scaled metadata

mr.tempshock.ko.all<-merge(meta.all.scaled,mr.tempshock.ko_table,by="SampleID")
mr.tempshock.ko.all[1:5,1:20]

#### Temp Shock Fxns & Env Vars - Scatterplots ####
## heatmaps of traits of interest

# prep env data for downstream visualizations/analyses
env.vars<-mr.tempshock.ko.all[,c(6,9:10,12:14,36:45)] # subset env var data
env.vars.names<-names(mr.tempshock.ko.all[,c(6,9:10,12:14,36:45)]) # create vector list of variable names

## create function name variable
tempshock.fxn.names<-names(mr.tempshock.ko.all)[names(mr.tempshock.ko.all) %in% c("cspA; cold shock protein",
                                                                   "htpX; heat shock protein HtpX [EC:3.4.24.-]")] # pull out names of columns in df that contain "Axis" in name


## Loop to Generate Scatter Plots
### comparing y ~ x where y = tempshocke relative coverage and x = whatever environmental variable
head(mr.tempshock.ko.all)
tempshock.scatter.plot.list<-list() # create empty list for each plot to be stored in
tempshock.fxn.names
env.vars.names # sanity check
plotnum <- 1 # counting our model numbers for indexes purposes in the loop

for (i in tempshock.fxn.names) {
  for (j in env.vars.names){
    tempshock.fxn.scat.plot=scat.plot.fxn(mr.tempshock.ko.all,mr.tempshock.ko.all[,j],mr.tempshock.ko.all[,i])
    plot.titled = tempshock.fxn.scat.plot + ggtitle(paste(as.character(i),"Relative Coverage ~",as.character(j),sep=" "),subtitle="Using Scaled Climate Data")
    tempshock.scatter.plot.list[[plotnum]]=plot(plot.titled)
    ggsave(plot.titled,filename = paste("figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/TempShock/SSD_",as.character(i),"~",as.character(j),"_scatterplot.png",sep=""), width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).
    #print(names(clsOI.meta.all[i]))
    #print(clsOI.meta.all[j])
    plotnum<-1+plotnum
  }
}


# check if our loop + functions worked!
tempshock.scatter.plot.list[[1]]

ggplot(data=mr.tempshock.ko.all,aes(x=ave.air_temp, y=`cspA; cold shock protein`, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Ave Air Temp", y="spmA KO Coverage", title="spmA Coverage x Average Air Temp",subtitle="Using Scaled Climate Data")

tempshock.scatter.plot.list[[15]]

ggplot(data=mr.tempshock.ko.all,aes(x=ave.air_temp, y=`htpX; heat shock protein HtpX [EC:3.4.24.-]`, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Ave Air Temp", y="wecP KO Coverage", title="wecP Coverage x Average Air Temp",subtitle="Using Scaled Climate Data")


# now visualize genes by sampdate
tempshock.date.plot.list<-list() # create empty list for each plot to be stored in

tempshock.fxn.names
env.vars.names
plotnum <- 1 # counting our model numbers for indexes purposes in the loop

for (i in tempshock.fxn.names) {
  tempshock.date.scat.plot=date.mgm.scat.plot.fxn(mr.tempshock.ko.all,mr.tempshock.ko.all[,i])
  date.plot.titled = tempshock.date.scat.plot + ggtitle(paste(as.character(i),"Relative Coverage by Collection Date",sep=" "))
  tempshock.date.plot.list[[plotnum]]=plot(date.plot.titled)
  ggsave(date.plot.titled,filename = paste("figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/TempShock/SSD_",as.character(i),"_by_SampDate_scatterplot.png",sep=""), width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

  #print(names(clsOI.meta.all[i]))
  #print(clsOI.meta.all[j])
  plotnum<-1+plotnum

}

tempshock.date.plot.list[[1]]
ggplot(mr.tempshock.ko.all,aes(x=SampDate, y=`cspA; cold shock protein`, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x = "",
       y = "",
       shape="",
       color = "R")

tempshock.date.plot.list[[2]]
ggplot(mr.tempshock.ko.all,aes(x=SampDate, y=`htpX; heat shock protein HtpX [EC:3.4.24.-]`, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x = "",
       y = "",
       shape="",
       color = "R")


#### What Env Vars Predict Distribution of TempShock Genes ####

# first with cspA
# is cspA normally distributed?
shapiro.test(mr.tempshock.ko.all$`cspA; cold shock protein`) # what is the p-value?
# W = 0.96169, p-value = 0.4733, suggests that this function is normally distributed...
hist(mr.tempshock.ko.all$`cspA; cold shock protein`)

# visualize Q-Q plot for alpha div
# The Q-Q plot, or quantile-quantile plot, is a graphical tool to help us assess if a set of data plausibly came from some theoretical distribution such as a normal or exponential.
# For example, if we run a statistical analysis that assumes our residuals are normally distributed, we can use a normal Q-Q plot to check that assumption
# more on Q-Q plots here: https://data.library.virginia.edu/understanding-q-q-plots/
# more here too: https://grodri.github.io/glms/notes/c2s9#:~:text=8%20The%20Q%2DQ%20Plot,versus%20quantiles%20of%20a%20distribution.
qqnorm(mr.tempshock.ko.all$`cspA; cold shock protein`, pch = 1, frame = FALSE)
qqline(mr.tempshock.ko.all$`cspA; cold shock protein`, col = "red", lwd = 2)

grep("cspA;", colnames(mr.tempshock.ko.all)) # find index for column we want

tempshock.glm.step1<-step(lm(formula = round(`cspA; cold shock protein`) ~ ., data=mr.tempshock.ko.all[,c(6,9:10,12:14,36:45,49)]))
summary(tempshock.glm.step1)
tempshock.glm.step1$anova

# final model
tempshock.cspA.glm1<-lm(formula =`cspA; cold shock protein` ~ ave.v_N.S.wind, data=mr.tempshock.ko.all)
summary(tempshock.cspA.glm1, corr = TRUE)
summary(tempshock.cspA.glm1, corr = TRUE)$correlation # true correlation between the predictor variables

# confirm that this is the best model by plotting the residuals and seeing their distribution
par(mar=c(1,1,1,1)) # change plot margins
plot(fitted(tempshock.cspA.glm1), tempshock.cspA.glm1$residuals) # produce residual vs. fitted plot
abline(0,0) # add a horizontal line at 0

# Q-Q norm plot of residuals
qqnorm(tempshock.cspA.glm1$residuals)
qqline(tempshock.cspA.glm1$residuals)

# density of residuals (how many points have certain residuals)
plot(density(tempshock.cspA.glm1$residuals))

# a better way to plot the residuals and their densities
ggplot(data=mr.tempshock.ko.all,aes(tempshock.cspA.glm1$residuals)) +
  geom_histogram(aes(y=..density..),color='black',fill='grey') +
  geom_density(alpha=.2,fill='lightpink')

#plot(tempshock.cspA.glm1) # another way to plot residuals

# now let's plot the model
cspA.plot1<-ggplot(data=mr.tempshock.ko.all,aes(x=ave.v_N.S.wind, y=`cspA; cold shock protein`))+
  stat_smooth(aes(x=ave.v_N.S.wind, y=round(`cspA; cold shock protein`)), method=stats::lm) +
  geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  labs(x="North to South Wind Vector Component (v)", y="cspA KO Coverage") +
  stat_poly_eq(use_label(c("adj.R2", "P")), formula=y~x,small.p = TRUE,inherit.aes = TRUE) +
  theme(text=element_text(size=16, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),
        axis.text = element_text(size=16),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=16),legend.text = element_text(size=15))

ggsave(cspA.plot1,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/TempShock/SSD_cspA_SouthtoNorth_WindVector_v_scatterplot.png", width=18, height=13, dpi=300,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

cspA.plot2<-cspA.plot1 + geom_vline(xintercept=-0.3, linetype="dashed",color = "black", size=1)

ggsave(cspA.plot2,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/TempShock/SSD_cspA_SouthtoNorth_WindVector_v_scatterplot_NvSLine.png", width=18, height=13, dpi=300,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

# next with htpX
# is htpX normally distributed?
shapiro.test(mr.tempshock.ko.all$`htpX; heat shock protein HtpX [EC:3.4.24.-]`) # what is the p-value?
# W = 0.94895, p-value = 0.2571, suggests that this function is normally distributed...
hist(mr.tempshock.ko.all$`htpX; heat shock protein HtpX [EC:3.4.24.-]`)

# visualize Q-Q plot for alpha div
# The Q-Q plot, or quantile-quantile plot, is a graphical tool to help us assess if a set of data plausibly came from some theoretical distribution such as a normal or exponential.
# For example, if we run a statistical analysis that assumes our residuals are normally distributed, we can use a normal Q-Q plot to check that assumption
# more on Q-Q plots here: https://data.library.virginia.edu/understanding-q-q-plots/
# more here too: https://grodri.github.io/glms/notes/c2s9#:~:text=8%20The%20Q%2DQ%20Plot,versus%20quantiles%20of%20a%20distribution.
qqnorm(mr.tempshock.ko.all$`htpX; heat shock protein HtpX [EC:3.4.24.-]`, pch = 1, frame = FALSE)
qqline(mr.tempshock.ko.all$`htpX; heat shock protein HtpX [EC:3.4.24.-]`, col = "red", lwd = 2)

grep("htpX;", colnames(mr.tempshock.ko.all)) # find index for column we want

tempshock.glm.step2<-step(lm(formula = `htpX; heat shock protein HtpX [EC:3.4.24.-]` ~ ., data=mr.tempshock.ko.all[,c(6,9:10,12:14,36:45,57)]))
summary(tempshock.glm.step2)
tempshock.glm.step2$anova

summary(lm(formula = `htpX; heat shock protein HtpX [EC:3.4.24.-]` ~ BarrenLand+OpenWater, data=mr.tempshock.ko.all))
summary(lm(formula = `htpX; heat shock protein HtpX [EC:3.4.24.-]` ~ ave.u_E.W.wind*OpenWater, data=mr.tempshock.ko.all))
summary(lm(formula = `htpX; heat shock protein HtpX [EC:3.4.24.-]` ~ OpenWater+SaltonSea, data=mr.tempshock.ko.all))

tempshock.glm2a<-lm(formula = `htpX; heat shock protein HtpX [EC:3.4.24.-]` ~ precip_24hr_accum+BarrenLand+OpenWater, data=mr.tempshock.ko.all)
tempshock.glm2b<-lm(formula = `htpX; heat shock protein HtpX [EC:3.4.24.-]` ~ BarrenLand+OpenWater, data=mr.tempshock.ko.all)
tempshock.glm2c<-lm(formula = `htpX; heat shock protein HtpX [EC:3.4.24.-]` ~ OpenWater+SaltonSea, data=mr.tempshock.ko.all)
tempshock.glm2d<-lm(formula = `htpX; heat shock protein HtpX [EC:3.4.24.-]` ~ ave.u_E.W.wind*OpenWater, data=mr.tempshock.ko.all)

tempshock.glm2e<-lm(formula = `htpX; heat shock protein HtpX [EC:3.4.24.-]` ~ BarrenLand+OpenWater+SaltonSea, data=mr.tempshock.ko.all)

# compare models via likehood ratio tests
# In summary, the LRT tells us if it is beneficial to add parameters to our model, or if we should stick with our simpler model.
# In their most basic form, the hypotheses for the LRT are:
# H0: You should use the nested model.
# Ha: You should use the complex model.
# Thus, if you reject the H0, you can conclude that the complex model is significantly more accurate than the nested model, and you would choose to use the complex model.
# If you fail to reject the H0, you can conclude that the complex model is NOT significantly more accurate than the nested model, so you would choose to use the nested model instead.
# ^ here, complex model is more variables than nested model

lrtest(tempshock.glm2a,tempshock.glm2b,tempshock.glm2c)
lrtest(tempshock.glm2e,tempshock.glm2b,tempshock.glm2c)

# compare AICs
AIC(tempshock.glm2a,tempshock.glm2b,tempshock.glm2c,tempshock.glm2e)

# compare the analysis of deviance between models with a Chisquared-based estimate
# which variables do we add or remove to improve the model
anova(tempshock.glm2a,tempshock.glm2b,tempshock.glm2c,test="Chisq")
anova(tempshock.glm2a,test="Chisq")
anova(tempshock.glm2e,test="Chisq")

# overall seems like BarrenLand+OpenWater (tempshock.glm2b) is the best model?

# final model:
summary(lm(formula = `htpX; heat shock protein HtpX [EC:3.4.24.-]` ~ BarrenLand+OpenWater, data=mr.tempshock.ko.all), corr = TRUE)

htpX.lm<-lm(formula = `htpX; heat shock protein HtpX [EC:3.4.24.-]` ~ BarrenLand+OpenWater, data=mr.tempshock.ko.all)
summary(htpX.lm, corr = TRUE)$correlation # true correlation between the predictor variables

# confirm that this is the best model by plotting the residuals and seeing their distribution
par(mar=c(1,1,1,1)) # change plot tempshock.glm2b
plot(fitted(tempshock.glm2b), tempshock.glm2b$residuals) # produce residual vs. fitted plot
abline(0,0) # add a horizontal line at 0

# Q-Q norm plot of residuals
qqnorm(tempshock.glm2b$residuals)
qqline(tempshock.glm2b$residuals)

# density of residuals (how many points have certain residuals)
plot(density(tempshock.glm2b$residuals))

# a better way to plot the residuals and their densities
ggplot(data=mr.tempshock.ko.all,aes(tempshock.glm2b$residuals)) +
  geom_histogram(aes(y=..density..),color='black',fill='grey') +
  geom_density(alpha=.2,fill='lightpink')

#plot(tempshock.glm2) # another way to plot residuals

# now let's plot the model!
htpX.plot1<-ggplot(data=mr.tempshock.ko.all,aes(x=BarrenLand, y=`htpX; heat shock protein HtpX [EC:3.4.24.-]`))+
  stat_smooth(aes(x=BarrenLand, y=round(`htpX; heat shock protein HtpX [EC:3.4.24.-]`)), method=stats::lm) +
  geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Barren Land STF", y="htpX KO Coverage", title="htpX Coverage x Barren Land STF",subtitle="Using Scaled Climate Data") +
  stat_poly_eq(use_label(c("adj.R2", "P")), formula=y~x,small.p = TRUE)

ggsave(htpX.plot1,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/TempShock/SSD_htpX_BarrenLandSTF_scatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

htpX.plot2<-ggplot(data=mr.tempshock.ko.all,aes(x=OpenWater, y=`htpX; heat shock protein HtpX [EC:3.4.24.-]`))+
  stat_smooth(aes(x=OpenWater, y=round(`htpX; heat shock protein HtpX [EC:3.4.24.-]`)), method=stats::lm) +
  geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Open Water STF", y="htpX KO Coverage", title="htpX Coverage x Open Water STF",subtitle="Using Scaled Climate Data") +
  stat_poly_eq(use_label(c("adj.R2", "P")), formula=y~x,small.p = TRUE)

ggsave(htpX.plot2,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/TempShock/SSD_htpX_OpenWaterSTF_scatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

#### Compare Variance of Temp Shock Genes Across Sites ####

# cspA is normally distributed; used rarefied counts to calculate ShanDiv
# use the following statisitcal tests for variance comparisons
## ANOVA: are variances significantly different between groups
## Tukey test: which groups' variances are significant different from one another
## Levene's test: is variance homogenous aka equal across samples?

head(mr.tempshock.ko.all)

fit1<-aov(`cspA; cold shock protein` ~ Site, data=mr.tempshock.ko.all)
# ANOVA is basically a regression but w/ categorical variables more info here https://www.statology.org/anova-vs-regression/
#pairwise.adonis(mr.tempshock.ko.all$`cspA; cold shock protein`, mr.tempshock.ko.all$SampDate, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(fit1)
# p = 0.82

# Tukey test - tells us which groups are significantly different from each other (more here: https://www.r-bloggers.com/2013/06/anova-and-tukeys-test-on-r/)
Tuk1<-TukeyHSD(fit1)
Tuk1$Site
# no sig differences

# fit.0<-aov(DustComplexity ~ as.factor(Elevation), data=mr.tempshock.ko.all)
# summary(fit.0)
# TukeyHSD(fit.0)

# Levene's test - test for homogeneity of variance
## Levene’s test is an inferential statistic used to check if the variances of a variable obtained for two or more groups are equal or not when data comes from a non-normal distribution
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## t test assumes that variances are the same, so Levene's test needs to be non significant
car::leveneTest(`cspA; cold shock protein` ~ Site, data = mr.tempshock.ko.all, center=mean)
# Levene's Test for Homogeneity of Variance (center = mean)
#       Df F value Pr(>F)
# group  3  0.7507 0.5347
# 20
# ^ p value is > 0.05 so we cannot reject the Null hypothesis -- variances are equal

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
# fligner.test(`cspA; cold shock protein` ~ SampDate, data = mr.tempshock.ko.all)
# Fligner-Killeen:med chi-squared = 1.9504, df = 2, p-value = 0.3771
# Which shows that the data do not deviate significantly from homogeneity.

#compare_means(`cspA; cold shock protein` ~ Site, data=mr.tempshock.ko.all, method="anova",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input


#### Do UV-Damaged DNA Repair Fxns of Interest Correlate? ####

mr.cov.sum.UVrep.ko[1:4,]
grep("lexA;", colnames(mr.cov.sum.UVrep.ko)) # find index for column we want
grep("uvrA;", colnames(mr.cov.sum.UVrep.ko)) # find index for column we want
grep("uvrB;", colnames(mr.cov.sum.UVrep.ko)) # find index for column we want
grep("uvrC;", colnames(mr.cov.sum.UVrep.ko)) # find index for column we want

UVrepOI.cor<-cor(mr.cov.sum.UVrep.ko[,c(2,5:7)])
UVrep.cor.p <- cor_pmat(mr.cov.sum.UVrep.ko[,c(2,5:7)])

ggcorrplot(UVrepOI.cor,method="square",lab = T,p.mat=UVrep.cor.p,sig.level = 0.05,
           hc.order=TRUE,outline.color="white",type = "lower")

par(mar=c(1,1,1,1))
plot(mr.cov.sum.UVrep.ko$`lexA; repressor LexA [EC:3.4.21.88]`,mr.cov.sum.UVrep.ko$`uvrA; excinuclease ABC subunit A`)
# NOTE: not all of the genes we looked at are required for UVrepulation, so it's not worrisome necessarily that these genes do not correlate in mgms

#### Create UV-Damage DNA Repair Coverage Table, then Merge with Metadata ####

head(mr.UVrep.ko)

## create table of Sample x KO_ID
mr.UVrep.ko_table<-as.data.frame(dcast(mr.UVrep.ko, SampleID~KO_Function.KEGG, value.var="MR_SumCovPerKO"))
mr.UVrep.ko_table[1:5,1:5]

## merge UVrep functions with scaled metadata

mr.UVrep.ko.all<-merge(meta.all.scaled,mr.UVrep.ko_table,by="SampleID")
mr.UVrep.ko.all[1:5,1:20]

#### UV Damaged DNA Repair Fxns & Env Vars - Scatterplots ####
## heatmaps of traits of interest


# prep env data for downstream visualizations/analyses
env.vars<-mr.UVrep.ko.all[,c(6,9:10,12:14,36:45)] # subset env var data
env.vars.names<-names(mr.UVrep.ko.all[,c(6,9:10,12:14,36:45)]) # create vector list of variable names

## create function name variable
UVrep.fxn.names<-names(mr.UVrep.ko.all)[names(mr.UVrep.ko.all) %in% c("uvrA; excinuclease ABC subunit A",
                                                                   "uvrB; excinuclease ABC subunit B",
                                                                   "uvrC; excinuclease ABC subunit C",
                                                                   "lexA; repressor LexA [EC:3.4.21.88]")] # pull out names of columns in df that contain "Axis" in name


## Loop to Generate Scatter Plots
### comparing y ~ x where y = UVrepe relative coverage and x = whatever environmental variable
head(mr.UVrep.ko.all)
UVrep.scatter.plot.list<-list() # create empty list for each plot to be stored in
UVrep.fxn.names
env.vars.names # sanity check
plotnum <- 1 # counting our model numbers for indexes purposes in the loop

for (i in UVrep.fxn.names) {
  for (j in env.vars.names){
    UVrep.fxn.scat.plot=scat.plot.fxn(mr.UVrep.ko.all,mr.UVrep.ko.all[,j],mr.UVrep.ko.all[,i])
    plot.titled = UVrep.fxn.scat.plot + ggtitle(paste(as.character(i),"Relative Coverage ~",as.character(j),sep=" "),subtitle="Using Scaled Climate Data")
    UVrep.scatter.plot.list[[plotnum]]=plot(plot.titled)
    ggsave(plot.titled,filename = paste("figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/UV_DNA_Repair/SSD_",as.character(i),"~",as.character(j),"_scatterplot.png",sep=""), width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).
    #print(names(clsOI.meta.all[i]))
    #print(clsOI.meta.all[j])
    plotnum<-1+plotnum
  }
}


# check if our loop + functions worked!
UVrep.scatter.plot.list[[1]]

ggplot(data=mr.UVrep.ko.all,aes(x=ave.air_temp, y=`lexA; repressor LexA [EC:3.4.21.88]`, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Ave Air Temp", y="spmA KO Coverage", title="spmA Coverage x Average Air Temp",subtitle="Using Scaled Climate Data")

UVrep.scatter.plot.list[[15]]

ggplot(data=mr.UVrep.ko.all,aes(x=Shrub, y=`lexA; repressor LexA [EC:3.4.21.88]`, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Shrub STF", y="lexA KO Coverage", title="lexA Coverage x Shrub STF",subtitle="Using Scaled Climate Data")


# now visualize genes by sampdate
UVrep.date.plot.list<-list() # create empty list for each plot to be stored in

UVrep.fxn.names
env.vars.names
plotnum <- 1 # counting our model numbers for indexes purposes in the loop

for (i in UVrep.fxn.names) {
  UVrep.date.scat.plot=date.mgm.scat.plot.fxn(mr.UVrep.ko.all,mr.UVrep.ko.all[,i])
  date.plot.titled = UVrep.date.scat.plot + ggtitle(paste(as.character(i),"Relative Coverage by Collection Date",sep=" "))
  UVrep.date.plot.list[[plotnum]]=plot(date.plot.titled)
  ggsave(date.plot.titled,filename = paste("figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/UV_DNA_Repair/SSD_",as.character(i),"_by_SampDate_scatterplot.png",sep=""), width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

  #print(names(clsOI.meta.all[i]))
  #print(clsOI.meta.all[j])
  plotnum<-1+plotnum

}

UVrep.date.plot.list[[1]]
ggplot(mr.UVrep.ko.all,aes(x=SampDate, y=`lexA; repressor LexA [EC:3.4.21.88]`, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x = "",
       y = "",
       shape="",
       color = "R")

UVrep.date.plot.list[[2]]
ggplot(mr.UVrep.ko.all,aes(x=SampDate, y=`uvrA; excinuclease ABC subunit A`, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x = "",
       y = "",
       shape="",
       color = "R")


#### What Env Vars Predict Distribution of UV-DNA Damage Repair Genes ####
# first with lexA
# is lexA normally distributed?
shapiro.test(mr.UVrep.ko.all$`lexA; repressor LexA [EC:3.4.21.88]`) # what is the p-value?
# W = 0.90529, p-value = 0.02791, suggests that this function is NOT normally distributed...
hist(mr.UVrep.ko.all$`lexA; repressor LexA [EC:3.4.21.88]`)

# visualize Q-Q plot for alpha div
# The Q-Q plot, or quantile-quantile plot, is a graphical tool to help us assess if a set of data plausibly came from some theoretical distribution such as a normal or exponential.
# For example, if we run a statistical analysis that assumes our residuals are normally distributed, we can use a normal Q-Q plot to check that assumption
# more on Q-Q plots here: https://data.library.virginia.edu/understanding-q-q-plots/
# more here too: https://grodri.github.io/glms/notes/c2s9#:~:text=8%20The%20Q%2DQ%20Plot,versus%20quantiles%20of%20a%20distribution.
qqnorm(mr.UVrep.ko.all$`lexA; repressor LexA [EC:3.4.21.88]`, pch = 1, frame = FALSE)
qqline(mr.UVrep.ko.all$`lexA; repressor LexA [EC:3.4.21.88]`, col = "red", lwd = 2)

grep("lexA;", colnames(mr.UVrep.ko.all)) # find index for column we want

# compare Poisson glm with glm.nb
# first Poisson
UVrep.glm.step1<-step(glm(formula = round(`lexA; repressor LexA [EC:3.4.21.88]`) ~ ., data=mr.UVrep.ko.all[,c(6,9:10,12:14,36:45,49)],family="poisson"))
summary(UVrep.glm.step1)
UVrep.glm.step1$anova

# now NB GlM
UVrep.glm.step2<-step(glm.nb(formula = round(`lexA; repressor LexA [EC:3.4.21.88]`) ~ ., data=mr.UVrep.ko.all[,c(6,9:10,12:14,36:45,49)]))
summary(UVrep.glm.step2)

# the smaller the AIC the better it's predictive power is
AIC(UVrep.glm.step1,UVrep.glm.step2) # Poisson is better, lower AIC

# lastly let's try a likelihood ratio test to compare
lrtest(UVrep.glm.step1,UVrep.glm.step2) #

anova(UVrep.glm.step1,UVrep.glm.step2,test="Chisq") # models seem to be the same

# results from Poisson GLM and NB GLM are the same -- go with Poisson

# final model:
summary(glm(formula = round(`lexA; repressor LexA [EC:3.4.21.88]`) ~ ave.v_N.S.wind, data=mr.UVrep.ko.all,family="poisson"), corr = TRUE)

UVrep.glm1<-glm(formula = round(`lexA; repressor LexA [EC:3.4.21.88]`) ~ ave.v_N.S.wind, data=mr.UVrep.ko.all,family="poisson")
summary(UVrep.glm1, corr = TRUE)$correlation # true correlation between the predictor variables

# calculate McFadden's Pseudo R-squared for model for GLM aka aka Analysis of Deviance
# more here: https://bookdown.org/egarpor/PM-UC3M/glm-deviance.html
with(summary(UVrep.glm1), 1 - deviance/null.deviance)

# confirm that this is the best model by plotting the residuals and seeing their distribution
par(mar=c(1,1,1,1)) # change plot margins
plot(fitted(UVrep.glm1), UVrep.glm1$residuals) # produce residual vs. fitted plot
abline(0,0) # add a horizontal line at 0

# Q-Q norm plot of residuals
qqnorm(UVrep.glm1$residuals)
qqline(UVrep.glm1$residuals)

# density of residuals (how many points have certain residuals)
plot(density(UVrep.glm1$residuals))

# a better way to plot the residuals and their densities
ggplot(data=mr.UVrep.ko.all,aes(UVrep.glm1$residuals)) +
  geom_histogram(aes(y=..density..),bindwidth=30,color='black',fill='grey') +
  geom_density(alpha=.2,fill='lightpink')

#plot(UVrep.glm1) # another way to plot residuals

# now let's plot the model!

lexA.plot1<-ggplot(data=mr.UVrep.ko.all,aes(x=ave.v_N.S.wind, y=round(`lexA; repressor LexA [EC:3.4.21.88]`)))+
  stat_smooth(aes(x=ave.v_N.S.wind, y=round(`lexA; repressor LexA [EC:3.4.21.88]`)), method=stats::glm, method.args=list(family="poisson")) +
  geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021"))+
  scale_shape_manual(name="Site",values = c(0,1,16,15)) + theme_classic() +
  theme(text=element_text(size=16, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),
        axis.text = element_text(size=16),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=16),legend.text = element_text(size=15)) +
  labs(x="North to South Wind Vector Component (v)", y="lexA KO Coverage") +
  stat_fit_tidy(method = "glm",method.args = list(family=poisson,formula=y~x),
                mapping = aes(label = sprintf('p = %.3f',after_stat(x_p.value))),label.y=10)

ggsave(lexA.plot1,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/UV_DNA_Repair/SSD_lexA_SouthtoNorth_WindVector_v_scatterplot.png", width=18, height=13, dpi=300,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

lexA.plot2<-lexA.plot1 + geom_vline(xintercept=-0.29, linetype="dashed",color = "black", size=1)

ggsave(lexA.plot2,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/UV_DNA_Repair/SSD_lexA_SouthtoNorth_WindVector_v_scatterplot_NvSLine.png", width=18, height=13, dpi=300,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

# # next with recA
# shapiro.test(mr.UVrep.ko.all$`recA; recombination protein RecA`) # what is the p-value?
# # W = 0.93157, p-value = 0.1057, suggests that this function is normally distributed...
# hist(mr.UVrep.ko.all$`recA; recombination protein RecA`)
#
# # visualize Q-Q plot for alpha div
# # The Q-Q plot, or quantile-quantile plot, is a graphical tool to help us assess if a set of data plausibly came from some theoretical distribution such as a normal or exponential.
# # For example, if we run a statistical analysis that assumes our residuals are normally distributed, we can use a normal Q-Q plot to check that assumption
# # more on Q-Q plots here: https://data.library.virginia.edu/understanding-q-q-plots/
# # more here too: https://grodri.github.io/glms/notes/c2s9#:~:text=8%20The%20Q%2DQ%20Plot,versus%20quantiles%20of%20a%20distribution.
# qqnorm(mr.UVrep.ko.all$`recA; recombination protein RecA`, pch = 1, frame = FALSE)
# qqline(mr.UVrep.ko.all$`recA; recombination protein RecA`, col = "red", lwd = 2)
#
# grep("recA;", colnames(mr.UVrep.ko.all)) # find index for column we want
#
# UVrep.glm.step2<-step(lm(formula = `recA; recombination protein RecA` ~ ., data=mr.UVrep.ko.all[,c(6,9:10,12:14,36:45,49)]))
# summary(UVrep.glm.step2)
# UVrep.glm.step2$anova
#
# summary(lm(formula = `recA; recombination protein RecA` ~ precip_24hr_accum, data=mr.UVrep.ko.all))
# # ^ seems to be the best model
#
# recA.plot1<-ggplot(data=mr.UVrep.ko.all,aes(x=precip_24hr_accum, y=`recA; recombination protein RecA`))+
#   stat_smooth(aes(x=precip_24hr_accum, y=`recA; recombination protein RecA`), method=stats::lm) +
#   geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
#                      labels=c("July 2020","August 2020","November 2020",
#                               "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
#   scale_shape_manual(name="Site",values = c(0,1,16,15)) +
#   theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
#         axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Accum Precip (24 hr)", y="recA KO Coverage", title="recA Coverage x Accum. Precip",subtitle="Using Scaled Climate Data") +
#   stat_poly_eq(use_label(c("adj.R2","P")), formula=y~x)
#
# ggsave(recA.plot1,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/UV_DNA_Repair/SSD_recA_PrecipAccum24hr_scatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

# next with uvrA
shapiro.test(mr.UVrep.ko.all$`uvrA; excinuclease ABC subunit A`) # what is the p-value?
# W = 0.94923, p-value = 0.2608, suggests that this function is normally distributed...
hist(mr.UVrep.ko.all$`uvrA; excinuclease ABC subunit A`)

# visualize Q-Q plot for alpha div
# The Q-Q plot, or quantile-quantile plot, is a graphical tool to help us assess if a set of data plausibly came from some theoretical distribution such as a normal or exponential.
# For example, if we run a statistical analysis that assumes our residuals are normally distributed, we can use a normal Q-Q plot to check that assumption
# more on Q-Q plots here: https://data.library.virginia.edu/understanding-q-q-plots/
# more here too: https://grodri.github.io/glms/notes/c2s9#:~:text=8%20The%20Q%2DQ%20Plot,versus%20quantiles%20of%20a%20distribution.
qqnorm(mr.UVrep.ko.all$`uvrA; excinuclease ABC subunit A`, pch = 1, frame = FALSE)
qqline(mr.UVrep.ko.all$`uvrA; excinuclease ABC subunit A`, col = "red", lwd = 2)

grep("uvrA;", colnames(mr.UVrep.ko.all)) # find index for column we want

UVrep.glm.step3<-step(lm(formula = `uvrA; excinuclease ABC subunit A` ~ ., data=mr.UVrep.ko.all[,c(6,9:10,12:14,36:45,52)]))
summary(UVrep.glm.step3)
UVrep.glm.step3$anova

summary(lm(formula = `uvrA; excinuclease ABC subunit A` ~ precip_24hr_accum+ave.v_N.S.wind+OpenWater, data=mr.UVrep.ko.all))
summary(lm(formula = `uvrA; excinuclease ABC subunit A` ~ ave.v_N.S.wind+OpenWater, data=mr.UVrep.ko.all))
# ^ best model
summary(lm(formula = `uvrA; excinuclease ABC subunit A` ~ ave.v_N.S.wind, data=mr.UVrep.ko.all))
summary(lm(formula = `uvrA; excinuclease ABC subunit A` ~ OpenWater, data=mr.UVrep.ko.all))
# OpenWater was near sig as well but not informative - really small contribution

# compare model with and without OpenWater var
UVrep.glm2a<-lm(formula = `uvrA; excinuclease ABC subunit A` ~ precip_24hr_accum+ave.v_N.S.wind+OpenWater, data=mr.UVrep.ko.all)
UVrep.glm2b<-lm(formula = `uvrA; excinuclease ABC subunit A` ~ ave.v_N.S.wind+OpenWater, data=mr.UVrep.ko.all)
UVrep.glm2c<-lm(formula = `uvrA; excinuclease ABC subunit A` ~ OpenWater, data=mr.UVrep.ko.all)
UVrep.glm2d<-lm(formula = `uvrA; excinuclease ABC subunit A` ~ ave.v_N.S.wind, data=mr.UVrep.ko.all)

# compare models via likehood ratio tests
lrtest(UVrep.glm2a,UVrep.glm2b)
lrtest(UVrep.glm2a,UVrep.glm2c)
lrtest(UVrep.glm2a,UVrep.glm2d)

# compare AICs
AIC(UVrep.glm2a,UVrep.glm2b,UVrep.glm2c,UVrep.glm2d)

# compare the analysis of deviance between models with a Chisquared-based estimate
# which variables do we add or remove to improve the model
anova(UVrep.glm2a,UVrep.glm2b,test="Chisq")
anova(UVrep.glm2a,UVrep.glm2c,test="Chisq")
anova(UVrep.glm2a,UVrep.glm2d,test="Chisq")

# final model
UVrep.glm2<-lm(formula = `uvrA; excinuclease ABC subunit A` ~ ave.v_N.S.wind+OpenWater, data=mr.UVrep.ko.all)
summary(UVrep.glm2, corr = TRUE)$correlation

# confirm that this is the best model by plotting the residuals and seeing their distribution
par(mar=c(1,1,1,1)) # change plot margins
plot(fitted(UVrep.glm2), UVrep.glm2$residuals) # produce residual vs. fitted plot
abline(0,0) # add a horizontal line at 0

# Q-Q norm plot of residuals
qqnorm(UVrep.glm2$residuals)
qqline(UVrep.glm2$residuals)

# density of residuals (how many points have certain residuals)
plot(density(UVrep.glm2$residuals))

# a better way to plot the residuals and their densities
ggplot(data=mr.UVrep.ko.all,aes(UVrep.glm2$residuals)) +
  geom_histogram(aes(y=..density..),color='black',fill='grey') +
  geom_density(alpha=.2,fill='lightpink')

#plot(UVrep.glm2a) # another way to plot residuals

# now let's plot the model!

uvrA.plot1<-ggplot(data=mr.UVrep.ko.all,aes(x=ave.v_N.S.wind, y=`uvrA; excinuclease ABC subunit A`))+
  stat_smooth(aes(x=ave.v_N.S.wind, y=`uvrA; excinuclease ABC subunit A`), method=stats::lm) +
  geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="North to South Wind Vector Component (v)", y="uvrA KO Coverage", title="uvrA Coverage x North to South Wind",subtitle="Using Scaled Climate Data") +
  stat_poly_eq(use_label(c("adj.R2","P")), formula=y~x)

ggsave(uvrA.plot1,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/UV_DNA_Repair/SSD_uvrA_SouthtoNorth_WindVector_v_scatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

uvrA.plot2<-ggplot(data=mr.UVrep.ko.all,aes(x=OpenWater, y=`uvrA; excinuclease ABC subunit A`))+
  stat_smooth(aes(x=OpenWater, y=`uvrA; excinuclease ABC subunit A`), method=stats::lm) +
  geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Open Water STF", y="uvrA KO Coverage", title="uvrA Coverage x Open Water STF",subtitle="Using Scaled Climate Data") +
  stat_poly_eq(use_label(c("adj.R2","P")), formula=y~x)

ggsave(uvrA.plot2,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/UV_DNA_Repair/SSD_uvrA_SouthtoNorth_WindVector_v_scatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

# NOTE: when graphing uvrA with OpenWater STF, you can see that higher coverage of uvrA is with samples that have lowest OpenWater STF, and as Open Water STF increases, uvrA coverage decreases
# inverse trend observed with Wind Vector Component v...
# could indicate that higher coverage of this gene appears to be at sites closer to Sea?
# also appears to have a seasonal trend...

# next with uvrB
shapiro.test(mr.UVrep.ko.all$`uvrB; excinuclease ABC subunit B`) # what is the p-value?
# W = 0.9334, p-value = 0.1161, suggests that this function is normally distributed...
hist(mr.UVrep.ko.all$`uvrB; excinuclease ABC subunit B`)

# visualize Q-Q plot for alpha div
# The Q-Q plot, or quantile-quantile plot, is a graphical tool to help us assess if a set of data plausibly came from some theoretical distribution such as a normal or exponential.
# For example, if we run a statistical analysis that assumes our residuals are normally distributed, we can use a normal Q-Q plot to check that assumption
# more on Q-Q plots here: https://data.library.virginia.edu/understanding-q-q-plots/
# more here too: https://grodri.github.io/glms/notes/c2s9#:~:text=8%20The%20Q%2DQ%20Plot,versus%20quantiles%20of%20a%20distribution.
qqnorm(mr.UVrep.ko.all$`uvrB; excinuclease ABC subunit B`, pch = 1, frame = FALSE)
qqline(mr.UVrep.ko.all$`uvrB; excinuclease ABC subunit B`, col = "red", lwd = 2)

grep("uvrB;", colnames(mr.UVrep.ko.all)) # find index for column we want

UVrep.glm.step3<-step(lm(formula = `uvrB; excinuclease ABC subunit B` ~ ., data=mr.UVrep.ko.all[,c(6,9:10,12:14,36:45,53)]))
summary(UVrep.glm.step3)
UVrep.glm.step3$anova

summary(lm(formula = `uvrB; excinuclease ABC subunit B` ~ ave.v_N.S.wind, data=mr.UVrep.ko.all))
summary(lm(formula = `uvrB; excinuclease ABC subunit B` ~ precip_24hr_accum+ave.v_N.S.wind+OpenWater, data=mr.UVrep.ko.all))
summary(lm(formula = `uvrB; excinuclease ABC subunit B` ~ ave.v_N.S.wind+OpenWater, data=mr.UVrep.ko.all))

# final model
summary(lm(formula = `uvrB; excinuclease ABC subunit B` ~ ave.v_N.S.wind+OpenWater, data=mr.UVrep.ko.all), corr = TRUE)
# ^ best model

UVrep.glm3<-lm(formula = `uvrB; excinuclease ABC subunit B` ~ ave.v_N.S.wind+OpenWater, data=mr.UVrep.ko.all)
summary(UVrep.glm1, corr = TRUE)$correlation # true correlation between the predictor variables

# confirm that this is the best model by plotting the residuals and seeing their distribution
par(mar=c(1,1,1,1)) # change plot margins
plot(fitted(UVrep.glm3), UVrep.glm3$residuals) # produce residual vs. fitted plot
abline(0,0) # add a horizontal line at 0

# Q-Q norm plot of residuals
qqnorm(UVrep.glm3$residuals)
qqline(UVrep.glm3$residuals)

# density of residuals (how many points have certain residuals)
plot(density(UVrep.glm3$residuals))

# a better way to plot the residuals and their densities
ggplot(data=mr.UVrep.ko.all,aes(UVrep.glm3$residuals)) +
  geom_histogram(aes(y=..density..),color='black',fill='grey') +
  geom_density(alpha=.2,fill='lightpink')

#plot(UVrep.glm3) # another way to plot residuals

# now let's plot the model!
# uvrB.plot1<-ggplot(data=mr.UVrep.ko.all,aes(x=precip_24hr_accum, y=`uvrB; excinuclease ABC subunit B`))+
#   stat_smooth(aes(x=precip_24hr_accum, y=`uvrB; excinuclease ABC subunit B`), method=stats::lm) +
#   geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
#                      labels=c("July 2020","August 2020","November 2020",
#                               "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
#   scale_shape_manual(name="Site",values = c(0,1,16,15)) +
#   theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
#         axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Accum. Precip (24 hrs)", y="uvrB KO Coverage", title="uvrB Coverage x Accum. Precipitation",subtitle="Using Scaled Climate Data") +
#   stat_poly_eq(use_label(c("adj.R2","P")), formula=y~x)
#
# ggsave(uvrB.plot1,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/UV_DNA_Repair/SSD_uvrB_AccumPrecip_scatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

uvrB.plot2<-ggplot(data=mr.UVrep.ko.all,aes(x=ave.v_N.S.wind, y=`uvrB; excinuclease ABC subunit B`))+
  stat_smooth(aes(x=ave.v_N.S.wind, y=`uvrB; excinuclease ABC subunit B`), method=stats::lm) +
  geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="North to South Wind Vector Component (v)", y="uvrB KO Coverage", title="uvrB Coverage x North to South Wind",subtitle="Using Scaled Climate Data") +
  stat_poly_eq(use_label(c("adj.R2","P")), formula=y~x)

ggsave(uvrB.plot2,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/UV_DNA_Repair/SSD_uvrB_SouthtoNorth_WindVector_vscatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

uvrB.plot3<-ggplot(data=mr.UVrep.ko.all,aes(x=OpenWater, y=`uvrB; excinuclease ABC subunit B`))+
  stat_smooth(aes(x=OpenWater, y=`uvrB; excinuclease ABC subunit B`), method=stats::lm) +
  geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Open Water STF", y="uvrB KO Coverage", title="uvrB Coverage x Open Water STF",subtitle="Using Scaled Climate Data") +
  stat_poly_eq(use_label(c("adj.R2","P")), formula=y~x)

ggsave(uvrB.plot3,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/UV_DNA_Repair/SSD_uvrB_OpenWater_scatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

# next with uvrC
shapiro.test(mr.UVrep.ko.all$`uvrC; excinuclease ABC subunit C`) # what is the p-value?
# W = 0.85326, p-value = 0.002514, suggests that this function is not normally distributed...
hist(mr.UVrep.ko.all$`uvrC; excinuclease ABC subunit C`)

# visualize Q-Q plot for alpha div
# The Q-Q plot, or quantile-quantile plot, is a graphical tool to help us assess if a set of data plausibly came from some theoretical distribution such as a normal or exponential.
# For example, if we run a statistical analysis that assumes our residuals are normally distributed, we can use a normal Q-Q plot to check that assumption
# more on Q-Q plots here: https://data.library.virginia.edu/understanding-q-q-plots/
# more here too: https://grodri.github.io/glms/notes/c2s9#:~:text=8%20The%20Q%2DQ%20Plot,versus%20quantiles%20of%20a%20distribution.
qqnorm(mr.UVrep.ko.all$`uvrC; excinuclease ABC subunit C`, pch = 1, frame = FALSE)
qqline(mr.UVrep.ko.all$`uvrC; excinuclease ABC subunit C`, col = "red", lwd = 2)

grep("uvrC;", colnames(mr.UVrep.ko.all)) # find index for column we want

# compare Poisson glm with glm.nb
# first Poisson
UVrep.glm.step3<-step(glm(formula = round(`uvrC; excinuclease ABC subunit C`) ~ ., data=mr.UVrep.ko.all[,c(6,9:10,12:14,36:45,54)],family="poisson"))
summary(UVrep.glm.step3)
UVrep.glm.step3$anova

# then NB GLM
UVrep.glm.step4<-step(glm.nb(formula = round(`uvrC; excinuclease ABC subunit C`) ~ ., data=mr.UVrep.ko.all[,c(6,9:10,12:14,36:45,54)]))
summary(UVrep.glm.step4)

# the smaller the AIC the better it's predictive power is
AIC(UVrep.glm.step3,UVrep.glm.step4) # Poisson is better

# lastly let's try a likelihood ratio test to compare
lrtest(UVrep.glm.step3,UVrep.glm.step4) # pois is significantly better model based on likelihood ratio test...

# now let's compare which variables to include or exclude...
summary(glm(formula = round(`uvrC; excinuclease ABC subunit C`) ~ precip_24hr_accum+ave.v_N.S.wind+OpenWater, data=mr.UVrep.ko.all,family="poisson"))
summary(glm(formula = round(`uvrC; excinuclease ABC subunit C`) ~ ave.v_N.S.wind+OpenWater, data=mr.UVrep.ko.all,family="poisson"))
summary(glm(formula = round(`uvrC; excinuclease ABC subunit C`) ~ ave.v_N.S.wind, data=mr.UVrep.ko.all,family="poisson"))
summary(glm(formula = round(`uvrC; excinuclease ABC subunit C`) ~ OpenWater, data=mr.UVrep.ko.all,family="poisson"))

uvrC.mod1<-glm(formula = round(`uvrC; excinuclease ABC subunit C`) ~ ave.v_N.S.wind+OpenWater, data=mr.UVrep.ko.all,family="poisson")
uvrC.mod2<-glm(formula = round(`uvrC; excinuclease ABC subunit C`) ~ ave.v_N.S.wind, data=mr.UVrep.ko.all,family="poisson")
uvrC.mod3<-glm(formula = round(`uvrC; excinuclease ABC subunit C`) ~ OpenWater, data=mr.UVrep.ko.all,family="poisson")

# compare models via likehood ratio tests
lrtest(uvrC.mod1,uvrC.mod2) # model with just Wind Vector Component v is significantly better model based on likelihood ratio test...
lrtest(uvrC.mod1,uvrC.mod3) # model with just OpenWater STF is significantly better model based on likelihood ratio test...

# compare the analysis of deviance between models with a Chisquared-based estimate
# which variables do we add or remove to improve the model
anova(uvrC.mod1,uvrC.mod2,test="Chisq") # model with just Wind Vector Component v is significantly better
anova(uvrC.mod1,uvrC.mod3,test="Chisq") # model with just oPen Water is significantly better

# compare models' AIC
AIC(uvrC.mod1,uvrC.mod2,uvrC.mod3)

# while ave.wind_direction+OpenWater model has lowest AIC, the model with only ave.wind_direction was significantly better
## the model with just OpenWater was marginally significantly better

# final model:
summary(glm(formula = round(`uvrC; excinuclease ABC subunit C`) ~ ave.v_N.S.wind+OpenWater, data=mr.UVrep.ko.all,family="poisson"), corr = TRUE)

uvrC.glm<-glm(formula = round(`uvrC; excinuclease ABC subunit C`) ~ ave.v_N.S.wind+OpenWater, data=mr.UVrep.ko.all,family="poisson")
summary(uvrC.glm, corr = TRUE)$correlation # true correlation between the predictor variables

1 - pchisq(deviance(uvrC.glm), df.residual(uvrC.glm)) # chi-square test p value; more here: https://www.r-bloggers.com/2022/05/calculate-the-p-value-from-chi-square-statistic-in-r/
# p = 0.4070469, cannot reject null hypothesis that distribution is different from Poisson

# calculate McFadden's Pseudo R-squared for model for GLM aka aka Analysis of Deviance
# more here: https://bookdown.org/egarpor/PM-UC3M/glm-deviance.html
with(summary(uvrC.glm), 1 - deviance/null.deviance)

# confirm that this is the best model by plotting the residuals and seeing their distribution
par(mar=c(1,1,1,1)) # change plot margins
plot(fitted(uvrC.glm), uvrC.glm$residuals) # produce residual vs. fitted plot
abline(0,0) # add a horizontal line at 0

# Q-Q norm plot of residuals
qqnorm(uvrC.glm$residuals)
qqline(uvrC.glm$residuals)

# density of residuals (how many points have certain residuals)
plot(density(uvrC.glm$residuals))

# a better way to plot the residuals and their densities
ggplot(data=mr.UVrep.ko.all,aes(uvrC.glm$residuals)) +
  geom_histogram(aes(y=..density..),color='black',fill='grey') +
  geom_density(alpha=.2,fill='lightpink')

#plot(uvrC.glm) # another way to plot residuals

# save summary of model for plottin the p value -- stat_fit_tidy is changing the p val for some reason?
uvrC.glm.sum<-summary(glm(formula = round(`uvrC; excinuclease ABC subunit C`) ~ ave.v_N.S.wind+OpenWater, data=mr.UVrep.ko.all,family="poisson"))
uvrC.glm.sum$coefficients[11] # the p value for Wind Vector Component v
uvrC.glm.sum$coefficients[12] # the p value for OpenWater

# now let's plot the model!
uvrC.plot1<-ggplot(data=mr.UVrep.ko.all,aes(x=ave.v_N.S.wind, y=`uvrC; excinuclease ABC subunit C`))+
  stat_smooth(aes(x=ave.v_N.S.wind, y=`uvrC; excinuclease ABC subunit C`), method=stats::glm, method.args=list(family="poisson")) +
  geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="North to South Wind Vector Component (v)", y="uvrC KO Coverage", title="uvrC Coverage x North to South Wind",subtitle="Using Scaled Climate Data") +
  stat_fit_tidy(method = "glm",method.args = list(family=poisson,formula=y~x),
                mapping = aes(label = sprintf("p = %.3g",uvrC.glm.sum$coefficients[11])),label.x=0.02,label.y=26)

ggsave(uvrC.plot1,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/UV_DNA_Repair/SSD_uvrC_SouthtoNorth_WindVector_v_scatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

uvrC.plot2<-ggplot(data=mr.UVrep.ko.all,aes(x=OpenWater, y=`uvrC; excinuclease ABC subunit C`))+
  stat_smooth(aes(x=OpenWater, y=`uvrC; excinuclease ABC subunit C`), method=stats::glm, method.args=list(family="poisson")) +
  geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Open Water STF", y="uvrC KO Coverage", title="uvrC Coverage x Open Water STF",subtitle="Using Scaled Climate Data") +
  stat_fit_tidy(method = "glm",method.args = list(family=poisson,formula=y~x),
                mapping = aes(label = sprintf("p = %.3g",uvrC.glm.sum$coefficients[12])),label.x=0.02,label.y=26)

ggsave(uvrC.plot2,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/UV_DNA_Repair/SSD_uvrC_OpenWaterSTF_scatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

#### Compare Variance of UV-Damage DNA Repair Genes Across Sites ####
# uvrA is normally distributed; used rarefied counts to calculate SR
fit1<-aov(`uvrA; excinuclease ABC subunit A` ~ Site, data=mr.UVrep.ko.all)
# ANOVA is basically a regression but w/ categorical variables more info here https://www.statology.org/anova-vs-regression/
#pairwise.adonis(bac.div.metadat.rar$AveShanDiv, bac.div.metadat.rar$SampDate, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(fit1)
#Df           Sum Sq Mean Sq    F value   Pr(>F)
#Site         3   99.5   33.18   0.412  0.746
#Residuals   20 1611.0   80.55

#p.adjust(summary(fit1)[[1]][["Pr(>F)"]][1],method="bonferroni")

# Tukey test - tells us which groups are significantly different from each other (more here: https://www.r-bloggers.com/2013/06/anova-and-tukeys-test-on-r/)
Tuk1<-TukeyHSD(fit1)
Tuk1

# Levene's test - test for homogeneity of variance
## Levene’s test is an inferential statistic used to check if the variances of a variable obtained for two or more groups are equal or not when data comes from a non-normal distribution
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## t test assumes that variances are the same, so Levene's test needs to be non significant
car::leveneTest(`uvrA; excinuclease ABC subunit A` ~ Site, data=mr.UVrep.ko.all, center=mean)
# Levene's Test for Homogeneity of Variance (center = mean)
#       Df F value Pr(>F)
# group  3  4.4429 0.01508 *
# 20
# ^ p value is < 0.05 so we CAN reject the Null hypothesis -- variances are NOT equal

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
# fligner.test(AveShanDiv ~ SampDate, data = bac.div.metadat.rar)
# Fligner-Killeen:med chi-squared = 1.9504, df = 2, p-value = 0.3771
# Which shows that the data do not deviate significantly from homogeneity.

#compare_means(`uvrA; excinuclease ABC subunit A` ~ Site, data=mr.UVrep.ko.all, method="anova",p.adjust.method = "bonferroni")

# lexA next, which is NOT normally distributed
# use the following statisitcal tests for variance comparisons
## Kruskal: are variances significantly different between groups
## Dunn test: which groups' variances are significant different from one another
## Fligner test: is variance homogenous aka equal across samples?

# Kruskal-Wallis test is an ANOVA for non-normal data
fit2<-kruskal.test(`lexA; repressor LexA [EC:3.4.21.88]` ~ Site, data=mr.UVrep.ko.all)

fit2
# Kruskal-Wallis chi-squared = 2.66, df = 3, p-value = 0.4471

# Instead of using Tukey test, we can use Dunn's test to see which groups significantly vary if Kruskal-Wallis test is significant
# ANOVA + Tukey for normally distributed data, Kruskal-Wallis + Dunn's test for non-normal data
rstatix::dunn_test(mr.UVrep.ko.all, `lexA; repressor LexA [EC:3.4.21.88]` ~ Site, p.adjust.method = "bonferroni", detailed = TRUE)

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(`lexA; repressor LexA [EC:3.4.21.88]` ~ Site, data = mr.UVrep.ko.all)
# Fligner-Killeen:med chi-squared = 1.1557, df = 3, p-value = 0.7636
# Which shows that the data do NOT deviate significantly from homogeneity.

#compare_means(`lexA; repressor LexA [EC:3.4.21.88]` ~ Site, data=mr.UVrep.ko.all, method="wilcox.test",p.adjust.method = "bonferroni")

#compare_means(`lexA; repressor LexA [EC:3.4.21.88]` ~ Site, data=mr.UVrep.ko.all, method="kruskal.test",p.adjust.method = "bonferroni")


#### Do Quorum Sensing Fxns of Interest Correlate? ####

mr.cov.sum.quorsens.ko[1:4,]
grep("bisR;", colnames(mr.cov.sum.quorsens.ko)) # find index for column we want
grep("prgX;", colnames(mr.cov.sum.quorsens.ko)) # find index for column we want

quorsensOI.cor<-cor(mr.cov.sum.quorsens.ko[,c(2,8)])
quorsens.cor.p <- cor_pmat(mr.cov.sum.quorsens.ko[,c(2,8)])

ggcorrplot(quorsensOI.cor,method="square",lab = T,p.mat=quorsens.cor.p,sig.level = 0.05,
           hc.order=TRUE,outline.color="white",type = "lower")

par(mar=c(1,1,1,1))
plot(mr.cov.sum.quorsens.ko$`bisR; LuxR family transcriptional regulator, quorum-sensing system regulator BisR`,mr.cov.sum.quorsens.ko$`prgX; HTH-type transcriptional regulator, major conjugation operon repressor`)
# NOTE: not all of the genes we looked at are required for quorsensulation, so it's not worrisome necessarily that these genes do not correlate in mgms

#### Create Quorum Sensing Coverage Table, & Merge with Metadata ####

head(mr.quorsens.ko)

## create table of Sample x KO_ID
mr.quorsens.ko_table<-as.data.frame(dcast(mr.quorsens.ko, SampleID~KO_Function.KEGG, value.var="MR_SumCovPerKO"))
mr.quorsens.ko_table[1:5,1:5]

## merge Quorum Sensing functions with scaled metadata

mr.quorsens.ko.all<-merge(meta.all.scaled,mr.quorsens.ko_table,by="SampleID")
mr.quorsens.ko.all[1:5,1:20]

#### Quorum Sensing & Env Vars - Scatterplots ####

# prep env data for downstream visualizations/analyses
env.vars<-mr.quorsens.ko.all[,c(6,9:10,12:14,36:45)] # subset env var data
env.vars.names<-names(mr.quorsens.ko.all[,c(6,9:10,12:14,36:45)]) # create vector list of variable names

## create function name variable
quorsens.fxn.names<-names(mr.quorsens.ko.all)[names(mr.quorsens.ko.all) %in% c("bisR; LuxR family transcriptional regulator, quorum-sensing system regulator BisR",
                                                                "prgX; HTH-type transcriptional regulator, major conjugation operon repressor")] # pull out names of columns in df that contain "Axis" in name

# first let's do some basic plotting and see what looks important...
# reminder:

ggplot(data=mr.quorsens.ko.all,aes(x=ave.air_temp, y=`bisR; LuxR family transcriptional regulator, quorum-sensing system regulator BisR`, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=3, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Ave Air Temp", y="bisR KO Coverage", title="bisR Coverage x Average Air Temp",subtitle="Using Scaled Climate Data")

ggplot(data=mr.quorsens.ko.all,aes(x=ave.air_temp, y=`prgX; HTH-type transcriptional regulator, major conjugation operon repressor`, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=3, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Ave Air Temp", y="prgX KO Coverage", title="prgX Coverage x Average Air Temp",subtitle="Using Scaled Climate Data")

## Loop to Generate Scatter Plots
### comparing y ~ x where y = quorsens.fxnilia relative coverage and x = whatever environmental variable
head(mr.quorsens.ko.all)
quorsens.scatter.plot.list<-list() # create empty list for each plot to be stored in
quorsens.fxn.names
env.vars.names # sanity check
plotnum <- 1 # counting our model numbers for indexes purposes in the loop

for (i in quorsens.fxn.names) {
  for (j in env.vars.names){
    quorsens.fxn.scat.plot=scat.plot.fxn(mr.quorsens.ko.all,mr.quorsens.ko.all[,j],mr.quorsens.ko.all[,i])
    plot.titled = quorsens.fxn.scat.plot + ggtitle(paste(as.character(i),"Relative Coverage ~",as.character(j),sep=" "),subtitle="Using Scaled Climate Data")
    quorsens.scatter.plot.list[[plotnum]]=plot(plot.titled)
    ggsave(plot.titled,filename = paste("figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/QuorumSensing/SSD_",as.character(i),"~",as.character(j),"_scatterplot.png",sep=""), width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).
    #print(names(clsOI.meta.all[i]))
    #print(clsOI.meta.all[j])
    plotnum<-1+plotnum
  }
}


# check if our loop + functions worked!
quorsens.scatter.plot.list[[16]]

ggplot(data=mr.quorsens.ko.all,aes(x=precip_24hr_accum, y=`prgX; HTH-type transcriptional regulator, major conjugation operon repressor`, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Precip 24 Hrs", y="prgX KO Coverage", title="prgX Coverage x Precip 24 Hrs",subtitle="Using Scaled Climate Data")

quorsens.scatter.plot.list[[2]]

ggplot(data=mr.quorsens.ko.all,aes(x=ave.air_temp, y=`bisR; LuxR family transcriptional regulator, quorum-sensing system regulator BisR`, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Ave Air Temp", y="bisR KO Coverage", title="bisR Coverage x Average Air Temp",subtitle="Using Scaled Climate Data")


# now visualize genes by sampdate
quorsens.date.plot.list<-list() # create empty list for each plot to be stored in

quorsens.fxn.names
env.vars.names
plotnum <- 1 # counting our model numbers for indexes purposes in the loop

for (i in quorsens.fxn.names) {
  quorsens.date.scat.plot=date.mgm.scat.plot.fxn(mr.quorsens.ko.all,mr.quorsens.ko.all[,i])
  date.plot.titled = quorsens.date.scat.plot + ggtitle(paste(as.character(i),"Relative Coverage by Collection Date",sep=" "))
  quorsens.date.plot.list[[plotnum]]=plot(date.plot.titled)
  ggsave(date.plot.titled,filename = paste("figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/QuorumSensing/SSD_",as.character(i),"_by_SampDate_scatterplot.png",sep=""), width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

  #print(names(clsOI.meta.all[i]))
  #print(clsOI.meta.all[j])
  plotnum<-1+plotnum

}

quorsens.date.plot.list[[1]]
ggplot(mr.quorsens.ko.all,aes(x=SampDate, y=`bisR; LuxR family transcriptional regulator, quorum-sensing system regulator BisR`, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x = "",
       y = "",
       shape="",
       color = "R")

quorsens.date.plot.list[[2]]
ggplot(mr.quorsens.ko.all,aes(x=SampDate, y=`prgX; HTH-type transcriptional regulator, major conjugation operon repressor`, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x = "",
       y = "",
       shape="",
       color = "R")


#### What Env Vars Predict Distribution of Quorum Sensing Genes ####
head(mr.quorsens.ko.all)

# first with bisR
# is bisR normally distributed?
shapiro.test(mr.quorsens.ko.all$`bisR; LuxR family transcriptional regulator, quorum-sensing system regulator BisR`) # what is the p-value?
# W = 0.94279, p-value = 0.1882, suggests that this function is normally distributed...
hist(mr.quorsens.ko.all$`bisR; LuxR family transcriptional regulator, quorum-sensing system regulator BisR`)

# visualize Q-Q plot for alpha div
# The Q-Q plot, or quantile-quantile plot, is a graphical tool to help us assess if a set of data plausibly came from some theoretical distribution such as a normal or exponential.
# For example, if we run a statistical analysis that assumes our residuals are normally distributed, we can use a normal Q-Q plot to check that assumption
# more on Q-Q plots here: https://data.library.virginia.edu/understanding-q-q-plots/
# more here too: https://grodri.github.io/glms/notes/c2s9#:~:text=8%20The%20Q%2DQ%20Plot,versus%20quantiles%20of%20a%20distribution.
qqnorm(mr.quorsens.ko.all$`bisR; LuxR family transcriptional regulator, quorum-sensing system regulator BisR`, pch = 1, frame = FALSE)
qqline(mr.quorsens.ko.all$`bisR; LuxR family transcriptional regulator, quorum-sensing system regulator BisR`, col = "red", lwd = 2)

# next let's run some models
grep("bisR;", colnames(mr.quorsens.ko.all)) # find index for column we want

quorsens.glm.step1<-step(glm(formula = round(`bisR; LuxR family transcriptional regulator, quorum-sensing system regulator BisR`) ~ ., data=mr.quorsens.ko.all[,c(6,9:10,12:14,36:45,49)],family="gaussian"))
summary(quorsens.glm.step1)

# now for variable selection...
summary(lm(formula = `bisR; LuxR family transcriptional regulator, quorum-sensing system regulator BisR` ~ precip_24hr_accum+CropLand+Shrub, data=mr.quorsens.ko.all))
summary(lm(formula = `bisR; LuxR family transcriptional regulator, quorum-sensing system regulator BisR` ~ precip_24hr_accum+Shrub, data=mr.quorsens.ko.all))
summary(lm(formula = `bisR; LuxR family transcriptional regulator, quorum-sensing system regulator BisR` ~ precip_24hr_accum+CropLand, data=mr.quorsens.ko.all))

bisR.mod1<-lm(formula = `bisR; LuxR family transcriptional regulator, quorum-sensing system regulator BisR` ~ precip_24hr_accum+CropLand+Shrub, data=mr.quorsens.ko.all)
# compare the analysis of deviance between models with a Chisquared-based estimate
# which variables do we add or remove to improve the model
anova(bisR.mod1,test="Chisq")

# seems like bisR ~ precip_24hr_accum+CropLand is the better model

# final model
quorsens.glm1<-lm(formula = `bisR; LuxR family transcriptional regulator, quorum-sensing system regulator BisR` ~ precip_24hr_accum+CropLand, data=mr.quorsens.ko.all)
quorsens.glm1.sum<-summary(quorsens.glm1, corr = TRUE)
quorsens.glm1.sum
summary(quorsens.glm1, corr = TRUE)$correlation # true correlation between the predictor variables

# confirm that this is the best model by plotting the residuals and seeing their distribution
par(mar=c(1,1,1,1)) # change plot margins
plot(fitted(quorsens.glm1), quorsens.glm1$residuals) # produce residual vs. fitted plot
abline(0,0) # add a horizontal line at 0

# Q-Q norm plot of residuals
qqnorm(quorsens.glm1$residuals)
qqline(quorsens.glm1$residuals)

# density of residuals (how many points have certain residuals)
plot(density(quorsens.glm1$residuals))

# a better way to plot the residuals and their densities
ggplot(data=mr.quorsens.ko.all,aes(quorsens.glm1$residuals)) +
  geom_histogram(aes(y=..density..),color='black',fill='grey') +
  geom_density(alpha=.2,fill='lightpink')

#plot(quorsens.glm1) # another way to plot residuals

# now let's plot the model!
bisR.plot1<-ggplot(data=mr.quorsens.ko.all,aes(x=precip_24hr_accum, y=`bisR; LuxR family transcriptional regulator, quorum-sensing system regulator BisR`))+
  stat_smooth(aes(x=precip_24hr_accum, y=`bisR; LuxR family transcriptional regulator, quorum-sensing system regulator BisR`), method=stats::glm,method.args = list(family = "gaussian")) +
  geom_jitter(aes(shape=Site,color=SampDate),size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Precip. Accum (24 hrs)", y="bisR KO Coverage", title="bisR Coverage x Accumulated Precipitation",subtitle="Using Scaled Climate Data") +
  stat_poly_eq(use_label(c("adj.R2", "P")), formula=y~x,label.y=100)

ggsave(bisR.plot1,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/QuorumSensing/SSD_bisR_PrecipAccum24_scatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

bisR.plot2<-ggplot(data=mr.quorsens.ko.all,aes(x=CropLand, y=`bisR; LuxR family transcriptional regulator, quorum-sensing system regulator BisR`))+
  stat_smooth(aes(x=CropLand, y=`bisR; LuxR family transcriptional regulator, quorum-sensing system regulator BisR`), method=stats::glm,method.args = list(family = "gaussian")) +
  geom_jitter(aes(shape=Site,color=SampDate),size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Crop Land STF", y="bisR KO Coverage", title="bisR Coverage x Crop Land STF",subtitle="Using Scaled Climate Data") +
  stat_poly_eq(use_label(c("adj.R2", "P")), formula=y~x,label.y=100)

ggsave(bisR.plot2,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/QuorumSensing/SSD_bisR_CropLandSTF_scatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

# next with prgX
# is prgX normally distributed?
shapiro.test(mr.quorsens.ko.all$`prgX; HTH-type transcriptional regulator, major conjugation operon repressor`) # what is the p-value?
# W = 0.79123, p-value = 0.0002125, suggests that this function is NOT normally distributed...
hist(mr.quorsens.ko.all$`prgX; HTH-type transcriptional regulator, major conjugation operon repressor`)

# visualize Q-Q plot for alpha div
# The Q-Q plot, or quantile-quantile plot, is a graphical tool to help us assess if a set of data plausibly came from some theoretical distribution such as a normal or exponential.
# For example, if we run a statistical analysis that assumes our residuals are normally distributed, we can use a normal Q-Q plot to check that assumption
# more on Q-Q plots here: https://data.library.virginia.edu/understanding-q-q-plots/
# more here too: https://grodri.github.io/glms/notes/c2s9#:~:text=8%20The%20Q%2DQ%20Plot,versus%20quantiles%20of%20a%20distribution.
qqnorm(mr.quorsens.ko.all$`prgX; HTH-type transcriptional regulator, major conjugation operon repressor`, pch = 1, frame = FALSE)
qqline(mr.quorsens.ko.all$`prgX; HTH-type transcriptional regulator, major conjugation operon repressor`, col = "red", lwd = 2)

grep("prgX;", colnames(mr.quorsens.ko.all)) # find index for column we want

# not sure if we should use a Poisson distribution or a negative binomial distribution for our GLM...
# so we are going to run both models, compare the results, and then calculate the log likelihood to see which is the better model
# for more info on the log likelihood, watch this: https://www.youtube.com/watch?v=8nogLkirA3I
## basically we are trying to find out which distribution fits the data (poisson or negative binomial in this case)
# more on why we use log likelihood and not just the likelihood here: https://math.stackexchange.com/questions/892832/why-we-consider-log-likelihood-instead-of-likelihood-in-gaussian-distribution
# and here: https://www.reddit.com/r/learnmachinelearning/comments/vewczb/why_do_we_maximize_loglikelihood_instead_of/
## likelihoods are a product, & usually you're mulitplying small values together which makes them even smaller and harder to measure
## log likelihoods turn these products into sums, and help to scale our small values without losing any information

# first a glm w/ Poisson distribution
quorsens.glm.step1<-step(glm(formula = round(`prgX; HTH-type transcriptional regulator, major conjugation operon repressor`) ~ ., data=mr.quorsens.ko.all[,c(6,9:10,12:14,36:45,55)],family="poisson"))
summary(quorsens.glm.step1)
quorsens.glm.step1$anova

# then a  NB glm...
quorsens.glm.step2<-step(glm.nb(formula = round(`prgX; HTH-type transcriptional regulator, major conjugation operon repressor`) ~ ., data=mr.quorsens.ko.all[,c(6,9:10,12:14,36:45,55)]))
summary(quorsens.glm.step2)
quorsens.glm.step2$anova

# Poisson and NB results are the same - go with Poisson

# so now let's re run the models and do some variable comparisons
summary(glm(formula = round(`prgX; HTH-type transcriptional regulator, major conjugation operon repressor`) ~ ave.wind_speed+ave.relative_humidity+Developed, data=mr.quorsens.ko.all, family="poisson"))
summary(glm(formula = round(`prgX; HTH-type transcriptional regulator, major conjugation operon repressor`) ~ ave.wind_speed+Developed, data=mr.quorsens.ko.all, family="poisson"))
summary(glm(formula = round(`prgX; HTH-type transcriptional regulator, major conjugation operon repressor`) ~ ave.wind_speed*ave.relative_humidity, data=mr.quorsens.ko.all, family="poisson"))
summary(glm(formula = round(`prgX; HTH-type transcriptional regulator, major conjugation operon repressor`) ~ ave.wind_speed+ave.relative_humidity, data=mr.quorsens.ko.all, family="poisson"))

quorsens.pois1<-glm(formula = round(`prgX; HTH-type transcriptional regulator, major conjugation operon repressor`) ~ ave.wind_speed+Developed+ave.relative_humidity, data=mr.quorsens.ko.all, family="poisson")
quorsens.pois2<-glm(formula = round(`prgX; HTH-type transcriptional regulator, major conjugation operon repressor`) ~ ave.wind_speed+Developed, data=mr.quorsens.ko.all, family="poisson")
quorsens.pois3a<-glm(formula = round(`prgX; HTH-type transcriptional regulator, major conjugation operon repressor`) ~ ave.wind_speed*ave.relative_humidity, data=mr.quorsens.ko.all, family="poisson")
quorsens.pois3b<-glm(formula = round(`prgX; HTH-type transcriptional regulator, major conjugation operon repressor`) ~ ave.wind_speed+ave.relative_humidity, data=mr.quorsens.ko.all, family="poisson")

# first let's see if this can tell us what variables to drop
anova(quorsens.pois1,test="Chisq")

# now we compare the loglikehoods * LRT only works for nested models
# In summary, the LRT tells us if it is beneficial to add parameters to our model, or if we should stick with our simpler model.
# In their most basic form, the hypotheses for the LRT are:
# H0: You should use the nested model.
# Ha: You should use the complex model.
# Thus, if you reject the H0, you can conclude that the complex model is significantly more accurate than the nested model, and you would choose to use the complex model.
# If you fail to reject the H0, you can conclude that the complex model is NOT significantly more accurate than the nested model, so you would choose to use the nested model instead.

## the higher the loglikelihood value, the better the model fits the data
# the smaller the AIC the better it's predictive power is

AIC(quorsens.pois1,quorsens.pois2,quorsens.pois3a,quorsens.pois3b)
# poisson is better... quorsens.pois

# lastly let's try a likelihood ratio test to compare
lrtest(quorsens.pois3b,quorsens.pois2,quorsens.pois1) #

anova(quorsens.pois1,quorsens.pois2,quorsens.pois3a,quorsens.pois3b,test="Chisq") #

# final model for prgX (w/ Poisson):
summary(glm(formula = round(`prgX; HTH-type transcriptional regulator, major conjugation operon repressor`) ~ ave.wind_speed+Developed, data=mr.quorsens.ko.all, family="poisson"), corr = TRUE)

quorsens.glm2<-glm(formula = round(`prgX; HTH-type transcriptional regulator, major conjugation operon repressor`) ~ ave.wind_speed+Developed, data=mr.quorsens.ko.all, family="poisson")
summary(quorsens.glm2, corr = TRUE)$correlation # true correlation between the predictor variables

# calculate McFadden's Pseudo R-squared for model for GLM aka aka Analysis of Deviance
# more here: https://bookdown.org/egarpor/PM-UC3M/glm-deviance.html
with(summary(quorsens.glm2), 1 - deviance/null.deviance)

# confirm that this is the best model by plotting the residuals and seeing their distribution
par(mar=c(1,1,1,1)) # change plot margins
plot(fitted(quorsens.glm2), quorsens.glm2$residuals) # produce residual vs. fitted plot
abline(0,0) # add a horizontal line at 0

# Q-Q norm plot of residuals
qqnorm(quorsens.glm2$residuals) # ave.wind_speed+ave.relative_humidity+Developed had the most normal residuals
qqline(quorsens.glm2$residuals)

# density of residuals (how many points have certain residuals)
plot(density(quorsens.glm2$residuals))

# a better way to plot the residuals and their densities
ggplot(data=mr.quorsens.ko.all,aes(quorsens.glm2$residuals)) +
  geom_histogram(aes(y=..density..),color='black',fill='grey') +
  geom_density(alpha=.2,fill='lightpink')

#plot(quorsens.glm2) # another way to plot residuals

# now let's plot the model!
prgX.plot1<-ggplot(data=mr.quorsens.ko.all,aes(x=ave.wind_speed, y=`prgX; HTH-type transcriptional regulator, major conjugation operon repressor`))+
  geom_smooth(aes(x=ave.wind_speed, y=`prgX; HTH-type transcriptional regulator, major conjugation operon repressor`), method=stats::glm,method.args = list(family = "poisson")) +
  geom_jitter(aes(shape=Site,color=SampDate),size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Ave Wind Speed", y="prgX KO Coverage", title="prgX Coverage x Ave Wind Speed",subtitle="Using Scaled Climate Data") +
  stat_fit_tidy(method = "glm",method.args = list(family=poisson,formula=y~x),
                mapping = aes(label = sprintf("p = %.3g",after_stat(x_p.value))))

ggsave(prgX.plot1,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/QuorumSensing/SSD_prgX_AveWindSpeed_scatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

# prgX.plot2<-ggplot(data=mr.quorsens.ko.all,aes(x=ave.relative_humidity, y=`prgX; HTH-type transcriptional regulator, major conjugation operon repressor`))+
#   geom_smooth(aes(x=ave.relative_humidity, y=`prgX; HTH-type transcriptional regulator, major conjugation operon repressor`), method=stats::glm,method.args = list(family = "poisson")) +
#   geom_jitter(aes(shape=Site,color=SampDate),size=4, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
#                      labels=c("July 2020","August 2020","November 2020",
#                               "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
#   scale_shape_manual(name="Site",values = c(0,1,16,15)) +
#   theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
#         axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Ave Relative Humidity", y="prgX KO Coverage", title="prgX Coverage x Ave Relative Humidity",subtitle="Using Scaled Climate Data") +
#   stat_fit_tidy(method = "glm",method.args = list(family=poisson,formula=y~x),
#                 mapping = aes(label = sprintf("p = %.3g",after_stat(x_p.value))))
#
# ggsave(prgX.plot2,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/QuorumSensing/SSD_prgX_AveRelHumidity_scatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

prgX.plot3<-ggplot(data=mr.quorsens.ko.all,aes(x=Developed, y=`prgX; HTH-type transcriptional regulator, major conjugation operon repressor`))+
  geom_smooth(aes(x=Developed, y=`prgX; HTH-type transcriptional regulator, major conjugation operon repressor`), method=stats::glm,method.args = list(family = "poisson")) +
  geom_jitter(aes(shape=Site,color=SampDate),size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Developed STF", y="prgX KO Coverage", title="prgX Coverage x Developed STF",subtitle="Using Scaled Climate Data") +
  stat_fit_tidy(method = "glm",method.args = list(family=poisson,formula=y~x),
                mapping = aes(label = sprintf("p = %.3g",after_stat(x_p.value))))

ggsave(prgX.plot3,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/QuorumSensing/SSD_prgX_Developed_scatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

#### Compare Variance of Quorum Sensing Genes Across Sites ####

# bisR is normally distributed; used rarefied counts to calculate ShanDiv
# use the following statisitcal tests for variance comparisons
## ANOVA: are variances significantly different between groups
## Tukey test: which groups' variances are significant different from one another
## Levene's test: is variance homogenous aka equal across samples?

head(mr.quorsens.ko.all)

fit1<-aov(`bisR; LuxR family transcriptional regulator, quorum-sensing system regulator BisR`~ Site, data=mr.quorsens.ko.all)
# ANOVA is basically a regression but w/ categorical variables more info here https://www.statology.org/anova-vs-regression/
#pairwise.adonis(mr.quorsens.ko.all$`spoIVCA; site-specific DNA recombinase`, mr.quorsens.ko.all$SampDate, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(fit1)
# p = 0.508

# Tukey test - tells us which groups are significantly different from each other (more here: https://www.r-bloggers.com/2013/06/anova-and-tukeys-test-on-r/)
Tuk1<-TukeyHSD(fit1)
Tuk1$Site
# no sig differences

# fit.0<-aov(DustComplexity ~ as.factor(Elevation), data=mr.quorsens.ko.all)
# summary(fit.0)
# TukeyHSD(fit.0)

# Levene's test - test for homogeneity of variance
## Levene’s test is an inferential statistic used to check if the variances of a variable obtained for two or more groups are equal or not when data comes from a non-normal distribution
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## t test assumes that variances are the same, so Levene's test needs to be non significant
car::leveneTest(`bisR; LuxR family transcriptional regulator, quorum-sensing system regulator BisR`~ Site, data = mr.quorsens.ko.all, center=mean)
# Levene's Test for Homogeneity of Variance (center = mean)
#       Df F value Pr(>F)
# group  3  1.1249 0.3627
# 20
# ^ p value is > 0.05 so we cannot reject the Null hypothesis -- variances are equal

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
# fligner.test(`bisR; LuxR family transcriptional regulator, quorum-sensing system regulator BisR`~ SampDate, data = mr.quorsens.ko.all)
# Fligner-Killeen:med chi-squared = 1.9504, df = 2, p-value = 0.3771
# Which shows that the data do not deviate significantly from homogeneity.

#compare_means(`bisR; LuxR family transcriptional regulator, quorum-sensing system regulator BisR`~ Site, data=mr.quorsens.ko.all, method="anova",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input


#### Do Osmoprotectant Fxns of Interest Correlate? ####

mr.cov.sum.osmo.ko[1:4,]
grep("osmY;", colnames(mr.cov.sum.osmo.ko)) # find index for column we want
grep("opuC;", colnames(mr.cov.sum.osmo.ko)) # find index for column we want
grep("osmB;", colnames(mr.cov.sum.osmo.ko)) # find index for column we want

osmoOI.cor1<-cor(mr.cov.sum.osmo.ko[,c(4:5,8)])
osmo.cor.p <- cor_pmat(mr.cov.sum.osmo.ko[,c(4:5,8)])

ggcorrplot(osmoOI.cor1,method="square",lab = T,p.mat=osmo.cor.p,sig.level = 0.05,
           hc.order=TRUE,outline.color="white",type = "lower")

par(mar=c(1,1,1,1))
plot(mr.cov.sum.osmo.ko$`osmY; hyperosmotically inducible periplasmic protein`,mr.cov.sum.osmo.ko$`opuC; osmoprotectant transport system substrate-binding protein`)
plot(mr.cov.sum.osmo.ko$`osmY; hyperosmotically inducible periplasmic protein`,mr.cov.sum.osmo.ko$`osmB; osmotically inducible lipoprotein OsmB`)
plot(mr.cov.sum.osmo.ko$`osmB; osmotically inducible lipoprotein OsmB`,mr.cov.sum.osmo.ko$`opuC; osmoprotectant transport system substrate-binding protein`)

#### Create Osmoprotectant Coverage Table, then Merge with Metadata ####

head(mr.osmo.ko)

## create table of Sample x KO_ID
mr.osmo.ko_table<-as.data.frame(dcast(mr.osmo.ko, SampleID~KO_Function.KEGG, value.var="MR_SumCovPerKO"))
mr.osmo.ko_table[1:5,1:5]
rownames(mr.osmo.ko_table)<-mr.osmo.ko_table$SampleID

## merge Osmoprotectant functions with scaled metadata

mr.osmo.ko.all<-merge(meta.all.scaled,mr.osmo.ko_table,by="SampleID")
mr.osmo.ko.all[1:5,]

#### Osmoprotectant Fxns & Env Vars - Scatterplots ####
## heatmaps of traits of interest

# prep env data for downstream visualizations/analyses
env.vars<-mr.osmo.ko.all[,c(6,9:10,12:14,36:45)] # subset env var data
env.vars.names<-names(mr.osmo.ko.all[,c(6,9:10,12:14,36:45)]) # create vector list of variable names

## create function name variable
osmo.fxn.names<-names(mr.osmo.ko.all)[names(mr.osmo.ko.all) %in% c("osmY; hyperosmotically inducible periplasmic protein",
                                                                "osmB; osmotically inducible lipoprotein OsmB",
                                                                "opuC; osmoprotectant transport system substrate-binding protein")] # pull out names of columns in df that contain "Axis" in name

## Loop to Generate Scatter Plots
### comparing y ~ x where y = osmoe relative coverage and x = whatever environmental variable
head(mr.osmo.ko.all)
osmo.scatter.plot.list<-list() # create empty list for each plot to be stored in
osmo.fxn.names
env.vars.names # sanity check
plotnum <- 1 # counting our model numbers for indexes purposes in the loop

for (i in osmo.fxn.names) {
  for (j in env.vars.names){
    osmo.fxn.scat.plot=scat.plot.fxn(mr.osmo.ko.all,mr.osmo.ko.all[,j],mr.osmo.ko.all[,i])
    plot.titled = osmo.fxn.scat.plot + ggtitle(paste(as.character(i),"Relative Coverage ~",as.character(j),sep=" "),subtitle="Using Scaled Climate Data")
    osmo.scatter.plot.list[[plotnum]]=plot(plot.titled)
    ggsave(plot.titled,filename = paste("figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/Osmoprotectant/SSD_",as.character(i),"~",as.character(j),"_scatterplot.png",sep=""), width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).
    #print(names(clsOI.meta.all[i]))
    #print(clsOI.meta.all[j])
    plotnum<-1+plotnum
  }
}


# check if our loop + functions worked!
osmo.scatter.plot.list[[2]]

ggplot(data=mr.osmo.ko.all,aes(x=ave.air_temp, y=`osmB; osmotically inducible lipoprotein OsmB`, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Ave Air Temp", y="osmB KO Coverage", title="osmB Coverage x Average Air Temp",subtitle="Using Scaled Climate Data")

osmo.scatter.plot.list[[18]]

ggplot(data=mr.osmo.ko.all,aes(x=ave.wind_speed, y=`opuC; osmoprotectant transport system substrate-binding protein`, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Ave Wind Speed", y="opuC KO Coverage", title="opuC Coverage x Ave Wind Speed",subtitle="Using Scaled Climate Data")


# now visualize genes by sampdate
osmo.date.plot.list<-list() # create empty list for each plot to be stored in

osmo.fxn.names
env.vars.names
plotnum <- 1 # counting our model numbers for indexes purposes in the loop

for (i in osmo.fxn.names) {
  osmo.date.scat.plot=date.mgm.scat.plot.fxn(mr.osmo.ko.all,mr.osmo.ko.all[,i])
  date.plot.titled = osmo.date.scat.plot + ggtitle(paste(as.character(i),"Relative Coverage by Collection Date",sep=" "))
  osmo.date.plot.list[[plotnum]]=plot(date.plot.titled)
  ggsave(date.plot.titled,filename = paste("figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/Osmoprotectant/SSD_",as.character(i),"_by_SampDate_scatterplot.png",sep=""), width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

  #print(names(clsOI.meta.all[i]))
  #print(clsOI.meta.all[j])
  plotnum<-1+plotnum

}

osmo.date.plot.list[[3]]
ggplot(mr.osmo.ko.all,aes(x=SampDate, y=`osmY; hyperosmotically inducible periplasmic protein`, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x = "",
       y = "",
       shape="",
       color = "R")

osmo.date.plot.list[[2]]
ggplot(mr.osmo.ko.all,aes(x=SampDate, y=`opuC; osmoprotectant transport system substrate-binding protein`, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x = "",
       y = "",
       shape="",
       color = "R")

#### What Env Vars Predict Distribution of Osmoprotectant Genes ####
# first with osmY
# is osmY normally distributed?
shapiro.test(mr.osmo.ko.all$`osmY; hyperosmotically inducible periplasmic protein`) # what is the p-value?
# W = 0.80367, p-value = 0.0003384, suggests that this function is NOT normally distributed...
hist(mr.osmo.ko.all$`osmY; hyperosmotically inducible periplasmic protein`)

# visualize Q-Q plot for alpha div
# The Q-Q plot, or quantile-quantile plot, is a graphical tool to help us assess if a set of data plausibly came from some theoretical distribution such as a normal or exponential.
# For example, if we run a statistical analysis that assumes our residuals are normally distributed, we can use a normal Q-Q plot to check that assumption
# more on Q-Q plots here: https://data.library.virginia.edu/understanding-q-q-plots/
# more here too: https://grodri.github.io/glms/notes/c2s9#:~:text=8%20The%20Q%2DQ%20Plot,versus%20quantiles%20of%20a%20distribution.
qqnorm(mr.osmo.ko.all$`osmY; hyperosmotically inducible periplasmic protein`, pch = 1, frame = FALSE)
qqline(mr.osmo.ko.all$`osmY; hyperosmotically inducible periplasmic protein`, col = "red", lwd = 2)

grep("osmY;", colnames(mr.osmo.ko.all)) # find index for column we want

# first a glm w/ Poisson distribution
osmo.glm.step1<-step(glm(formula = round(`osmY; hyperosmotically inducible periplasmic protein`) ~ ., data=mr.osmo.ko.all[,c(6,9:10,12:14,36:45,55)],family="poisson"))
summary(osmo.glm.step1)
osmo.glm.step1$anova

# then a NB glm...
osmo.glm.step2<-step(glm.nb(formula = round(`osmY; hyperosmotically inducible periplasmic protein`) ~ ., data=mr.osmo.ko.all[,c(6,9:10,12:14,36:45,55)]))
summary(osmo.glm.step2)
osmo.glm.step2$anova

# compare which is better, Poisson or NB, with likelihood ratio test and Chisq test (basically same thing, different ways of running it and different information)
# more here: https://stats.stackexchange.com/questions/59879/logistic-regression-anova-chi-square-test-vs-significance-of-coefficients-ano
lrtest(osmo.glm.step1,osmo.glm.step2) # NB does not seem significantly better -- will go with Poisson
anova(osmo.glm.step1,osmo.glm.step2,test="Chisq")

# best model is below, but what happens if we remove precip_24hr_accum -- does our model get better?
summary(glm(formula = round(`osmY; hyperosmotically inducible periplasmic protein`) ~ precip_24hr_accum+ave.wind_speed+ave.v_N.S.wind, data=mr.osmo.ko.all, family="poisson"))

# let's compare the model with and without precipitation, then do likelihood ratio test
osmo.test1<-glm(formula = round(`osmY; hyperosmotically inducible periplasmic protein`) ~ precip_24hr_accum+ave.wind_speed+ave.u_E.W.wind*ave.v_N.S.wind, data=mr.osmo.ko.all, family="poisson")
osmo.test2<-glm(formula = round(`osmY; hyperosmotically inducible periplasmic protein`) ~ precip_24hr_accum+ave.wind_speed+ave.v_N.S.wind, data=mr.osmo.ko.all, family="poisson")
osmo.test3<-glm(formula = round(`osmY; hyperosmotically inducible periplasmic protein`) ~ precip_24hr_accum+ave.wind_speed, data=mr.osmo.ko.all, family="poisson")

# In summary, the LRT tells us if it is beneficial to add parameters to our model, or if we should stick with our simpler model.
# In their most basic form, the hypotheses for the LRT are:
# H0: You should use the nested model.
# Ha: You should use the complex model.
# Thus, if you reject the H0, you can conclude that the complex model is significantly more accurate than the nested model, and you would choose to use the complex model.
# If you fail to reject the H0, you can conclude that the complex model is NOT significantly more accurate than the nested model, so you would choose to use the nested model instead.
# ^ here, complex model is more variables than nested model
lrtest(osmo.test1,osmo.test2,osmo.test3) # put nested model second in lrtest()
AIC(osmo.test1,osmo.test2,osmo.test3)
anova(osmo.test1,osmo.test2,osmo.test3,test="Chisq")

anova(osmo.test2,osmo.test3,test="Chisq") # put nested model second in lrtest()

# final model:
osmo.pois<-glm(formula = round(`osmY; hyperosmotically inducible periplasmic protein`) ~ precip_24hr_accum+ave.wind_speed+ave.v_N.S.wind, data=mr.osmo.ko.all, family="poisson")
osmo.pois.sum<-summary(osmo.pois, corr = TRUE)
osmo.pois.sum
summary(osmo.pois, corr = TRUE)$correlation # true correlation between the predictor variables

# calculate McFadden's Pseudo R-squared for model for GLM aka aka Analysis of Deviance
# more here: https://bookdown.org/egarpor/PM-UC3M/glm-deviance.html
with(summary(osmo.pois), 1 - deviance/null.deviance)

# confirm that this is the best model by plotting the residuals and seeing their distribution
par(mar=c(1,1,1,1)) # change plot margins
plot(fitted(osmo.pois), osmo.pois$residuals) # produce residual vs. fitted plot
abline(0,0) # add a horizontal line at 0

# Q-Q norm plot of residuals
qqnorm(osmo.pois$residuals)
qqline(osmo.pois$residuals)

# density of residuals (how many points have certain residuals)
plot(density(osmo.pois$residuals))

# a better way to plot the residuals and their densities
ggplot(data=mr.osmo.ko.all,aes(osmo.pois$residuals)) +
  geom_histogram(aes(y=..density..),color='black',fill='grey') +
  geom_density(alpha=.2,fill='lightpink')

#plot(osmo.pois) # another way to plot residuals

# now let's plot the model!
osmY.plot1<-ggplot(data=mr.osmo.ko.all,aes(x=precip_24hr_accum, y=`osmY; hyperosmotically inducible periplasmic protein`))+
  stat_smooth(aes(x=precip_24hr_accum, y=round(`osmY; hyperosmotically inducible periplasmic protein`)), method=stats::glm,method.args = list(family = "poisson")) +
  geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Accum Precip (24 hrs)", y="osmY KO Coverage", title="osmY Coverage x Accum. Precip",subtitle="Using Scaled Climate Data") +
  stat_fit_tidy(method = "glm",method.args = list(family=poisson,formula=y~x),
                mapping = aes(label = sprintf("p = %.3g",after_stat(x_p.value))))

ggsave(osmY.plot1,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/Osmoprotectant/SSD_osmY_AccumPrecip_scatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

osmY.plot2<-ggplot(data=mr.osmo.ko.all,aes(x=ave.wind_speed, y=`osmY; hyperosmotically inducible periplasmic protein`))+
  stat_smooth(aes(x=ave.wind_speed, y=round(`osmY; hyperosmotically inducible periplasmic protein`)), method=stats::glm,method.args = list(family = "poisson")) +
  geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Ave Wind Speed", y="osmY KO Coverage", title="osmY Coverage x Ave Wind Speed",subtitle="Using Scaled Climate Data") +
  stat_fit_tidy(method = "glm",method.args = list(family=poisson,formula=y~x),
                mapping = aes(label = sprintf("p = %.3g",after_stat(x_p.value))))

ggsave(osmY.plot2,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/Osmoprotectant/SSD_osmY_AveWindSpeed_scatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

osmY.plot3<-ggplot(data=mr.osmo.ko.all,aes(x=ave.v_N.S.wind, y=`osmY; hyperosmotically inducible periplasmic protein`))+
  stat_smooth(aes(x=ave.v_N.S.wind, y=round(`osmY; hyperosmotically inducible periplasmic protein`)), method=stats::glm,method.args = list(family = "poisson")) +
  geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="North to South Wind Vector Component (v)", y="osmY KO Coverage", title="osmY Coverage x North to South Wind",subtitle="Using Scaled Climate Data") +
  stat_fit_tidy(method = "glm",method.args = list(family=poisson,formula=y~x),
                mapping = aes(label = sprintf("p = %.3g",after_stat(x_p.value))))

ggsave(osmY.plot3,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/Osmoprotectant/SSD_osmY_SouthtoNorth_WindVector_v_scatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

# next with opuC
# is opuC normally distributed?
shapiro.test(mr.osmo.ko.all$`opuC; osmoprotectant transport system substrate-binding protein`) # what is the p-value?
# W = 0.89909, p-value = 0.02059, suggests that this function is NOT normally distributed...
hist(mr.osmo.ko.all$`opuC; osmoprotectant transport system substrate-binding protein`)

# visualize Q-Q plot for alpha div
# The Q-Q plot, or quantile-quantile plot, is a graphical tool to help us assess if a set of data plausibly came from some theoretical distribution such as a normal or exponential.
# For example, if we run a statistical analysis that assumes our residuals are normally distributed, we can use a normal Q-Q plot to check that assumption
# more on Q-Q plots here: https://data.library.virginia.edu/understanding-q-q-plots/
# more here too: https://grodri.github.io/glms/notes/c2s9#:~:text=8%20The%20Q%2DQ%20Plot,versus%20quantiles%20of%20a%20distribution.
qqnorm(mr.osmo.ko.all$`opuC; osmoprotectant transport system substrate-binding protein`, pch = 1, frame = FALSE)
qqline(mr.osmo.ko.all$`opuC; osmoprotectant transport system substrate-binding protein`, col = "red", lwd = 2)

grep("opuC;", colnames(mr.osmo.ko.all)) # find index for column we want

# not sure if we should use a Poisson distribution or a negative binomial distribution for our GLM...
# so we are going to run both models, compare the results, and then calculate the log likelihood to see which is the better model
# for more info on the log likelihood, watch this: https://www.youtube.com/watch?v=8nogLkirA3I
## basically we are trying to find out which distribution fits the data (poisson or negative binomial in this case)
# more on why we use log likelihood and not just the likelihood here: https://math.stackexchange.com/questions/892832/why-we-consider-log-likelihood-instead-of-likelihood-in-gaussian-distribution
# and here: https://www.reddit.com/r/learnmachinelearning/comments/vewczb/why_do_we_maximize_loglikelihood_instead_of/
## likelihoods are a product, & usually you're mulitplying small values together which makes them even smaller and harder to measure
## log likelihoods turn these products into sums, and help to scale our small values without losing any information

# then a glm w/ Poisson distribution
osmo.glm.step3<-step(glm(formula = round(`opuC; osmoprotectant transport system substrate-binding protein`) ~ ., data=mr.osmo.ko.all[,c(6,9:10,12:14,36:45,52)],family="poisson"))
summary(osmo.glm.step3)
osmo.glm.step3$anova

# then a NB glm...
osmo.glm.step4<-step(glm.nb(formula = round(`opuC; osmoprotectant transport system substrate-binding protein`) ~ ., data=mr.osmo.ko.all[,c(6,9:10,12:14,36:45,52)]))
summary(osmo.glm.step4)
osmo.glm.step4$anova

# the smaller the AIC the better it's predictive power is
AIC(osmo.glm.step3,osmo.glm.step4) # Poisson seems to be better since there really isn't a difference

# now we compare the loglikehoods * LRT only works for nested models
# In summary, the LRT tells us if it is beneficial to add parameters to our model, or if we should stick with our simpler model.
# In their most basic form, the hypotheses for the LRT are:
# H0: You should use the nested model.
# Ha: You should use the complex model.
# Thus, if you reject the H0, you can conclude that the complex model is significantly more accurate than the nested model, and you would choose to use the complex model.
# If you fail to reject the H0, you can conclude that the complex model is NOT significantly more accurate than the nested model, so you would choose to use the nested model instead.

## the higher the value, the better the model fits the data

lrtest(osmo.glm.step3,osmo.glm.step4) # put nested model second in lrtest()
# ^ cannot reject Null that Poisson is worse than NB
anova(osmo.glm.step3,osmo.glm.step4,test="Chisq")

# compare the final model between NB and Poisson distributions...
summary(glm.nb(formula = round(`opuC; osmoprotectant transport system substrate-binding protein`) ~ precip_24hr_accum*SaltonSea, data=mr.osmo.ko.all))
summary(glm(formula = round(`opuC; osmoprotectant transport system substrate-binding protein`) ~ precip_24hr_accum*SaltonSea, data=mr.osmo.ko.all,family="poisson"))

# NB model is better -- far less deviance than Poisson model

summary(glm.nb(formula = round(`opuC; osmoprotectant transport system substrate-binding protein`) ~ precip_24hr_accum*ave.air_temp*ave.wind_speed, data=mr.osmo.ko.all))
summary(glm.nb(formula = round(`opuC; osmoprotectant transport system substrate-binding protein`) ~ precip_24hr_accum*ave.wind_speed, data=mr.osmo.ko.all))
summary(glm.nb(formula = round(`opuC; osmoprotectant transport system substrate-binding protein`) ~ precip_24hr_accum*SaltonSea, data=mr.osmo.ko.all))

# precip_24hr_accum*ave.wind_speed has the lowest residual deviance - probably best best

# let's compare the model with and without precipitation, then do likelihood ratio test
osmo.test3<-glm.nb(formula = round(`opuC; osmoprotectant transport system substrate-binding protein`) ~ precip_24hr_accum*ave.air_temp*ave.wind_speed, data=mr.osmo.ko.all)
osmo.test4<-glm.nb(formula = round(`opuC; osmoprotectant transport system substrate-binding protein`) ~ precip_24hr_accum*ave.wind_speed, data=mr.osmo.ko.all)
osmo.test5<-glm.nb(formula = round(`opuC; osmoprotectant transport system substrate-binding protein`) ~ precip_24hr_accum*SaltonSea, data=mr.osmo.ko.all)

lrtest(osmo.test3,osmo.test4,osmo.test5) # put nested model second in lrtest()
AIC(osmo.test3,osmo.test4,osmo.test5)
anova(osmo.test3,osmo.test4,osmo.test5,test="Chisq") #

anova(osmo.test3,test="Chisq") #

# NOTES:
# seems like precip_24hr_accum*ave.wind_speed is the best model due to its lowest residual deviance, second lowest AIC (and not far behind more complex model) and (next line)
## the anova(osmo.test3,test="Chisq) results show that p value for ave wind speed is lower than ave air temp, and the interaction with precip has a slightly lower p val than the interaction between precip * air temp

# final model:
osmo.pois2<-glm.nb(formula = round(`opuC; osmoprotectant transport system substrate-binding protein`) ~ precip_24hr_accum*ave.wind_speed, data=mr.osmo.ko.all)
osmo.pois2.sum<-summary(osmo.pois2, corr = TRUE)
osmo.pois2.sum
summary(osmo.pois2, corr = TRUE)$correlation # true correlation between the predictor variables

# calculate McFadden's Pseudo R-squared for model for GLM aka aka Analysis of Deviance
# more here: https://bookdown.org/egarpor/PM-UC3M/glm-deviance.html
with(summary(osmo.pois2), 1 - deviance/null.deviance)

# now we compar the loglikehoods and see which one is greater
## the higher the value, the better the model fits the data
logLik(osmo.pois2)

# the smaller the AIC the better it's predictive power is
AIC(osmo.pois2)

# confirm that this is the best model by plotting the residuals and seeing their distribution
par(mar=c(1,1,1,1)) # change plot margins
plot(fitted(osmo.pois2), osmo.pois2$residuals) # produce residual vs. fitted plot
abline(0,0) # add a horizontal line at 0

# Q-Q norm plot of residuals
qqnorm(osmo.pois2$residuals)
qqline(osmo.pois2$residuals)

# density of residuals (how many points have certain residuals)
plot(density(osmo.pois2$residuals))

# a better way to plot the residuals and their densities
ggplot(data=mr.osmo.ko.all,aes(osmo.pois2$residuals)) +
  geom_histogram(aes(y=..density..),color='black',fill='grey') +
  geom_density(alpha=.2,fill='lightpink')

#plot(osmo.pois2) # another way to plot residuals

# now let's plot the model!
opuC.plot1<-ggplot(data=mr.osmo.ko.all,aes(x=precip_24hr_accum, y=`opuC; osmoprotectant transport system substrate-binding protein`))+
  stat_smooth(aes(x=precip_24hr_accum, y=round(`opuC; osmoprotectant transport system substrate-binding protein`)), method=stats::glm) +
  geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Precip. Accum (24 hrs)", y="opuC KO Coverage", title="opuC Coverage x Precipitation Accumulated",subtitle="Using Scaled Climate Data") +
  stat_fit_tidy(method = "glm",method.args = list(formula=y~x),
                mapping = aes(label = sprintf("p = %.3g",after_stat(x_p.value))),label.y=26)

ggsave(opuC.plot1,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/Osmoprotectant/SSD_opuC_PrecipAccum24hr_scatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

opuC.plot2<-ggplot(data=mr.osmo.ko.all,aes(x=ave.wind_speed, y=`opuC; osmoprotectant transport system substrate-binding protein`))+
  stat_smooth(aes(x=ave.wind_speed, y=round(`opuC; osmoprotectant transport system substrate-binding protein`)), method=stats::glm) +
  geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Ave Wind Speed", y="opuC KO Coverage", title="opuC Coverage x Salton Sea STF",subtitle="Using Scaled Climate Data") +
  stat_fit_tidy(method = "glm",method.args = list(formula=y~x),
                mapping = aes(label = sprintf("p = %.3g",after_stat(x_p.value))),label.y=26)

ggsave(opuC.plot2,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/Osmoprotectant/SSD_opuC_AveWindSpeed_scatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

# next with osmB
# is osmB normally distributed?
shapiro.test(mr.osmo.ko.all$`osmB; osmotically inducible lipoprotein OsmB`) # what is the p-value?
# W = 0.72296, p-value = 2.057e-05, suggests that this function is NOT normally distributed...
hist(mr.osmo.ko.all$`osmB; osmotically inducible lipoprotein OsmB`)

# visualize Q-Q plot for alpha div
# The Q-Q plot, or quantile-quantile plot, is a graphical tool to help us assess if a set of data plausibly came from some theoretical distribution such as a normal or exponential.
# For example, if we run a statistical analysis that assumes our residuals are normally distributed, we can use a normal Q-Q plot to check that assumption
# more on Q-Q plots here: https://data.library.virginia.edu/understanding-q-q-plots/
# more here too: https://grodri.github.io/glms/notes/c2s9#:~:text=8%20The%20Q%2DQ%20Plot,versus%20quantiles%20of%20a%20distribution.
qqnorm(mr.osmo.ko.all$`osmB; osmotically inducible lipoprotein OsmB`, pch = 1, frame = FALSE)
qqline(mr.osmo.ko.all$`osmB; osmotically inducible lipoprotein OsmB`, col = "red", lwd = 2)

grep("osmB;", colnames(mr.osmo.ko.all)) # find index for column we want

# not sure if we should use a Poisson distribution or a negative binomial distribution for our GLM...
# so we are going to run both models, compare the results, and then calculate the log likelihood to see which is the better model
# for more info on the log likelihood, watch this: https://www.youtube.com/watch?v=8nogLkirA3I
## basically we are trying to find out which distribution fits the data (poisson or negative binomial in this case)
# more on why we use log likelihood and not just the likelihood here: https://math.stackexchange.com/questions/892832/why-we-consider-log-likelihood-instead-of-likelihood-in-gaussian-distribution
# and here: https://www.reddit.com/r/learnmachinelearning/comments/vewczb/why_do_we_maximize_loglikelihood_instead_of/
## likelihoods are a product, & usually you're mulitplying small values together which makes them even smaller and harder to measure
## log likelihoods turn these products into sums, and help to scale our small values without losing any information

# first a glm w/ Poisson distribution
osmo.glm.step5<-step(glm(formula = round(`osmB; osmotically inducible lipoprotein OsmB`) ~ ., data=mr.osmo.ko.all[,c(6,9:10,12:14,36:45,53)],family="poisson"))
summary(osmo.glm.step5)
osmo.glm.step5$anova

# then a NB glm...
osmo.glm.step6<-step(glm.nb(formula = round(`osmB; osmotically inducible lipoprotein OsmB`) ~ ., data=mr.osmo.ko.all[,c(6,9:10,12:14,36:45,53)]))
summary(osmo.glm.step6)
osmo.glm.step6$anova

# Poisson and NB results are the same - go with Poisson

summary(glm(formula = round(`osmB; osmotically inducible lipoprotein OsmB`) ~ precip_24hr_accum+CropLand+SaltonSea, data=mr.osmo.ko.all,family="poisson"))
summary(glm(formula = round(`osmB; osmotically inducible lipoprotein OsmB`) ~ CropLand+SaltonSea, data=mr.osmo.ko.all,family="poisson"))

# let's compare the model with and without precipitation, then do likelihood ratio test
osmo.test6<-glm(formula = round(`osmB; osmotically inducible lipoprotein OsmB`) ~ precip_24hr_accum+CropLand+SaltonSea, data=mr.osmo.ko.all,family="poisson")
osmo.test7<-glm(formula = round(`osmB; osmotically inducible lipoprotein OsmB`) ~ CropLand+SaltonSea, data=mr.osmo.ko.all,family="poisson")

# now we compare the loglikehoods * LRT only works for nested models!
# In summary, the LRT tells us if it is beneficial to add parameters to our model, or if we should stick with our simpler model.
# In their most basic form, the hypotheses for the LRT are:
# H0: You should use the nested model.
# Ha: You should use the complex model.
# Thus, if you reject the H0, you can conclude that the complex model is significantly more accurate than the nested model, and you would choose to use the complex model.
# If you fail to reject the H0, you can conclude that the complex model is NOT significantly more accurate than the nested model, so you would choose to use the nested model instead.

lrtest(osmo.test6,osmo.test7) # put nested model second in lrtest()
AIC(osmo.test6,osmo.test7)
anova(osmo.test6,osmo.test7,test="Chisq")

# p value for lrtest() is < 0.05 -- reject H0, use more complex model!

# final model below appears to be the best model
osmB.pois<-glm(formula = round(`osmB; osmotically inducible lipoprotein OsmB`) ~ precip_24hr_accum+CropLand+SaltonSea, data=mr.osmo.ko.all,family="poisson")
osmB.pois.sum<-summary(osmB.pois,corr=TRUE)
osmB.pois.sum
summary(osmB.pois, corr = TRUE)$correlation # true correlation between the predictor variables

# calculate McFadden's Pseudo R-squared for model for GLM aka aka Analysis of Deviance
# more here: https://bookdown.org/egarpor/PM-UC3M/glm-deviance.html
with(summary(osmB.pois), 1 - deviance/null.deviance)

# now we compar the loglikehoods and see which one is greater
## the higher the value, the better the model fits the data
logLik(osmB.pois) #
AIC(osmB.pois)

osmB.plot1<-ggplot(data=mr.osmo.ko.all,aes(x=precip_24hr_accum, y=`osmB; osmotically inducible lipoprotein OsmB`))+
  stat_smooth(aes(x=precip_24hr_accum, y=round(`osmB; osmotically inducible lipoprotein OsmB`)), method=stats::glm) +
  geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Accum Precip (24 hr)", y="osmB KO Coverage", title="osmB Coverage x Accum. Precipitation (24 hrs)",subtitle="Using Scaled Climate Data") +
  stat_fit_tidy(method = "glm",method.args = list(formula=y~x),
                mapping = aes(label = sprintf("p = %.3g",after_stat(x_p.value))),label.y=26)

ggsave(osmB.plot1,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/Osmoprotectant/SSD_osmB_AccumPrecip_scatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

osmB.plot2<-ggplot(data=mr.osmo.ko.all,aes(x=CropLand, y=`osmB; osmotically inducible lipoprotein OsmB`))+
  stat_smooth(aes(x=CropLand, y=round(`osmB; osmotically inducible lipoprotein OsmB`)), method=stats::glm) +
  geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Crop Land STF", y="osmB KO Coverage", title="osmB Coverage x Crop Land STF",subtitle="Using Scaled Climate Data") +
  stat_fit_tidy(method = "glm",method.args = list(formula=y~x),
                mapping = aes(label = sprintf("p = %.3g",after_stat(x_p.value))),label.y=26)

ggsave(osmB.plot2,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/Osmoprotectant/SSD_osmB_CropLand_scatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

osmB.plot3<-ggplot(data=mr.osmo.ko.all,aes(x=SaltonSea, y=`osmB; osmotically inducible lipoprotein OsmB`))+
  stat_smooth(aes(x=SaltonSea, y=round(`osmB; osmotically inducible lipoprotein OsmB`)), method=stats::glm) +
  geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Salton Sea STF", y="osmB KO Coverage", title="osmB Coverage x Salton Sea STF",subtitle="Using Scaled Climate Data") +
  stat_fit_tidy(method = "glm",method.args = list(formula=y~x),
                mapping = aes(label = sprintf("p = %.3g",after_stat(x_p.value))),label.y=26)

ggsave(osmB.plot3,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/Osmoprotectant/SSD_osmB_SaltonSea_scatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

#### Compare Variance of Osmoprotectant Genes Across Sites ####

# osmY is normally distributed; used rarefied counts to calculate ShanDiv
# use the following statisitcal tests for variance comparisons
## ANOVA: are variances significantly different between groups
## Tukey test: which groups' variances are significant different from one another
## Levene's test: is variance homogenous aka equal across samples?

head(mr.osmo.ko.all)

fit1<-aov(`osmY; hyperosmotically inducible periplasmic protein` ~ Site, data=mr.osmo.ko.all)
# ANOVA is basically a regression but w/ categorical variables more info here https://www.statology.org/anova-vs-regression/
#pairwise.adonis(mr.osmo.ko.all$`osmY; hyperosmotically inducible periplasmic protein`, mr.osmo.ko.all$SampDate, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(fit1)
# p = 0.662

# Tukey test - tells us which groups are significantly different from each other (more here: https://www.r-bloggers.com/2013/06/anova-and-tukeys-test-on-r/)
Tuk1<-TukeyHSD(fit1)
Tuk1$Site
# no sig differences

# fit.0<-aov(DustComplexity ~ as.factor(Elevation), data=mr.osmo.ko.all)
# summary(fit.0)
# TukeyHSD(fit.0)

# Levene's test - test for homogeneity of variance
## Levene’s test is an inferential statistic used to check if the variances of a variable obtained for two or more groups are equal or not when data comes from a non-normal distribution
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## t test assumes that variances are the same, so Levene's test needs to be non significant
car::leveneTest(`osmY; hyperosmotically inducible periplasmic protein` ~ Site, data = mr.osmo.ko.all, center=mean)
# Levene's Test for Homogeneity of Variance (center = mean)
#       Df F value Pr(>F)
# group  3  4.3757 0.01597 *
# 20
# ^ p value is < 0.05 so we can reject the Null hypothesis -- variances are NOT equal

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
# fligner.test(`osmY; hyperosmotically inducible periplasmic protein` ~ SampDate, data = mr.osmo.ko.all)
# Fligner-Killeen:med chi-squared = 1.9504, df = 2, p-value = 0.3771
# Which shows that the data do not deviate significantly from homogeneity.

#compare_means(`osmY; hyperosmotically inducible periplasmic protein` ~ Site, data=mr.osmo.ko.all, method="anova",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input

#### Do Antibiotic Resistance Genes Fxns of Interest Correlate? ####

mr.cov.sum.ARG.ko[1:4,]
#grep("osmY;", colnames(mr.cov.sum.ARG.ko)) # find index for column we want
#grep("opuC;", colnames(mr.cov.sum.ARG.ko)) # find index for column we want
#grep("osmB;", colnames(mr.cov.sum.ARG.ko)) # find index for column we want

ARGOI.cor1<-cor(mr.cov.sum.ARG.ko[,-1])
ARG.cor.p <- cor_pmat(mr.cov.sum.ARG.ko[,-1])

ggcorrplot(ARGOI.cor1,method="square",lab = T,p.mat=ARG.cor.p,sig.level = 0.05,
           hc.order=TRUE,outline.color="white",type = "lower")

par(mar=c(1,1,1,1))
plot(mr.cov.sum.ARG.ko$`macrolide phosphotransferase`, mr.cov.sum.ARG.ko$`beta-lactamase class D OXA-9`)

#### Create Antibiotic Resistance Genes Coverage Table, then Merge with Metadata ####

head(mr.ARG.ko)

## create table of Sample x KO_ID
mr.ARG.ko_table<-as.data.frame(dcast(mr.ARG.ko, SampleID~KO_Function.KEGG, value.var="MR_SumCovPerKO"))
mr.ARG.ko_table[1:5,1:5]
rownames(mr.ARG.ko_table)<-mr.ARG.ko_table$SampleID

## merge Antibiotic Resistance Genes functions with scaled metadata

mr.ARG.ko.all<-merge(meta.all.scaled,mr.ARG.ko_table,by="SampleID")
mr.ARG.ko.all[1:5,]

#### Antibiotic Resistance Genes Fxns & Env Vars - Scatterplots ####
## heatmaps of traits of interest

# prep env data for downstream visualizations/analyses
env.vars<-mr.ARG.ko.all[,c(6,9:10,12:14,36:45)] # subset env var data
env.vars.names<-names(mr.ARG.ko.all[,c(6,9:10,12:14,36:45)]) # create vector list of variable names

## create function name variable
ARG.fxn.names<-names(mr.ARG.ko.all)[names(mr.ARG.ko.all) %in% c("aminoglycoside 6'-N-acetyltransferase I","beta-lactamase class D OXA-9",
                                                                   "dihydrofolate reductase DfrA","macrolide phosphotransferase",
                                                                   "ribosomal protection tetracycline resistance protein")] # pull out names of columns in df that contain "Axis" in name

## Loop to Generate Scatter Plots
### comparing y ~ x where y = ARGe relative coverage and x = whatever environmental variable
head(mr.ARG.ko.all)
ARG.scatter.plot.list<-list() # create empty list for each plot to be stored in
ARG.fxn.names
env.vars.names # sanity check
plotnum <- 1 # counting our model numbers for indexes purposes in the loop

for (i in ARG.fxn.names) {
  for (j in env.vars.names){
    ARG.fxn.scat.plot=scat.plot.fxn(mr.ARG.ko.all,mr.ARG.ko.all[,j],mr.ARG.ko.all[,i])
    plot.titled = ARG.fxn.scat.plot + ggtitle(paste(as.character(i),"Relative Coverage ~",as.character(j),sep=" "),subtitle="Using Scaled Climate Data")
    ARG.scatter.plot.list[[plotnum]]=plot(plot.titled)
    ggsave(plot.titled,filename = paste("figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/AntibioticResistanceGenes/SSD_",as.character(i),"~",as.character(j),"_scatterplot.png",sep=""), width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).
    #print(names(clsOI.meta.all[i]))
    #print(clsOI.meta.all[j])
    plotnum<-1+plotnum
  }
}


# check if our loop + functions worked!
ARG.scatter.plot.list[[2]]

ggplot(data=mr.ARG.ko.all,aes(x=ave.air_temp, y=`aminoglycoside 6'-N-acetyltransferase I`, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Ave Air Temp", y="aminoglycoside KO Coverage", title=" Coverage x Average Air Temp",subtitle="Using Scaled Climate Data")

ARG.scatter.plot.list[[18]]

ggplot(data=mr.ARG.ko.all,aes(x=ave.air_temp, y=`beta-lactamase class D OXA-9`, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Ave Air Temp", y="beta-lactamase KO Coverage", title=" Coverage x Ave Wind Speed",subtitle="Using Scaled Climate Data")


# now visualize genes by sampdate
ARG.date.plot.list<-list() # create empty list for each plot to be stored in

ARG.fxn.names
env.vars.names
plotnum <- 1 # counting our model numbers for indexes purposes in the loop

for (i in ARG.fxn.names) {
  ARG.date.scat.plot=date.mgm.scat.plot.fxn(mr.ARG.ko.all,mr.ARG.ko.all[,i])
  date.plot.titled = ARG.date.scat.plot + ggtitle(paste(as.character(i),"Relative Coverage by Collection Date",sep=" "))
  ARG.date.plot.list[[plotnum]]=plot(date.plot.titled)
  ggsave(date.plot.titled,filename = paste("figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/AntibioticResistanceGenes/SSD_",as.character(i),"_by_SampDate_scatterplot.png",sep=""), width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

  #print(names(clsOI.meta.all[i]))
  #print(clsOI.meta.all[j])
  plotnum<-1+plotnum

}

ARG.date.plot.list[[3]]
ggplot(mr.ARG.ko.all,aes(x=SampDate, y=mr.ARG.ko.all$`beta-lactamase class D OXA-9`, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x = "",
       y = "",
       shape="",
       color = "R")

#### What Env Vars Predict Distribution of Antibiotic Resistance Genes Genes ####
# first with beta-lactamase
# is beta-lactamase normally distributed?
shapiro.test(mr.ARG.ko.all$`beta-lactamase class D OXA-9`) # what is the p-value?
# W = 0.61628, p-value = 9.437e-07, suggests that this function is NOT normally distributed...
hist(mr.ARG.ko.all$`beta-lactamase class D OXA-9`)

# visualize Q-Q plot for alpha div
# The Q-Q plot, or quantile-quantile plot, is a graphical tool to help us assess if a set of data plausibly came from some theoretical distribution such as a normal or exponential.
# For example, if we run a statistical analysis that assumes our residuals are normally distributed, we can use a normal Q-Q plot to check that assumption
# more on Q-Q plots here: https://data.library.virginia.edu/understanding-q-q-plots/
# more here too: https://grodri.github.io/glms/notes/c2s9#:~:text=8%20The%20Q%2DQ%20Plot,versus%20quantiles%20of%20a%20distribution.
qqnorm(mr.ARG.ko.all$`beta-lactamase class D OXA-9`, pch = 1, frame = FALSE)
qqline(mr.ARG.ko.all$`beta-lactamase class D OXA-9`, col = "red", lwd = 2)

grep("beta-lactamase", colnames(mr.ARG.ko.all)) # find index for column we want

# first a glm w/ Poisson distribution
ARG.glm.step1<-step(glm(formula = round(mr.ARG.ko.all$`beta-lactamase class D OXA-9`) ~ ., data=mr.ARG.ko.all[,c(6,9:10,12:14,36:45,50)],family="poisson"))
summary(ARG.glm.step1)
ARG.glm.step1$anova

# then a NB glm...
ARG.glm.step2<-step(glm.nb(formula = round(mr.ARG.ko.all$`beta-lactamase class D OXA-9`) ~ ., data=mr.ARG.ko.all[,c(6,9:10,12:14,36:45,50)]))
summary(ARG.glm.step2)
ARG.glm.step2$anova

# compare which is better, Poisson or NB, with likelihood ratio test and Chisq test (basically same thing, different ways of running it and different information)
# more here: https://stats.stackexchange.com/questions/59879/logistic-regression-anova-chi-square-test-vs-significance-of-coefficients-ano
lrtest(ARG.glm.step1,ARG.glm.step2) # Poisson seems to be nearly same as NB - go with Poisson
anova(ARG.glm.step1,ARG.glm.step2,test="Chisq")

# best model is below, but what happens if we play with the variables -- does our model get better?
summary(glm(formula = round(mr.ARG.ko.all$`beta-lactamase class D OXA-9`) ~ ave.air_temp+ave.relative_humidity+OpenWater+Shrub, data=mr.ARG.ko.all, family="poisson"))

# let's compare the model with and without precipitation, then do likelihood ratio test
ARG.test1<-glm(formula = round(mr.ARG.ko.all$`beta-lactamase class D OXA-9`) ~ ave.air_temp+ave.relative_humidity+OpenWater+Shrub, data=mr.ARG.ko.all, family="poisson")
ARG.test2<-glm(formula = round(mr.ARG.ko.all$`beta-lactamase class D OXA-9`) ~ ave.air_temp+ave.relative_humidity+Shrub, data=mr.ARG.ko.all, family="poisson")
ARG.test3<-glm(formula = round(mr.ARG.ko.all$`beta-lactamase class D OXA-9`) ~ ave.air_temp+OpenWater+Shrub, data=mr.ARG.ko.all, family="poisson")
ARG.test4<-glm(formula = round(mr.ARG.ko.all$`beta-lactamase class D OXA-9`) ~ ave.air_temp+Shrub, data=mr.ARG.ko.all, family="poisson")

# In summary, the LRT tells us if it is beneficial to add parameters to our model, or if we should stick with our simpler model.
# In their most basic form, the hypotheses for the LRT are:
# H0: You should use the nested model.
# Ha: You should use the complex model.
# Thus, if you reject the H0, you can conclude that the complex model is significantly more accurate than the nested model, and you would choose to use the complex model.
# If you fail to reject the H0, you can conclude that the complex model is NOT significantly more accurate than the nested model, so you would choose to use the nested model instead.
# ^ here, complex model is more variables than nested model
lrtest(ARG.test1,ARG.test2,ARG.test3,ARG.test4) # put nested model second in lrtest()
AIC(ARG.test1,ARG.test2,ARG.test3,ARG.test4)
anova(ARG.test1,ARG.test2,ARG.test3,ARG.test4,test="Chisq")

anova(ARG.test1,test="Chisq") # put nested model second in lrtest()
ARG.test1a<-glm(formula = round(mr.ARG.ko.all$`beta-lactamase class D OXA-9`) ~ ave.air_temp+OpenWater+ave.relative_humidity+Shrub, data=mr.ARG.ko.all, family="poisson")
anova(ARG.test1a,test="Chisq") # relative humidity should probably be removed from final model

# final model:
ARG.pois<-glm(formula = round(mr.ARG.ko.all$`beta-lactamase class D OXA-9`) ~ ave.air_temp+OpenWater+Shrub, data=mr.ARG.ko.all, family="poisson")
ARG.pois.sum<-summary(ARG.pois, corr = TRUE)
ARG.pois.sum
summary(ARG.pois, corr = TRUE)$correlation # true correlation between the predictor variables

# calculate McFadden's Pseudo R-squared for model for GLM aka aka Analysis of Deviance
# more here: https://bookdown.org/egarpor/PM-UC3M/glm-deviance.html
with(summary(ARG.pois), 1 - deviance/null.deviance)

# confirm that this is the best model by plotting the residuals and seeing their distribution
par(mar=c(1,1,1,1)) # change plot margins
plot(fitted(ARG.pois), ARG.pois$residuals) # produce residual vs. fitted plot
abline(0,0) # add a horizontal line at 0

# Q-Q norm plot of residuals
qqnorm(ARG.pois$residuals)
qqline(ARG.pois$residuals)

# density of residuals (how many points have certain residuals)
plot(density(ARG.pois$residuals))

# a better way to plot the residuals and their densities
ggplot(data=mr.ARG.ko.all,aes(ARG.pois$residuals)) +
  geom_histogram(aes(y=..density..),color='black',fill='grey') +
  geom_density(alpha=.2,fill='lightpink')

#plot(ARG.pois) # another way to plot residuals

# now let's plot the model!
betalac.plot1<-ggplot(data=mr.ARG.ko.all,aes(x=ave.air_temp, y=mr.ARG.ko.all$`beta-lactamase class D OXA-9`))+
  stat_smooth(aes(x=ave.air_temp, y=round(mr.ARG.ko.all$`beta-lactamase class D OXA-9`)), method=stats::glm,method.args = list(family = "poisson")) +
  geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Ave Air Temp", y="beta-lactamase KO Coverage", title="beta-lactamase Coverage x Ave Air Temp",subtitle="Using Scaled Climate Data") +
  stat_fit_tidy(method = "glm",method.args = list(family=poisson,formula=y~x),
                mapping = aes(label = sprintf("p = %.3g",after_stat(x_p.value))))

ggsave(betalac.plot1,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/AntibioticResistanceGenes/SSD_BetaLactamase_AveAirTemp_scatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

betalac.plot2<-ggplot(data=mr.ARG.ko.all,aes(x=OpenWater, y=mr.ARG.ko.all$`beta-lactamase class D OXA-9`))+
  stat_smooth(aes(x=OpenWater, y=round(mr.ARG.ko.all$`beta-lactamase class D OXA-9`)), method=stats::glm,method.args = list(family = "poisson")) +
  geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Open Water STF", y="beta-lactamase KO Coverage", title="beta-lactamase Coverage x Open Water STF",subtitle="Using Scaled Climate Data") +
  stat_fit_tidy(method = "glm",method.args = list(family=poisson,formula=y~x),
                mapping = aes(label = sprintf("p = %.3g",after_stat(x_p.value))))

ggsave(betalac.plot2,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/AntibioticResistanceGenes/SSD_BetaLactamase_OpenWater_scatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

betalac.plot3<-ggplot(data=mr.ARG.ko.all,aes(x=Shrub, y=mr.ARG.ko.all$`beta-lactamase class D OXA-9`))+
  stat_smooth(aes(x=Shrub, y=round(mr.ARG.ko.all$`beta-lactamase class D OXA-9`)), method=stats::glm,method.args = list(family = "poisson")) +
  geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Shrub STF", y="beta-lactamase KO Coverage", title="beta-lactamase Coverage x Shrub STF",subtitle="Using Scaled Climate Data") +
  stat_fit_tidy(method = "glm",method.args = list(family=poisson,formula=y~x),
                mapping = aes(label = sprintf("p = %.3g",after_stat(x_p.value))))

ggsave(betalac.plot3,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/AntibioticResistanceGenes/SSD_BetaLactamase_Shrub_scatterplot.png", width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

#### Compare Variance of Antibiotic Resistance Genes Genes Across Sites ####

# beta-lactamase is NOT normally distributed; used rarefied counts to calculate ShanDiv
# use the following statisitcal tests for variance comparisons
## ANOVA: are variances significantly different between groups
## Tukey test: which groups' variances are significant different from one another
## Levene's test: is variance homogenous aka equal across samples?

head(mr.ARG.ko.all)

fit1<-aov(mr.ARG.ko.all$`beta-lactamase class D OXA-9` ~ Site, data=mr.ARG.ko.all)
# ANOVA is basically a regression but w/ categorical variables more info here https://www.statology.org/anova-vs-regression/
#pairwise.adonis(mr.ARG.ko.all$`beta-lactamase class D OXA-9`, mr.ARG.ko.all$SampDate, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(fit1)
# p = 0.593

# Tukey test - tells us which groups are significantly different from each other (more here: https://www.r-bloggers.com/2013/06/anova-and-tukeys-test-on-r/)
Tuk1<-TukeyHSD(fit1)
Tuk1$Site
# no sig differences

# fit.0<-aov(DustComplexity ~ as.factor(Elevation), data=mr.ARG.ko.all)
# summary(fit.0)
# TukeyHSD(fit.0)

# Levene's test - test for homogeneity of variance
## Levene’s test is an inferential statistic used to check if the variances of a variable obtained for two or more groups are equal or not when data comes from a non-normal distribution
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## t test assumes that variances are the same, so Levene's test needs to be non significant
car::leveneTest(mr.ARG.ko.all$`beta-lactamase class D OXA-9` ~ Site, data = mr.ARG.ko.all, center=mean)
# Levene's Test for Homogeneity of Variance (center = mean)
#       Df F value Pr(>F)
# group  3  4.3757 0.03319 *
# 20
# ^ p value is < 0.05 so we can reject the Null hypothesis -- variances are NOT equal

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
# fligner.test(mr.ARG.ko.all$`beta-lactamase class D OXA-9` ~ SampDate, data = mr.ARG.ko.all)
# Fligner-Killeen:med chi-squared = 1.9504, df = 2, p-value = 0.3771
# Which shows that the data do not deviate significantly from homogeneity.

#compare_means(mr.ARG.ko.all$`beta-lactamase class D OXA-9` ~ Site, data=mr.ARG.ko.all, method="anova",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input

#### Does the Coverage of Specific Traits Correlate with Each Other ####

## first create house keeping name object
housekeep.fxn.name<-"recA; recombination protein RecA"

# then let's review which functions we want to compare
lps.fxn.names
spor.fxn.names
tempshock.fxn.names
UVrep.fxn.names
quorsens.fxn.names
osmo.fxn.names<-names(mr.osmo.ko.all)[names(mr.osmo.ko.all) %in% c("osmY; hyperosmotically inducible periplasmic protein",
                                                                   "opuC; osmoprotectant transport system substrate-binding protein")] # pull out names of columns in df that contain "Axis" in name
# ^ recreate osmo.fxn.names since we do not discuss osmB in manuscript ("osmB; osmotically inducible lipoprotein OsmB")

# check if functions are in the same order by SampleID
mr.lps.ko_table$SampleID %in% mr.spor.ko.all$SampleID
mr.spor.ko.all$SampleID %in% mr.tempshock.ko.all$SampleID
mr.tempshock.ko.all$SampleID %in% mr.UVrep.ko.all$SampleID
mr.UVrep.ko.all$SampleID %in% meta.all.scaled$SampleID
mr.quorsens.ko.all$SampleID %in% meta.all.scaled$SampleID
mr.osmo.ko.all$SampleID %in% meta.all.scaled$SampleID
mr.HK.ko.all$SampleID %in% meta.all.scaled$SampleID

# pull out functions & coverage by category of functions
lpsOI<-mr.lps.ko.all[names(mr.lps.ko.all) %in% lps.fxn.names]
#lpsOI$Category<-"LPS"

sporOI<-mr.spor.ko.all[names(mr.spor.ko.all) %in% spor.fxn.names]
#sporOI$Category<-"Spore"

TSOI<-mr.tempshock.ko.all[names(mr.tempshock.ko.all) %in% tempshock.fxn.names]
#TSOI$Category<-"Temp Shock"

UVrepOI<-mr.UVrep.ko.all[names(mr.UVrep.ko.all) %in% UVrep.fxn.names]
#UVrepOI$Category<-"UV Damage Repair"

QSOI<-mr.quorsens.ko.all[names(mr.quorsens.ko.all) %in% quorsens.fxn.names]
#QSOI$Category<-"Quorum Sensing"

osmoOI<-mr.osmo.ko.all[names(mr.osmo.ko.all) %in% osmo.fxn.names]
#osmoOI$Category<-"Osmoprotectant"

# create dataframe of functions of interest
FxnsOI<-cbind(SampleID=meta.all.scaled$SampleID,lpsOI,sporOI,TSOI,UVrepOI,QSOI,osmoOI)
head(FxnsOI)
rownames(FxnsOI)<-FxnsOI$SampleID

# Visualize with a corrplot
cor.fxns1 <- cor(FxnsOI[,-1], method='pearson')
cor.fxns1

p.fxns <- cor_pmat(FxnsOI[,-1])

symnum(cor.fxns1)

ggcorrplot(cor.fxns1,method="square",lab = T,p.mat=p.fxns,
           hc.order=TRUE,outline.color="white",type = "lower")

png('figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/SSD_MGMs_FxnsOfInterest_CorrPlot.png',width = 1100, height = 1100)
ggcorrplot(cor.fxns1,p.mat=p.fxns,hc.order=TRUE,
           method="square",outline.color="white",type = "lower")

# env variables with a correlation of <|0.7| is a good threshold for determining if predictors correlate
dev.off()

png('figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/SSD_MGMs_FxnsOfInterest_CorrPlot_labeled.png',width = 1100, height = 1100)
ggcorrplot(cor.fxns1,method="square",lab = T,p.mat=p.fxns,sig.level = 0.05,insig="blank",
           hc.order=TRUE,outline.color="white",type = "lower")

# env variables with a correlation of <|0.7| is a good threshold for determining if predictors correlate
dev.off()

# save corrplots into 1 plot
fxn.corrplot1<-ggcorrplot(cor.fxns1,p.mat=p.fxns,hc.order=TRUE,
           method="square",outline.color="white",type = "lower")

fxn.corrplot2<-ggcorrplot(cor.fxns1,method="square",lab = T,p.mat=p.fxns,sig.level = 0.05,insig="blank",
                          hc.order=TRUE,outline.color="white",type = "lower")

fxn.corrplot.together<-ggarrange(fxn.corrplot1,fxn.corrplot2, ncol = 1, nrow = 2,common.legend=TRUE)

ggsave(fxn.corrplot.together,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/SSD_MGMs_FxnsOfInterest_CorrPlot_Combined.png", width=25, height=35, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).


# # subset out correlations with recA, house keeping gene
# cor.fxns.df<-as.data.frame(cor.fxns1)
# p.fxns.df<-as.data.frame(p.fxns)
#
# recA.corr<-data.frame(SampleID=rownames(cor.fxns.df),FxnCorr.recA=cor.fxns.df$`recA; recombination protein RecA`)
# recA.pval<-data.frame(SampleID=rownames(p.fxns.df),FxnCorrPval.recA=p.fxns.df$`recA; recombination protein RecA`)
#
# recA.corr.df<-merge(recA.corr,recA.pval,by="SampleID")
# recA.corr.df<-recA.corr.df[order(recA.corr.df$FxnCorrPval.recA),]

#### Does Coverage of GOI Differ from House Keeping Gene Coverage ####

HKOI<-mr.HK.ko.all[names(mr.HK.ko.all) %in% housekeep.fxn.name]
#HKOI$Category<-"House Keeping Gene"

# re-create dataframe of functions of interest
FxnsOI_HK<-cbind(SampleID=meta.all.scaled$SampleID,lpsOI,sporOI,TSOI,UVrepOI,QSOI,osmoOI,HKOI)
head(FxnsOI_HK)
rownames(FxnsOI_HK)<-FxnsOI_HK$SampleID

head(FxnsOI_HK)

kruskal.test(`uvrA; excinuclease ABC subunit A` ~ `recA; recombination protein RecA`, data = FxnsOI_HK)
wilcox.test(FxnsOI_HK$`uvrA; excinuclease ABC subunit A`, FxnsOI_HK$`recA; recombination protein RecA`, p.adjust.method = "bonf") # returns p values

kw.HK_<- vector('list', ncol(FxnsOI_HK[,-ncol(FxnsOI_HK)]) * 1) # create empty list where the corr output is stored
wcx.HK_<- vector('list', ncol(FxnsOI_HK[,-ncol(FxnsOI_HK)]) * 1) # create empty list where the corr output is stored
results_<- vector('list', ncol(FxnsOI_HK[,-ncol(FxnsOI_HK)]) * 1) # create an empty list where the corr summaries are stored
sig.results<-vector('list', ncol(FxnsOI_HK[,-ncol(FxnsOI_HK)]) * 1)
compnum <- 1 # counting our model numbers for indexes purposes in the loop

for (j in 1:ncol(FxnsOI_HK[,-c(1,ncol(FxnsOI_HK))])){ # for each column in indep.var.df
  # cor.test(x, y)
  kw.HK_[[compnum]] <- kruskal.test(FxnsOI_HK[,j] ~ `recA; recombination protein RecA`, data = FxnsOI_HK)
  #wcx.HK_[[compnum]] <- wilcox.test(FxnsOI_HK[,j], FxnsOI_HK$`recA; recombination protein RecA`, p.adjust.method = "bonf") # returns p values
  results_[[compnum]] <- kw.HK_[[compnum]] # save results of corr into list called results
  names(results_)[compnum]<-paste(names(FxnsOI_HK)[j],"~ recA") # rename list element to contain the name of the columns used in the model

  # save only significant corrs to another list called sig.results
  ## if p-value < 0.05, save to sig.results list
  ifelse(results_[[compnum]]$p.value < 0.05, sig.results[[compnum]]<-results_[[compnum]], "Not Sig")
  names(sig.results)[compnum]<-paste(names(FxnsOI_HK)[j],"~ recA")
  compnum <- compnum + 1 # add 1 to modelnumber so we keep track of # of models (for indexing purposes in list)

}

for (j in 2:ncol(FxnsOI_HK)){ # for each column in indep.var.df
  #print(FxnsOI_HK[,j])
  print(wilcox.test(FxnsOI_HK[,j], FxnsOI_HK[,ncol(FxnsOI_HK)], p.adjust.method = "bonf"))
  #print(kruskal.test(FxnsOI_HK[,j] ~ FxnsOI_HK[,ncol(FxnsOI_HK)], data = FxnsOI_HK))
  print(paste(names(FxnsOI_HK[j]),"~ recA"))
}

#### Create Multi-Panel Plot for Manuscript ####

# sig
# lexA.plotA<-ggplot(data=mr.UVrep.ko.all,aes(x=ave.v_N.S.wind, y=round(`lexA; repressor LexA [EC:3.4.21.88]`)))+
#   stat_smooth(aes(x=ave.v_N.S.wind, y=round(`lexA; repressor LexA [EC:3.4.21.88]`)), method=stats::glm, method.args=list(family="poisson")) +
#   geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
#                      labels=c("July 2020","August 2020","November 2020",
#                               "July 2021","August 2021","September 2021","December 2021"))+
#   scale_shape_manual(name="Site",values = c(0,1,16,15)) + theme_classic() +
#   stat_poly_eq(use_label(c("adj.R2", "P")), formula=y~x,small.p = TRUE) +
#   theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=18),axis.title.y = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=18),legend.text = element_text(size=15),plot.title = element_text(size=15)) +
#   labs(x="North to South Wind Vector Component (v)", y="lexA KO Coverage")

# not sig
uvrA.plotA<-ggplot(data=mr.UVrep.ko.all,aes(x=ave.v_N.S.wind, y=`uvrA; excinuclease ABC subunit A`))+
  stat_smooth(aes(x=ave.v_N.S.wind, y=`uvrA; excinuclease ABC subunit A`), method=stats::lm) +
  geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  stat_poly_eq(use_label(c("adj.R2", "P")), formula=y~x,small.p = TRUE,inherit.aes = TRUE) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=18),axis.title.y = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=18),legend.text = element_text(size=15),plot.title = element_text(size=15)) +
  labs(x="North to South Wind Vector Component (v)", y="uvrA KO Coverage")

# sig
uvrB.plotA<-ggplot(data=mr.UVrep.ko.all,aes(x=ave.v_N.S.wind, y=`uvrB; excinuclease ABC subunit B`))+
  stat_smooth(aes(x=ave.v_N.S.wind, y=`uvrB; excinuclease ABC subunit B`), method=stats::lm) +
  geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  stat_poly_eq(use_label(c("adj.R2", "P")), formula=y~x,small.p = TRUE,inherit.aes = TRUE) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=18),axis.title.y = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=18),legend.text = element_text(size=15),plot.title = element_text(size=15)) +
  labs(x="North to South Wind Vector Component (v)", y="uvrB KO Coverage")

# not sig
uvrC.plotA<-ggplot(data=mr.UVrep.ko.all,aes(x=ave.v_N.S.wind, y=`uvrC; excinuclease ABC subunit C`))+
  stat_smooth(aes(x=ave.v_N.S.wind, y=`uvrC; excinuclease ABC subunit C`), method=stats::glm, method.args=list(family="poisson")) +
  geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  stat_poly_eq(use_label(c("adj.R2", "P")), formula=y~x,small.p = TRUE,inherit.aes = TRUE) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=18),axis.title.y = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=18),legend.text = element_text(size=15),plot.title = element_text(size=15)) +
  labs(x="North to South Wind Vector Component (v)", y="uvrC KO Coverage")

# not sig
spmA.plotA<-ggplot(data=mr.spor.ko.all,aes(x=ave.v_N.S.wind, y=`spmA; spore maturation protein A`))+
  stat_smooth(aes(x=ave.v_N.S.wind, y=round(`spmA; spore maturation protein A`)), method=stats::glm,method.args = list(family = "poisson")) +
  geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  stat_poly_eq(use_label(c("adj.R2", "P")), formula=y~x,small.p = TRUE) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=18),axis.title.y = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=18),legend.text = element_text(size=15),plot.title = element_text(size=15)) +
  labs(x="North to South Wind Vector Component (v)", y="spmA KO Coverage")

# not sig
osmY.plotA<-ggplot(data=mr.osmo.ko.all,aes(x=ave.v_N.S.wind, y=`osmY; hyperosmotically inducible periplasmic protein`))+
  stat_smooth(aes(x=ave.v_N.S.wind, y=round(`osmY; hyperosmotically inducible periplasmic protein`)), method=stats::glm,method.args = list(family = "poisson")) +
  geom_jitter(aes(shape=Site,color=SampDate),size=6, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("seagreen3","orange","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  stat_poly_eq(use_label(c("adj.R2", "P")), formula=y~x,label.x = "left",
               label.y = "top",small.p = TRUE) +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=18),axis.title.y = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=18),legend.text = element_text(size=15),plot.title = element_text(size=15)) +
  labs(x="North to South Wind Vector Component (v)", y="osmY KO Coverage")


sigwind.comp<-ggarrange(cspA.plot2,lexA.plot2,ncol = 2, nrow = 1,legend="right",common.legend = TRUE)
ggsave(sigwind.comp,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/SSD_MGMs_SigFxns_WindCompV_Comparisons_Combined_updated.png", width=35, height=15, dpi=300,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

windcomp.v.compare.plots<-ggarrange(uvrA.plotA,uvrB.plotA,
                                    ncol = 2, nrow = 1,legend="right",common.legend = TRUE)
windcomp.v.compare.plots

ggsave(windcomp.v.compare.plots,filename = "figures/MGM_Figs/Contigs/MGM_Fxn_Corrs/SSD_MGMs_WindCompV_Comparisons_Combined_NotAllSig_updated.png", width=35, height=15, dpi=300,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

#### Save Progress ####

save.image("data/SSD_mgm_Fxn_GLMs_Ready.Rdata") # load Rdata to global env
