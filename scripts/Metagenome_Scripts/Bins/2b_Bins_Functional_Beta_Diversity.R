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
  #library(apeglm)
  library(plyr)
  library(dplyr)
  library(viridis)
  library(readxl)
  library(metagenomeSeq)
  #library(heatmaply)
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
  library(shades)
  #library(ALDEx2)
  library(rstatix)
  library(devtools)
  library(pairwiseAdonis)
  library(extrafont)
})

#### Load Data & See Info About Data ####
#load("data/SSD_MGM_Bins_analysis_data.Rdata") # load Rdata to global env
load("data/Metagenomes/SSD_MGM_Bins_Fxn_BetaDiv.Rdata")

head(all_bin_meta)
arsen.bin.fxns[1:4,]
head(mgm.bin_mr_melt)

# Before your start!:
## create column for SampDates that removes the period between the month and year (for plotting)
all_bin_meta$SampDate_Plot<-gsub("\\."," ",all_bin_meta$SampDate)
all_bin_meta$SampDate_Plot[(all_bin_meta$SampDate_Plot) == "August 2020"] <- "Aug 2020"
all_bin_meta$SampDate_Plot[(all_bin_meta$SampDate_Plot) == "November 2020"] <- "Nov 2020"
all_bin_meta$SampDate_Plot[(all_bin_meta$SampDate_Plot) == "August 2021"] <- "Aug 2021"
all_bin_meta$SampDate_Plot[(all_bin_meta$SampDate_Plot) == "September 2021"] <- "Sept 2021"
all_bin_meta$SampDate_Plot[(all_bin_meta$SampDate_Plot) == "December 2021"] <- "Dec 2021"
all_bin_meta$SampDate_Plot<-factor(all_bin_meta$SampDate_Plot,levels=c("July 2020", "Aug 2020", "Nov 2020", "July 2021",
                                                                             "Aug 2021", "Sept 2021", "Dec 2021"))
all_bin_meta$SampDate_Plot # sanity check

# ABOUT THE DATA:
# Before transformations (i.e., VST, CLR, etc) were done, the following was performed:
# contigs were co-assembled, so genes are found in same assembly (same set of contigs)
# non-normalized reads were mapped to genes in contigs (reads mapped by metagenome aka sample)
# featureCounts counted reads that mapped to genes on contigs
# Reads mapped to genes were divided by gene length for all genes across all samples because same KO can be assigned to multiple genes
# Gene coverage was then added together for each KO ID, since multiple genes were assigned the same KO ID
# Summed coverage per KO was then transformed and/or normalized via median-ratio normalization (mr), vst, and clr

## For pathway analyses -- after gene coverage was calculated and added together per KO ID, they were added together for each pathway
## summed coverages per KO ID, then per pathway were transformed by CLR

# For bins only: pseudocount of 1 was added before MR transformation so that log ratio of geometric means could be calculated
# NOTE about CLR transformation:
## uses a pseudocount of 1 to replace 0s, which is why not all 0s are treated equally
## need to look into robustCLR, which uses CLR transformation without 0s. Need more info on this methodology...

# #### Functional Beta Diversity - MRN data ####
# # MR = median-ratio normalization
# mgm.bin.mr[1:4,1:4] # sample IDs are rows, genes are columns
# mgm.bin_fxn.cov_table[1:4,1:4] # sanity check --> KOs with low coverage are still in this df
#
# # check rownames of MR & Mr transformed feature count data & metadata
# rownames(mgm.bin.mr) %in% rownames(all_bin_meta)
#
# ## PCOA with VST transformed data first
# # calculate our Euclidean distance matrix using VST data
# mgm.bin.euc_dist.mr <- dist(mgm.bin.mr, method = "euclidean")
#
# # creating our hierarcical clustering dendrogram
# mgm.bin.euc.mr_clust <- hclust(mgm.bin.euc_dist.mr, method="ward.D2")
#
# # let's make it a little nicer...
# mgm.bin.euc.mr_dend <- as.dendrogram(mgm.bin.euc.mr_clust, hang=0.2)
# mgm.bin.dend_cols <- as.character(all_bin_meta$SampDate_Color[order.dendrogram(mgm.bin.euc.mr_dend)])
# labels_colors(mgm.bin.euc.mr_dend) <- mgm.bin.dend_cols
#
# #png(filename="figures/BetaDiversity/SSD_16S_CLR_EucDist_Dendrogram1.png",width = 7, height = 7, units = "in",res = 800)
# plot(mgm.bin.euc.mr_dend, ylab="MR Euclidean Distance", cex = 0.1,horiz=TRUE) + title(main = "MAG Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 2)
# #legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c( "#ef781c","#03045e","#059c3f"),pch = 15, bty = "n")
# # Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
# dev.off()
#
# # let's use our Euclidean distance matrix from before
# mgm.bin.pcoa.mr <- pcoa(mgm.bin.euc_dist.mr) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
# ##save.image("data/ssd_mr.euc.dist1_3.7.23.Rdata")
#
# # The proportion of variances explained is in its element values$Relative_eig
# mgm.bin.pcoa.mr$values
# #     Eigenvalues Relative_eig Broken_stick Cumul_eig Cumul_br_stick
# # 1  4.985687e+04 0.5479262 0.0675112521 0.5479262     0.06751125
# # 2  1.317962e+04 0.1448439 0.0536223632 0.6927701     0.12113362
#
# # extract principal coordinates
# mgm.bin.pcoa.mr.vectors<-data.frame(mgm.bin.pcoa.mr$vectors)
# mgm.bin.pcoa.mr.vectors$Bin_ID<-rownames(mgm.bin.pcoa.mr$vectors)
#
# # merge pcoa coordinates w/ metadata
# mgm.bin.pcoa.mr.meta<-merge(mgm.bin.pcoa.mr.vectors, all_bin_meta, by.x="Bin_ID", by.y="Bin_ID")
# mgm.bin.pcoa.mr.meta$SampleMonth
# mgm.bin.pcoa.mr.meta$SampDate
#
# head(mgm.bin.pcoa.mr.meta)
#
# head(mgm.bin.pcoa.mr$values) # pull out Relative (Relative_eig) variation % to add to axes labels
#
# # create PCoA ggplot fig
# pcoa1<-ggplot(mgm.bin.pcoa.mr.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(SampDate),shape=Site), size=4)+theme_bw()+
#   labs(title="PCoA: Bacterial Functions in Salton Sea Dust",subtitle="Using Median-Ratio Normalized Feature Data",xlab="PC1 [47.13%]", ylab="PC2 [25.79%]",color="Sample Date")+theme_classic()+ theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+xlab("PC1 [47.13%]") + ylab("PC2 [25.79%]") + scale_shape_manual(values = c(7,10, 15,16))
#
# ggsave(pcoa1,filename = "figures/MGM_Figs/Bins/FxnDiv/PCoAs/MedianRatioNormalization/SSD_MGM_Bins_pcoa_MR_sampdate.png", width=12, height=10, dpi=600,create.dir = TRUE)
#
# pcoa1a<-ggplot(mgm.bin.pcoa.mr.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(SampDate),shape=Site,size=CollectionYear))+theme_bw()+
#   labs(title="PCoA: Bacterial Functions in Salton Sea Dust",subtitle="Using Median-Ratio Normalized Feature Data",xlab="PC1 [47.13%]", ylab="PC2 [25.79%]",color="Sample Date")+theme_classic()+ theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_size_manual(values = c("2020" = 7, "2021"=4),labels=c("2020","2021")) +
#   scale_color_manual(name ="Sample Date",values=unique(mgm.bin.pcoa.mr.meta$SampDate_Color[order(mgm.bin.pcoa.mr.meta$SampDate)]),labels=c("July 2020", "August 2020", "November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
#   xlab("PC1 [47.13%]") + ylab("PC2 [25.79%]") + scale_shape_manual(values = c(7,10, 15,16))
#
# ggsave(pcoa1a,filename = "figures/MGM_Figs/Bins/FxnDiv/PCoAs/MedianRatioNormalization/SSD_MGM_Bins_pcoa_MR_sampdate_v2.png", width=12, height=10, dpi=600,create.dir = TRUE)
#
# ## betadisper to look at homogeneity of group dispersions (aka variance) when considering multiple variables
# # multivariate analogue to Levene's test of homogeneity of variances
# # program finds spatial median or centroid of the group, & compare distances of group to centroid/spatial median via ANOVA
#
# #While PERMANOVA tests differences in group means (analogous to MANOVA),
# ## a related test called PERMDISP can be used to evaluate homogeneity of group dispersion
# #(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
# ## * need a distance matrix!
#
# # check rownames of MR & Mr transformed feature count data & metadata
# # rownames(mgm.bin.mr) %in% rownames(all_bin_meta)
# # mgm.bin.euc_dist.mr came from mgm.bin.mr (pre-Euclidean distance calculation)
#
# # first by Site
# mgm.bin.disper1<-betadisper(mgm.bin.euc_dist.mr, all_bin_meta$Site)
# mgm.bin.disper1
#
# permutest(mgm.bin.disper1, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
# #Pairwise comparisons:
# #  (Observed p-value below diagonal, permuted p-value above diagonal)
# #         PD     BDC      DP    WI
# # PD          0.27100 0.31900 0.421
# # BDC 0.27470         0.77000 0.767
# # DP  0.32008 0.74119         0.969
# # WI  0.40856 0.75608 0.95284
#
# anova(mgm.bin.disper1) # p = 0.5682 --> accept the Null H, spatial medians are NOT significantly difference across sample dates
#
# TukeyHSD(mgm.bin.disper1) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other
#
# #                             diff        lwr       upr     p adj
# # BDC-PD -174.453220 -540.4040 191.4975 0.5528869
# # DP-PD  -140.272644 -506.2234 225.6781 0.7094941
# # WI-PD  -133.797926 -499.7487 232.1528 0.7379750
# # DP-BDC   34.180576 -331.7702 400.1313 0.9935378
# # WI-BDC   40.655294 -325.2955 406.6060 0.9892566
# # WI-DP     6.474718 -359.4760 372.4255 0.9999548
#
# # Visualize dispersions
# png('figures/MGM_Figs/Bins/FxnDiv/PCoAs/MedianRatioNormalization/SSD_MGM_Bins_pcoa_MR_betadispersion_sampledate.png',width = 700, height = 600, res=100)
# plot(mgm.bin.disper1,main = "Centroids and Dispersion (Median-Ratio Data)", col=colorset6$Site_Color)
# dev.off()
#
# png('figures/MGM_Figs/Bins/FxnDiv/PCoAs/MedianRatioNormalization/SSD_MGM_Bins_boxplot_MR_centroid_distance_sampledate.png',width = 700, height = 600, res=100)
# boxplot(mgm.bin.disper1,xlab="Site", main = "Distance to Centroid by Category (Median-Ratio Data)", sub="Euclidean Distance of Median-Ratio Normalized Data", col=colorset6$Site_Color)
# dev.off()
#
# # next by sampdate
# mgm.bin.disper2<-betadisper(mgm.bin.euc_dist.mr, all_bin_meta$SampDate)
# mgm.bin.disper2
#
# permutest(mgm.bin.disper2, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
# #Pairwise comparisons:
# #  (Observed p-value below diagonal, permuted p-value above diagonal)
# #               July.2020 August.2020 November.2020 July.2021 August.2021 September.2021 December.2021
# # July.2020                    0.51500       0.39900   0.52900                    0.22500         0.867
# # August.2020      0.48699                   0.82000   0.95700                    0.41000         0.314
# # November.2020    0.38643     0.79863                 0.87800                    0.50900         0.254
# # July.2021        0.48033     0.93268       0.84588                              0.49200         0.343
# # August.2021
# # September.2021   0.24016     0.38975       0.45927   0.46416                                    0.175
# # December.2021    0.84776     0.32450       0.24531   0.31308                    0.18285
#
# anova(mgm.bin.disper2) # p = 0.3704 --> accept the Null H, spatial medians are NOT significantly difference across sample dates
#
# TukeyHSD(mgm.bin.disper2) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other
#
# #                             diff        lwr       upr     p adj
# # August.2020-July.2020          99.805636  -405.8284  605.4397 0.9934205
# # November.2020-July.2020       125.834303  -379.7998  631.4684 0.9783936
# # July.2021-July.2020           107.687701  -438.4595  653.8349 0.9934575
# # August.2021-July.2020        -204.827073 -1004.3047  594.6506 0.9750820
# # September.2021-July.2020      276.692803  -228.9413  782.3269 0.5529359
# # December.2021-July.2020       -30.673435  -536.3075  474.9606 0.9999923
# # November.2020-August.2020      26.028667  -479.6054  531.6627 0.9999971
# # July.2021-August.2020           7.882065  -538.2651  554.0293 1.0000000
# # August.2021-August.2020      -304.632709 -1104.1104  494.8449 0.8571585
# # September.2021-August.2020    176.887167  -328.7469  682.5212 0.8982181
# # December.2021-August.2020    -130.479071  -636.1131  375.1550 0.9741934
# # July.2021-November.2020       -18.146602  -564.2938  528.0006 0.9999998
# # August.2021-November.2020    -330.661376 -1130.1390  468.8163 0.8067577
# # September.2021-November.2020  150.858500  -354.7756  656.4926 0.9488910
# # December.2021-November.2020  -156.507738  -662.1418  349.1263 0.9397101
# # August.2021-July.2021        -312.514774 -1138.2117  513.1822 0.8607982
# # September.2021-July.2021      169.005102  -377.1421  715.1523 0.9397775
# # December.2021-July.2021      -138.361136  -684.5083  407.7861 0.9764135
# # September.2021-August.2021    481.519876  -317.9578 1280.9975 0.4470247
# # December.2021-August.2021     174.153638  -625.3240  973.6313 0.9889772
# # December.2021-September.2021 -307.366238  -813.0003  198.2678 0.4367483
#
# # Visualize dispersions
# png('figures/MGM_Figs/Bins/FxnDiv/PCoAs/MedianRatioNormalization/SSD_MGM_Bins_pcoa_MR_betadispersion_sampledate.png',width = 700, height = 600, res=100)
# plot(mgm.bin.disper2,main = "Centroids and Dispersion (Median-Ratio Data)", col=colorset8$SampDate_Color)
# dev.off()
#
# png('figures/MGM_Figs/Bins/FxnDiv/PCoAs/MedianRatioNormalization/SSD_MGM_Bins_boxplot_MR_centroid_distance_sampledate.png',width = 700, height = 600, res=100)
# boxplot(mgm.bin.disper2,xlab="Collection Date", main = "Distance to Centroid by Category (Median-Ratio Data)", sub="Euclidean Distance of Median-Ratio Normalized Data", col=colorset8$SampDate_Color)
# dev.off()
#
#
#

#### Pull Out LPS Metabolic Fxns from MR data ####
## heatmaps of traits of interest
mgm.bin.mr.lps<-merge(lps.bin.fxns, mgm.bin_mr_melt, by=c("KO_ID"),allow.cartesian = TRUE)
head(mgm.bin.mr.lps)

mgm.bin.mr[1:4,1:4]

# pull out LPS functions from MR normalized, summed coverages (summed gene coverage per KO)
lps.bins.ko<-mgm.bin.mr[,which(colnames(mgm.bin.mr) %in% lps.bin.fxns$KO_ID)] # merge MR data w/ S fxns found in contigs from KOFamScan
lps.bins.ko$Bin_ID<-rownames(lps.bins.ko)
lps.bins.ko.melt<-melt(lps.bins.ko, by="Bin_ID")
colnames(lps.bins.ko.melt)[which(names(lps.bins.ko.melt) == "variable")] <- "KO_ID"
colnames(lps.bins.ko.melt)[which(names(lps.bins.ko.melt) == "value")] <- "MR_SumCovPerKO"
head(lps.bins.ko.melt) #sanity check

mr.lps.bins.ko<-merge(lps.bins.ko.melt,lps.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(mr.lps.bins.ko)
colnames(mr.lps.bins.ko)[which(names(mr.lps.bins.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
mr.cov.sum.lps.bins.ko<-as.data.frame(dcast(mr.lps.bins.ko, Bin_ID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(mr.cov.sum.lps.bins.ko)<-mr.cov.sum.lps.bins.ko$Bin_ID
mr.cov.sum.lps.bins.ko[1:4,]

# sanity check
mr.cov.sum.lps.bins.ko$`lpxL, htrB; Kdo2-lipid IVA lauroyltransferase/acyltransferase [EC:2.3.1.241 2.3.1.-]`[1:4]
head(mr.lps.bins.ko)

#### LPS Heat Maps ####
# see max & mean of summed
max(mr.cov.sum.lps.bins.ko[,-1])
min(mr.cov.sum.lps.bins.ko[,-1])
max(mr.cov.sum.lps.bins.ko[,-1])/2

# first heat map of sulfur KOs
#heatmap(as.matrix(mr.cov.sum.lps.bins.ko[,-1]), scale = "none")

colSums(mr.cov.sum.lps.bins.ko[,-1])
#mr.cov.sum.lps.bins.ko2 <- mr.cov.sum.lps.bins.ko[,which(colSums(mr.cov.sum.lps.bins.ko[,-1])>10)]

#heatmap(as.matrix(mr.cov.sum.lps.bins.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap -- using MR data that includes NAs so they are blocked out on heatmaps
mr.lps.bins.ko[1:4,]
mr.lps.bins.all<-merge(mr.lps.bins.ko,all_bin_meta,by="Bin_ID")
head(mr.lps.bins.all)
mr.lps.bins.all$PlotID = factor(mr.lps.bins.all$Bin_ID, levels=unique(mr.lps.bins.all$Bin_ID[order(mr.lps.bins.all$Site,mr.lps.bins.all$SampDate)]), ordered=TRUE)

unique(mr.lps.bins.all$LPS_Structure)
# mr.lps.bins.all$LPS_Structure<-mr.lps.bins.all$LPS_Structure
# mr.lps.bins.all$LPS_Structure[(mr.lps.bins.all$LPS_Structure) == "KDO2-lipid A biosynthesis, Raetz LPS_Structure"] <- "Lipid A"
# mr.lps.bins.all$LPS_Structure[(mr.lps.bins.all$LPS_Structure) == "CMP-KDO biosynthesis"] <- "CMP-KDO"
# mr.lps.bins.all$LPS_Structure[(mr.lps.bins.all$LPS_Structure) == "KDO2-lipid A biosynthesis, Raetz LPS_Structure, non-LpxL-LpxM type"] <- "Lipid A, non-LpxL-LpxM"
# mr.lps.bins.all$LPS_Structure[(mr.lps.bins.all$LPS_Structure) == "ADP-L-glycero-D-manno-heptose biosynthesis"] <- "ALgDmh"
# mr.lps.bins.all$LPS_Structure[(mr.lps.bins.all$LPS_Structure) == "KDO2-lipid A modification LPS_Structure"] <- "Lipid A Mod"
#
# mr.lps.bins.all$LPS_Structure<-factor(mr.lps.bins.all$LPS_Structure,levels=c("KDO2-lipid A biosynthesis, Raetz LPS_Structure","CMP-KDO biosynthesis","KDO2-lipid A biosynthesis, Raetz LPS_Structure, non-LpxL-LpxM type","ADP-L-glycero-D-manno-heptose biosynthesis","KDO2-lipid A modification LPS_Structure"))
# mr.lps.bins.all$LPS_Structure<-factor(mr.lps.bins.all$LPS_Structure,levels=c("Lipid A","CMP-KDO","Lipid A, non-LpxL-LpxM","ALgDmh","Lipid A Mod"))

#mr.lps.bins.all$KO_Function.KEGG = factor(mr.lps.bins.all$KO_Function.KEGG, levels=unique(mr.lps.bins.all$KO_Function.KEGG[order(mr.lps.bins.all$LPS_Structure)]), ordered=TRUE)

head(mr.lps.bins.all)

# convert all 0s to NAs so they appear gray on the heatmap
mr.lps.bins.all$MR_SumCovPerKO_NA<-ifelse(mr.lps.bins.all$MR_SumCovPerKO==0,NA,mr.lps.bins.all$MR_SumCovPerKO)
head(mr.lps.bins.all)

# For heatmap color gradient
max(mr.lps.bins.all$MR_SumCovPerKO, na.rm=TRUE)
max(mr.lps.bins.all$MR_SumCovPerKO, na.rm=TRUE)/2
min(mr.lps.bins.all$MR_SumCovPerKO, na.rm=TRUE)

# merge LPS KO coverages with MAG taxonomic annotations
mr.lps.taxa.all<-merge(mr.lps.bins.all,mgm.bin.tax,by=c("Bin_ID","SampleID"))
head(mr.lps.taxa.all)

# Figures below
# by Bin_ID

lps.hm1a<-ggplot(mr.lps.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="LPS Biosynthesis in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(lps.hm1a,filename = "figures/MGM_Figs/Bins/FxnDiv/LPS_Biosynthesis/LPS_KOFxns_MGMs_Bins_BinID_by_Function_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

lps.hm1a2<-ggplot(mr.lps.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="LPS Biosynthesis in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(lps.hm1a2,filename = "figures/MGM_Figs/Bins/FxnDiv/LPS_Biosynthesis/LPS_KOFxns_MGMs_Bins_BinID_by_Function_Site_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)

lps.hm1a3<-ggplot(mr.lps.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="LPS Biosynthesis in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(LPS_Structure~Site,scales="free", space = "free")

ggsave(lps.hm1a3,filename = "figures/MGM_Figs/Bins/FxnDiv/LPS_Biosynthesis/LPS_KOFxns_MGMs_Bins_by_BinID_by_Function_LPS_Structure_Site_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

lps.hm1a3a<-ggplot(mr.lps.bins.all[mr.lps.bins.all$LPS_Structure=="O-antigen Repeat Unit",], aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="O-Antigen (LPS) Biosynthesis in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(LPS_Structure~Site,scales="free", space = "free")

ggsave(lps.hm1a3a,filename = "figures/MGM_Figs/Bins/FxnDiv/LPS_Biosynthesis/LPS_KOFxns_MGMs_Bins_by_BinID_by_Function_OAntigen_Only_Site_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

lps.hm1a3b<-ggplot(mr.lps.bins.all[mr.lps.bins.all$LPS_Structure=="Lipid A",], aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Lipid A (LPS) Biosynthesis in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(LPS_Structure~Site,scales="free", space = "free")

ggsave(lps.hm1a3b,filename = "figures/MGM_Figs/Bins/FxnDiv/LPS_Biosynthesis/LPS_KOFxns_MGMs_Bins_by_BinID_by_Function_LipidA_Only_Site_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

lps.hm1a3c<-ggplot(mr.lps.bins.all[mr.lps.bins.all$LPS_Structure=="Core Region",], aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="LPS Core Region Biosynthesis in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(LPS_Structure~Site,scales="free", space = "free")

ggsave(lps.hm1a3c,filename = "figures/MGM_Figs/Bins/FxnDiv/LPS_Biosynthesis/LPS_KOFxns_MGMs_Bins_by_BinID_by_Function_CoreRegion_Only_Site_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

# lps.hm1a4<-ggplot(mr.lps.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.15) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="LPS Biosynthesis in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
#   theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(LPS_Structure~.,scales="free_y", space = "free")
#
# ggsave(lps.hm1a4,filename = "figures/MGM_Figs/Bins/FxnDiv/LPS_Biosynthesis/LPS_KOFxns_MGMs_Bins_BinID_by_Function_LPS_Structure_heatmap2.png", width=10, height=15, dpi=600,create.dir = TRUE)

# lps.hm1a5<-ggplot(mr.lps.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="LPS Biosynthesis in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
#   theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(LPS_Structure~Site+SampDate, scales="free", space = "free")
#
# ggsave(lps.hm1a5,filename = "figures/MGM_Figs/Bins/FxnDiv/LPS_Biosynthesis/LPS_KOFxns_MGMs_Bins_by_BinID_by_Function_LPS_Structure_Site_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)

lps.hm1a6<-ggplot(mr.lps.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="LPS Biosynthesis in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(LPS_Structure~SampDate, scales="free", space = "free")

ggsave(lps.hm1a6,filename = "figures/MGM_Figs/Bins/FxnDiv/LPS_Biosynthesis/LPS_KOFxns_MGMs_Bins_by_BinID_by_Function_LPS_Structure_SampDate_heatmap.png", width=49.5, height=15, dpi=600,create.dir = TRUE)

lps.hm1a7<-ggplot(mr.lps.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="LPS Biosynthesis in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(LPS_Structure~SampDate+Site, scales="free", space = "free")

ggsave(lps.hm1a7,filename = "figures/MGM_Figs/Bins/FxnDiv/LPS_Biosynthesis/LPS_KOFxns_MGMs_Bins_by_BinID_by_Function_LPS_Structure_SampDate_Site_heatmap.png", width=49.5, height=15, dpi=600,create.dir = TRUE)

## by MAG Taxonomic ID, not Bin ID

lps.hm1b<-ggplot(mr.lps.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="LPS Biosynthesis in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(lps.hm1b,filename = "figures/MGM_Figs/Bins/FxnDiv/LPS_Biosynthesis/LPS_KOFxns_MGMs_Bins_Genus_by_Function_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

lps.hm1b2<-ggplot(mr.lps.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="LPS Biosynthesis in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(lps.hm1b2,filename = "figures/MGM_Figs/Bins/FxnDiv/LPS_Biosynthesis/LPS_KOFxns_MGMs_Bins_Genus_by_Function_Site_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)

lps.hm1b3<-ggplot(mr.lps.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="LPS Biosynthesis in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(LPS_Structure~Site,scales="free", space = "free")

ggsave(lps.hm1b3,filename = "figures/MGM_Figs/Bins/FxnDiv/LPS_Biosynthesis/LPS_KOFxns_MGMs_Bins_by_Genus_by_Function_LPS_Structure_Site_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

lps.hm1b3a<-ggplot(mr.lps.taxa.all[mr.lps.taxa.all$LPS_Structure=="O-antigen Repeat Unit",], aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="O-Antigen (LPS) Biosynthesis in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(LPS_Structure~Site,scales="free", space = "free")

ggsave(lps.hm1b3a,filename = "figures/MGM_Figs/Bins/FxnDiv/LPS_Biosynthesis/LPS_KOFxns_MGMs_Bins_by_BinID_by_Function_OAntigen_Only_Site_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

lps.hm1b3b<-ggplot(mr.lps.taxa.all[mr.lps.taxa.all$LPS_Structure=="Lipid A",], aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Lipid A (LPS) Biosynthesis in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(LPS_Structure~Site,scales="free", space = "free")

ggsave(lps.hm1b3b,filename = "figures/MGM_Figs/Bins/FxnDiv/LPS_Biosynthesis/LPS_KOFxns_MGMs_Bins_by_BinID_by_Function_LipidA_Only_Site_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

lps.hm1b3c<-ggplot(mr.lps.taxa.all[mr.lps.taxa.all$LPS_Structure=="Core Region",], aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="LPS Core Region Biosynthesis in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(LPS_Structure~Site,scales="free", space = "free")

ggsave(lps.hm1b3c,filename = "figures/MGM_Figs/Bins/FxnDiv/LPS_Biosynthesis/LPS_KOFxns_MGMs_Bins_by_BinID_by_Function_CoreRegion_Only_Site_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

lps.hm1b4<-ggplot(mr.lps.taxa.all[mr.lps.taxa.all$MR_SumCovPerKO>=3,], aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="LPS Biosynthesis in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(LPS_Structure~Site,scales="free_y", space = "free")

ggsave(lps.hm1b4,filename = "figures/MGM_Figs/Bins/FxnDiv/LPS_Biosynthesis/LPS_KOFxns_MGMs_Bins_Genus_by_Function_LPS_Structure_heatmap2.png", width=10, height=15, dpi=600,create.dir = TRUE)

# lps.hm1b5<-ggplot(mr.lps.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="LPS Biosynthesis in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
#   theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(LPS_Structure~SampDate+Site, scales="free", space = "free")
#
# ggsave(lps.hm1b5,filename = "figures/MGM_Figs/Bins/FxnDiv/LPS_Biosynthesis/LPS_KOFxns_MGMs_Bins_by_Genus_by_Function_LPS_Structure_SampDate_Site_heatmap.png", width=49.5, height=15, dpi=600,create.dir = TRUE)

lps.hm1b6<-ggplot(mr.lps.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="LPS Biosynthesis in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(LPS_Structure~SampDate, scales="free", space = "free")

ggsave(lps.hm1b6,filename = "figures/MGM_Figs/Bins/FxnDiv/LPS_Biosynthesis/LPS_KOFxns_MGMs_Bins_by_Genus_by_Function_LPS_Structure_SampDate_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)

lps.hm1b7<-ggplot(mr.lps.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="LPS Biosynthesis in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(LPS_Structure~SampDate+Site, scales="free", space = "free")

ggsave(lps.hm1b7,filename = "figures/MGM_Figs/Bins/FxnDiv/LPS_Biosynthesis/LPS_KOFxns_MGMs_Bins_by_Genus_by_Function_LPS_Structure_SampDate_Site_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)

#### Pull Out Quorum Sensing Metabolic Fxns from MR data ####
## heatmaps of traits of interest
mgm.bin.mr.quorsens<-merge(QuorSens.bin.fxns, mgm.bin_mr_melt, by=c("KO_ID"),allow.cartesian = TRUE)
head(mgm.bin.mr.quorsens)

mgm.bin.mr[1:4,1:4]

# pull out LPS functions from MR normalized, summed coverages (summed gene coverage per KO)
quorsens.bins.ko<-mgm.bin.mr[,which(colnames(mgm.bin.mr) %in% QuorSens.bin.fxns$KO_ID)] # merge MR data w/ S fxns found in contigs from KOFamScan
quorsens.bins.ko$Bin_ID<-rownames(quorsens.bins.ko)
quorsens.bins.ko.melt<-melt(quorsens.bins.ko, by="Bin_ID")
colnames(quorsens.bins.ko.melt)[which(names(quorsens.bins.ko.melt) == "variable")] <- "KO_ID"
colnames(quorsens.bins.ko.melt)[which(names(quorsens.bins.ko.melt) == "value")] <- "MR_SumCovPerKO"
head(quorsens.bins.ko.melt) #sanity check

mr.quorsens.bins.ko<-merge(quorsens.bins.ko.melt,QuorSens.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(mr.quorsens.bins.ko)
colnames(mr.quorsens.bins.ko)[which(names(mr.quorsens.bins.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
mr.cov.sum.quorsens.bins.ko<-as.data.frame(dcast(mr.quorsens.bins.ko, Bin_ID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(mr.cov.sum.quorsens.bins.ko)<-mr.cov.sum.quorsens.bins.ko$Bin_ID
mr.cov.sum.quorsens.bins.ko[1:4,]

# sanity check
mr.cov.sum.quorsens.bins.ko$`spoVS; stage V QuorumSensing protein S`[1:4]
head(mr.quorsens.bins.ko)

#### Quorum Sensing Heat Maps ####
# see max & mean of summed
max(mr.cov.sum.quorsens.bins.ko[,-1])
min(mr.cov.sum.quorsens.bins.ko[,-1])
max(mr.cov.sum.quorsens.bins.ko[,-1])/2

# first heat map of  KOs
#heatmap(as.matrix(mr.cov.sum.quorsens.bins.ko[,-1]), scale = "none")

colSums(mr.cov.sum.quorsens.bins.ko[,-1])
#mr.cov.sum.quorsens.bins.ko2 <- mr.cov.sum.quorsens.bins.ko[,which(colSums(mr.cov.sum.quorsens.bins.ko[,-1])>10)]

#heatmap(as.matrix(mr.cov.sum.quorsens.bins.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap -- using MR data that includes NAs so they are blocked out on heatmaps
mr.quorsens.bins.ko[1:4,]
mr.quorsens.bins.all<-merge(mr.quorsens.bins.ko,all_bin_meta,by="Bin_ID")
head(mr.quorsens.bins.all)
#mr.quorsens.bins.all$PlotID<-mr.quorsens.bins.all$Bin_ID
mr.quorsens.bins.all$PlotID = factor(mr.quorsens.bins.all$Bin_ID, levels=unique(mr.quorsens.bins.all$Bin_ID[order(mr.quorsens.bins.all$Site,mr.quorsens.bins.all$SampDate)]), ordered=TRUE)

unique(mr.quorsens.bins.all$Pathway)

head(mr.quorsens.bins.all)

# convert all 0s to NAs so they appear gray on the heatmap
mr.quorsens.bins.all$MR_SumCovPerKO_NA<-ifelse(mr.quorsens.bins.all$MR_SumCovPerKO==0,NA,mr.quorsens.bins.all$MR_SumCovPerKO)
head(mr.quorsens.bins.all)

# For heatmap color gradient
max(mr.quorsens.bins.all$MR_SumCovPerKO, na.rm=TRUE)
max(mr.quorsens.bins.all$MR_SumCovPerKO, na.rm=TRUE)/2
min(mr.quorsens.bins.all$MR_SumCovPerKO, na.rm=TRUE)

# merge QuorumSensing KO coverages with MAG taxonomic annotations
mr.quorsens.taxa.all<-merge(mr.quorsens.bins.all,mgm.bin.tax,by=c("Bin_ID","SampleID"))
head(mr.quorsens.taxa.all)

# Figures below
# by Bin_ID

quorsens.hm1a<-ggplot(mr.quorsens.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Quorum Sensing in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(quorsens.hm1a,filename = "figures/MGM_Figs/Bins/FxnDiv/QuorumSensing/QuorumSensing_KOFxns_MGMs_Bins_BinID_by_Function_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

quorsens.hm1a2<-ggplot(mr.quorsens.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Quorum Sensing in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(quorsens.hm1a2,filename = "figures/MGM_Figs/Bins/FxnDiv/QuorumSensing/QuorumSensing_KOFxns_MGMs_Bins_BinID_by_Function_Site_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)

quorsens.hm1a3<-ggplot(mr.quorsens.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Quorum Sensing in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate_Plot,scales="free_x", space = "free")

ggsave(quorsens.hm1a3,filename = "figures/MGM_Figs/Bins/FxnDiv/QuorumSensing/QuorumSensing_KOFxns_MGMs_Bins_BinID_by_Function_SampDate_heatmap.png", width=49.5, height=15, dpi=600,create.dir = TRUE)

quorsens.hm1a4<-ggplot(mr.quorsens.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Quorum Sensing in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate_Plot+Site,scales="free_x", space = "free")

ggsave(quorsens.hm1a4,filename = "figures/MGM_Figs/Bins/FxnDiv/QuorumSensing/QuorumSensing_KOFxns_MGMs_Bins_BinID_by_Function_SampDate_Site_heatmap.png", width=49.5, height=15, dpi=600,create.dir = TRUE)

## by MAG Taxonomic ID, not Bin ID

quorsens.hm1b<-ggplot(mr.quorsens.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Quorum Sensing in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(quorsens.hm1b,filename = "figures/MGM_Figs/Bins/FxnDiv/QuorumSensing/QuorumSensing_KOFxns_MGMs_Bins_Genus_by_Function_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

quorsens.hm1b2<-ggplot(mr.quorsens.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Quorum Sensing in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(quorsens.hm1b2,filename = "figures/MGM_Figs/Bins/FxnDiv/QuorumSensing/QuorumSensing_KOFxns_MGMs_Bins_Genus_by_Function_Site_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)

quorsens.hm1b3<-ggplot(mr.quorsens.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Quorum Sensing in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate_Plot,scales="free_x", space = "free")

ggsave(quorsens.hm1b3,filename = "figures/MGM_Figs/Bins/FxnDiv/QuorumSensing/QuorumSensing_KOFxns_MGMs_Bins_BinID_by_Function_SampDate_heatmap.png", width=49.5, height=15, dpi=600,create.dir = TRUE)

quorsens.hm1b4<-ggplot(mr.quorsens.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Quorum Sensing in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate_Plot+Site,scales="free_x", space = "free")

ggsave(quorsens.hm1b4,filename = "figures/MGM_Figs/Bins/FxnDiv/QuorumSensing/QuorumSensing_KOFxns_MGMs_Bins_BinID_by_Function_SampDate_Site_heatmap.png", width=49.5, height=15, dpi=600,create.dir = TRUE)

#### Pull Out Sporulation Metabolic Fxns from MR data ####
## heatmaps of traits of interest
mgm.bin.mr.spor<-merge(spor.bin.fxns, mgm.bin_mr_melt, by=c("KO_ID"),allow.cartesian = TRUE)
head(mgm.bin.mr.spor)

mgm.bin.mr[1:4,1:4]

# pull out LPS functions from MR normalized, summed coverages (summed gene coverage per KO)
spor.bins.ko<-mgm.bin.mr[,which(colnames(mgm.bin.mr) %in% spor.bin.fxns$KO_ID)] # merge MR data w/ S fxns found in contigs from KOFamScan
spor.bins.ko$Bin_ID<-rownames(spor.bins.ko)
spor.bins.ko.melt<-melt(spor.bins.ko, by="Bin_ID")
colnames(spor.bins.ko.melt)[which(names(spor.bins.ko.melt) == "variable")] <- "KO_ID"
colnames(spor.bins.ko.melt)[which(names(spor.bins.ko.melt) == "value")] <- "MR_SumCovPerKO"
head(spor.bins.ko.melt) #sanity check

mr.spor.bins.ko<-merge(spor.bins.ko.melt,spor.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(mr.spor.bins.ko)
colnames(mr.spor.bins.ko)[which(names(mr.spor.bins.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
mr.cov.sum.spor.bins.ko<-as.data.frame(dcast(mr.spor.bins.ko, Bin_ID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(mr.cov.sum.spor.bins.ko)<-mr.cov.sum.spor.bins.ko$Bin_ID
mr.cov.sum.spor.bins.ko[1:4,]

# sanity check
mr.cov.sum.spor.bins.ko$`spoVS; stage V sporulation protein S`[1:4]
head(mr.spor.bins.ko)

#### Sporulation Heat Maps ####
# see max & mean of summed
max(mr.cov.sum.spor.bins.ko[,-1])
min(mr.cov.sum.spor.bins.ko[,-1])
max(mr.cov.sum.spor.bins.ko[,-1])/2

# first heat map of sulfur KOs
#heatmap(as.matrix(mr.cov.sum.spor.bins.ko[,-1]), scale = "none")

colSums(mr.cov.sum.spor.bins.ko[,-1])
#mr.cov.sum.spor.bins.ko2 <- mr.cov.sum.spor.bins.ko[,which(colSums(mr.cov.sum.spor.bins.ko[,-1])>10)]

#heatmap(as.matrix(mr.cov.sum.spor.bins.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap -- using MR data that includes NAs so they are blocked out on heatmaps
mr.spor.bins.ko[1:4,]
mr.spor.bins.all<-merge(mr.spor.bins.ko,all_bin_meta,by="Bin_ID")
head(mr.spor.bins.all)
#mr.spor.bins.all$PlotID<-mr.spor.bins.all$Bin_ID
mr.spor.bins.all$PlotID = factor(mr.spor.bins.all$Bin_ID, levels=unique(mr.spor.bins.all$Bin_ID[order(mr.spor.bins.all$Site,mr.spor.bins.all$SampDate)]), ordered=TRUE)

unique(mr.spor.bins.all$Pathway)

head(mr.spor.bins.all)

# convert all 0s to NAs so they appear gray on the heatmap
mr.spor.bins.all$MR_SumCovPerKO_NA<-ifelse(mr.spor.bins.all$MR_SumCovPerKO==0,NA,mr.spor.bins.all$MR_SumCovPerKO)
head(mr.spor.bins.all)

# For heatmap color gradient
max(mr.spor.bins.all$MR_SumCovPerKO, na.rm=TRUE)
max(mr.spor.bins.all$MR_SumCovPerKO, na.rm=TRUE)/2
min(mr.spor.bins.all$MR_SumCovPerKO, na.rm=TRUE)

# merge Sporulation KO coverages with MAG taxonomic annotations
mr.spor.taxa.all<-merge(mr.spor.bins.all,mgm.bin.tax,by=c("Bin_ID","SampleID"))
head(mr.spor.taxa.all)

# Figures below
# by Bin_ID

spor.hm1a<-ggplot(mr.spor.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sporulation Functions in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(spor.hm1a,filename = "figures/MGM_Figs/Bins/FxnDiv/Sporulation/Sporulation_KOFxns_MGMs_Bins_BinID_by_Function_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

spor.hm1a2<-ggplot(mr.spor.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sporulation Functions in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(spor.hm1a2,filename = "figures/MGM_Figs/Bins/FxnDiv/Sporulation/Sporulation_KOFxns_MGMs_Bins_BinID_by_Function_Site_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)

spor.hm1a3<-ggplot(mr.spor.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sporulation Functions in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate,scales="free_x", space = "free")

ggsave(spor.hm1a3,filename = "figures/MGM_Figs/Bins/FxnDiv/Sporulation/Sporulation_KOFxns_MGMs_Bins_BinID_by_Function_SampDate_heatmap.png", width=49.5, height=15, dpi=600,create.dir = TRUE)

spor.hm1a4<-ggplot(mr.spor.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sporulation Functions in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate+Site,scales="free_x", space = "free")

ggsave(spor.hm1a4,filename = "figures/MGM_Figs/Bins/FxnDiv/Sporulation/Sporulation_KOFxns_MGMs_Bins_BinID_by_Function_SampDate_Site_heatmap.png", width=49.5, height=15, dpi=600,create.dir = TRUE)

## by MAG Taxonomic ID, not Bin ID

spor.hm1b<-ggplot(mr.spor.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sporulation Functions in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(spor.hm1b,filename = "figures/MGM_Figs/Bins/FxnDiv/Sporulation/Sporulation_KOFxns_MGMs_Bins_Genus_by_Function_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

spor.hm1b2<-ggplot(mr.spor.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sporulation Functions in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(spor.hm1b2,filename = "figures/MGM_Figs/Bins/FxnDiv/Sporulation/Sporulation_KOFxns_MGMs_Bins_Genus_by_Function_Site_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)

spor.hm1b3<-ggplot(mr.spor.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sporulation Functions in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate,scales="free_x", space = "free")

ggsave(spor.hm1b3,filename = "figures/MGM_Figs/Bins/FxnDiv/Sporulation/Sporulation_KOFxns_MGMs_Bins_Genus_by_Function_SampDate_heatmap.png", width=49.5, height=15, dpi=600,create.dir = TRUE)

spor.hm1b4<-ggplot(mr.spor.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sporulation Functions in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate+Site,scales="free_x", space = "free")

ggsave(spor.hm1b4,filename = "figures/MGM_Figs/Bins/FxnDiv/Sporulation/Sporulation_KOFxns_MGMs_Bins_Genus_by_Function_SampDate_Site_heatmap.png", width=49.5, height=15, dpi=600,create.dir = TRUE)

#### Pull Out Temp Shock Metabolic Fxns from MRN data ####
## heatmaps of traits of interest

mgm.bin.mr[1:4,1:4]

# pull out heat shock functions from MR normalized, summed coverages (summed gene coverage per KO)
tempshock.bins.ko<-mgm.bin.mr[,which(colnames(mgm.bin.mr) %in% tempshock.bin.fxns$KO_ID)] # merge MRN data w/ S fxns found in contigs from KOFamScan
tempshock.bins.ko$Bin_ID<-rownames(tempshock.bins.ko)
tempshock.bins.ko.melt<-melt(tempshock.bins.ko, by="Bin_ID")
colnames(tempshock.bins.ko.melt)[which(names(tempshock.bins.ko.melt) == "variable")] <- "KO_ID"
colnames(tempshock.bins.ko.melt)[which(names(tempshock.bins.ko.melt) == "value")] <- "MR_SumCovPerKO"
head(tempshock.bins.ko.melt) #sanity check

mr.tempshock.bins.ko<-merge(tempshock.bins.ko.melt,tempshock.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(mr.tempshock.bins.ko)
colnames(mr.tempshock.bins.ko)[which(names(mr.tempshock.bins.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
mr.cov.sum.tempshock.bins.ko<-as.data.frame(dcast(mr.tempshock.bins.ko, Bin_ID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(mr.cov.sum.tempshock.bins.ko)<-mr.cov.sum.tempshock.bins.ko$Bin_ID
mr.cov.sum.tempshock.bins.ko[1:4,]

# sanity check
mr.cov.sum.tempshock.bins.ko$`hslR; ribosome-associated heat shock protein Hsp15`[1:4]
head(mr.tempshock.bins.ko)

#### Temp Shock KO Heat Maps ####
# see max & mean of summed
max(mr.cov.sum.tempshock.bins.ko[,-1])
min(mr.cov.sum.tempshock.bins.ko[,-1])
max(mr.cov.sum.tempshock.bins.ko[,-1])/2

# first heat map of tempshock KOs
#heatmap(as.matrix(mr.cov.sum.tempshock.bins.ko[,-1]), scale = "none")

colSums(mr.cov.sum.tempshock.bins.ko[,-1])
#mr.cov.sum.tempshock.bins.ko2 <- mr.cov.sum.tempshock.bins.ko[,which(colSums(mr.cov.sum.tempshock.bins.ko[,-1])>10)]

##heatmap(as.matrix(mr.cov.sum.tempshock.bins.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap -- using MR data that includes NAs so they are blocked out on heatmaps
mr.tempshock.bins.ko[1:4,]
mr.tempshock.bins.all<-merge(mr.tempshock.bins.ko,all_bin_meta,by="Bin_ID")
head(mr.tempshock.bins.all)
#mr.tempshock.bins.all$PlotID<-mr.tempshock.bins.all$Bin_ID
mr.tempshock.bins.all$PlotID = factor(mr.tempshock.bins.all$Bin_ID, levels=unique(mr.tempshock.bins.all$Bin_ID[order(mr.tempshock.bins.all$Site,mr.tempshock.bins.all$SampDate)]), ordered=TRUE)

#mr.tempshock.bins.all$PathShort<-mr.tempshock.bins.all$Pathway
# mr.tempshock.bins.all$PathShort[(mr.tempshock.bins.all$PathShort) == "Dissimilatory Sulfate Redox"] <- "D.SO4 RedOx"
# mr.tempshock.bins.all$PathShort[(mr.tempshock.bins.all$PathShort) == "Assimilatory Sulfate Reduction"] <- "A.SO4 Red"
# mr.tempshock.bins.all$PathShort[(mr.tempshock.bins.all$PathShort) == "Multiple Pathways"] <- "Multi Paths"
# mr.tempshock.bins.all$PathShort[(mr.tempshock.bins.all$PathShort) == "S Disproportionation"] <- "S Disprop."

#mr.tempshock.bins.all$Pathway<-factor(mr.tempshock.bins.all$Pathway,levels=c("Assimilatory Sulfate Reduction","Dissimilatory Sulfate Redox","Multiple Pathways","SOX","S Disproportionation"))
#mr.tempshock.bins.all$PathShort<-factor(mr.tempshock.bins.all$PathShort,levels=c("A.SO4 Red","D.SO4 RedOx","Multi Paths","SOX","S Disprop."))

head(mr.tempshock.bins.all)

# convert all 0s to NAs so they appear gray on the heatmap
mr.tempshock.bins.all$MR_SumCovPerKO_NA<-ifelse(mr.tempshock.bins.all$MR_SumCovPerKO==0,NA,mr.tempshock.bins.all$MR_SumCovPerKO)
head(mr.tempshock.bins.all)

# For heatmap color gradient
max(mr.tempshock.bins.all$MR_SumCovPerKO, na.rm=TRUE)
max(mr.tempshock.bins.all$MR_SumCovPerKO, na.rm=TRUE)/2
min(mr.tempshock.bins.all$MR_SumCovPerKO, na.rm=TRUE)

# merge temp shock KO coverages with MAG taxonomic annotations
mr.tempshock.taxa.all<-merge(mr.tempshock.bins.all,mgm.bin.tax,by=c("Bin_ID","SampleID"))
head(mr.tempshock.taxa.all)

# Figures below
# by Bin_ID

tempshock.hm1a<-ggplot(mr.tempshock.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Temperature Shock Proteins in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(tempshock.hm1a,filename = "figures/MGM_Figs/Bins/FxnDiv/TempShock/TempShock_KOFxns_MGMs_Bins_BinID_by_Function_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

tempshock.hm1a2<-ggplot(mr.tempshock.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Temperature Shock Proteins in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(tempshock.hm1a2,filename = "figures/MGM_Figs/Bins/FxnDiv/TempShock/TempShock_KOFxns_MGMs_Bins_BinID_by_Function_Site_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

tempshock.hm1a3<-ggplot(mr.tempshock.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Temperature Shock Proteins in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate,scales="free_x", space = "free")

ggsave(tempshock.hm1a3,filename = "figures/MGM_Figs/Bins/FxnDiv/TempShock/TempShock_KOFxns_MGMs_Bins_by_BinID_by_Function_SampDate_heatmap.png", width=49.5, height=13, dpi=600,create.dir = TRUE)

tempshock.hm1a4<-ggplot(mr.tempshock.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Temperature Shock Proteins in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate+Site,scales="free_x", space = "free")

ggsave(tempshock.hm1a4,filename = "figures/MGM_Figs/Bins/FxnDiv/TempShock/TempShock_KOFxns_MGMs_Bins_by_BinID_by_Function_SampDate_Site_heatmap.png", width=49.5, height=13, dpi=600,create.dir = TRUE)

#  by MAG Taxonomic ID, not Bin ID

tempshock.hm1b<-ggplot(mr.tempshock.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Temperature Shock Proteins in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(tempshock.hm1b,filename = "figures/MGM_Figs/Bins/FxnDiv/TempShock/TempShock_KOFxns_MGMs_Bins_Genus_by_Function_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

tempshock.hm1b2<-ggplot(mr.tempshock.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Temperature Shock Proteins in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(tempshock.hm1b2,filename = "figures/MGM_Figs/Bins/FxnDiv/TempShock/TempShock_KOFxns_MGMs_Bins_Genus_by_Function_Site_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

tempshock.hm1b3<-ggplot(mr.tempshock.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Temperature Shock Proteins in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate,scales="free_x", space = "free")

ggsave(tempshock.hm1b3,filename = "figures/MGM_Figs/Bins/FxnDiv/TempShock/TempShock_KOFxns_MGMs_Bins_by_Genus_by_Function_SampDate_heatmap.png", width=49.5, height=13, dpi=600,create.dir = TRUE)

tempshock.hm1b4<-ggplot(mr.tempshock.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Temperature Shock Proteins in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate+Site,scales="free_x", space = "free")

ggsave(tempshock.hm1b4,filename = "figures/MGM_Figs/Bins/FxnDiv/TempShock/TempShock_KOFxns_MGMs_Bins_by_Genus_by_Function_SampDate_Site_heatmap.png", width=49.5, height=13, dpi=600,create.dir = TRUE)

#### Pull Out Aerotaxis Fxns from MRN data ####
## heatmaps of traits of interest

mgm.bin.mr[1:4,1:4]

# check if aerotaxis functions found in MAGs
mgm.bin.mr[,which(colnames(mgm.bin.mr) %in% "K03776")]

# pull out aerotaxis functions from MR normalized, summed coverages (summed gene coverage per KO)
# aerotaxis.bins.ko<-data.frame(K03776=mgm.bin.mr[,which(colnames(mgm.bin.mr) %in% "K03776")],Bin_ID=rownames(mgm.bin.mr)) # merge MRN data w/ S fxns found in contigs from KOFamScan
# aerotaxis.bins.ko.melt<-melt(aerotaxis.bins.ko, by="Bin_ID")
# colnames(aerotaxis.bins.ko.melt)[which(names(aerotaxis.bins.ko.melt) == "variable")] <- "KO_ID"
# colnames(aerotaxis.bins.ko.melt)[which(names(aerotaxis.bins.ko.melt) == "value")] <- "MR_SumCovPerKO"
# head(aerotaxis.bins.ko.melt) #sanity check
#
# #mr.aerotaxis.bins.ko<-merge(aerotaxis.bins.ko.melt,heatshock.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
# #colnames(mr.aerotaxis.bins.ko)[which(names(mr.aerotaxis.bins.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
# aerotaxis.bins.ko.melt$KO_Function.KEGG<-"aer; aerotaxis receptor"
# mr.cov.sum.aerotaxis.bins.ko<-as.data.frame(dcast(aerotaxis.bins.ko.melt, Bin_ID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###just dcast, nothing is being added here!
# rownames(mr.cov.sum.aerotaxis.bins.ko)<-mr.cov.sum.aerotaxis.bins.ko$Bin_ID
# mr.cov.sum.aerotaxis.bins.ko[1:4,]
#
# # sanity check
# mr.cov.sum.aerotaxis.bins.ko$`aer; aerotaxis receptor`[1:4]
# head(mr.cov.sum.aerotaxis.bins.ko)
#
# #### Aerotaxis KO Heat Maps ####
# # see max & mean of summed
# max(mr.cov.sum.aerotaxis.bins.ko[,-1])
# min(mr.cov.sum.aerotaxis.bins.ko[,-1])
# max(mr.cov.sum.aerotaxis.bins.ko[,-1])/2
#
# # first heat map of sulfur KOs
# # #heatmap(as.matrix(mr.cov.sum.aerotaxis.bins.ko[,-1]), scale = "none")
# #
# # colSums(mr.cov.sum.aerotaxis.bins.ko[,-1])
# # #mr.cov.sum.aerotaxis.bins.ko2 <- mr.cov.sum.aerotaxis.bins.ko[,which(colSums(mr.cov.sum.aerotaxis.bins.ko[,-1])>10)]
#
# #heatmap(as.matrix(mr.cov.sum.aerotaxis.bins.ko[,-1]), scale = "none")
#
# # prep for ggplot2 heatmap -- using MR data that includes NAs so they are blocked out on heatmaps
# aerotaxis.bins.ko.melt[1:4,]
# mr.aerotaxis.bins.all<-merge(aerotaxis.bins.ko.melt,all_bin_meta,by="Bin_ID")
# head(mr.aerotaxis.bins.all)
# mr.aerotaxis.bins.all$PlotID<-mr.aerotaxis.bins.all$Bin_ID
# mr.aerotaxis.bins.all$PlotID = factor(mr.aerotaxis.bins.all$PlotID, levels=unique(mr.aerotaxis.bins.all$PlotID[order(mr.aerotaxis.bins.all$Site,mr.aerotaxis.bins.all$SampDate)]), ordered=TRUE)
#
# #mr.aerotaxis.bins.all$PathShort<-mr.aerotaxis.bins.all$Pathway
# # mr.aerotaxis.bins.all$PathShort[(mr.aerotaxis.bins.all$PathShort) == "Dissimilatory Sulfate Redox"] <- "D.SO4 RedOx"
# # mr.aerotaxis.bins.all$PathShort[(mr.aerotaxis.bins.all$PathShort) == "Assimilatory Sulfate Reduction"] <- "A.SO4 Red"
# # mr.aerotaxis.bins.all$PathShort[(mr.aerotaxis.bins.all$PathShort) == "Multiple Pathways"] <- "Multi Paths"
# # mr.aerotaxis.bins.all$PathShort[(mr.aerotaxis.bins.all$PathShort) == "S Disproportionation"] <- "S Disprop."
#
# #mr.aerotaxis.bins.all$Pathway<-factor(mr.aerotaxis.bins.all$Pathway,levels=c("Assimilatory Sulfate Reduction","Dissimilatory Sulfate Redox","Multiple Pathways","SOX","S Disproportionation"))
# #mr.aerotaxis.bins.all$PathShort<-factor(mr.aerotaxis.bins.all$PathShort,levels=c("A.SO4 Red","D.SO4 RedOx","Multi Paths","SOX","S Disprop."))
#
# head(mr.aerotaxis.bins.all)
#
# # convert all 0s to NAs so they appear gray on the heatmap
# mr.aerotaxis.bins.all$MR_SumCovPerKO_NA<-ifelse(mr.aerotaxis.bins.all$MR_SumCovPerKO==0,NA,mr.aerotaxis.bins.all$MR_SumCovPerKO)
#
# # For heatmap color gradient
# max(aerotaxis.bins.ko.melt$MR_SumCovPerKO, na.rm=TRUE)
# max(aerotaxis.bins.ko.melt$MR_SumCovPerKO, na.rm=TRUE)/2
# min(aerotaxis.bins.ko.melt$MR_SumCovPerKO, na.rm=TRUE)
#
# # merge aerotaxis KO coverages with MAG taxonomic annotations
# mr.aerotaxis.taxa.all<-merge(mr.aerotaxis.bins.all,mgm.bin.tax,by=c("Bin_ID","SampleID"))
# head(mr.aerotaxis.taxa.all)
#
# # Figures below
# # by Bin_ID
#
# aerotaxis.hm1a<-ggplot(mr.aerotaxis.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Aerotaxis Receptor in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
#   theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))
#
# ggsave(aerotaxis.hm1a,filename = "figures/MGM_Figs/Bins/FxnDiv/Aerotaxis/Aerotaxis_KOFxns_MGMs_Bins_BinID_by_Function_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)
#
# aerotaxis.hm1a2<-ggplot(mr.aerotaxis.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.15) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Aerotaxis Receptors in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
#   theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")
#
# ggsave(aerotaxis.hm1a2,filename = "figures/MGM_Figs/Bins/FxnDiv/Aerotaxis/Aerotaxis_KOFxns_MGMs_Bins_BinID_by_Function_Site_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)
#
# ## by MAG Taxonomic ID, not Bin ID
#
# aerotaxis.hm1b<-ggplot(mr.aerotaxis.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Aerotaxis Functions in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
#   theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))
#
# ggsave(aerotaxis.hm1b,filename = "figures/MGM_Figs/Bins/FxnDiv/Aerotaxis/Aerotaxis_KOFxns_MGMs_Bins_Genus_by_Function_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)
#
# aerotaxis.hm1b2<-ggplot(mr.aerotaxis.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.15) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Aerotaxis Functions in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
#   theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")
#
# ggsave(aerotaxis.hm1b2,filename = "figures/MGM_Figs/Bins/FxnDiv/Aerotaxis/Aerotaxis_KOFxns_MGMs_Bins_Genus_by_Function_Site_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)

#### Pull Out UV-Damage DNA Repair Fxns from MRN Data ####
## heatmaps of traits of interest

mgm.bin.mr[1:4,1:4]

# pull out heat shock functions from MR normalized, summed coverages (summed gene coverage per KO)
UVrep.bins.ko<-mgm.bin.mr[,which(colnames(mgm.bin.mr) %in% UVrep.bin.fxns$KO_ID)] # merge MRN data w/ S fxns found in contigs from KOFamScan
UVrep.bins.ko$Bin_ID<-rownames(UVrep.bins.ko)
UVrep.bins.ko.melt<-melt(UVrep.bins.ko, by="Bin_ID")
colnames(UVrep.bins.ko.melt)[which(names(UVrep.bins.ko.melt) == "variable")] <- "KO_ID"
colnames(UVrep.bins.ko.melt)[which(names(UVrep.bins.ko.melt) == "value")] <- "MR_SumCovPerKO"
head(UVrep.bins.ko.melt) #sanity check

mr.UVrep.bins.ko<-merge(UVrep.bins.ko.melt,UV.DNA.rep.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(mr.UVrep.bins.ko)
colnames(mr.UVrep.bins.ko)[which(names(mr.UVrep.bins.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
mr.cov.sum.UVrep.bins.ko<-as.data.frame(dcast(mr.UVrep.bins.ko, Bin_ID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(mr.cov.sum.UVrep.bins.ko)<-mr.cov.sum.UVrep.bins.ko$Bin_ID
mr.cov.sum.UVrep.bins.ko[1:4,]

# sanity check
mr.cov.sum.UVrep.bins.ko$`uvrC; excinuclease ABC subunit C`[1:4]
head(mr.UVrep.bins.ko)

#### UV-Damage DNA Repair KO Heat Maps ####
# see max & mean of summed
max(mr.cov.sum.UVrep.bins.ko[,-1])
min(mr.cov.sum.UVrep.bins.ko[,-1])
max(mr.cov.sum.UVrep.bins.ko[,-1])/2

# first heat map of UVrep KOs
##heatmap(as.matrix(mr.cov.sum.UVrep.bins.ko[,-1]), scale = "none")

colSums(mr.cov.sum.UVrep.bins.ko[,-1])
#mr.cov.sum.UVrep.bins.ko2 <- mr.cov.sum.UVrep.bins.ko[,which(colSums(mr.cov.sum.UVrep.bins.ko[,-1])>10)]

##heatmap(as.matrix(mr.cov.sum.UVrep.bins.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap -- using MR data that includes NAs so they are blocked out on heatmaps
mr.UVrep.bins.ko[1:4,]
mr.UVrep.bins.all<-merge(mr.UVrep.bins.ko,all_bin_meta,by="Bin_ID")
head(mr.UVrep.bins.all)
#mr.UVrep.bins.all$PlotID<-mr.UVrep.bins.all$Bin_ID
mr.UVrep.bins.all$PlotID = factor(mr.UVrep.bins.all$Bin_ID, levels=unique(mr.UVrep.bins.all$Bin_ID[order(mr.UVrep.bins.all$Site,mr.UVrep.bins.all$SampDate)]), ordered=TRUE)

#mr.UVrep.bins.all$PathShort<-mr.UVrep.bins.all$Pathway
# mr.UVrep.bins.all$PathShort[(mr.UVrep.bins.all$PathShort) == "Dissimilatory Sulfate Redox"] <- "D.SO4 RedOx"

head(mr.UVrep.bins.all)

# convert all 0s to NAs so they appear gray in heat map
mr.UVrep.bins.all$MR_SumCovPerKO_NA<-ifelse(mr.UVrep.bins.all$MR_SumCovPerKO==0,NA,mr.UVrep.bins.all$MR_SumCovPerKO)
head(mr.UVrep.bins.all)

# For heatmap color gradient
max(mr.UVrep.bins.all$MR_SumCovPerKO, na.rm=TRUE)
max(mr.UVrep.bins.all$MR_SumCovPerKO, na.rm=TRUE)/2
min(mr.UVrep.bins.all$MR_SumCovPerKO, na.rm=TRUE)

# shorten KO fxn names for plotting
mr.UVrep.bins.all$KOFxn_Short<-gsub("; .*","",mr.UVrep.bins.all$KO_Function.KEGG)

unique(mr.UVrep.bins.all$SampDate)
mr.UVrep.bins.all$SampDatePlot<-as.character(mr.UVrep.bins.all$SampDate)
mr.UVrep.bins.all$SampDatePlot[which(mr.UVrep.bins.all$SampDatePlot == "August.2020")] <- "Aug 2020"
mr.UVrep.bins.all$SampDatePlot[which(mr.UVrep.bins.all$SampDatePlot == "November.2020")] <- "Nov 2020"
mr.UVrep.bins.all$SampDatePlot[which(mr.UVrep.bins.all$SampDatePlot == "August.2021")] <- "Aug 2021"
mr.UVrep.bins.all$SampDatePlot[which(mr.UVrep.bins.all$SampDatePlot == "September.2021")] <- "Sep 2021"
mr.UVrep.bins.all$SampDatePlot[which(mr.UVrep.bins.all$SampDatePlot == "December.2021")] <- "Dec 2021"
mr.UVrep.bins.all$SampDatePlot<-gsub("\\."," ",mr.UVrep.bins.all$SampDatePlot)

unique(mr.UVrep.bins.all$SampDatePlot)
mr.UVrep.bins.all$SampDatePlot = factor(mr.UVrep.bins.all$SampDatePlot, levels=unique(mr.UVrep.bins.all$SampDatePlot[order(mr.UVrep.bins.all$SampDate)]), ordered=TRUE)
unique(mr.UVrep.bins.all$SampDatePlot)

# merge UV repair KO coverages with MAG taxonomic annotations
mr.UVrep.taxa.all<-merge(mr.UVrep.bins.all,mgm.bin.tax,by=c("Bin_ID","SampleID"))
head(mr.UVrep.taxa.all)
mr.UVrep.taxa.all$Genus = factor(mr.UVrep.taxa.all$Genus, levels=unique(mr.UVrep.taxa.all$Genus[order(mr.UVrep.taxa.all$Site,mr.UVrep.taxa.all$SampDate)]), ordered=TRUE)

# Figures below
# by Bin_ID

UVrep.hm1a<-ggplot(mr.UVrep.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="UV-Damage DNA Repair Functions in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(UVrep.hm1a,filename = "figures/MGM_Figs/Bins/FxnDiv/UV_DNA_Repair/UV.DNA.Repair_KOFxns_MGMs_Bins_BinID_by_Function_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

UVrep.hm1a2<-ggplot(mr.UVrep.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=21),
        axis.text.y = element_text(size=22),axis.text.x = element_text(angle=45, hjust=1,size=19),legend.text = element_text(size=19),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),strip.text.x = element_text(size = 26,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(UVrep.hm1a2,filename = "figures/MGM_Figs/Bins/FxnDiv/UV_DNA_Repair/UV.DNA.Repair_KOFxns_MGMs_Bins_BinID_by_Function_Site_heatmap.png", width=47, height=15, dpi=300,create.dir = TRUE)

UVrep.hm1a3<-ggplot(mr.UVrep.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=21),
        axis.text.y = element_text(size=22),axis.text.x = element_text(angle=45, hjust=1,size=19),legend.text = element_text(size=19),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),strip.text.x = element_text(size = 26,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDatePlot,scales="free_x", space = "free")

ggsave(UVrep.hm1a3,filename = "figures/MGM_Figs/Bins/FxnDiv/UV_DNA_Repair/UV.DNA.Repair_KOFxns_MGMs_Bins_by_BinID_by_Function_SampDate_heatmap.png", width=49.5, height=13, dpi=300,create.dir = TRUE)

UVrep.hm1a4<-ggplot(mr.UVrep.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="UV-Damage DNA Repair Functions in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate_Plot+Site,scales="free_x", space = "free")

ggsave(UVrep.hm1a4,filename = "figures/MGM_Figs/Bins/FxnDiv/UV_DNA_Repair/UV.DNA.Repair_KOFxns_MGMs_Bins_by_BinID_by_Function_SampDate_Site_heatmap.png", width=49.5, height=13, dpi=600,create.dir = TRUE)

## by MAG Taxonomic ID, not Bin ID

UVrep.hm1b<-ggplot(mr.UVrep.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="UV-Damage DNA RepairFunctions in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(UVrep.hm1b,filename = "figures/MGM_Figs/Bins/FxnDiv/UV_DNA_Repair/UV.DNA.Repair_KOFxns_MGMs_Bins_Genus_by_Function_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

UVrep.hm1b2<-ggplot(mr.UVrep.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="UV-Damage DNA RepairFunctions in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=17),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=17),plot.title = element_blank(),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 17,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(UVrep.hm1b2,filename = "figures/MGM_Figs/Bins/FxnDiv/UV_DNA_Repair/UV.DNA.Repair_KOFxns_MGMs_Bins_Genus_by_Function_Site_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)

UVrep.hm1b3<-ggplot(mr.UVrep.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="UV-Damage DNA RepairFunctions in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=16, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=17),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate_Plot,scales="free_x", space = "free")

ggsave(UVrep.hm1b3,filename = "figures/MGM_Figs/Bins/FxnDiv/UV_DNA_Repair/UV.DNA.Repair_KOFxns_MGMs_Bins_Genus_by_Function_SampDate_heatmap.png", width=49.5, height=15, dpi=300,create.dir = TRUE)

UVrep.hm1b4<-ggplot(mr.UVrep.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="UV-Damage DNA RepairFunctions in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=16, family="Arial"), axis.title.x = element_text(size=25),axis.title.y = element_text(size=25),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=16),plot.title = element_blank(),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate_Plot+Site,scales="free_x", space = "free")

ggsave(UVrep.hm1b4,filename = "figures/MGM_Figs/Bins/FxnDiv/UV_DNA_Repair/UV.DNA.Repair_KOFxns_MGMs_Bins_Genus_by_Function_SampDate_Site_heatmap.png", width=49.5, height=15, dpi=600,create.dir = TRUE)

#### Pull out House Keeping Gene Coverage ####
# will use this later!
mgm.bin.mr$K03553 # recA is a house keeping gene

# pull out HK functions from MR normalized, summed coverages (summed gene coverage per KO)
HK.bin.ko<-data.frame(KO_ID="K03553",MR_SumCovPerKO=mgm.bin.mr$K03553,Bin_ID=rownames(mgm.bin.mr)) # merge MR data w/ S fxns found in contigs from KOFamScan
head(HK.bin.ko) #sanity check

mr.HK.bin.ko<-merge(HK.bin.ko,Housekeep.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(mr.HK.bin.ko)
colnames(mr.HK.bin.ko)[which(names(mr.HK.bin.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website

head(mr.HK.bin.ko)

## create table of Sample x KO_ID
mr.HK.bin.ko_table<-as.data.frame(dcast(mr.HK.bin.ko, Bin_ID~KO_Function.KEGG, value.var="MR_SumCovPerKO"))
mr.HK.bin.ko_table[1:5,]

## merge HK functions with scaled metadata
HK.bin.ko.all<-merge(meta.all.scaled,mr.HK.bin.ko_table,by="Bin_ID")
HK.bin.ko.all[1:5,]

# Find mean coverage of recA
HK.bin.mgm.means<-data.frame(MeanScaledCov=mean(HK.bin.ko.all$`recA; recombination protein RecA`),KO_ID=names(HK.bin.ko.all)[ncol(HK.bin.ko.all)])
rownames(HK.bin.mgm.means)<-HK.bin.mgm.means$KO_ID
HK.bin.mgm.means$Category<-"House Keeping"
HK.bin.mgm.means

#### Pull Out Metal Resistance/Tolerance Fxns from MR data ####
met.bins.ko<-mgm.bin.mr[,which(colnames(mgm.bin.mr) %in% metal.bin.fxns$KO_ID)]
met.bins.ko$Bin_ID<-rownames(met.bins.ko)
met.bins.ko.melt<-melt(met.bins.ko, by="Bin_ID")
colnames(met.bins.ko.melt)[which(names(met.bins.ko.melt) == "variable")] <- "KO_ID"
colnames(met.bins.ko.melt)[which(names(met.bins.ko.melt) == "value")] <- "MR_SumCovPerKO"
met.bins.ko.melt #sanity check

mr.met.bins.ko<-merge(met.bins.ko.melt,metal.re.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(mr.met.bins.ko)
colnames(mr.met.bins.ko)[which(names(mr.met.bins.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
mr.cov.sum.met.bins.ko<-as.data.frame(dcast(mr.met.bins.ko, Bin_ID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(mr.cov.sum.met.bins.ko)<-mr.cov.sum.met.bins.ko$Bin_ID
mr.cov.sum.met.bins.ko[1:4,]

#### Metal Resistance KO Heat Maps ####
# see max & mean of summed
max(mr.cov.sum.met.bins.ko[,-1])
min(mr.cov.sum.met.bins.ko[,-1])
max(mr.cov.sum.met.bins.ko[,-1])/2

# first heat map of met KOs
##heatmap(as.matrix(mr.cov.sum.met.bins.ko[,-1]), scale = "none")

colSums(mr.cov.sum.met.bins.ko[,-1])
#mr.cov.sum.met.bins.ko2 <- mr.cov.sum.met.bins.ko[,which(colSums(mr.cov.sum.met.bins.ko[,-1])>10)]

##heatmap(as.matrix(mr.cov.sum.met.bins.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap -- using MR data that includes NAs so they are blocked out on heatmaps
mr.met.bins.ko[1:4,]
mr.met.bins.all<-merge(mr.met.bins.ko,all_bin_meta,by="Bin_ID")
head(mr.met.bins.all)
#mr.met.bins.all$PlotID<-mr.met.bins.all$Bin_ID
mr.met.bins.all$PlotID = factor(mr.met.bins.all$Bin_ID, levels=unique(mr.met.bins.all$Bin_ID[order(mr.met.bins.all$Site,mr.met.bins.all$SampDate)]), ordered=TRUE)

#mr.met.bins.all$PathShort<-mr.met.bins.all$Pathway
# mr.met.bins.all$PathShort[(mr.met.bins.all$PathShort) == "Dissimilatory Sulfate Redox"] <- "D.SO4 RedOx"

head(mr.met.bins.all)

# convert all 0s to NAs so they appear gray in heat map
mr.met.bins.all$MR_SumCovPerKO_NA<-ifelse(mr.met.bins.all$MR_SumCovPerKO==0,NA,mr.met.bins.all$MR_SumCovPerKO)
head(mr.met.bins.all)

# For heatmap color gradient
max(mr.met.bins.all$MR_SumCovPerKO, na.rm=TRUE)
max(mr.met.bins.all$MR_SumCovPerKO, na.rm=TRUE)/2
min(mr.met.bins.all$MR_SumCovPerKO, na.rm=TRUE)

# merge metal resistance coverages with MAG taxonomic annotations
mr.met.taxa.all<-merge(mr.met.bins.all,mgm.bin.tax,by=c("Bin_ID","SampleID"))
head(mr.met.taxa.all)
mr.met.taxa.all$Genus = factor(mr.met.taxa.all$Genus, levels=unique(mr.met.taxa.all$Genus[order(mr.met.taxa.all$Site,mr.met.taxa.all$SampDate)]), ordered=TRUE)

# Figures below
# by Bin_ID

met.hm1a<-ggplot(mr.met.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Metal Resistance in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(met.hm1a,filename = "figures/MGM_Figs/Bins/FxnDiv/MetalResistance/MetalResistance_KOFxns_MGMs_Bins_BinID_by_Function_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

met.hm1a2<-ggplot(mr.met.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Metal Resistance in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(met.hm1a2,filename = "figures/MGM_Figs/Bins/FxnDiv/MetalResistance/MetalResistance_KOFxns_MGMs_Bins_BinID_by_Function_Site_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)

met.hm1a3<-ggplot(mr.met.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Metal Resistance in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate,scales="free_x", space = "free")

ggsave(met.hm1a3,filename = "figures/MGM_Figs/Bins/FxnDiv/MetalResistance/MetalResistance_KOFxns_MGMs_Bins_by_BinID_by_Function_SampDate_heatmap.png", width=49.5, height=13, dpi=600,create.dir = TRUE)

met.hm1a4<-ggplot(mr.met.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Metal Resistance in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate+Site,scales="free_x", space = "free")

ggsave(met.hm1a4,filename = "figures/MGM_Figs/Bins/FxnDiv/MetalResistance/MetalResistance_KOFxns_MGMs_Bins_by_BinID_by_Function_SampDate_Site_heatmap.png", width=49.5, height=13, dpi=600,create.dir = TRUE)

## by MAG Taxonomic ID, not Bin ID

met.hm1b<-ggplot(mr.met.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="UV-Damage DNA RepairFunctions in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(met.hm1b,filename = "figures/MGM_Figs/Bins/FxnDiv/MetalResistance/MetalResistance_KOFxns_MGMs_Bins_Genus_by_Function_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

met.hm1b2<-ggplot(mr.met.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="UV-Damage DNA RepairFunctions in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(met.hm1b2,filename = "figures/MGM_Figs/Bins/FxnDiv/MetalResistance/MetalResistance_KOFxns_MGMs_Bins_Genus_by_Function_Site_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)

met.hm1b3<-ggplot(mr.met.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="UV-Damage DNA RepairFunctions in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate,scales="free_x", space = "free")

ggsave(met.hm1b3,filename = "figures/MGM_Figs/Bins/FxnDiv/MetalResistance/MetalResistance_KOFxns_MGMs_Bins_Genus_by_Function_SampDate_heatmap.png", width=49.5, height=15, dpi=600,create.dir = TRUE)

met.hm1b4<-ggplot(mr.met.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="UV-Damage DNA RepairFunctions in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate+Site,scales="free_x", space = "free")

ggsave(met.hm1b4,filename = "figures/MGM_Figs/Bins/FxnDiv/MetalResistance/MetalResistance_KOFxns_MGMs_Bins_Genus_by_Function_SampDate_Site_heatmap.png", width=49.5, height=15, dpi=600,create.dir = TRUE)

#### Pull Out Osmoprotectant Production Fxns from MR data ####
osmo.bins.ko<-mgm.bin.mr[,which(colnames(mgm.bin.mr) %in% osmo.bin.fxns$KO_ID)]
osmo.bins.ko$Bin_ID<-rownames(osmo.bins.ko)
osmo.bins.ko.melt<-melt(osmo.bins.ko, by="Bin_ID")
colnames(osmo.bins.ko.melt)[which(names(osmo.bins.ko.melt) == "variable")] <- "KO_ID"
colnames(osmo.bins.ko.melt)[which(names(osmo.bins.ko.melt) == "value")] <- "MR_SumCovPerKO"
osmo.bins.ko.melt #sanity check

mr.osmo.bins.ko<-merge(osmo.bins.ko.melt,osmo.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(mr.osmo.bins.ko)
colnames(mr.osmo.bins.ko)[which(names(mr.osmo.bins.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
mr.cov.sum.osmo.bins.ko<-as.data.frame(dcast(mr.osmo.bins.ko, Bin_ID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(mr.cov.sum.osmo.bins.ko)<-mr.cov.sum.osmo.bins.ko$Bin_ID
mr.cov.sum.osmo.bins.ko[1:4,]


#### Osmoprotectant Production KO Heat Maps ####
# see max & mean of summed
max(mr.cov.sum.osmo.bins.ko[,-1])
min(mr.cov.sum.osmo.bins.ko[,-1])
max(mr.cov.sum.osmo.bins.ko[,-1])/2

# first heat map of osmo KOs
##heatmap(as.matrix(mr.cov.sum.osmo.bins.ko[,-1]), scale = "none")

colSums(mr.cov.sum.osmo.bins.ko[,-1])
#mr.cov.sum.osmo.bins.ko2 <- mr.cov.sum.osmo.bins.ko[,which(colSums(mr.cov.sum.osmo.bins.ko[,-1])>10)]

##heatmap(as.matrix(mr.cov.sum.osmo.bins.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap -- using MR data that includes NAs so they are blocked out on heatmaps
mr.osmo.bins.ko[1:4,]
mr.osmo.bins.all<-merge(mr.osmo.bins.ko,all_bin_meta,by="Bin_ID")
head(mr.osmo.bins.all)
#mr.osmo.bins.all$PlotID<-mr.osmo.bins.all$Bin_ID
mr.osmo.bins.all$PlotID = factor(mr.osmo.bins.all$Bin_ID, levels=unique(mr.osmo.bins.all$Bin_ID[order(mr.osmo.bins.all$Site,mr.osmo.bins.all$SampDate)]), ordered=TRUE)

#mr.osmo.bins.all$PathShort<-mr.osmo.bins.all$Pathway
# mr.osmo.bins.all$PathShort[(mr.osmo.bins.all$PathShort) == "Dissimilatory Sulfate Redox"] <- "D.SO4 RedOx"

head(mr.osmo.bins.all)

# convert all 0s to NAs so they appear gray in heat map
mr.osmo.bins.all$MR_SumCovPerKO_NA<-ifelse(mr.osmo.bins.all$MR_SumCovPerKO==0,NA,mr.osmo.bins.all$MR_SumCovPerKO)
head(mr.osmo.bins.all)

# For heatmap color gradient
max(mr.osmo.bins.all$MR_SumCovPerKO, na.rm=TRUE)
max(mr.osmo.bins.all$MR_SumCovPerKO, na.rm=TRUE)/2
min(mr.osmo.bins.all$MR_SumCovPerKO, na.rm=TRUE)

# merge osmoprotectant production coverages with MAG taxonomic annotations
mr.osmo.taxa.all<-merge(mr.osmo.bins.all,mgm.bin.tax,by=c("Bin_ID","SampleID"))
head(mr.osmo.taxa.all)
mr.osmo.taxa.all$Genus = factor(mr.osmo.taxa.all$Genus, levels=unique(mr.osmo.taxa.all$Genus[order(mr.osmo.taxa.all$Site,mr.osmo.taxa.all$SampDate)]), ordered=TRUE)

# Figures below
# by Bin_ID

osmo.hm1a<-ggplot(mr.osmo.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Osmoprotectant Production in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(osmo.hm1a,filename = "figures/MGM_Figs/Bins/FxnDiv/Osmoprotectant/Osmoprotectant_KOFxns_MGMs_Bins_BinID_by_Function_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

osmo.hm1a2<-ggplot(mr.osmo.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Osmoprotectant Production in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(osmo.hm1a2,filename = "figures/MGM_Figs/Bins/FxnDiv/Osmoprotectant/Osmoprotectant_KOFxns_MGMs_Bins_BinID_by_Function_Site_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)

osmo.hm1a3<-ggplot(mr.osmo.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Osmoprotectant Production in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate,scales="free_x", space = "free")

ggsave(osmo.hm1a3,filename = "figures/MGM_Figs/Bins/FxnDiv/Osmoprotectant/Osmoprotectant_KOFxns_MGMs_Bins_by_BinID_by_Function_SampDate_heatmap.png", width=49.5, height=13, dpi=600,create.dir = TRUE)

osmo.hm1a4<-ggplot(mr.osmo.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Osmoprotectant Production in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate+Site,scales="free_x", space = "free")

ggsave(osmo.hm1a4,filename = "figures/MGM_Figs/Bins/FxnDiv/Osmoprotectant/Osmoprotectant_KOFxns_MGMs_Bins_by_BinID_by_Function_SampDate_Site_heatmap.png", width=49.5, height=13, dpi=600,create.dir = TRUE)

## by MAG Taxonomic ID, not Bin ID

osmo.hm1b<-ggplot(mr.osmo.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Osmoprotectant Production in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(osmo.hm1b,filename = "figures/MGM_Figs/Bins/FxnDiv/Osmoprotectant/Osmoprotectant_KOFxns_MGMs_Bins_Genus_by_Function_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

osmo.hm1b2<-ggplot(mr.osmo.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Osmoprotectant Production  in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(osmo.hm1b2,filename = "figures/MGM_Figs/Bins/FxnDiv/Osmoprotectant/Osmoprotectant_KOFxns_MGMs_Bins_Genus_by_Function_Site_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)

osmo.hm1b3<-ggplot(mr.osmo.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Osmoprotectant Production in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate,scales="free_x", space = "free")

ggsave(osmo.hm1b3,filename = "figures/MGM_Figs/Bins/FxnDiv/Osmoprotectant/Osmoprotectant_KOFxns_MGMs_Bins_Genus_by_Function_SampDate_heatmap.png", width=49.5, height=15, dpi=600,create.dir = TRUE)

osmo.hm1b4<-ggplot(mr.osmo.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Osmoprotectant Production in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate+Site,scales="free_x", space = "free")

ggsave(osmo.hm1b4,filename = "figures/MGM_Figs/Bins/FxnDiv/Osmoprotectant/Osmoprotectant_KOFxns_MGMs_Bins_Genus_by_Function_SampDate_Site_heatmap.png", width=49.5, height=15, dpi=600,create.dir = TRUE)

#### Pull Out ALL Genes of Interest per Pathway/Cycle from MR data ####
## heatmaps of traits of interest

mgm.bin.mr[1:4,1:4]

# pull out sulfur functions from MR transformed, summed coverages (summed coverage per KO)
All_GOI.bins.ko<-mgm.bin.mr[,which(colnames(mgm.bin.mr) %in% All_GOI.bin.fxns$KO_ID)] # merge MR data w/ All_GOI-related fxns found in contigs from KOFamScan
All_GOI.bins.ko$Bin_ID<-rownames(All_GOI.bins.ko)
NA %in% All_GOI.bins.ko
grep("Bin_ID", names(All_GOI.bins.ko)) # find index for column we want aka "SampleID" column

mean(as.matrix(All_GOI.bins.ko[,-151])) # mean KO coverage
min(as.matrix(All_GOI.bins.ko[,-151])) # mean KO coverage

All_GOI.bins.ko.melt<-melt(All_GOI.bins.ko, by="Bin_ID")
colnames(All_GOI.bins.ko.melt)[which(names(All_GOI.bins.ko.melt) == "variable")] <- "KO_ID"
colnames(All_GOI.bins.ko.melt)[which(names(All_GOI.bins.ko.melt) == "value")] <- "MR_SumCovPerKO"
head(All_GOI.bins.ko.melt) #sanity check

mr.bins.all_GOI.bins.ko<-merge(All_GOI.bins.ko.melt,all_goi.kegg,by.x=c("KO_ID"),by.y=c("KO_ID"))
head(mr.bins.all_GOI.bins.ko)
colnames(mr.bins.all_GOI.bins.ko)[which(names(mr.bins.all_GOI.bins.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
mr.cov.sum.bins.all_GOI.bins.ko<-as.data.frame(dcast(mr.bins.all_GOI.bins.ko, Bin_ID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###
rownames(mr.cov.sum.bins.all_GOI.bins.ko)<-mr.cov.sum.bins.all_GOI.bins.ko$Bin_ID
mr.cov.sum.bins.all_GOI.bins.ko[1:4,1:4]

# how many MGMs have certain shared functions of interest?
# created this section after reviewing heatmaps to see how many MAGs contain certain functions

## UV DNA repair
length(mr.cov.sum.bins.all_GOI.bins.ko$`uvrA; excinuclease ABC subunit A`[mr.cov.sum.bins.all_GOI.bins.ko$`uvrA; excinuclease ABC subunit A`>0])
# 83 MAGs have uvrA
length(mr.cov.sum.bins.all_GOI.bins.ko$`uvrB; excinuclease ABC subunit B`[mr.cov.sum.bins.all_GOI.bins.ko$`uvrB; excinuclease ABC subunit B`>0])
# 85 MAGs have uvrB
length(mr.cov.sum.bins.all_GOI.bins.ko$`uvrC; excinuclease ABC subunit C`[mr.cov.sum.bins.all_GOI.bins.ko$`uvrC; excinuclease ABC subunit C`>0])
# 87 MAGs have uvrC
length(mr.cov.sum.bins.all_GOI.bins.ko$`lexA; repressor LexA [EC:3.4.21.88]`[mr.cov.sum.bins.all_GOI.bins.ko$`lexA; repressor LexA [EC:3.4.21.88]`>0])
# 60 MAGs have lexA

# Tempshock
length(mr.cov.sum.bins.all_GOI.bins.ko$`cspA; cold shock protein`[mr.cov.sum.bins.all_GOI.bins.ko$`cspA; cold shock protein`>0])
# 88 MAGs have cspA
length(mr.cov.sum.bins.all_GOI.bins.ko$`hslR; ribosome-associated heat shock protein Hsp15`[mr.cov.sum.bins.all_GOI.bins.ko$`hslR; ribosome-associated heat shock protein Hsp15`>0])
# 67 MAGs have hsIR
length(mr.cov.sum.bins.all_GOI.bins.ko$`htpX; heat shock protein HtpX [EC:3.4.24.-]`[mr.cov.sum.bins.all_GOI.bins.ko$`htpX; heat shock protein HtpX [EC:3.4.24.-]`>0])
# 39 MAGs have htpX

## Quorum sensing
length(mr.cov.sum.bins.all_GOI.bins.ko$`bisR; LuxR family transcriptional regulator, quorum-sensing system regulator BisR`[mr.cov.sum.bins.all_GOI.bins.ko$`bisR; LuxR family transcriptional regulator, quorum-sensing system regulator BisR`>0])
# 53 MAGs have bisR

## Osmoprotectant transport
length(mr.cov.sum.bins.all_GOI.bins.ko$`opuC; osmoprotectant transport system substrate-binding protein`[mr.cov.sum.bins.all_GOI.bins.ko$`opuC; osmoprotectant transport system substrate-binding protein`>0])
# 35 MAGs have opuC
length(mr.cov.sum.bins.all_GOI.bins.ko$`osmY; hyperosmotically inducible periplasmic protein`[mr.cov.sum.bins.all_GOI.bins.ko$`osmY; hyperosmotically inducible periplasmic protein`>0])
# 28 MAGs have osmY

## Sporulation
length(mr.cov.sum.bins.all_GOI.bins.ko$`spmA; spore maturation protein A`[mr.cov.sum.bins.all_GOI.bins.ko$`spmA; spore maturation protein A`>0])
# 24 MAGs have spmA
length(mr.cov.sum.bins.all_GOI.bins.ko$`spoVR; stage V sporulation protein R`[mr.cov.sum.bins.all_GOI.bins.ko$`spoVR; stage V sporulation protein R`>0])
# 25 MAGs have spoVR

## LPS
length(mr.cov.sum.bins.all_GOI.bins.ko$`wecP; UDP-GalNAc:undecaprenyl-phosphate GalNAc-1-phosphate transferase [EC:2.7.8.40]`[mr.cov.sum.bins.all_GOI.bins.ko$`wecP; UDP-GalNAc:undecaprenyl-phosphate GalNAc-1-phosphate transferase [EC:2.7.8.40]`>0])
# 41 MAGs have wecP
length(mr.cov.sum.bins.all_GOI.bins.ko$`wbgP; O55-antigen biosynthesis glycosyltransferase`[mr.cov.sum.bins.all_GOI.bins.ko$`wbgP; O55-antigen biosynthesis glycosyltransferase`>0])
# 51 MAGs have wbgP
length(mr.cov.sum.bins.all_GOI.bins.ko$`wbbC; O7-antigen biosynthesis mannosyltransferase`[mr.cov.sum.bins.all_GOI.bins.ko$`wbbC; O7-antigen biosynthesis mannosyltransferase`>0])
# 20 MAGs have wbbC
length(mr.cov.sum.bins.all_GOI.bins.ko$`lpxD; UDP-3-O-[3-hydroxymyristoyl] glucosamine N-acyltransferase [EC:2.3.1.191]`[mr.cov.sum.bins.all_GOI.bins.ko$`lpxD; UDP-3-O-[3-hydroxymyristoyl] glucosamine N-acyltransferase [EC:2.3.1.191]`>0])
# 41 MAGs have lpxD
length(mr.cov.sum.bins.all_GOI.bins.ko$`gtrB; polyisoprenyl-phosphate glycosyltransferase [EC:2.4.-.-]`[mr.cov.sum.bins.all_GOI.bins.ko$`gtrB; polyisoprenyl-phosphate glycosyltransferase [EC:2.4.-.-]`>0])
# 43 MAGs hae gtrB

## House keeping gene recA
length(mr.cov.sum.bins.all_GOI.bins.ko$`recA; recombination protein RecA`[mr.cov.sum.bins.all_GOI.bins.ko$`recA; recombination protein RecA`>0])
# 90 MAGs have recA

# how many genes are we looking for per category of interest?
length(all_goi.kegg$Category[all_goi.kegg$Category=="LPS Biosynthesis"]) # 161 total genes we are looking for
length(all_goi.kegg$Category[all_goi.kegg$Category=="Sporulation"]) # 38 total genes we are looking for
length(all_goi.kegg$Category[all_goi.kegg$Category=="Temperature Shock Resistance"]) # 29 total genes we are looking for
length(all_goi.kegg$Category[all_goi.kegg$Category=="Osmoprotectant Production/Accumulation"]) # 11 total genes we are looking for
length(all_goi.kegg$Category[all_goi.kegg$Category=="Quorum Sensing"]) # 37 total genes we are looking for
length(all_goi.kegg$Category[all_goi.kegg$Category=="UV-Damange DNA Repair"]) # 6 total genes we are looking for

all_goi.kegg[all_goi.kegg$Category=="Temperature Shock Resistance",] # sanity check - just to review which functions these are

### ALL Genes of Interest per Pathway Heat Maps ####
# see max & mean of summed
max(mr.cov.sum.bins.all_GOI.bins.ko[,-1])
mean(as.matrix(mr.cov.sum.bins.all_GOI.bins.ko[,-1]))

# first heat map of sulfur KOs
#heatmap(as.matrix(mr.cov.sum.bins.all_GOI.bins.ko[,-1]), scale = "none")

colSums(mr.cov.sum.bins.all_GOI.bins.ko[,-1])

#heatmap(as.matrix(mr.cov.sum.bins.all_GOI.bins.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap
mr.bins.all_GOI.bins.ko[1:4,]
mr.bins.all_GOI.all<-merge(mr.bins.all_GOI.bins.ko,all_bin_meta,by="Bin_ID")
head(mr.bins.all_GOI.all)
mr.bins.all_GOI.all$PlotID = factor(mr.bins.all_GOI.all$Bin_ID, levels=unique(mr.bins.all_GOI.all$Bin_ID[order(mr.bins.all_GOI.all$Site,mr.bins.all_GOI.all$SampDate)]), ordered=TRUE)

unique(mr.bins.all_GOI.all$Category)
mr.bins.all_GOI.all$CatShort<-mr.bins.all_GOI.all$Category
mr.bins.all_GOI.all$CatShort[(mr.bins.all_GOI.all$CatShort) == "LPS Biosynthesis"] <- "LPS"
mr.bins.all_GOI.all$CatShort[(mr.bins.all_GOI.all$CatShort) == "Sporulation"] <- "Spore"
mr.bins.all_GOI.all$CatShort[(mr.bins.all_GOI.all$CatShort) == "Temperature Shock Resistance"] <- "Temp."
mr.bins.all_GOI.all$CatShort[(mr.bins.all_GOI.all$CatShort) == "Osmoprotectant Production/Accumulation"] <- "Osmo"
mr.bins.all_GOI.all$CatShort[(mr.bins.all_GOI.all$CatShort) == "Quorum Sensing"] <- "QS"
mr.bins.all_GOI.all$CatShort[(mr.bins.all_GOI.all$CatShort) == "UV-Damange DNA Repair"] <- "UV"
mr.bins.all_GOI.all$CatShort[(mr.bins.all_GOI.all$CatShort) == "House Keeping Gene"] <- "House"
unique(mr.bins.all_GOI.all$CatShort)

head(mr.bins.all_GOI.all)

# order PlotID by site & sampdate
mr.bins.all_GOI.all$PlotID = factor(mr.bins.all_GOI.all$Bin_ID, levels=unique(mr.bins.all_GOI.all$Bin_ID[order(mr.bins.all_GOI.all$Site,mr.bins.all_GOI.all$SampDate)]), ordered=TRUE)

# convert all 0s to NAs so they appear gray in heat map
mr.bins.all_GOI.all$MR_SumCovPerKO_NA<-ifelse(mr.bins.all_GOI.all$MR_SumCovPerKO==0,NA,mr.bins.all_GOI.all$MR_SumCovPerKO)
head(mr.bins.all_GOI.all)

# how many genes do we have across the MAGs per category? aka genes with a normalized coverage > 0
length(unique(mr.bins.all_GOI.all$KO_Function.KEGG[mr.bins.all_GOI.all$CatShort=="LPS" & mr.bins.all_GOI.all$MR_SumCovPerKO>0]))
# 82 out of 161 genes investigated here
length(unique(mr.bins.all_GOI.all$KO_Function.KEGG[mr.bins.all_GOI.all$CatShort=="Spore" & mr.bins.all_GOI.all$MR_SumCovPerKO>0]))
# 35 out of 38 genes investigated here
length(unique(mr.bins.all_GOI.all$KO_Function.KEGG[mr.bins.all_GOI.all$CatShort=="Temp." & mr.bins.all_GOI.all$MR_SumCovPerKO>0]))
# 8 out of 29 genes investigated here
length(unique(mr.bins.all_GOI.all$KO_Function.KEGG[mr.bins.all_GOI.all$CatShort=="Osmo" & mr.bins.all_GOI.all$MR_SumCovPerKO>0]))
# 7 out of 11 genes investigated here
length(unique(mr.bins.all_GOI.all$KO_Function.KEGG[mr.bins.all_GOI.all$CatShort=="QS" & mr.bins.all_GOI.all$MR_SumCovPerKO>0]))
# 11 out of the 37 genes investigated here
length(unique(mr.bins.all_GOI.all$KO_Function.KEGG[mr.bins.all_GOI.all$CatShort=="UV" & mr.bins.all_GOI.all$MR_SumCovPerKO>0]))
# 6 out of the 6 genes investigated here

# shorten KO fxn names for plotting
mr.bins.all_GOI.all$KOFxn_Short<-gsub("; .*","",mr.bins.all_GOI.all$KO_Function.KEGG)

# merge GOI KO coverages with MAG taxonomic annotations
mr.All_GOI.taxa.all<-merge(mr.bins.all_GOI.all,mgm.bin.tax,by=c("Bin_ID","SampleID"))
head(mr.All_GOI.taxa.all)
mr.All_GOI.taxa.all$Genus = factor(mr.All_GOI.taxa.all$Genus, levels=unique(mr.All_GOI.taxa.all$Genus[order(mr.All_GOI.taxa.all$Site,mr.All_GOI.taxa.all$SampDate)]), ordered=TRUE)

# Figures below
# By Bin ID

All_GOI.hm1a<-ggplot(mr.bins.all_GOI.all, aes(PlotID, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Dust Microbial Survival in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(All_GOI.hm1a,filename = "figures/MGM_Figs/Bins/FxnDiv/All_GOI/ALL_GOI_KOFxns_MGMs_Bins_by_Function_heatmap.png", width=45, height=30, dpi=600,create.dir = TRUE)

All_GOI.hm1a2<-ggplot(mr.bins.all_GOI.all, aes(PlotID, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Dust Microbial Survival in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(All_GOI.hm1a2,filename = "figures/MGM_Figs/Bins/FxnDiv/All_GOI/ALL_GOI_KOFxns_MGMs_Bins_by_Function_Site_heatmap.png", width=45, height=30, dpi=600,create.dir = TRUE)

All_GOI.hm1a2b<-ggplot(mr.bins.all_GOI.all[mr.bins.all_GOI.all$MR_SumCovPerKO>=3,], aes(PlotID, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Dust Microbial Survival in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Normalized Coverages >= 3)",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(All_GOI.hm1a2b,filename = "figures/MGM_Figs/Bins/FxnDiv/All_GOI/ALL_GOI_KOFxns_MGMs_Bins_by_Function_Site_HigherCov_heatmap.png", width=45, height=30, dpi=600,create.dir = TRUE)

All_GOI.hm1a2c<-ggplot(mr.bins.all_GOI.all[mr.bins.all_GOI.all$MR_SumCovPerKO>=3,], aes(PlotID, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Dust Microbial Survival in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Normalized Coverages >= 3)",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate_Plot,scales="free_x", space = "free")

ggsave(All_GOI.hm1a2c,filename = "figures/MGM_Figs/Bins/FxnDiv/All_GOI/ALL_GOI_KOFxns_MGMs_Bins_by_Function_SampDate_HigherCov_heatmap.png", width=45, height=30, dpi=600,create.dir = TRUE)

All_GOI.hm1a2d<-ggplot(mr.bins.all_GOI.all[mr.bins.all_GOI.all$MR_SumCovPerKO>=3,], aes(PlotID, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) + geom_text(label=round(mr.bins.all_GOI.all$MR_SumCovPerKO[mr.bins.all_GOI.all$MR_SumCovPerKO>=3],2)) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Dust Microbial Survival in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Normalized Coverages >= 3)",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate_Plot,scales="free_x", space = "free")

ggsave(All_GOI.hm1a2d,filename = "figures/MGM_Figs/Bins/FxnDiv/All_GOI/ALL_GOI_KOFxns_MGMs_Bins_by_Function_SampDate_HigherCov_heatmap_labeled.png", width=45, height=30, dpi=600,create.dir = TRUE)

All_GOI.hm1a3<-ggplot(mr.bins.all_GOI.all, aes(PlotID, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Dust Microbial Survival in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(CatShort~Site, scales="free", space = "free")

ggsave(All_GOI.hm1a3,filename = "figures/MGM_Figs/Bins/FxnDiv/All_GOI/ALL_GOI_KOFxns_MGMs_Bins_by_Function_KO_Category_Site_heatmap.png", width=45, height=30, dpi=600,create.dir = TRUE)

All_GOI.hm1a3a<-ggplot(mr.bins.all_GOI.all[mr.bins.all_GOI.all$MR_SumCovPerKO>=3,], aes(PlotID, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Dust Microbial Survival in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Normalized Coverages >= 3)",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(CatShort~Site, scales="free", space = "free")

ggsave(All_GOI.hm1a3a,filename = "figures/MGM_Figs/Bins/FxnDiv/All_GOI/ALL_GOI_KOFxns_MGMs_Bins_by_Function_KO_Category_Site_HigherCov_heatmap.png", width=45, height=30, dpi=600,create.dir = TRUE)

All_GOI.hm1a3b<-ggplot(mr.bins.all_GOI.all[mr.bins.all_GOI.all$MR_SumCovPerKO>=3,], aes(PlotID, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) + geom_text(label=round(mr.bins.all_GOI.all$MR_SumCovPerKO[mr.bins.all_GOI.all$MR_SumCovPerKO>=3],2)) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Dust Microbial Survival in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Normalized Coverages >= 3)",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(CatShort~Site, scales="free", space = "free")

ggsave(All_GOI.hm1a3b,filename = "figures/MGM_Figs/Bins/FxnDiv/All_GOI/ALL_GOI_KOFxns_MGMs_Bins_by_Function_KO_Category_Site_HigherCov_heatmap_labeled.png", width=45, height=30, dpi=600,create.dir = TRUE)

All_GOI.hm1a4<-ggplot(mr.bins.all_GOI.all, aes(PlotID, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(fill="MR Coverage Per KO") +
  theme(text=element_text(size=15, family="Arial"), legend.title.align=0.5, legend.title = element_text(size=22),
        axis.text = element_text(size=20),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 15,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(CatShort~SampDate_Plot, scales="free", space = "free")

ggsave(All_GOI.hm1a4,filename = "figures/MGM_Figs/Bins/FxnDiv/All_GOI/ALL_GOI_KOFxns_MGMs_Bins_by_Function_KO_Category_SampDate_heatmap.png", width=45, height=40, dpi=600,create.dir = TRUE)

All_GOI.hm1a4a<-ggplot(mr.bins.all_GOI.all[mr.bins.all_GOI.all$MR_SumCovPerKO>=3,], aes(PlotID, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Dust Microbial Survival in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Normalized Coverages >= 3)",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(CatShort~SampDate_Plot, scales="free", space = "free")

ggsave(All_GOI.hm1a4a,filename = "figures/MGM_Figs/Bins/FxnDiv/All_GOI/ALL_GOI_KOFxns_MGMs_Bins_by_Function_KO_Category_SampDate_HigherCov_heatmap.png", width=45, height=30, dpi=600,create.dir = TRUE)

All_GOI.hm1a4b<-ggplot(mr.bins.all_GOI.all[mr.bins.all_GOI.all$MR_SumCovPerKO>=3,], aes(PlotID, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) + geom_text(label=round(mr.bins.all_GOI.all$MR_SumCovPerKO[mr.bins.all_GOI.all$MR_SumCovPerKO>=3],2)) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Dust Microbial Survival in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Normalized Coverages >= 3)",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(CatShort~SampDate_Plot, scales="free", space = "free")

ggsave(All_GOI.hm1a4b,filename = "figures/MGM_Figs/Bins/FxnDiv/All_GOI/ALL_GOI_KOFxns_MGMs_Bins_by_Function_KO_Category_SampDate_HigherCov_heatmap_labeled.png", width=45, height=30, dpi=600,create.dir = TRUE)

# by Taxa ID

All_GOI.hm1b<-ggplot(mr.All_GOI.taxa.all, aes(Genus, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Dust Microbial Survival in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(All_GOI.hm1b,filename = "figures/MGM_Figs/Bins/FxnDiv/All_GOI/ALL_GOI_KOFxns_MGMs_Bins_by_Function_heatmap.png", width=30, height=30, dpi=600,create.dir = TRUE)

All_GOI.hm1b2<-ggplot(mr.All_GOI.taxa.all, aes(Genus, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Dust Microbial Survival in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(All_GOI.hm1b2,filename = "figures/MGM_Figs/Bins/FxnDiv/All_GOI/ALL_GOI_KOFxns_MGMs_Bins_by_Function_Site_heatmap.png", width=30, height=30, dpi=600,create.dir = TRUE)

All_GOI.hm1b2b<-ggplot(mr.All_GOI.taxa.all[mr.All_GOI.taxa.all$MR_SumCovPerKO>=3,], aes(Genus, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Dust Microbial Survival in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Normalized Coverages >= 3)",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(All_GOI.hm1b2b,filename = "figures/MGM_Figs/Bins/FxnDiv/All_GOI/ALL_GOI_KOFxns_MGMs_Bins_by_Function_Site_HigherCov_heatmap.png", width=30, height=30, dpi=600,create.dir = TRUE)

All_GOI.hm1b2c<-ggplot(mr.All_GOI.taxa.all[mr.All_GOI.taxa.all$MR_SumCovPerKO>=3,], aes(Genus, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Dust Microbial Survival in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Normalized Coverages >= 3)",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate_Plot,scales="free_x", space = "free")

ggsave(All_GOI.hm1b2c,filename = "figures/MGM_Figs/Bins/FxnDiv/All_GOI/ALL_GOI_KOFxns_MGMs_Bins_by_Function_SampDate_HigherCov_heatmap.png", width=30, height=30, dpi=600,create.dir = TRUE)

All_GOI.hm1b2d<-ggplot(mr.All_GOI.taxa.all[mr.All_GOI.taxa.all$MR_SumCovPerKO>=3,], aes(Genus, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) + geom_text(label=round(mr.All_GOI.taxa.all$MR_SumCovPerKO[mr.All_GOI.taxa.all$MR_SumCovPerKO>=3],2)) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Dust Microbial Survival in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Normalized Coverages >= 3)",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate_Plot,scales="free_x", space = "free")

ggsave(All_GOI.hm1b2d,filename = "figures/MGM_Figs/Bins/FxnDiv/All_GOI/ALL_GOI_KOFxns_MGMs_Bins_by_Function_SampDate_HigherCov_heatmap_labeled.png", width=35, height=30, dpi=600,create.dir = TRUE)

All_GOI.hm1b3<-ggplot(mr.All_GOI.taxa.all, aes(Genus, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Dust Microbial Survival in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(CatShort~Site, scales="free", space = "free")

ggsave(All_GOI.hm1b3,filename = "figures/MGM_Figs/Bins/FxnDiv/All_GOI/ALL_GOI_KOFxns_MGMs_Bins_Genus_by_Function_KO_Category_Site_heatmap.png", width=35, height=30, dpi=600,create.dir = TRUE)

All_GOI.hm1b3a<-ggplot(mr.All_GOI.taxa.all[mr.All_GOI.taxa.all$MR_SumCovPerKO>=3,], aes(Genus, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Dust Microbial Survival in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Normalized Coverages >= 3)",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(CatShort~Site, scales="free", space = "free")

ggsave(All_GOI.hm1b3a,filename = "figures/MGM_Figs/Bins/FxnDiv/All_GOI/ALL_GOI_KOFxns_MGMs_Bins_Genus_by_Function_KO_Category_Site_HigherCov_heatmap.png", width=35, height=30, dpi=600,create.dir = TRUE)

All_GOI.hm1b3b<-ggplot(mr.All_GOI.taxa.all[mr.All_GOI.taxa.all$MR_SumCovPerKO>=3,], aes(Genus, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) + geom_text(label=round(mr.All_GOI.taxa.all$MR_SumCovPerKO[mr.All_GOI.taxa.all$MR_SumCovPerKO>=3],2)) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Dust Microbial Survival in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Normalized Coverages >= 3)",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(CatShort~Site, scales="free", space = "free")

ggsave(All_GOI.hm1b3b,filename = "figures/MGM_Figs/Bins/FxnDiv/All_GOI/ALL_GOI_KOFxns_MGMs_Bins_Genus_by_Function_KO_Category_Site_HigherCov_heatmap_labeled.png", width=35, height=30, dpi=600,create.dir = TRUE)

All_GOI.hm1b4<-ggplot(mr.All_GOI.taxa.all, aes(Genus, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Dust Microbial Survival in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(CatShort~SampDate_Plot, scales="free", space = "free")

ggsave(All_GOI.hm1b4,filename = "figures/MGM_Figs/Bins/FxnDiv/All_GOI/ALL_GOI_KOFxns_MGMs_Bins_Genus_by_Function_KO_Category_SampDate_heatmap.png", width=35, height=30, dpi=600,create.dir = TRUE)

All_GOI.hm1b4a<-ggplot(mr.All_GOI.taxa.all[mr.All_GOI.taxa.all$MR_SumCovPerKO>=3,], aes(Genus, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Dust Microbial Survival in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Normalized Coverages >= 3)",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(CatShort~SampDate_Plot, scales="free", space = "free")

ggsave(All_GOI.hm1b4a,filename = "figures/MGM_Figs/Bins/FxnDiv/All_GOI/ALL_GOI_KOFxns_MGMs_Bins_Genus_by_Function_KO_Category_SampDate_HigherCov_heatmap.png", width=35, height=30, dpi=600,create.dir = TRUE)

All_GOI.hm1b4b<-ggplot(mr.All_GOI.taxa.all[mr.All_GOI.taxa.all$MR_SumCovPerKO>=3,], aes(Genus, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) + geom_text(label=round(mr.All_GOI.taxa.all$MR_SumCovPerKO[mr.All_GOI.taxa.all$MR_SumCovPerKO>=3],2)) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Dust Microbial Survival in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Normalized Coverages >= 3)",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(CatShort~SampDate_Plot, scales="free", space = "free")

ggsave(All_GOI.hm1b4b,filename = "figures/MGM_Figs/Bins/FxnDiv/All_GOI/ALL_GOI_KOFxns_MGMs_Bins_Genus_by_Function_KO_Category_SampDate_HigherCov_heatmap_labeled.png", width=35, height=30, dpi=600,create.dir = TRUE)


#### Pull Out Antibiotic Resistance Genes from MR data ####
## heatmaps of traits of interest

mgm.bin.mr[1:4,1:4]

# pull out sulfur functions from MR transformed, summed coverages (summed coverage per KO)
ARG.bins.ko<-mgm.bin.mr[,which(colnames(mgm.bin.mr) %in% ARG.bin.fxns$KO_ID)] # merge MR data w/ ARG-related fxns found in contigs from KOFamScan
ARG.bins.ko$Bin_ID<-rownames(ARG.bins.ko)
NA %in% ARG.bins.ko
grep("Bin_ID", names(ARG.bins.ko)) # find index for column we want aka "SampleID" column

mean(as.matrix(ARG.bins.ko[,-151])) # mean KO coverage
min(as.matrix(ARG.bins.ko[,-151])) # mean KO coverage

ARG.bins.ko.melt<-melt(ARG.bins.ko, by="Bin_ID")
colnames(ARG.bins.ko.melt)[which(names(ARG.bins.ko.melt) == "variable")] <- "KO_ID"
colnames(ARG.bins.ko.melt)[which(names(ARG.bins.ko.melt) == "value")] <- "MR_SumCovPerKO"
head(ARG.bins.ko.melt) #sanity check

mr.bins.ARG.bins.ko<-merge(ARG.bins.ko.melt,ARG.kegg,by.x=c("KO_ID"),by.y=c("KO_ID"))
head(mr.bins.ARG.bins.ko)
colnames(mr.bins.ARG.bins.ko)[which(names(mr.bins.ARG.bins.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
mr.cov.sum.bins.ARG.bins.ko<-as.data.frame(dcast(mr.bins.ARG.bins.ko, Bin_ID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###
rownames(mr.cov.sum.bins.ARG.bins.ko)<-mr.cov.sum.bins.ARG.bins.ko$Bin_ID
mr.cov.sum.bins.ARG.bins.ko[1:4,1:4]

### Antibiotic Resistance Genes Heat Maps ####
# see max & mean of summed
max(mr.cov.sum.bins.ARG.bins.ko[,-1])
mean(as.matrix(mr.cov.sum.bins.ARG.bins.ko[,-1]))

# first heat map of sulfur KOs
#heatmap(as.matrix(mr.cov.sum.bins.ARG.bins.ko[,-1]), scale = "none")

colSums(mr.cov.sum.bins.ARG.bins.ko[,-1])

#heatmap(as.matrix(mr.cov.sum.bins.ARG.bins.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap
mr.bins.ARG.bins.ko[1:4,]
mr.bins.ARG.all<-merge(mr.bins.ARG.bins.ko,all_bin_meta,by="Bin_ID")
head(mr.bins.ARG.all)
mr.bins.ARG.all$PlotID = factor(mr.bins.ARG.all$Bin_ID, levels=unique(mr.bins.ARG.all$Bin_ID[order(mr.bins.ARG.all$Site,mr.bins.ARG.all$SampDate)]), ordered=TRUE)

unique(mr.bins.ARG.all$DrugGroup)
mr.bins.ARG.all$DrugGroupShort<-mr.bins.ARG.all$DrugGroup
mr.bins.ARG.all$DrugGroupShort[(mr.bins.ARG.all$DrugGroupShort) == "Macrolide Antibiotic"] <- "Macrolide"
mr.bins.ARG.all$DrugGroupShort[(mr.bins.ARG.all$DrugGroupShort) == "Narrow-spectrum Penicillin"] <- "Penicillin"
unique(mr.bins.ARG.all$DrugGroupShort) # sanity check

mr.bins.ARG.all$DrugGroup<-factor(mr.bins.ARG.all$DrugGroup,levels=c("Macrolide Antibiotic","Narrow-spectrum Penicillin","Trimethoprim","Aminoglycoside","Tetracycline"))
mr.bins.ARG.all$DrugGroupShort<-factor(mr.bins.ARG.all$DrugGroupShort,levels=c("Macrolide","Penicillin","Trimethoprim","Aminoglycoside","Tetracycline"))

head(mr.bins.ARG.all)

# order PlotID by site & sampdate
mr.bins.ARG.all$PlotID = factor(mr.bins.ARG.all$Bin_ID, levels=unique(mr.bins.ARG.all$Bin_ID[order(mr.bins.ARG.all$Site,mr.bins.ARG.all$SampDate)]), ordered=TRUE)

# convert all 0s to NAs so they appear gray in heat map
mr.bins.ARG.all$MR_SumCovPerKO_NA<-ifelse(mr.bins.ARG.all$MR_SumCovPerKO==0,NA,mr.bins.ARG.all$MR_SumCovPerKO)
head(mr.bins.ARG.all)

unique(mr.bins.ARG.all$SampDate)
mr.bins.ARG.all$SampDatePlot<-as.character(mr.bins.ARG.all$SampDate)
mr.bins.ARG.all$SampDatePlot[which(mr.bins.ARG.all$SampDatePlot == "August.2020")] <- "Aug 2020"
mr.bins.ARG.all$SampDatePlot[which(mr.bins.ARG.all$SampDatePlot == "November.2020")] <- "Nov 2020"
mr.bins.ARG.all$SampDatePlot[which(mr.bins.ARG.all$SampDatePlot == "August.2021")] <- "Aug 2021"
mr.bins.ARG.all$SampDatePlot[which(mr.bins.ARG.all$SampDatePlot == "September.2021")] <- "Sep 2021"
mr.bins.ARG.all$SampDatePlot[which(mr.bins.ARG.all$SampDatePlot == "December.2021")] <- "Dec 2021"
mr.bins.ARG.all$SampDatePlot<-gsub("\\."," ",mr.bins.ARG.all$SampDatePlot)

unique(mr.bins.ARG.all$SampDatePlot)
mr.bins.ARG.all$SampDatePlot = factor(mr.bins.ARG.all$SampDatePlot, levels=unique(mr.bins.ARG.all$SampDatePlot[order(mr.bins.ARG.all$SampDate)]), ordered=TRUE)
unique(mr.bins.ARG.all$SampDatePlot)

# merge GOI KO coverages with MAG taxonomic annotations
mr.ARG.taxa.all<-merge(mr.bins.ARG.all,mgm.bin.tax,by=c("Bin_ID","SampleID"))
head(mr.ARG.taxa.all)
mr.ARG.taxa.all$Genus = factor(mr.ARG.taxa.all$Genus, levels=unique(mr.ARG.taxa.all$Genus[order(mr.ARG.taxa.all$Site,mr.ARG.taxa.all$SampDate)]), ordered=TRUE)

# Figures below
# By Bin ID

ARG.hm1a<-ggplot(mr.bins.ARG.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Antibiotic Resistance Genes in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(ARG.hm1a,filename = "figures/MGM_Figs/Bins/FxnDiv/AntibioticResistanceGenes/ARG_KOFxns_MGMs_Bins_by_Function_heatmap.png", width=45, height=30, dpi=600,create.dir = TRUE)

ARG.hm1a2<-ggplot(mr.bins.ARG.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Antibiotic Resistance Genes in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=17),axis.text.x = element_text(angle=45, hjust=1,size=17),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 15,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(ARG.hm1a2,filename = "figures/MGM_Figs/Bins/FxnDiv/AntibioticResistanceGenes/ARG_KOFxns_MGMs_Bins_by_Function_Site_heatmap.png", width=45, height=30, dpi=600,create.dir = TRUE)

ARG.hm1a2b<-ggplot(mr.bins.ARG.all[mr.bins.ARG.all$MR_SumCovPerKO>=1,], aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Antibiotic Resistance Genes in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Normalized Coverages >= 1)",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=17),axis.text.x = element_text(angle=45, hjust=1,size=17),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 15,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(ARG.hm1a2b,filename = "figures/MGM_Figs/Bins/FxnDiv/AntibioticResistanceGenes/ARG_KOFxns_MGMs_Bins_by_Function_Site_HigherCov_heatmap.png", width=45, height=30, dpi=600,create.dir = TRUE)

ARG.hm1a2c<-ggplot(mr.bins.ARG.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=21),
        axis.text.y = element_text(size=22),axis.text.x = element_text(angle=45, hjust=1,size=19),legend.text = element_text(size=19),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),strip.text.x = element_text(size = 26,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDatePlot,scales="free_x", space = "free")

ggsave(ARG.hm1a2c,filename = "figures/MGM_Figs/Bins/FxnDiv/AntibioticResistanceGenes/ARG_KOFxns_MGMs_Bins_by_Function_SampDate_heatmap.png", width=49.5, height=11, dpi=300,create.dir = TRUE)

ARG.hm1a2d<-ggplot(mr.bins.ARG.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Antibiotic Resistance Genes in Salton Sea Dust MAGs",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=17),axis.text.x = element_text(angle=45, hjust=1,size=17),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 15,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDatePlot,scales="free_x", space = "free")

ggsave(ARG.hm1a2d,filename = "figures/MGM_Figs/Bins/FxnDiv/AntibioticResistanceGenes/ARG_KOFxns_MGMs_Bins_by_Function_SampDate_heatmap.png", width=40, height=15, dpi=300,create.dir = TRUE)

ARG.hm1a3<-ggplot(mr.bins.ARG.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Antibiotic Resistance Genes in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=17),axis.text.x = element_text(angle=45, hjust=1,size=17),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 15,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(DrugGroupShort~Site, scales="free", space = "free")

ggsave(ARG.hm1a3,filename = "figures/MGM_Figs/Bins/FxnDiv/AntibioticResistanceGenes/ARG_KOFxns_MGMs_Bins_by_Function_KO_DrugGroup_Site_heatmap.png", width=45, height=30, dpi=600,create.dir = TRUE)

ARG.hm1a3a<-ggplot(mr.bins.ARG.all[mr.bins.ARG.all$MR_SumCovPerKO>=1,], aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Antibiotic Resistance Genes in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Normalized Coverages >= 1)",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=17),axis.text.x = element_text(angle=45, hjust=1,size=17),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 15,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(DrugGroupShort~Site, scales="free", space = "free")

ggsave(ARG.hm1a3a,filename = "figures/MGM_Figs/Bins/FxnDiv/AntibioticResistanceGenes/ARG_KOFxns_MGMs_Bins_by_Function_KO_DrugGroup_Site_HigherCov_heatmap.png", width=45, height=30, dpi=600,create.dir = TRUE)

ARG.hm1a3b<-ggplot(mr.bins.ARG.all[mr.bins.ARG.all$MR_SumCovPerKO>=1,], aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) + geom_text(label=round(mr.bins.ARG.all$MR_SumCovPerKO[mr.bins.ARG.all$MR_SumCovPerKO>=1],2)) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Antibiotic Resistance Genes in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Normalized Coverages >= 1)",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=17),axis.text.x = element_text(angle=45, hjust=1,size=17),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 15,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(DrugGroupShort~Site, scales="free", space = "free")

ggsave(ARG.hm1a3b,filename = "figures/MGM_Figs/Bins/FxnDiv/AntibioticResistanceGenes/ARG_KOFxns_MGMs_Bins_by_Function_KO_DrugGroup_Site_HigherCov_heatmap_labeled.png", width=45, height=30, dpi=600,create.dir = TRUE)

ARG.hm1a4<-ggplot(mr.bins.ARG.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Antibiotic Resistance Genes in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=17),axis.text.x = element_text(angle=45, hjust=1,size=17),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 15,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(DrugGroupShort~SampDate_Plot, scales="free", space = "free")

ggsave(ARG.hm1a4,filename = "figures/MGM_Figs/Bins/FxnDiv/AntibioticResistanceGenes/ARG_KOFxns_MGMs_Bins_by_Function_KO_DrugGroup_SampDate_heatmap.png", width=45, height=30, dpi=600,create.dir = TRUE)

ARG.hm1a4a<-ggplot(mr.bins.ARG.all[mr.bins.ARG.all$MR_SumCovPerKO>=1,], aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Antibiotic Resistance Genes in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Normalized Coverages >= 1)",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=17),axis.text.x = element_text(angle=45, hjust=1,size=17),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 15,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(DrugGroupShort~SampDate_Plot, scales="free", space = "free")

ggsave(ARG.hm1a4a,filename = "figures/MGM_Figs/Bins/FxnDiv/AntibioticResistanceGenes/ARG_KOFxns_MGMs_Bins_by_Function_KO_DrugGroup_SampDate_HigherCov_heatmap.png", width=45, height=30, dpi=600,create.dir = TRUE)

# ARG.hm1a4b<-ggplot(mr.bins.ARG.all[mr.bins.ARG.all$MR_SumCovPerKO>=1,], aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.15) + geom_text(label=round(mr.bins.ARG.all$MR_SumCovPerKO[mr.bins.ARG.all$MR_SumCovPerKO>=1],2)) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Antibiotic Resistance Genes in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Normalized Coverages >= 1)",fill="MR Coverage Per KO") +
#   theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(DrugGroupShort~SampDate_Plot, scales="free", space = "free")
#
# ggsave(ARG.hm1a4b,filename = "figures/MGM_Figs/Bins/FxnDiv/AntibioticResistanceGenes/ARG_KOFxns_MGMs_Bins_by_Function_KO_DrugGroup_SampDate_HigherCov_heatmap_labeled.png", width=45, height=30, dpi=600,create.dir = TRUE)

# by Taxa ID

ARG.hm1b<-ggplot(mr.ARG.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Antibiotic Resistance Genes in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=17),axis.text.x = element_text(angle=45, hjust=1,size=17),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 15,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(ARG.hm1b,filename = "figures/MGM_Figs/Bins/FxnDiv/AntibioticResistanceGenes/ARG_KOFxns_MGMs_Bins_by_Function_heatmap.png", width=30, height=30, dpi=600,create.dir = TRUE)

ARG.hm1b2<-ggplot(mr.ARG.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Antibiotic Resistance Genes in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=17),axis.text.x = element_text(angle=45, hjust=1,size=17),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 15,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(ARG.hm1b2,filename = "figures/MGM_Figs/Bins/FxnDiv/AntibioticResistanceGenes/ARG_KOFxns_MGMs_Bins_by_Function_Site_heatmap.png", width=30, height=30, dpi=600,create.dir = TRUE)

ARG.hm1b2b<-ggplot(mr.ARG.taxa.all[mr.ARG.taxa.all$MR_SumCovPerKO>=1,], aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Antibiotic Resistance Genes in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Normalized Coverages >= 1)",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=17),axis.text.x = element_text(angle=45, hjust=1,size=17),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 15,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(ARG.hm1b2b,filename = "figures/MGM_Figs/Bins/FxnDiv/AntibioticResistanceGenes/ARG_KOFxns_MGMs_Bins_by_Function_Site_HigherCov_heatmap.png", width=30, height=30, dpi=600,create.dir = TRUE)

ARG.hm1b2c<-ggplot(mr.ARG.taxa.all[mr.ARG.taxa.all$MR_SumCovPerKO>=1,], aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Antibiotic Resistance Genes in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Normalized Coverages >= 1)",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=17),axis.text.x = element_text(angle=45, hjust=1,size=17),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 15,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate_Plot,scales="free_x", space = "free")

ggsave(ARG.hm1b2c,filename = "figures/MGM_Figs/Bins/FxnDiv/AntibioticResistanceGenes/ARG_KOFxns_MGMs_Bins_by_Function_SampDate_HigherCov_heatmap.png", width=30, height=30, dpi=600,create.dir = TRUE)

ARG.hm1b2d<-ggplot(mr.ARG.taxa.all[mr.ARG.taxa.all$MR_SumCovPerKO>=1,], aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) + geom_text(label=round(mr.ARG.taxa.all$MR_SumCovPerKO[mr.ARG.taxa.all$MR_SumCovPerKO>=1],2)) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Antibiotic Resistance Genes in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Normalized Coverages >= 1)",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=17),axis.text.x = element_text(angle=45, hjust=1,size=17),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 15,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate_Plot,scales="free_x", space = "free")

ggsave(ARG.hm1b2d,filename = "figures/MGM_Figs/Bins/FxnDiv/AntibioticResistanceGenes/ARG_KOFxns_MGMs_Bins_by_Function_SampDate_HigherCov_heatmap_labeled.png", width=35, height=30, dpi=600,create.dir = TRUE)

ARG.hm1b3<-ggplot(mr.ARG.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Antibiotic Resistance Genes in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=17),axis.text.x = element_text(angle=45, hjust=1,size=17),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 15,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(DrugGroupShort~Site, scales="free", space = "free")

ggsave(ARG.hm1b3,filename = "figures/MGM_Figs/Bins/FxnDiv/AntibioticResistanceGenes/ARG_KOFxns_MGMs_Bins_Genus_by_Function_KO_DrugGroup_Site_heatmap.png", width=35, height=30, dpi=600,create.dir = TRUE)

ARG.hm1b3a<-ggplot(mr.ARG.taxa.all[mr.ARG.taxa.all$MR_SumCovPerKO>=1,], aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Antibiotic Resistance Genes in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Normalized Coverages >= 1)",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=17),axis.text.x = element_text(angle=45, hjust=1,size=17),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 15,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(DrugGroupShort~Site, scales="free", space = "free")

ggsave(ARG.hm1b3a,filename = "figures/MGM_Figs/Bins/FxnDiv/AntibioticResistanceGenes/ARG_KOFxns_MGMs_Bins_Genus_by_Function_KO_DrugGroup_Site_HigherCov_heatmap.png", width=35, height=30, dpi=600,create.dir = TRUE)

ARG.hm1b3b<-ggplot(mr.ARG.taxa.all[mr.ARG.taxa.all$MR_SumCovPerKO>=1,], aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) + geom_text(label=round(mr.ARG.taxa.all$MR_SumCovPerKO[mr.ARG.taxa.all$MR_SumCovPerKO>=1],2)) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Antibiotic Resistance Genes in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Normalized Coverages >= 1)",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=17),axis.text.x = element_text(angle=45, hjust=1,size=17),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 15,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(DrugGroupShort~Site, scales="free", space = "free")

ggsave(ARG.hm1b3b,filename = "figures/MGM_Figs/Bins/FxnDiv/AntibioticResistanceGenes/ARG_KOFxns_MGMs_Bins_Genus_by_Function_KO_DrugGroup_Site_HigherCov_heatmap_labeled.png", width=35, height=30, dpi=600,create.dir = TRUE)

ARG.hm1b4<-ggplot(mr.ARG.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Antibiotic Resistance Genes in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=17),axis.text.x = element_text(angle=45, hjust=1,size=17),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 15,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(DrugGroupShort~SampDate_Plot, scales="free", space = "free")

ggsave(ARG.hm1b4,filename = "figures/MGM_Figs/Bins/FxnDiv/AntibioticResistanceGenes/ARG_KOFxns_MGMs_Bins_Genus_by_Function_KO_DrugGroup_SampDate_heatmap.png", width=35, height=30, dpi=600,create.dir = TRUE)

ARG.hm1b4a<-ggplot(mr.ARG.taxa.all[mr.ARG.taxa.all$MR_SumCovPerKO>=1,], aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Antibiotic Resistance Genes in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Normalized Coverages >= 1)",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=17),axis.text.x = element_text(angle=45, hjust=1,size=17),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 15,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(DrugGroupShort~SampDate_Plot, scales="free", space = "free")

ggsave(ARG.hm1b4a,filename = "figures/MGM_Figs/Bins/FxnDiv/AntibioticResistanceGenes/ARG_KOFxns_MGMs_Bins_Genus_by_Function_KO_DrugGroup_SampDate_HigherCov_heatmap.png", width=35, height=30, dpi=600,create.dir = TRUE)

ARG.hm1b4b<-ggplot(mr.ARG.taxa.all[mr.ARG.taxa.all$MR_SumCovPerKO>=1,], aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) + geom_text(label=round(mr.ARG.taxa.all$MR_SumCovPerKO[mr.ARG.taxa.all$MR_SumCovPerKO>=1],2)) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Antibiotic Resistance Genes in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Normalized Coverages >= 1)",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=17),axis.text.x = element_text(angle=45, hjust=1,size=17),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 15,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(DrugGroupShort~SampDate_Plot, scales="free", space = "free")

ggsave(ARG.hm1b4b,filename = "figures/MGM_Figs/Bins/FxnDiv/AntibioticResistanceGenes/ARG_KOFxns_MGMs_Bins_Genus_by_Function_KO_DrugGroup_SampDate_HigherCov_heatmap_labeled.png", width=35, height=30, dpi=600,create.dir = TRUE)


#### Pull Out Arsenic Metabolic Fxns from MR data ####
ars.bins.ko<-mgm.bin.mr[,which(colnames(mgm.bin.mr) %in% arsen.bin.fxns$KO_ID)]
ars.bins.ko$Bin_ID<-rownames(ars.bins.ko)
ars.bins.ko.melt<-melt(ars.bins.ko, by="Bin_ID")
colnames(ars.bins.ko.melt)[which(names(ars.bins.ko.melt) == "variable")] <- "KO_ID"
colnames(ars.bins.ko.melt)[which(names(ars.bins.ko.melt) == "value")] <- "MR_SumCovPerKO"
ars.bins.ko.melt #sanity check

mr.ars.bins.ko<-merge(ars.bins.ko.melt,arsen.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(mr.ars.bins.ko)
colnames(mr.ars.bins.ko)[which(names(mr.ars.bins.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
mr.cov.sum.ars.bins.ko<-as.data.frame(dcast(mr.ars.bins.ko, Bin_ID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(mr.cov.sum.ars.bins.ko)<-mr.cov.sum.ars.bins.ko$Bin_ID
mr.cov.sum.ars.bins.ko[1:4,]

#### Arsenic Metabolism KO Heat Maps ####
# see max & mean of summed
max(mr.cov.sum.ars.bins.ko[,-1])
min(mr.cov.sum.ars.bins.ko[,-1])
max(mr.cov.sum.ars.bins.ko[,-1])/2

# first heat map of ars KOs
##heatmap(as.matrix(mr.cov.sum.ars.bins.ko[,-1]), scale = "none")

colSums(mr.cov.sum.ars.bins.ko[,-1])
#mr.cov.sum.ars.bins.ko2 <- mr.cov.sum.ars.bins.ko[,which(colSums(mr.cov.sum.ars.bins.ko[,-1])>10)]

##heatmap(as.matrix(mr.cov.sum.ars.bins.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap -- using MR data that includes NAs so they are blocked out on heatmaps
mr.ars.bins.ko[1:4,]
mr.ars.bins.all<-merge(mr.ars.bins.ko,all_bin_meta,by="Bin_ID")
head(mr.ars.bins.all)
#mr.ars.bins.all$PlotID<-mr.ars.bins.all$Bin_ID
mr.ars.bins.all$PlotID = factor(mr.ars.bins.all$Bin_ID, levels=unique(mr.ars.bins.all$Bin_ID[order(mr.ars.bins.all$Site,mr.ars.bins.all$SampDate)]), ordered=TRUE)

#mr.ars.bins.all$PathShort<-mr.ars.bins.all$Pathway
# mr.ars.bins.all$PathShort[(mr.ars.bins.all$PathShort) == "Dissimilatory Sulfate Redox"] <- "D.SO4 RedOx"

head(mr.ars.bins.all)

# convert all 0s to NAs so they appear gray in heat map
mr.ars.bins.all$MR_SumCovPerKO_NA<-ifelse(mr.ars.bins.all$MR_SumCovPerKO==0,NA,mr.ars.bins.all$MR_SumCovPerKO)
head(mr.ars.bins.all)

# For heatmap color gradient
max(mr.ars.bins.all$MR_SumCovPerKO, na.rm=TRUE)
max(mr.ars.bins.all$MR_SumCovPerKO, na.rm=TRUE)/2
min(mr.ars.bins.all$MR_SumCovPerKO, na.rm=TRUE)

# merge arsic metabolism coverages with MAG taxonomic annotations
mr.ars.taxa.all<-merge(mr.ars.bins.all,mgm.bin.tax,by=c("Bin_ID","SampleID"))
head(mr.ars.taxa.all)
mr.ars.taxa.all$Genus = factor(mr.ars.taxa.all$Genus, levels=unique(mr.ars.taxa.all$Genus[order(mr.ars.taxa.all$Site,mr.ars.taxa.all$SampDate)]), ordered=TRUE)

# Figures below
# by Bin_ID

ars.hm1a<-ggplot(mr.ars.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Arsenic Metabolism in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(ars.hm1a,filename = "figures/MGM_Figs/Bins/FxnDiv/Arsenic/Arsenic_KOFxns_MGMs_Bins_BinID_by_Function_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

ars.hm1a2<-ggplot(mr.ars.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Arsenic Metabolism in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(ars.hm1a2,filename = "figures/MGM_Figs/Bins/FxnDiv/Arsenic/Arsenic_KOFxns_MGMs_Bins_BinID_by_Function_Site_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)

ars.hm1a3<-ggplot(mr.ars.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Arsenic Metabolism in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate,scales="free_x", space = "free")

ggsave(ars.hm1a3,filename = "figures/MGM_Figs/Bins/FxnDiv/Arsenic/Arsenic_KOFxns_MGMs_Bins_by_BinID_by_Function_SampDate_heatmap.png", width=49.5, height=13, dpi=600,create.dir = TRUE)

ars.hm1a4<-ggplot(mr.ars.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Arsenic Metabolism in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate+Site,scales="free_x", space = "free")

ggsave(ars.hm1a4,filename = "figures/MGM_Figs/Bins/FxnDiv/Arsenic/Arsenic_KOFxns_MGMs_Bins_by_BinID_by_Function_SampDate_Site_heatmap.png", width=49.5, height=13, dpi=600,create.dir = TRUE)

## by MAG Taxonomic ID, not Bin ID

ars.hm1b<-ggplot(mr.ars.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Arsenic Metabolism in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(ars.hm1b,filename = "figures/MGM_Figs/Bins/FxnDiv/Arsenic/Arsenic_KOFxns_MGMs_Bins_Genus_by_Function_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

ars.hm1b2<-ggplot(mr.ars.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Arsenic Metabolism in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(ars.hm1b2,filename = "figures/MGM_Figs/Bins/FxnDiv/Arsenic/Arsenic_KOFxns_MGMs_Bins_Genus_by_Function_Site_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)

ars.hm1b3<-ggplot(mr.ars.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Arsenic Metabolism in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate,scales="free_x", space = "free")

ggsave(ars.hm1b3,filename = "figures/MGM_Figs/Bins/FxnDiv/Arsenic/Arsenic_KOFxns_MGMs_Bins_Genus_by_Function_SampDate_heatmap.png", width=49.5, height=15, dpi=600,create.dir = TRUE)

ars.hm1b4<-ggplot(mr.ars.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Arsenic Metabolism in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate+Site,scales="free_x", space = "free")

ggsave(ars.hm1b4,filename = "figures/MGM_Figs/Bins/FxnDiv/Arsenic/Arsenic_KOFxns_MGMs_Bins_Genus_by_Function_SampDate_Site_heatmap.png", width=49.5, height=15, dpi=600,create.dir = TRUE)

#### Pull Out Selenium Metabolism Fxns from MR data ####
sel.bins.ko<-mgm.bin.mr[,which(colnames(mgm.bin.mr) %in% selen.bin.fxns$KO_ID)]
sel.bins.ko$Bin_ID<-rownames(sel.bins.ko)
sel.bins.ko.melt<-melt(sel.bins.ko, by="Bin_ID")
colnames(sel.bins.ko.melt)[which(names(sel.bins.ko.melt) == "variable")] <- "KO_ID"
colnames(sel.bins.ko.melt)[which(names(sel.bins.ko.melt) == "value")] <- "MR_SumCovPerKO"
sel.bins.ko.melt #sanity check

mr.selen.bins.ko<-merge(sel.bins.ko.melt,selen.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(mr.selen.bins.ko)
colnames(mr.selen.bins.ko)[which(names(mr.selen.bins.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
mr.cov.sum.selen.bins.ko<-as.data.frame(dcast(mr.selen.bins.ko, Bin_ID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(mr.cov.sum.selen.bins.ko)<-mr.cov.sum.selen.bins.ko$Bin_ID
mr.cov.sum.selen.bins.ko[1:4,]

#### Selenium Metabolism KO Heat Maps ####
# see max & mean of summed
max(mr.cov.sum.selen.bins.ko[,-1])
min(mr.cov.sum.selen.bins.ko[,-1])
max(mr.cov.sum.selen.bins.ko[,-1])/2

# first heat map of selen KOs
##heatmap(as.matrix(mr.cov.sum.selen.bins.ko[,-1]), scale = "none")

colSums(mr.cov.sum.selen.bins.ko[,-1])
#mr.cov.sum.selen.bins.ko2 <- mr.cov.sum.selen.bins.ko[,which(colSums(mr.cov.sum.selen.bins.ko[,-1])>10)]

##heatmap(as.matrix(mr.cov.sum.selen.bins.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap -- using MR data that includes NAs so they are blocked out on heatmaps
mr.selen.bins.ko[1:4,]
mr.selen.bins.all<-merge(mr.selen.bins.ko,all_bin_meta,by="Bin_ID")
head(mr.selen.bins.all)
#mr.selen.bins.all$PlotID<-mr.selen.bins.all$Bin_ID
mr.selen.bins.all$PlotID = factor(mr.selen.bins.all$Bin_ID, levels=unique(mr.selen.bins.all$Bin_ID[order(mr.selen.bins.all$Site,mr.selen.bins.all$SampDate)]), ordered=TRUE)

#mr.selen.bins.all$PathShort<-mr.selen.bins.all$Pathway
# mr.selen.bins.all$PathShort[(mr.selen.bins.all$PathShort) == "Dissimilatory Sulfate Redox"] <- "D.SO4 RedOx"

head(mr.selen.bins.all)

# convert all 0s to NAs so they appear gray in heat map
mr.selen.bins.all$MR_SumCovPerKO_NA<-ifelse(mr.selen.bins.all$MR_SumCovPerKO==0,NA,mr.selen.bins.all$MR_SumCovPerKO)
head(mr.selen.bins.all)

# For heatmap color gradient
max(mr.selen.bins.all$MR_SumCovPerKO, na.rm=TRUE)
max(mr.selen.bins.all$MR_SumCovPerKO, na.rm=TRUE)/2
min(mr.selen.bins.all$MR_SumCovPerKO, na.rm=TRUE)

# merge selenium KO coverages with MAG taxonomic annotations
mr.selen.taxa.all<-merge(mr.selen.bins.all,mgm.bin.tax,by=c("Bin_ID","SampleID"))
head(mr.selen.taxa.all)
mr.selen.taxa.all$Genus = factor(mr.selen.taxa.all$Genus, levels=unique(mr.selen.taxa.all$Genus[order(mr.selen.taxa.all$Site,mr.selen.taxa.all$SampDate)]), ordered=TRUE)

# Figures below
# by Bin_ID

selen.hm1a<-ggplot(mr.selen.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Selenium Metabolism in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(selen.hm1a,filename = "figures/MGM_Figs/Bins/FxnDiv/Selenium/Selenium_KOFxns_MGMs_Bins_BinID_by_Function_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

selen.hm1a2<-ggplot(mr.selen.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Selenium Metabolism in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(selen.hm1a2,filename = "figures/MGM_Figs/Bins/FxnDiv/Selenium/Selenium_KOFxns_MGMs_Bins_BinID_by_Function_Site_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)

selen.hm1a3<-ggplot(mr.selen.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Selenium Metabolism in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate,scales="free_x", space = "free")

ggsave(selen.hm1a3,filename = "figures/MGM_Figs/Bins/FxnDiv/Selenium/Selenium_KOFxns_MGMs_Bins_by_BinID_by_Function_SampDate_heatmap.png", width=49.5, height=13, dpi=600,create.dir = TRUE)

selen.hm1a4<-ggplot(mr.selen.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Selenium Metabolism in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate+Site,scales="free_x", space = "free")

ggsave(selen.hm1a4,filename = "figures/MGM_Figs/Bins/FxnDiv/Selenium/Selenium_KOFxns_MGMs_Bins_by_BinID_by_Function_SampDate_Site_heatmap.png", width=49.5, height=13, dpi=600,create.dir = TRUE)

## by MAG Taxonomic ID, not Bin ID

selen.hm1b<-ggplot(mr.selen.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Selenium Metabolism in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(selen.hm1b,filename = "figures/MGM_Figs/Bins/FxnDiv/Selenium/Selenium_KOFxns_MGMs_Bins_Genus_by_Function_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

selen.hm1b2<-ggplot(mr.selen.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Selenium Metabolism in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(selen.hm1b2,filename = "figures/MGM_Figs/Bins/FxnDiv/Selenium/Selenium_KOFxns_MGMs_Bins_Genus_by_Function_Site_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)

selen.hm1b3<-ggplot(mr.selen.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Selenium Metabolism in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate,scales="free_x", space = "free")

ggsave(selen.hm1b3,filename = "figures/MGM_Figs/Bins/FxnDiv/Selenium/Selenium_KOFxns_MGMs_Bins_Genus_by_Function_SampDate_heatmap.png", width=49.5, height=15, dpi=600,create.dir = TRUE)

selen.hm1b4<-ggplot(mr.selen.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Selenium Metabolism in Salton Sea Dust MAGs",subtitle="Using Median-Ratio Normalized, Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate+Site,scales="free_x", space = "free")

ggsave(selen.hm1b4,filename = "figures/MGM_Figs/Bins/FxnDiv/Selenium/Selenium_KOFxns_MGMs_Bins_Genus_by_Function_SampDate_Site_heatmap.png", width=49.5, height=15, dpi=600,create.dir = TRUE)

#### Pull Out Phototrophy Fxns from MR data ####
## heatmaps of traits of interest

mgm.bin.mr[1:4,1:4]

# pull out Phototrophy functions from MR transformed, summed coverages (summed coverage per KO)
photo.bins.ko<-mgm.bin.mr[,which(colnames(mgm.bin.mr) %in% photo.bin.fxns$KO_ID)] # merge MR data w/ photoon-related fxns found in contigs from KOFamScan
photo.bins.ko$Bin_ID<-rownames(photo.bins.ko)
photo.bins.ko.melt<-melt(photo.bins.ko, by="Bin_ID")
colnames(photo.bins.ko.melt)[which(names(photo.bins.ko.melt) == "variable")] <- "KO_ID"
colnames(photo.bins.ko.melt)[which(names(photo.bins.ko.melt) == "value")] <- "MR_SumCovPerKO"
head(photo.bins.ko.melt) #sanity check

mr.photo.bins.ko<-merge(photo.bins.ko.melt,photo.kegg,by.x=c("KO_ID"),by.y=c("KO_ID"))
head(mr.photo.bins.ko)
colnames(mr.photo.bins.ko)[which(names(mr.photo.bins.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
mr.cov.sum.photo.bins.ko<-as.data.frame(dcast(mr.photo.bins.ko, Bin_ID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###
rownames(mr.cov.sum.photo.bins.ko)<-mr.cov.sum.photo.bins.ko$Bin_ID
mr.cov.sum.photo.bins.ko[1:4,]

### Phototrophy Heat Maps ####
# see max & mean of summed
max(mr.cov.sum.photo.bins.ko[,-1])
mean(as.matrix(mr.cov.sum.photo.bins.ko[,-1]))

# first heat map of sulfur KOs
#heatmap(as.matrix(mr.cov.sum.photo.bins.ko[,-1]), scale = "none")

colSums(mr.cov.sum.photo.bins.ko[,-1])

#heatmap(as.matrix(mr.cov.sum.photo.bins.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap
mr.photo.bins.ko[1:4,]
mr.photo.bins.all<-merge(mr.photo.bins.ko,all_bin_meta,by="Bin_ID")
head(mr.photo.bins.all)
mr.photo.bins.all$PlotID = factor(mr.photo.bins.all$Bin_ID, levels=unique(mr.photo.bins.all$Bin_ID[order(mr.photo.bins.all$SampDate,mr.photo.bins.all$Site)]), ordered=TRUE)
mr.photo.bins.all$SampDate<-gsub("\\."," ",mr.photo.bins.all$SampDate)
mr.photo.bins.all$SampDate<-factor(mr.photo.bins.all$SampDate, levels=c("August 2021","December 2021","April 2022"))

unique(mr.photo.bins.all$Pathway)
mr.photo.bins.all$PathShort<-mr.photo.bins.all$Pathway
mr.photo.bins.all$PathShort[(mr.photo.bins.all$PathShort) == "Proteorhodopsin"] <- "PR"
#mr.photo.bins.all$PathShort[(mr.photo.bins.all$PathShort) == "Bacteriorhodopsin"] <- "BR"
mr.photo.bins.all$PathShort[(mr.photo.bins.all$PathShort) == "Sensory Rhodopsin"] <- "SR"
#mr.photo.bins.all$PathShort[(mr.photo.bins.all$PathShort) == "Halorhodopsin"] <- "HR"
mr.photo.bins.all$PathShort[(mr.photo.bins.all$PathShort) == "Photosystem II"] <- "PS II"
mr.photo.bins.all$PathShort[(mr.photo.bins.all$PathShort) == "Photosystem I"] <- "PS I"
mr.photo.bins.all$PathShort[(mr.photo.bins.all$PathShort) == "Anoxygenic Photosystem II"] <- "AnOx PS"

mr.photo.bins.all$Pathway<-factor(mr.photo.bins.all$Pathway,levels=c("Proteorhodopsin","Sensory Rhodopsin","Photosystem II","Photosystem I","Anoxygenic Photosystem II"))
mr.photo.bins.all$PathShort<-factor(mr.photo.bins.all$PathShort,levels=c("PR","SR","PS II","PS I","AnOx PS"))

mr.photo.bins.all$MethShort<-mr.photo.bins.all$Method
mr.photo.bins.all$MethShort[(mr.photo.bins.all$MethShort) == "Bacterial Rhodopsin"] <- "Bac Rhod"
mr.photo.bins.all$MethShort[(mr.photo.bins.all$MethShort) == "Anoxygenic Photosynthesis"] <- "AnOx PS"
mr.photo.bins.all$MethShort[(mr.photo.bins.all$MethShort) == "Oxygenic Photosynthesis"] <- "Ox PS"

mr.photo.bins.all$Method<-factor(mr.photo.bins.all$Method,levels=c("Bacterial Rhodopsin","Oxygenic Photosynthesis","Anoxygenic Photosynthesis"))
mr.photo.bins.all$MethShort<-factor(mr.photo.bins.all$MethShort,levels=c("Bac Rhod","Ox PS","AnOx PS"))

unique(mr.photo.bins.all$Phototrophy)
mr.photo.bins.all$Phototrophy<-factor(mr.photo.bins.all$Phototrophy,levels=c("Hetero","Auto"))

mr.photo.bins.all$KO_Function.KEGG = factor(mr.photo.bins.all$KO_Function.KEGG, levels=unique(mr.photo.bins.all$KO_Function.KEGG[order(mr.photo.bins.all$Phototrophy)]), ordered=TRUE)

head(mr.photo.bins.all)

# convert all 0s to NAs so they appear gray in heat map
mr.photo.bins.all$MR_SumCovPerKO_NA<-ifelse(mr.photo.bins.all$MR_SumCovPerKO==0,NA,mr.photo.bins.all$MR_SumCovPerKO)
head(mr.photo.bins.all)

# For heatmap color gradient
max(mr.photo.bins.all$MR_SumCovPerKO, na.rm=TRUE)
median(mr.photo.bins.all$MR_SumCovPerKO, na.rm=TRUE)
min(mr.photo.bins.all$MR_SumCovPerKO, na.rm=TRUE)

# merge phototrophy KO coverages with MAG taxonomic annotations
mr.photo.taxa.all<-merge(mr.photo.bins.all,mgm.bin.tax,by=c("Bin_ID","SampleID"))
head(mr.photo.taxa.all)
mr.photo.taxa.all$Genus = factor(mr.photo.taxa.all$Genus, levels=unique(mr.photo.taxa.all$Genus[order(mr.photo.taxa.all$Site,mr.photo.taxa.all$SampDate)]), ordered=TRUE)

# Figures below
# by Bin_ID

photo.hm1a<-ggplot(mr.photo.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Phototrophy Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(photo.hm1a,filename = "figures/MGM_Figs/Bins/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_Bins_BinID_by_Function_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

photo.hm1b<-ggplot(mr.photo.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Phototrophy Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Phototrophy~.,scales="free_y", space = "free")

ggsave(photo.hm1b,filename = "figures/MGM_Figs/Bins/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_Bins_BinID_by_Function_Phototrophy_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)

photo.hm1c<-ggplot(mr.photo.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Phototrophy Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate+Site,scales="free_x", space = "free")

ggsave(photo.hm1c,filename = "figures/MGM_Figs/Bins/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_Bins_by_BinID_by_Function_SampDate_Site_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

photo.hm1d<-ggplot(mr.photo.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Phototrophy Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate,scales="free_y", space = "free")

ggsave(photo.hm1d,filename = "figures/MGM_Figs/Bins/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_Bins_BinID_by_Function_Phototrophy_System_SampDate_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)

photo.hm1e<-ggplot(mr.photo.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Phototrophy Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate+Site, scales="free", space = "free")

ggsave(photo.hm1e,filename = "figures/MGM_Figs/Bins/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_Bins_by_BinID_by_Function_Phototrophy_SampDate_Site_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)

photo.hm1f<-ggplot(mr.photo.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Phototrophy Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate+Site, scales="free", space = "free")

ggsave(photo.hm1f,filename = "figures/MGM_Figs/Bins/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_Bins_by_BinID_by_Function_Phototrophy_SampDate_Site_System_heatmap2.png", width=45, height=15, dpi=600,create.dir = TRUE)

photo.hm1g<-ggplot(mr.photo.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Phototrophy Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(MethShort~SampDate+Site, scales="free", space = "free")

ggsave(photo.hm1g,filename = "figures/MGM_Figs/Bins/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_Bins_by_BinID_by_Function_Phototrophy_SampDate_Site_Method_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)

## by MAG Taxonomic ID, not Bin ID

# photo.hm1b<-ggplot(mr.photo.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Phototrophy Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
#   theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))
#
# ggsave(photo.hm1b,filename = "figures/MGM_Figs/Bins/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_Bins_Genus_by_Function_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

photo.hm1b2<-ggplot(mr.photo.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Phototrophy Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Phototrophy~.,scales="free_y", space = "free")

ggsave(photo.hm1b2,filename = "figures/MGM_Figs/Bins/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_Bins_Genus_by_Function_Phototrophy_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)

photo.hm1c2<-ggplot(mr.photo.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Phototrophy Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate+Site,scales="free_x", space = "free")

ggsave(photo.hm1c2,filename = "figures/MGM_Figs/Bins/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_Bins_Genus_by_Function_SampDate_Site_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

photo.hm1d2<-ggplot(mr.photo.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Phototrophy Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate,scales="free_y", space = "free")

ggsave(photo.hm1d2,filename = "figures/MGM_Figs/Bins/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_Bins_Genus_by_Function_Phototrophy_System_SampDate_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)

photo.hm1e2<-ggplot(mr.photo.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Phototrophy Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Phototrophy~SampDate+Site, scales="free", space = "free")

ggsave(photo.hm1e2,filename = "figures/MGM_Figs/Bins/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_Bins_Genus_by_Function_Phototrophy_SampDate_Site_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)

photo.hm1f2<-ggplot(mr.photo.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Phototrophy Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate+Site, scales="free", space = "free")

ggsave(photo.hm1f2,filename = "figures/MGM_Figs/Bins/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_Bins_Genus_by_Function_Phototrophy_SampDate_Site_System_heatmap2.png", width=45, height=15, dpi=600,create.dir = TRUE)

photo.hm1e2<-ggplot(mr.photo.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Phototrophy Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(MethShort~SampDate+Site, scales="free", space = "free")

ggsave(photo.hm1e2,filename = "figures/MGM_Figs/Bins/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_Bins_Genus_by_Function_Phototrophy_SampDate_Site_Method_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)


#### Pull Out Sulfur Metabolic Fxns from MR data ####
## heatmaps of traits of interest

mgm.bin.mr[1:4,1:4]

# pull out sulfur functions from MR transformed, summed coverages (summed gene coverage per KO)
sulf.bins.ko<-mgm.bin.mr[,which(colnames(mgm.bin.mr) %in% sulfur.bin.fxns$KO_ID)] # merge MR data w/ S fxns found in contigs from KOFamScan
sulf.bins.ko$Bin_ID<-rownames(sulf.bins.ko)
sulf.bins.ko.melt<-melt(sulf.bins.ko, by="Bin_ID")
colnames(sulf.bins.ko.melt)[which(names(sulf.bins.ko.melt) == "variable")] <- "KO_ID"
colnames(sulf.bins.ko.melt)[which(names(sulf.bins.ko.melt) == "value")] <- "MR_SumCovPerKO"
head(sulf.bins.ko.melt) #sanity check

mr.sulf.bins.ko<-merge(sulf.bins.ko.melt,sulf.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(mr.sulf.bins.ko)
colnames(mr.sulf.bins.ko)[which(names(mr.sulf.bins.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
mr.cov.sum.sulf.bins.ko<-as.data.frame(dcast(mr.sulf.bins.ko, Bin_ID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(mr.cov.sum.sulf.bins.ko)<-mr.cov.sum.sulf.bins.ko$Bin_ID
mr.cov.sum.sulf.bins.ko[1:4,]

# sanity check
mr.cov.sum.sulf.bins.ko$`cysJ; sulfite reductase (NADPH) flavoprotein alpha-component [EC:1.8.1.2]`[1:4]
head(mr.sulf.bins.ko)

#### Sulfur Heat Maps ####
# see max & mean of summed
max(mr.cov.sum.sulf.bins.ko[,-1])
min(mr.cov.sum.sulf.bins.ko[,-1])
max(mr.cov.sum.sulf.bins.ko[,-1])/2

# first heat map of sulfur KOs
#heatmap(as.matrix(mr.cov.sum.sulf.bins.ko[,-1]), scale = "none")

colSums(mr.cov.sum.sulf.bins.ko[,-1])
#mr.cov.sum.sulf.bins.ko2 <- mr.cov.sum.sulf.bins.ko[,which(colSums(mr.cov.sum.sulf.bins.ko[,-1])>10)]

#heatmap(as.matrix(mr.cov.sum.sulf.bins.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap -- using MR data that includes NAs so they are blocked out on heatmaps
mr.sulf.bins.ko[1:4,]
mr.sulf.bins.all<-merge(mr.sulf.bins.ko,all_bin_meta,by="Bin_ID")
head(mr.sulf.bins.all)
mr.sulf.bins.all$PlotID = factor(mr.sulf.bins.all$Bin_ID, levels=unique(mr.sulf.bins.all$Bin_ID[order(mr.sulf.bins.all$Site,mr.sulf.bins.all$SampDate)]), ordered=TRUE)
# mr.sulf.bins.all$SampDate<-gsub("\\."," ",mr.sulf.bins.all$SampDate)

mr.sulf.bins.all$PathShort<-mr.sulf.bins.all$Pathway
mr.sulf.bins.all$PathShort[(mr.sulf.bins.all$PathShort) == "Dissimilatory Sulfate Redox"] <- "D.SO4 RedOx"
mr.sulf.bins.all$PathShort[(mr.sulf.bins.all$PathShort) == "Assimilatory Sulfate Reduction"] <- "A.SO4 Red"
mr.sulf.bins.all$PathShort[(mr.sulf.bins.all$PathShort) == "Multiple Pathways"] <- "Multi Paths"
mr.sulf.bins.all$PathShort[(mr.sulf.bins.all$PathShort) == "S Disproportionation"] <- "S Disprop."

mr.sulf.bins.all$Pathway<-factor(mr.sulf.bins.all$Pathway,levels=c("Assimilatory Sulfate Reduction","Dissimilatory Sulfate Redox","Multiple Pathways","SOX","S Disproportionation"))
mr.sulf.bins.all$PathShort<-factor(mr.sulf.bins.all$PathShort,levels=c("A.SO4 Red","D.SO4 RedOx","Multi Paths","SOX","S Disprop."))

mr.sulf.bins.all$PathSpecShort<-mr.sulf.bins.all$PathwaySpecific
mr.sulf.bins.all$PathSpecShort[(mr.sulf.bins.all$PathSpecShort) == "Dissimilatory Sulfate Redox"] <- "D.SO4 RedOx"
mr.sulf.bins.all$PathSpecShort[(mr.sulf.bins.all$PathSpecShort) == "Assimilatory Sulfate Reduction"] <- "A.SO4 Red"
mr.sulf.bins.all$PathSpecShort[(mr.sulf.bins.all$PathSpecShort) == "Multiple Pathways"] <- "Multi"
mr.sulf.bins.all$PathSpecShort[(mr.sulf.bins.all$PathSpecShort) == "Sulfur Disproportionation"] <- "Dispro"
mr.sulf.bins.all$PathSpecShort[(mr.sulf.bins.all$PathSpecShort) == "Sulfide Oxidation"] <- "H2S Ox"
mr.sulf.bins.all$PathSpecShort[(mr.sulf.bins.all$PathSpecShort) == "Sulfite Oxidation"] <- "SO3 Ox"
mr.sulf.bins.all$PathSpecShort[(mr.sulf.bins.all$PathSpecShort) == "Thiosulfate Oxidation"] <- "S2O3 Ox"

mr.sulf.bins.all$PathwaySpecific<-factor(mr.sulf.bins.all$PathwaySpecific,levels=c("Assimilatory Sulfate Reduction","Dissimilatory Sulfate Redox","Multiple Pathways","SOX","Sulfur Disproportionation","Sulfide Oxidation","Sulfite Oxidation","Thiosulfate Oxidation"))
mr.sulf.bins.all$PathSpecShort<-factor(mr.sulf.bins.all$PathSpecShort,levels=c("A.SO4 Red","D.SO4 RedOx","Multi","Dispro","H2S Ox","SO3 Ox","S2O3 Ox"))

mr.sulf.bins.all$KO_Function.KEGG = factor(mr.sulf.bins.all$KO_Function.KEGG, levels=unique(mr.sulf.bins.all$KO_Function.KEGG[order(mr.sulf.bins.all$Pathway)]), ordered=TRUE)

head(mr.sulf.bins.all)

# convert all 0s to NAs so they appear gray in heat map
mr.sulf.bins.all$MR_SumCovPerKO_NA<-ifelse(mr.sulf.bins.all$MR_SumCovPerKO==0,NA,mr.sulf.bins.all$MR_SumCovPerKO)
head(mr.sulf.bins.all)

# For heatmap color gradient
max(mr.sulf.bins.all$MR_SumCovPerKO, na.rm=TRUE)
max(mr.sulf.bins.all$MR_SumCovPerKO, na.rm=TRUE)/2
min(mr.sulf.bins.all$MR_SumCovPerKO, na.rm=TRUE)

# merge sulfur KO coverages with MAG taxonomic annotations
mr.sulf.taxa.all<-merge(mr.sulf.bins.all,mgm.bin.tax,by=c("Bin_ID","SampleID"))
head(mr.sulf.taxa.all)
mr.sulf.taxa.all$Genus = factor(mr.sulf.taxa.all$Genus, levels=unique(mr.sulf.taxa.all$Genus[order(mr.sulf.taxa.all$Site,mr.sulf.taxa.all$SampDate)]), ordered=TRUE)

# Figures below
# by Bin_ID

sulf.hm1a<-ggplot(mr.sulf.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sulfur Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(sulf.hm1a,filename = "figures/MGM_Figs/Bins/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_BinID_by_Function_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

sulf.hm1a2<-ggplot(mr.sulf.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sulfur Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Pathway~.,scales="free_y", space = "free")

ggsave(sulf.hm1a2,filename = "figures/MGM_Figs/Bins/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_BinID_by_Function_Pathway_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)

sulf.hm1a3<-ggplot(mr.sulf.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sulfur Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~Site,scales="free_x", space = "free")

ggsave(sulf.hm1a3,filename = "figures/MGM_Figs/Bins/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_by_BinID_by_Function_Site_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

# sulf.hm1a4<-ggplot(mr.sulf.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.15) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sulfur Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
#   theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")
#
# ggsave(sulf.hm1a4,filename = "figures/MGM_Figs/Bins/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_BinID_by_Function_Pathway_heatmap2.png", width=45, height=15, dpi=600,create.dir = TRUE)

sulf.hm1a5<-ggplot(mr.sulf.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sulfur Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Pathway~Site, scales="free", space = "free")

ggsave(sulf.hm1a5,filename = "figures/MGM_Figs/Bins/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_by_BinID_by_Function_Pathway_Site_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)

sulf.hm1a6<-ggplot(mr.sulf.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sulfur Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~Site, scales="free", space = "free")

ggsave(sulf.hm1a6,filename = "figures/MGM_Figs/Bins/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_by_BinID_by_Function_Pathway_Site_heatmap2.png", width=45, height=15, dpi=600,create.dir = TRUE)

sulf.hm1a7<-ggplot(mr.sulf.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sulfur Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathwaySpecific~Site, scales="free", space = "free")

ggsave(sulf.hm1a7,filename = "figures/MGM_Figs/Bins/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_by_BinID_by_Function_PathwaySpecific_Site_heatmap.png", width=45, height=20, dpi=600,create.dir = TRUE)

sulf.hm1a8<-ggplot(mr.sulf.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sulfur Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathSpecShort~Site, scales="free", space = "free")

ggsave(sulf.hm1a8,filename = "figures/MGM_Figs/Bins/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_by_BinID_by_Function_PathwaySpecific_Site_heatmap2.png", width=45, height=20, dpi=600,create.dir = TRUE)

sulf.hm1a9<-ggplot(mr.sulf.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sulfur Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathSpecShort~SampDate, scales="free", space = "free")

ggsave(sulf.hm1a9,filename = "figures/MGM_Figs/Bins/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_by_BinID_by_Function_PathwaySpecific_SampDate_heatmap.png", width=49.5, height=20, dpi=600,create.dir = TRUE)

sulf.hm1a10<-ggplot(mr.sulf.bins.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sulfur Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathSpecShort~SampDate+Site, scales="free", space = "free")

ggsave(sulf.hm1a10,filename = "figures/MGM_Figs/Bins/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_by_BinID_by_Function_PathwaySpecific_SampDate_Site_heatmap.png", width=49.5, height=20, dpi=600,create.dir = TRUE)

## by MAG Taxonomic ID, not Bin ID

sulf.hm1b<-ggplot(mr.sulf.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sulfur Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(sulf.hm1b,filename = "figures/MGM_Figs/Bins/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_Genus_by_Function_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

sulf.hm1b2<-ggplot(mr.sulf.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sulfur Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Pathway~.,scales="free_y", space = "free")

ggsave(sulf.hm1b2,filename = "figures/MGM_Figs/Bins/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_Genus_by_Function_Pathway_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)

sulf.hm1b3<-ggplot(mr.sulf.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sulfur Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~Site,scales="free_x", space = "free")

ggsave(sulf.hm1b3,filename = "figures/MGM_Figs/Bins/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_by_Genus_by_Function_Site_heatmap.png", width=45, height=13, dpi=600,create.dir = TRUE)

# sulf.hm1b4<-ggplot(mr.sulf.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.15) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sulfur Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
#   theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")
#
# ggsave(sulf.hm1b4,filename = "figures/MGM_Figs/Bins/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_Genus_by_Function_Pathway_heatmap2.png", width=45, height=15, dpi=600,create.dir = TRUE)

sulf.hm1b5<-ggplot(mr.sulf.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sulfur Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Pathway~Site, scales="free", space = "free")

ggsave(sulf.hm1b5,filename = "figures/MGM_Figs/Bins/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_by_Genus_by_Function_Pathway_Site_heatmap.png", width=45, height=15, dpi=600,create.dir = TRUE)

sulf.hm1b6<-ggplot(mr.sulf.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sulfur Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~Site, scales="free", space = "free")

ggsave(sulf.hm1b6,filename = "figures/MGM_Figs/Bins/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_by_Genus_by_Function_Pathway_Site_heatmap2.png", width=45, height=15, dpi=600,create.dir = TRUE)

sulf.hm1b7<-ggplot(mr.sulf.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sulfur Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathwaySpecific~Site, scales="free", space = "free")

ggsave(sulf.hm1b7,filename = "figures/MGM_Figs/Bins/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_by_Genus_by_Function_PathwaySpecific_Site_heatmap.png", width=45, height=20, dpi=600,create.dir = TRUE)

sulf.hm1b8<-ggplot(mr.sulf.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sulfur Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathSpecShort~Site, scales="free", space = "free")

ggsave(sulf.hm1b8,filename = "figures/MGM_Figs/Bins/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_by_Genus_by_Function_PathwaySpecific_Site_heatmap2.png", width=45, height=20, dpi=600,create.dir = TRUE)

sulf.hm1b9<-ggplot(mr.sulf.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sulfur Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathSpecShort~SampDate, scales="free", space = "free")

ggsave(sulf.hm1b9,filename = "figures/MGM_Figs/Bins/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_by_Genus_by_Function_PathwaySpecific_SampDate_heatmap.png", width=49.5, height=20, dpi=600,create.dir = TRUE)

sulf.hm1b10<-ggplot(mr.sulf.taxa.all, aes(Genus, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sulfur Metabolism in Salton Sea Dust MAGs",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
  theme(text=element_text(size=12, family="Arial"), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathSpecShort~SampDate+Site, scales="free", space = "free")

ggsave(sulf.hm1b10,filename = "figures/MGM_Figs/Bins/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_by_Genus_by_Function_PathwaySpecific_SampDate_Site_heatmap.png", width=49.5, height=20, dpi=600,create.dir = TRUE)

### Export Global Env for Other Scripts ####
save.image("data/SSD_MGM_Bins_Fxn_BetaDiv.Rdata")
# ^ includes all data combined in object bac.dat.bins.all, ASV table (samples are rows, ASVs are columns), all_bin_meta, and an ASV count table (where ASVs are rows, not columns)
# Version Information
sessionInfo()