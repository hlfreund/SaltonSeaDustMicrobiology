#### Set WD & Load Libraries ####
getwd() # use setwd("path/to/files") if you are not in the right directory
setwd("/Volumes/HLF_SSD/Aronson_Lab_Data/Salton_Sea/SaltonSeaWater")
suppressPackageStartupMessages({ # load packages quietly
  library(devtools)
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
})

#### Load Global Env to Import Count/ASV Tables ####
load("data/SSeawater_Data_Ready.Rdata") # save global env to Rdata file
load("data/SSeawater_BetaDiv_Data.Rdata")
#load("data/SSW_16S_CLR_EucDist_PCoA_Ready.Rdata")


#save.image("data/Env_Seqs_All/env.seq_analysis.Rdata") # save global env to Rdata file
bac.dat.all[1:6,1:6]
bac.ASV_table[1:4,1:4]
bac.ASV_table[(nrow(bac.ASV_table)-4):(nrow(bac.ASV_table)),(ncol(bac.ASV_table)-4):(ncol(bac.ASV_table))] # last 4 rows & cols
head(meta_scaled)

#### Beta Diversity ####
rownames(bac.ASV_table)
bac.ASV_table[1:4,1:4]

# CLR transformation of ASV table
# df must have rownames are SampleIDs, columns are ASV IDs for vegan functions below
b.clr<-decostand(bac.ASV_table[,-1],method = "clr", pseudocount = 1) #CLR transformation
b.clr[1:4,1:4]

# check rownames of CLR transformed ASV data & metadata
rownames(b.clr) %in% rownames(meta_scaled)
meta_scaled=meta_scaled[rownames(b.clr),] ## reorder metadata to match order of CLR data

# calculate our Euclidean distance matrix using CLR data
b.euc_dist <- dist(b.clr, method = "euclidean")

# creating our hierarcical clustering dendrogram
b.euc_clust <- hclust(b.euc_dist, method="ward.D2")

# let's make it a little nicer...
b.euc_dend <- as.dendrogram(b.euc_clust, hang=0.2)
b.dend_cols <- as.character(meta_scaled$SampDate_Color[order.dendrogram(b.euc_dend)])
labels_colors(b.euc_dend) <- b.dend_cols

## DO NOT RUN THIS LINE, THIS IS YOUR COLOR REFERENCE!!!!
#(August.2021="#ef781c",December.2021="#03045e",April.2022="#059c3f")

plot(b.euc_dend, ylab="CLR Euclidean Distance",cex = 0.5) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 2)
#legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c("#ef781c","#03045e","#059c3f"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# PCOA w/ Euclidean distance matrix (of CLR data)
b.pcoa <- pcoa(b.euc_dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
#save.image("data/SSW_16S_CLR_EucDist_PCoA_Ready.Rdata")

# The proportion of variances explained is in its element values$Relative_eig
b.pcoa$values

# extract principal coordinates
b.pcoa.vectors<-data.frame(b.pcoa$vectors)
b.pcoa.vectors$SampleID<-rownames(b.pcoa$vectors)

# merge pcoa coordinates w/ metadata
b.pcoa.meta<-merge(b.pcoa.vectors, meta_scaled, by.x="SampleID", by.y="SampleID")
b.pcoa.meta$SampleMonth
b.pcoa.meta$SampDate

head(b.pcoa.meta)

head(b.pcoa$values) # pull out Relative (Relative_eig) variation % to add to axes labels
save.image("data/SSW_16S_CLR_EucDist_PCoA_Ready.Rdata")

#### Visualize PCoAs ####
# create PCoA ggplot fig
pcoa1<-ggplot(b.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(SampDate)), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Seawater",subtitle="Using Centered-Log Ratio Data",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Type",values=unique(b.pcoa.meta$SampDate_Color[order(b.pcoa.meta$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("PC1 [31.99%]") + ylab("PC2 [27.38%]")

ggsave(pcoa1,filename = "figures/BetaDiversity/SSW_16S_pcoa_CLR_sampdate.png", width=12, height=10, dpi=600)

# sample month shape, depth color
pcoa2<-ggplot(b.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=as.numeric(as.character(Depth_m)),shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Seawater",subtitle="Using Centered-Log Ratio Data",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  xlab("PC1 [31.99%]") + ylab("PC2 [27.38%]")

ggsave(pcoa2,filename = "figures/BetaDiversity/SSW_16S_pcoa_CLR_depth_sampdate.png", width=12, height=10, dpi=600)

#### Homogeneity of Variance & PERMANOVA tests - Composition by Groups ####
## betadisper to look at homogeneity of group dispersions (aka variance) when considering multiple variables
# multivariate analogue to Levene's test of homogeneity of variances
# program finds spatial median or centroid of the group, & compare distances of group to centroid/spatial median via ANOVA

#While PERMANOVA tests differences in group means (analogous to MANOVA),
## a related test called PERMDISP can be used to evaluate homogeneity of group dispersion
#(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
## * need a distance matrix!

rownames(meta_scaled) %in% rownames(b.clr) #b.clr was used to make the distance matrix b.euc_dist

# first by compare dispersions by sampling date
b.disper1<-betadisper((vegdist(b.clr,method="euclidean")), meta_scaled$SampDate)
b.disper1

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:

permutest(b.disper1, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#                 August.2021 December.2021 April.2022
# August.2021                   0.0040000      0.001
# December.2021   0.0073988                    0.115
# April.2022      0.0012603     0.1201286

anova(b.disper1) # p = 0.0003451 --> reject the Null H, spatial medians (a measure of dispersion) are significantly difference across sample dates

TukeyHSD(b.disper1) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other
#                               diff       lwr       upr     p adj
# December.2021-August.2021 -4.556897 -7.631983 -1.481811 0.0033548
# April.2022-August.2021    -5.605532 -8.680617 -2.530446 0.0004428
# April.2022-December.2021  -1.048634 -4.123720  2.026451 0.6710007

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

pnova1<-adonis2(b.clr ~ SampDate,data=meta_scaled,method = "euclidean",by="terms",permutations=1000)
pnova1 # p-value = 0.000999

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

b.clr.dist = (vegdist(b.clr, "euclidean", na.rm = TRUE)) #distance matrix using Bray's dissimilarity index for trait distribution (traits of interest only)
pair.mod1<-pairwise.adonis(b.clr.dist,meta_scaled$SampDate, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
pair.mod1
#                           pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
# 1  December.2021 vs April.2022  1  7006.287 17.13323 0.5503196   0.001      0.003   *
# 2 December.2021 vs August.2021  1  7416.531 13.48028 0.4905437   0.001      0.003   *
# 3    April.2022 vs August.2021  1  7987.429 15.15788 0.5198554   0.001      0.003   *

# Visualize dispersions
png('figures/BetaDiversity/pcoa_betadispersion_sampledate.png',width = 700, height = 600, res=100)
plot(b.disper1,main = "Centroids and Dispersion based on Aitchison Distance", col=colorset1$SampDate_Color)
dev.off()

png('boxplot_centroid_distance_sampledate.png',width = 700, height = 600, res=100)
boxplot(b.disper1,xlab="Sample Collection Date", main = "Distance to Centroid by Category", sub="Based on Aitchison Distance", col=colorset1$SampDate_Color)
dev.off()

# Next compare dispersions by depth
b.disper2<-betadisper((vegdist(b.clr,method="euclidean")), meta_scaled$Depth_m)
b.disper2

permutest(b.disper2, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons

anova(b.disper2) # p = 0.6277 --> accept the Null H, spatial medians are NOT significantly difference across sample dates

TukeyHSD(b.disper2) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other
# no sig results

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

pnova2<-adonis2(b.clr ~ Depth_m,data=meta_scaled,method = "euclidean",by="terms",permutations=1000)
pnova2 # p-value = 1

#b.clr.dist = (vegdist(b.clr, "euclidean", na.rm = TRUE)) #distance matrix using Bray's dissimilarity index for trait distribution (traits of interest only)
pair.mod2<-pairwise.adonis(b.clr.dist,meta_scaled$Depth_m, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
pair.mod2
# none are significantly different

col.depth <- colorRampPalette(c("red", "blue"))
col.depth(9)

# Visualize dispersions
png('figures/BetaDiversity/pcoa_betadispersion_depth.png',width = 700, height = 600, res=100)
plot(b.disper2,main = "Centroids and Dispersion based on Aitchison Distance", col=col.depth(9))
dev.off()

png('figures/BetaDiversity/boxplot_centroid_distance_depth.png',width = 700, height = 600, res=100)
boxplot(b.disper2,xlab="Sample Collection Depth", main = "Distance to Centroid by Category", sub="Based on Aitchison Distance", col=col.depth(9))
dev.off()

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

# create column for Depth that is a numeric version of this variable, rather than a factor
meta_scaled$Depth.num<-as.numeric(as.character(meta_scaled$Depth_m))

# now make sure your data frames you're comparing are in the same exact order!!
rownames(b.clr) %in% rownames(meta_scaled)
meta_scaled=meta_scaled[rownames(b.clr),] ## reorder metadata to match order of CLR data
perm <- with(meta_scaled, how(nperm = 1000)) # using SampDate as block because there is a significant difference between sample dates, trying to remove this effect when looking at permanovas

pnova1<-adonis2(b.clr ~ DO_Percent_Local*ORP_mV*Temp_DegC*Dissolved_OrganicMatter_RFU*Depth.num*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova1
# nothing

adonis2(b.clr ~ DO_Percent_Local*ORP_mV*Temp_DegC*Dissolved_OrganicMatter_RFU*Depth.num*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)
#         Df SumOfSqs     R2    F Pr(>F)
#Model    23    25343  1
#Residual  0        0  0
#Total    23    25343  1

pnova2<-adonis2(b.clr ~ DO_Percent_Local*ORP_mV*Temp_DegC*Dissolved_OrganicMatter_RFU*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova2
# nothing significant

adonis2(b.clr ~ DO_Percent_Local*ORP_mV*Temp_DegC*Dissolved_OrganicMatter_RFU*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)
#         Df SumOfSqs      R2      F   Pr(>F)
#Model    23    34412 0.73114 1.8918 0.4615
#Residual 16    12654 0.26886
#Total    39    47066 1.00000

pnova3<-adonis2(b.clr ~ DO_Percent_Local*ORP_mV*Temp_DegC*Dissolved_OrganicMatter_RFU*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova3

adonis2(b.clr ~ ORP_mV*DO_Percent_Local*Dissolved_OrganicMatter_RFU*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)

pnova4<-adonis2(b.clr ~ ORP_mV*DO_Percent_Local*Dissolved_OrganicMatter_RFU*Sulfate_milliM,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova4
#                                         Df SumOfSqs      R2       F   Pr(>F)
# ORP_mV                                                              1   3435.1 0.13555  8.0320 0.000999 ***
# DO_Percent_Local                                                    1   4009.5 0.15821  9.3749 0.000999 ***
# Dissolved_OrganicMatter_RFU                                         1   5843.0 0.23056 13.6620 0.000999 ***
# Sulfate_milliM                                                      1    961.4 0.03794  2.2480 0.042957 *
# ORP_mV:DO_Percent_Local                                             1   1322.9 0.05220  3.0932 0.002997 **
# DO_Percent_Local:Dissolved_OrganicMatter_RFU                        1   1413.3 0.05577  3.3045 0.005994 **
# DO_Percent_Local:Sulfate_milliM                                     1    794.0 0.03133  1.8566 0.072927 .

adonis2(b.clr ~ ORP_mV*DO_Percent_Local*Dissolved_OrganicMatter_RFU*Sulfate_milliM,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)
#         Df SumOfSqs      R2      F   Pr(>F)
#Model    15  21921.1 0.86499 3.417 0.000999 ***
#Residual  8   3421.5 0.13501
#Total    23  25342.6 1.00000

pnova4b<-adonis2(b.clr ~ DO_Percent_Local*Dissolved_OrganicMatter_RFU*ORP_mV,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova4b
#                                   Df SumOfSqs      R2       F   Pr(>F)
# DO_Percent_Local                                     1   5830.3 0.23006 12.2377 0.000999 ***
#Dissolved_OrganicMatter_RFU                          1   6344.6 0.25036 13.3172 0.000999 ***
# ORP_mV                                               1   1112.7 0.04391  2.3355 0.027972 *
# DO_Percent_Local:Dissolved_OrganicMatter_RFU         1   1268.0 0.05004  2.6616 0.010989 *
# DO_Percent_Local:ORP_mV                              1   1792.5 0.07073  3.7624 0.000999 ***
# Dissolved_OrganicMatter_RFU:ORP_mV                   1    652.0 0.02573  1.3686 0.173826
# DO_Percent_Local:Dissolved_OrganicMatter_RFU:ORP_mV  1    719.5 0.02839  1.5103 0.159840
# Residual                                            16   7622.8 0.30079
# Total                                               23  25342.6 1.00000

pnova5<-adonis2(b.clr ~ DO_Percent_Local*Dissolved_OrganicMatter_RFU,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova5
#                                               Df SumOfSqs      R2       F   Pr(>F)
# DO_Percent_Local                              1   5830.3 0.23006  9.5297 0.000999 ***
#   Dissolved_OrganicMatter_RFU                   1   6344.6 0.25036 10.3704 0.000999 ***
#   DO_Percent_Local:Dissolved_OrganicMatter_RFU  1    931.5 0.03676  1.5226 0.134865
# Residual                                     20  12236.1 0.48283
# Total                                        23  25342.6 1.00000

adonis2(b.clr ~ DO_Percent_Local*Dissolved_OrganicMatter_RFU,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)
#         Df SumOfSqs      R2      F   Pr(>F)
#Model     3    13106 0.51717 7.1409 0.000999 ***
#Residual 20    12236 0.48283
#Total    23    25343 1.00000

### SELF REMINDER FOR R^2
### Coefficient of Determination, denoted R2 or r2
### is the proportion of the variance in the dependent variable that is predictable from the independent variable(s)

### Pseudo F stat for PERMANOVA
### pseudo F-ratio: It compares the total sum of squared dissimilarities (or ranked dissimilarities) among objects belonging to different groups to that of objects belonging to the same group.
### Larger F-ratios indicate more pronounced group separation, however, the significance of this ratio is usually of more interest than its magnitude.

#### Save Everything ####
save.image("data/SSeawater_BetaDiv_Data.Rdata")
