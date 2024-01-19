#### Set WD & Load Libraries ####
getwd() # use setwd("path/to/files") if you are not in the right directory
setwd("/Volumes/HLF_SSD/Aronson_Lab_Data/Salton_Sea/SaltonSeaWater")
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
})

#### Load Global Env to Import Count/ASV Tables ####
load("data/SSeawater_Data_Ready.Rdata") # save global env to Rdata file
load("data/SSeawater_AlphaDiv_Data.Rdata")
#load("data/ssw_clr.euc.dist_2.21.23.Rdata")

#save.image("data/Env_Seqs_All/env.seq_analysis.Rdata") # save global env to Rdata file
bac.dat.all[1:6,1:6]
bac.ASV_table[1:4,1:4]
bac.ASV_table[(nrow(bac.ASV_table)-4):(nrow(bac.ASV_table)),(ncol(bac.ASV_table)-4):(ncol(bac.ASV_table))] # last 4 rows & cols
head(meta_scaled)

## DO NOT RUN THIS LINE, THIS IS YOUR COLOR REFERENCE!!!!
#(August.2021="#ef781c",December.2021="#03045e",April.2022="#059c3f")

#### Rarefaction & Species Accumulation Curves ####
# bacteria/archaea

# Species Accumulation Curve
sc2<-specaccum(bac.ASV_table[,-1],"random")
plot(sc2, ci.type="poly", col="darkgreen", lwd=2, ci.lty=0, ci.col="lightgreen")
boxplot(sc2, col="yellow", add=TRUE, pch=20)

# Prep for Rarefaction Curve
rowSums(bac.ASV_table[,-1]) # total # ASVs per sample, excluding SampleID from calculation
sort(colSums(bac.ASV_table[,-1]))

# Create Rarefaction curve
png('figures/AlphaDiversity/SSW_16S_rarecurve.png')
rarecurve(as.matrix(bac.ASV_table[,-1]),col=metadata$SampDate_Color, step=1000, label=F,ylab="ASVs")
# to show sampel labels per curve, change label=T
dev.off()

# ASVs per Sample
total_asvs<-data.frame(ASV_Total=rowSums(bac.ASV_table[,-1]),metadata)

ggplot(data=total_asvs, aes(x=SampleID, y=ASV_Total,fill=Sample_Type)) +
  geom_bar(stat="identity",colour="black")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## average ASV per sample month & depth
aggregate(bac.ASV_all$Count, list(bac.ASV_all$SampleMonth), FUN=mean)
aggregate(bac.ASV_all$Count, list(bac.ASV_all$Depth_m), FUN=mean)

#### Alpha Diversity & Species Richness ####

## Calculate Shannon Diversity (abundance + richness considered in diversity calculation)
# if you have another package loaded that has a diversity function, you can specify that you want to use vegan's diversity function as shown below
Shan_ent.16s<-vegan::diversity(bac.ASV_table[,-1], index="shannon") # Shannon entropy
Shan_div.16s<- exp(Shan_ent.16s) # Shannon Diversity aka Hill number 1

# create data frame with Shannon entropy and Shannon diversity values
div_16s<-data.frame(Bac_Shannon_Entropy=Shan_ent.16s,Bac_Shannon_Diversity=Shan_div.16s)
class(div_16s)
div_16s$SampleID<-rownames(div_16s)
head(div_16s)

# Calculate species richness (number of species per sample)
specnumber(bac.ASV_table[,-1])

# Create a DF with Species Richness
S_16s<-data.frame(Bac_Species_Richness=specnumber(bac.ASV_table[,-1]), SampleID=rownames(bac.ASV_table)) # finds # of species per sample using RAW count data; if MARGIN = 2 it finds frequencies of species

# merge richness and diversity dataframes together
d.r_16s<-merge(div_16s, S_16s, by.x="SampleID", by.y="SampleID")

# merge w/ metadata
bac.div.metadat <- merge(d.r_16s,meta_scaled, by.x="SampleID", by.y="SampleID")
head(bac.div.metadat)
class(bac.div.metadat) # want data frame

unique(bac.div.metadat$SampleMonth) # see how many elements there are in the Group variable
unique(bac.div.metadat$Depth_m) # see how many elements there are in the Group variable
bac.div.metadat$Depth_m<-factor(bac.div.metadat$Depth_m, levels=c("0","3","4","5","7","9","10","11"))

# drop the outliers
#bac.div.metadat<-bac.div.metadat[bac.div.metadat$Bac_Shannon_Diversity<300 & bac.div.metadat$Bac_Species_Richness>100,]

# create numeric variable for depth to be used for models later
bac.div.metadat$Depth.num<-as.numeric(as.character(bac.div.metadat$Depth_m))

# save diversity data
save.image("data/SSeawater_AlphaDiv_Data.Rdata")

#### Visualize Alpha Diversity & Species Richness ####
## Shannon Diversity by Sample Month & Depth
bac.a.div<-ggplot(bac.div.metadat, aes(x=SampDate, y=Bac_Shannon_Diversity,fill=SampDate)) +geom_boxplot(color="black")+scale_x_discrete(labels=c("August 2021","December 2021","April 2022"))+theme_bw()+
  scale_fill_manual(values=unique(bac.div.metadat$SampDate_Color[order(bac.div.metadat$SampDate)]), name ="Sample Date",labels=c("August 2021","December 2021","April 2022"))+theme_classic()+
  labs(title = "Bacterial Shannon Diversity by Sample Date", x="Sample Date", y="Shannon Diversity", fill="Sample Month")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3)), method="t.test", hide.ns = TRUE,label = "p.format")

ggsave(bac.a.div,filename = "figures/AlphaDiversity/SSW_Bacterial_alpha_diversity_sampledate.png", width=13, height=10, dpi=600)

bac.div.metadat$Depth_m=as.numeric(levels(bac.div.metadat$Depth_m))[bac.div.metadat$Depth_m]
# ^ note: cannot turn numbers that are factors in R into numeric values...
## have to convert factor levels into numeric, then use the numeric "levels" to pull out numbers from Depth_m column in df to make sure the Depth_m columns is now numeric, not a factor

bac.a.div2<-ggplot(bac.div.metadat, aes(x=as.factor(Depth_m), y=Bac_Shannon_Diversity)) +geom_boxplot(aes(fill=Depth_m),color="black")+
  labs(title = "Bacterial Shannon Diversity by Sampling Depth", x="Depth (m)", y="Shannon Diversity", fill="Depth (m)")+
  scale_fill_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  coord_flip() + scale_x_discrete(limits=rev)

ggsave(bac.a.div2,filename = "figures/AlphaDiversity/SSW_Bacterial_alpha_diversity_depth_v1.png", width=13, height=10, dpi=600)

bac.a.div3<-ggplot(bac.div.metadat, aes(x=as.factor(Depth_m), y=Bac_Shannon_Diversity)) +geom_boxplot(aes(fill=Depth_m),color="black")+
  labs(title = "Bacterial Shannon Diversity by Sampling Depth", x="Depth (m)", y="Shannon Diversity", fill="Depth (m)")+
  scale_fill_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))

ggsave(bac.a.div3,filename = "figures/AlphaDiversity/SSW_Bacterial_alpha_diversity_depth_v2.png", width=13, height=10, dpi=600)

## Species Richness by Sample Type
bac.a.sr<-ggplot(bac.div.metadat, aes(x=SampDate, y=Bac_Species_Richness,fill=SampDate)) +geom_boxplot(color="black")+scale_x_discrete(labels=c("August 2021","December 2021","April 2022"))+theme_bw()+
  scale_fill_manual(values=unique(bac.div.metadat$SampDate_Color[order(bac.div.metadat$SampDate)]), name ="Sample Date",labels=c("August 2021","December 2021","April 2022"))+theme_classic()+
  labs(title = "Bacterial Species Richness by Sample Date", x="Sample Date", y="Species Richness")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3)), method="t.test", hide.ns = TRUE,label = "p.format")

ggsave(bac.a.sr,filename = "figures/AlphaDiversity/SSW_Bacterial_species_richness_samplemonth.png", width=13, height=10, dpi=600)

bac.a.sr2<-ggplot(bac.div.metadat, aes(x=as.factor(Depth_m), y=Bac_Species_Richness,fill=bac.div.metadat$Depth_m)) +geom_boxplot(aes(fill=as.numeric(bac.div.metadat$Depth_m)),color="black")+
  labs(title = "Bacterial Species Richness by Sampling Depth", x="Depth (m)", y="Species Richness", fill="Depth (m)")+
  scale_fill_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  coord_flip() + scale_x_discrete(limits=rev)

ggsave(bac.a.sr2,filename = "figures/AlphaDiversity/SSW_Bacterial_species_richness_depth_v1.png", width=13, height=10, dpi=600)

bac.a.sr3<-ggplot(bac.div.metadat, aes(x=as.factor(Depth_m), y=Bac_Species_Richness,fill=bac.div.metadat$Depth_m)) +geom_boxplot(aes(fill=as.numeric(bac.div.metadat$Depth_m)),color="black")+
  labs(title = "Bacterial Species Richness by Sampling Depth", x="Depth (m)", y="Species Richness", fill="Depth (m)")+
  scale_fill_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))

ggsave(bac.a.sr3,filename = "figures/AlphaDiversity/SSW_Bacterial_species_richness_depth_v2.png", width=13, height=10, dpi=600)

#### Using Shapiro-Wilk test for Normality ####
shapiro.test(bac.div.metadat$Bac_Shannon_Diversity) # what is the p-value?
# p-value = 0.5816
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(bac.div.metadat$Bac_Shannon_Diversity, col="blue") # with outliars

# visualize Q-Q plot for alpha div
# The Q-Q plot, or quantile-quantile plot, is a graphical tool to help us assess if a set of data plausibly came from some theoretical distribution such as a normal or exponential.
# For example, if we run a statistical analysis that assumes our residuals are normally distributed, we can use a normal Q-Q plot to check that assumption
# more on Q-Q plots here: https://data.library.virginia.edu/understanding-q-q-plots/
# more here too: https://grodri.github.io/glms/notes/c2s9#:~:text=8%20The%20Q%2DQ%20Plot,versus%20quantiles%20of%20a%20distribution.
qqnorm(bac.div.metadat$Bac_Shannon_Diversity, pch = 1, frame = FALSE)
qqline(bac.div.metadat$Bac_Shannon_Diversity, col = "red", lwd = 2)

shapiro.test(bac.div.metadat$Bac_Species_Richness) # what is the p-value?
# p-value = 0.6063
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(bac.div.metadat$Bac_Species_Richness, col="blue")

# visualize Q-Q plot for species richness
qqnorm(bac.div.metadat$Bac_Species_Richness, pch = 1, frame = FALSE) # with outliars
qqline(bac.div.metadat$Bac_Species_Richness, col = "red", lwd = 2)

qqnorm(bac.div.metadat$Bac_Species_Richness, pch = 1, frame = FALSE) # without outliars
qqline(bac.div.metadat$Bac_Species_Richness, col = "red", lwd = 2)

shapiro.test(bac.div.metadat$DO_Percent_Local) # p-value = 0.02586
hist(bac.div.metadat$DO_Percent_Local, col="blue")
# visualize Q-Q plot for species richness
qqnorm(bac.div.metadat$DO_Percent_Local, pch = 1, frame = FALSE) # with outliars
qqline(bac.div.metadat$DO_Percent_Local, col = "red", lwd = 2)

shapiro.test(bac.div.metadat$ORP_mV) # p-value = 1.731e-08
hist(bac.div.metadat$ORP_mV, col="blue")
# visualize Q-Q plot for species richness
qqnorm(bac.div.metadat$ORP_mV, pch = 1, frame = FALSE) # with outliars
qqline(bac.div.metadat$ORP_mV, col = "red", lwd = 2)

shapiro.test(bac.div.metadat$Temp_DegC) # p-value = 0.0002829
hist(bac.div.metadat$Temp_DegC, col="blue")
# visualize Q-Q plot for species richness
qqnorm(bac.div.metadat$Temp_DegC, pch = 1, frame = FALSE) # with outliars
qqline(bac.div.metadat$Temp_DegC, col = "red", lwd = 2)

shapiro.test(bac.div.metadat$Dissolved_OrganicMatter_RFU) #  p-value = 0.05411
hist(bac.div.metadat$Dissolved_OrganicMatter_RFU, col="blue")
# visualize Q-Q plot for species richness
qqnorm(bac.div.metadat$Dissolved_OrganicMatter_RFU, pch = 1, frame = FALSE) # with outliars
qqline(bac.div.metadat$Dissolved_OrganicMatter_RFU, col = "red", lwd = 2)

shapiro.test(bac.div.metadat$Sulfate_milliM) # p-value = 0.1912
hist(bac.div.metadat$Sulfate_milliM, col="blue")
# visualize Q-Q plot for species richness
qqnorm(bac.div.metadat$Sulfate_milliM, pch = 1, frame = FALSE) # with outliars
qqline(bac.div.metadat$Sulfate_milliM, col = "red", lwd = 2)

shapiro.test(bac.div.metadat$Sulfide_microM) # p-value = 3.813e-08
hist(bac.div.metadat$Sulfide_microM, col="blue")
# visualize Q-Q plot for species richness
qqnorm(bac.div.metadat$Sulfide_microM, pch = 1, frame = FALSE) # with outliars
qqline(bac.div.metadat$Sulfide_microM, col = "red", lwd = 2)

#### Linear Regression Comparisons - Shannon Diversity ####
## here the focus is comparing dust complexity to alpha diversity, species richness, & elevation
head(bac.div.metadat)

# just look at everything at once in step-wise fashion
step1<-step(glm(formula = Bac_Shannon_Diversity ~ ., data=bac.div.metadat[,c(3,11,13:14,18:20)]))
#                 Estimate Std. Error t value Pr(>|t|)
# ORP_mV           -7.049      3.865  -1.824 0.083128 .
# Temp_DegC       -20.691      4.716  -4.387 0.000285 ***
# Sulfate_milliM  -10.211      4.099  -2.491 0.021651 *

div.glm.all<-glm(formula = Bac_Shannon_Diversity ~ ., data=bac.div.metadat[,c(3,11,13:14,18:20)])
summary(div.glm.all)

div.glm.p<-coef(summary(div.glm.all))[,4] # p-values
Div.GLM.Pval<-data.frame(Div.GLM.AdjPval=p.adjust(div.glm.p, method="bonferroni",n=length(div.glm.p)),Div.GLM.Pval=div.glm.p)
#                               Div.GLM.AdjPval Div.GLM.Pval
# (Intercept)                    1.209565e-14 1.727951e-15
# DO_Percent_Local               1.000000e+00 5.749260e-01
# ORP_mV                         1.000000e+00 3.469644e-01
# Temp_DegC                      7.803731e-02 1.114819e-02 -- near sig after adjustment
# Dissolved_OrganicMatter_RFU    1.000000e+00 9.064130e-01
# Sulfate_milliM                 3.569710e-01 5.099586e-02
# Sulfide_microM                 1.000000e+00 6.472779e-01

# [Example Code for Different Models]

s.div.glm.fit1<-glm(formula = Bac_Shannon_Diversity ~ DO_Percent_Local, data=bac.div.metadat)%>%
  adjust_pvalue(method="bonferroni")

# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.div.glm.fit1)

# sanity check that lm() vs glm(familiy=Gaussian) is the same thing - and it is!
summary(lm(formula = Bac_Shannon_Diversity ~ DO_Percent_Local, data=bac.div.metadat)%>%
          adjust_pvalue(method="bonferroni"))

# code for mixed effects model
mixed1 = lmer(Bac_Shannon_Diversity ~ DO_Percent_Local+ (1 | Depth_m), data = bac.div.metadat)
summary(mixed1)

# Shan Div ~ DO%
plot(Bac_Shannon_Diversity ~ DO_Percent_Local, data=bac.div.metadat,col=SampDate_Color)

# Shan Div ~ ORP

plot(Bac_Shannon_Diversity ~ ORP_mV, data=bac.div.metadat,col=SampDate_Color)

# Shan Div ~ Temp (C)

plot(Bac_Shannon_Diversity ~ Temp_DegC, data=bac.div.metadat,col=SampDate_Color)

# Shan Div ~ DOM

plot(Bac_Shannon_Diversity ~ Dissolved_OrganicMatter_RFU, data=bac.div.metadat,col=SampDate_Color)

# Shan Div ~ Sulfate

plot(Bac_Shannon_Diversity ~ Sulfate_milliM, data=bac.div.metadat,col=SampDate_Color)

# Shan Div ~ Sulfide

plot(Bac_Shannon_Diversity ~ Sulfide_microM, data=bac.div.metadat,col=SampDate_Color)

# Shan Div ~ Depth
plot(Bac_Shannon_Diversity ~ Depth.num, data=bac.div.metadat,col=SampDate_Color)

fit1<-aov(Bac_Shannon_Diversity ~ SampDate, data=bac.div.metadat)
#pairwise.adonis(bac.div.metadat$Bac_Shannon_Diversity, bac.div.metadat$SampDate, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(fit1)
#Df           Sum Sq Mean Sq    F value   Pr(>F)
#SampDate     2   3214  1607.1   5.317 0.0135 *
#Residuals   21   6347   302.2

# Tukey test - tells us which groups are significantly different from each other (more here: https://www.r-bloggers.com/2013/06/anova-and-tukeys-test-on-r/)
Tuk1<-TukeyHSD(fit1)
Tuk1$SampDate
#                             diff        lwr      upr      p adj
# December.2021-August.2021 24.7566611   2.846671 46.66665 0.02503300 *
# April.2022-August.2021    24.3365005   2.426510 46.24649 0.02778694 *
# April.2022-December.2021  -0.4201606 -22.330151 21.48983 0.99871280

# fit.0<-aov(DustComplexity ~ as.factor(Elevation), data=bac.div.metadat)
# summary(fit.0)
# TukeyHSD(fit.0)
# Levene's test with one independent variable
## Levene's tests whether variances of 2 samples are equal
## we want variances to be the same -- want NON SIGNIFICANCE!
## t test assumes that variances are the same, so Levene's test needs to be non significant
## Fligner's test is a Levene's test for data that are not normally distributed
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(Bac_Shannon_Diversity ~ SampDate, data = bac.div.metadat)
# Fligner-Killeen:med chi-squared = 1.1963, df = 7, p-value = 0.991
# Which shows that the data do not deviate significantly from homogeneity.
compare_means(Bac_Shannon_Diversity ~ SampDate, data=bac.div.metadat, method="anova",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input

#### Linear Regression Comparisons - Species Richness ####
## here the focus is comparing dust complexity to alpha diversity, species richness, & elevation
head(bac.div.metadat) # bac.div.metadat - excludes outliar with very high Shannon diversity

# just look at everything at once in step-wise fashion
step2<-step(glm(formula = Bac_Species_Richness ~ ., data=bac.div.metadat[,c(4,11,13:14,18:20)]))
summary(step2)
#                               Estimate Std. Error t value Pr(>|t|)
# Temp_DegC                     114.82      62.20   1.846   0.0798 .
# Dissolved_OrganicMatter_RFU   128.54      52.69   2.439   0.0241 *

sr.glm.all<-glm(formula = Bac_Species_Richness ~ ., data=bac.div.metadat[,c(4,11,13:14,18:20)])
summary(sr.glm.all)

sr.glm.p<-coef(summary(sr.glm.all))[,4] # p-values
SR.GLM.Pval<-data.frame(SR.GLM.AdjPval=p.adjust(sr.glm.p, method="bonferroni",n=length(sr.glm.p)),SR.GLM.Pval=sr.glm.p)
#                               SR.GLM.AdjPval  SR.GLM.Pval
# (Intercept)                   1.388731e-14 1.983901e-15
# DO_Percent_Local              1.000000e+00 1.836612e-01
# ORP_mV                        1.000000e+00 2.813919e-01
# Temp_DegC                     1.000000e+00 1.847486e-01
# Dissolved_OrganicMatter_RFU   5.309113e-01 7.584447e-02
# Sulfate_milliM                1.000000e+00 4.890501e-01
# Sulfide_microM                1.000000e+00 3.932963e-01

# Species Richness ~ DO%

plot(Bac_Species_Richness ~ DO_Percent_Local, data=bac.div.metadat,col=SampDate_Color)

# Species Richness ~ ORP

plot(Bac_Species_Richness ~ ORP_mV, data=bac.div.metadat,col=SampDate_Color)

# Species Richness ~ Temp

plot(Bac_Species_Richness ~ Temp_DegC, data=bac.div.metadat,col=SampDate_Color)

# Species Richness ~ DOM

plot(Bac_Species_Richness ~ Dissolved_OrganicMatter_RFU, data=bac.div.metadat,col=SampDate_Color)

# Species Richness ~ Sulfate

plot(Bac_Species_Richness ~ Sulfate_milliM, data=bac.div.metadat,col=SampDate_Color)

# Species Richness ~ Sulfide

plot(Bac_Species_Richness ~ Sulfide_microM, data=bac.div.metadat,col=SampDate_Color)

# Species Richness ~ Depth

plot(Bac_Species_Richness ~ Depth.num, data=bac.div.metadat,col=SampDate_Color)

fit2<-aov(Bac_Species_Richness ~ SampDate, data=bac.div.metadat)
#pairwise.adonis(bac.div.metadat$Bac_Species_Richness, bac.div.metadat$Depth_m, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(fit2)
#Df           Sum Sq Mean Sq    F value   Pr(>F)
#            Df Sum Sq Mean Sq F value  Pr(>F)
#SampDate     2 378899  189450   6.465 0.00649 **
#Residuals   21 615429   29306

# Tukey test - tells us which groups are significantly different from each other (more here: https://www.r-bloggers.com/2013/06/anova-and-tukeys-test-on-r/)
Tuk2<-TukeyHSD(fit2)
Tuk2$SampDate
#                               diff       lwr       upr       p adj
# December.2021-August.2021 -182.500 -398.2486  33.24863 0.107382247
# April.2022-August.2021    -305.875 -521.6236 -90.12637 0.004886271 **
# April.2022-December.2021  -123.375 -339.1236  92.37363 0.338647759

# Levene's test with one independent variable
## Levene's tests whether variances of 2 samples are equal
## we want variances to be the same -- want NON SIGNIFICANCE!
## t test assumes that variances are the same, so Levene's test needs to be non significant
## Fligner's test is a Levene's test for data that are not normally distributed
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(Bac_Species_Richness ~ SampDate, data = bac.div.metadat)
# Fligner-Killeen:med chi-squared = 0.45636, df = 2, p-value = 0.796
# Which shows that the data do not deviate significantly from homogeneity.
compare_means(Bac_Species_Richness ~ SampDate, data=bac.div.metadat, method="anova",p.adjust.method = "bonferroni")

#### Prep Data for Linear Regressions within Timepoints ####
## here the focus is comparing dust complexity to alpha diversity, species richness, & elevation
head(bac.div.metadat)

# create the dataframes
aug21.div<-subset(bac.div.metadat, bac.div.metadat$SampDate=="August.2021")
dec21.div<-subset(bac.div.metadat, bac.div.metadat$SampDate=="December.2021")
apr22.div<-subset(bac.div.metadat, bac.div.metadat$SampDate=="April.2022")

#### August - Shannon Diversity ####
# August 2021
aug21.div.glm.fit1<-glm(formula = Bac_Shannon_Diversity ~ DO_Percent_Local, data=aug21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(aug21.div.glm.fit1)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)       0.0112706  0.0004664  24.166   <2e-16 ***
#DO_Percent_Local -0.0009547  0.0005092  -1.875   0.0687 .

aug21.div.glm.fit2<-glm(formula = Bac_Shannon_Diversity ~ ORP_mV, data=aug21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(aug21.div.glm.fit2)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)   0.0112343  0.0004731  23.745   <2e-16 ***
#ORP_mV      -0.0001512  0.0004968  -0.304    0.763

aug21.div.glm.fit3<-glm(formula = Bac_Shannon_Diversity ~ Temp_DegC, data=aug21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(aug21.div.glm.fit3)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept) 0.0113225  0.0004451  25.439   <2e-16 ***
#Temp_DegC   0.0011947  0.0004861   2.458   0.0188 *

aug21.div.glm.fit4<-glm(formula = Bac_Shannon_Diversity ~ Dissolved_OrganicMatter_RFU, data=aug21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(aug21.div.glm.fit4)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)                 0.0112493  0.0004708  23.896   <2e-16 ***
#Dissolved_OrganicMatter_RFU 0.0004269  0.0004659   0.916    0.365

aug21.div.glm.fit5<-glm(formula = Bac_Shannon_Diversity ~ Sulfate_milliM, data=aug21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(aug21.div.glm.fit5)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)     0.0112325  0.0004772  23.536   <2e-16 ***
#Sulfate_milliM -0.0002664  0.0004893  -0.545    0.589

aug21.div.glm.fit6<-glm(formula = Bac_Shannon_Diversity ~ Sulfide_microM, data=aug21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(aug21.div.glm.fit6)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)    0.0112379  0.0004723   23.80   <2e-16 ***
#Sulfide_microM 0.0002944  0.0005160    0.57    0.572

fit1<-aov(Bac_Shannon_Diversity ~ as.factor(Depth_m), data=aug21.div)
#pairwise.adonis(aug21.div$Bac_Shannon_Diversity, aug21.div$Depth_m, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(fit1)
#Df           Sum Sq Mean Sq    F value   Pr(>F)
#Depth_m      7   4097   585.3   1.114   0.38
#Residuals   31  16294   525.6
Tuk1<-TukeyHSD(fit1)
Tuk1$Depth_m

#plot(Bac_Shannon_Diversity ~ Depth_m, data=aug21.div)
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
fligner.test(Bac_Shannon_Diversity ~ Depth_m, data = aug21.div)
# Fligner-Killeen:med chi-squared = 4.091, df = 7, p-value = 0.7692
# Which shows that the data do not deviate significantly from homogeneity.
compare_means(Bac_Shannon_Diversity ~ Depth_m, data=aug21.div, method="anova",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input


#### August - Species Richness ####
# August 2021
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
# December 2021
dec21.div.glm.fit1<-glm(formula = Bac_Shannon_Diversity ~ DO_Percent_Local, data=dec21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(dec21.div.glm.fit1)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)       0.0112706  0.0004664  24.166   <2e-16 ***
#DO_Percent_Local -0.0009547  0.0005092  -1.875   0.0687 .

dec21.div.glm.fit2<-glm(formula = Bac_Shannon_Diversity ~ ORP_mV, data=dec21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(dec21.div.glm.fit2)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)   0.0112343  0.0004731  23.745   <2e-16 ***
#ORP_mV      -0.0001512  0.0004968  -0.304    0.763

dec21.div.glm.fit3<-glm(formula = Bac_Shannon_Diversity ~ Temp_DegC, data=dec21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(dec21.div.glm.fit3)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept) 0.0113225  0.0004451  25.439   <2e-16 ***
#Temp_DegC   0.0011947  0.0004861   2.458   0.0188 *

dec21.div.glm.fit4<-glm(formula = Bac_Shannon_Diversity ~ Dissolved_OrganicMatter_RFU, data=dec21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(dec21.div.glm.fit4)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)                 0.0112493  0.0004708  23.896   <2e-16 ***
#Dissolved_OrganicMatter_RFU 0.0004269  0.0004659   0.916    0.365

dec21.div.glm.fit5<-glm(formula = Bac_Shannon_Diversity ~ Sulfate_milliM, data=dec21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(dec21.div.glm.fit5)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)     0.0112325  0.0004772  23.536   <2e-16 ***
#Sulfate_milliM -0.0002664  0.0004893  -0.545    0.589

dec21.div.glm.fit6<-glm(formula = Bac_Shannon_Diversity ~ Sulfide_microM, data=dec21.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(dec21.div.glm.fit6)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)    0.0112379  0.0004723   23.80   <2e-16 ***
#Sulfide_microM 0.0002944  0.0005160    0.57    0.572

fit1<-aov(Bac_Shannon_Diversity ~ as.factor(Depth_m), data=dec21.div)
#pairwise.adonis(dec21.div$Bac_Shannon_Diversity, dec21.div$Depth_m, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(fit1)
#Df           Sum Sq Mean Sq    F value   Pr(>F)
#Depth_m      7   4097   585.3   1.114   0.38
#Residuals   31  16294   525.6
Tuk1<-TukeyHSD(fit1)
Tuk1$Depth_m

#plot(Bac_Shannon_Diversity ~ Depth_m, data=dec21.div)
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
fligner.test(Bac_Shannon_Diversity ~ Depth_m, data = dec21.div)
# Fligner-Killeen:med chi-squared = 4.091, df = 7, p-value = 0.7692
# Which shows that the data do not deviate significantly from homogeneity.
compare_means(Bac_Shannon_Diversity ~ Depth_m, data=dec21.div, method="anova",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input

#### December - Species Richness ####
# December 2021
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
# April 2022
apr22.div.glm.fit1<-glm(formula = Bac_Shannon_Diversity ~ DO_Percent_Local, data=apr22.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(apr22.div.glm.fit1)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)       0.0112706  0.0004664  24.166   <2e-16 ***
#DO_Percent_Local -0.0009547  0.0005092  -1.875   0.0687 .

apr22.div.glm.fit2<-glm(formula = Bac_Shannon_Diversity ~ ORP_mV, data=apr22.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(apr22.div.glm.fit2)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)   0.0112343  0.0004731  23.745   <2e-16 ***
#ORP_mV      -0.0001512  0.0004968  -0.304    0.763

apr22.div.glm.fit3<-glm(formula = Bac_Shannon_Diversity ~ Temp_DegC, data=apr22.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(apr22.div.glm.fit3)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept) 0.0113225  0.0004451  25.439   <2e-16 ***
#Temp_DegC   0.0011947  0.0004861   2.458   0.0188 *

apr22.div.glm.fit4<-glm(formula = Bac_Shannon_Diversity ~ Dissolved_OrganicMatter_RFU, data=apr22.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(apr22.div.glm.fit4)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)                 0.0112493  0.0004708  23.896   <2e-16 ***
#Dissolved_OrganicMatter_RFU 0.0004269  0.0004659   0.916    0.365

apr22.div.glm.fit5<-glm(formula = Bac_Shannon_Diversity ~ Sulfate_milliM, data=apr22.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(apr22.div.glm.fit5)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)     0.0112325  0.0004772  23.536   <2e-16 ***
#Sulfate_milliM -0.0002664  0.0004893  -0.545    0.589

apr22.div.glm.fit6<-glm(formula = Bac_Shannon_Diversity ~ Sulfide_microM, data=apr22.div)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(apr22.div.glm.fit6)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)    0.0112379  0.0004723   23.80   <2e-16 ***
#Sulfide_microM 0.0002944  0.0005160    0.57    0.572

fit1<-aov(Bac_Shannon_Diversity ~ as.factor(Depth_m), data=apr22.div)
#pairwise.adonis(apr22.div$Bac_Shannon_Diversity, apr22.div$Depth_m, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(fit1)
#Df           Sum Sq Mean Sq    F value   Pr(>F)
#Depth_m      7   4097   585.3   1.114   0.38
#Residuals   31  16294   525.6
Tuk1<-TukeyHSD(fit1)
Tuk1$Depth_m

#plot(Bac_Shannon_Diversity ~ Depth_m, data=apr22.div)
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
fligner.test(Bac_Shannon_Diversity ~ Depth_m, data = apr22.div)
# Fligner-Killeen:med chi-squared = 4.091, df = 7, p-value = 0.7692
# Which shows that the data do not deviate significantly from homogeneity.
compare_means(Bac_Shannon_Diversity ~ Depth_m, data=apr22.div, method="anova",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input



#### April - Species Richness ####
# April 2022
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

ggplot(bac.div.metadat, aes(x = DO_Percent_Local, y = Bac_Shannon_Diversity)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampDate), size=3) + theme_classic() +
  stat_smooth(method = "glm", col = "black", se=FALSE, size=1)+ labs(title="Dissolved Oxygen x 16S Shannon Diversity", color="Depth (m)")+ylab("Shannon Diversity")+xlab("Dissolved Oxygen (%)")+
  scale_colour_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) +
  scale_shape_discrete(labels=c("June 2021","August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  stat_cor(label.y = 150, label.x=3) +
  stat_regline_equation(aes(label=paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),label.y = 160,label.x=3)

ggplot(bac.div.metadat, aes(x = DO_Percent_Local, y = Bac_Shannon_Diversity)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampDate), size=3) + theme_classic() +
  labs(title="Dissolved Oxygen x 16S Shannon Diversity", color="Depth (m)")+ylab("Shannon Diversity")+xlab("Dissolved Oxygen (%)")+
  scale_colour_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) +
  scale_shape_discrete(labels=c("June 2021","August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))

ggplot(bac.div.metadat, aes(x = ORP_mV, y = Bac_Shannon_Diversity)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampDate), size=3) + theme_classic() +
  stat_smooth(method = "glm", col = "black", se=FALSE, size=1)+ labs(title="Oxidation-Reduction Potential x 16S Shannon Diversity", color="Depth (m)")+ylab("Shannon Diversity")+xlab("Redox Potential (mV)")+
  scale_colour_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) +
  scale_shape_discrete(labels=c("June 2021","August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  stat_cor(label.y = 3, label.x=1) +
  stat_regline_equation(aes(label=paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),label.y = 3.2,label.x=1)

ggplot(bac.div.metadat, aes(x = Temp_DegC, y = Bac_Shannon_Diversity)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampDate), size=3) + theme_classic() + labs(title="Temperature x 16S Shannon Diversity", color="Depth (m)")+ylab("16S Shannon Diversity")+xlab("Temperature (C)")+
  scale_colour_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) +
  scale_shape_discrete(labels=c("June 2021","August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))

## Species Richness & Environmental Variables
# note: R (correlation coefficient) vs R^2 (coefficient of determination): https://towardsdatascience.com/r%C2%B2-or-r%C2%B2-when-to-use-what-4968eee68ed3

ggplot(bac.div.metadat, aes(x = DO_Percent_Local, y = Bac_Species_Richness)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampDate), size=3) + theme_classic() +
  stat_smooth(method = "glm", col = "black", se=FALSE, size=1)+ labs(title="Dissolved Oxygen x 16S Species Richness", color="Depth (m)")+ylab("Species Richness")+xlab("Dissolved Oxygen (%)")+
  scale_colour_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) +
  scale_shape_discrete(labels=c("June 2021","August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  stat_cor(label.y = 3, label.x=1) +
  stat_regline_equation(aes(label=paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),label.y = 3.1,label.x=1)

ggplot(bac.div.metadat, aes(x = ORP_mV, y = Bac_Species_Richness)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampDate), size=3) + theme_classic() +
  stat_smooth(method = "glm", col = "black", se=FALSE, size=1)+ labs(title="Oxidation-Reduction Potential x 16S Species Richness", color="Depth (m)")+ylab("Species Richness")+xlab("Redox Potential (mV)")+
  scale_colour_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) +
  scale_shape_discrete(labels=c("June 2021","August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  stat_cor(label.y = 3, label.x=1) +
  stat_regline_equation(aes(label=paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),label.y = 3.2,label.x=1)

ggplot(bac.div.metadat, aes(x = Temp_DegC, y = Bac_Species_Richness)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampDate), size=3) + theme_classic() + labs(title="Temperature x 16S Species Richness", color="Depth (m)")+ylab("16S Species Richness")+xlab("Temperature (C)")+
  scale_colour_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) +
  scale_shape_discrete(labels=c("June 2021","August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))


#### Save Everything ####
save.image("data/SSeawater_AlphaDiv_Data.Rdata")
