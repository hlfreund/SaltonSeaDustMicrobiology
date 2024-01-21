#### Set WD & Load Libraries ####
getwd() # use setwd("path/to/files") if you are not in the right directory
setwd("/Volumes/HLF_SSD/Aronson_Lab_Data/Salton_Sea/SaltonSeaDust")
suppressPackageStartupMessages({ # load packages quietly
  library(devtools)
  library(phyloseq)
  library(ggplot2)
  library(vegan)
  library(lme4)
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
  library(RColorBrewer)
})

#### Load Global Env to Import Count/ASV Tables ####
load("data/Amplicon/SSDust_16S.V3V4_W23_Data_Ready.Rdata") # save global env to Rdata file

head(b.dust.all)
bac.ASV_table[1:4,1:4]
bac.ASV_table[(nrow(bac.ASV_table)-4):(nrow(bac.ASV_table)),(ncol(bac.ASV_table)-4):(ncol(bac.ASV_table))] # last 4 rows & cols
head(dust_meta)
head(SurfTypFreq)

#### Using Shapiro-Wilk test for Normality ####

shapiro.test(SurfTypFreq$BarrenLand) # what is the p-value?
# my p-value was p-value =  0.0007935
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(SurfTypFreq$BarrenLand, col="blue")

# visualize Q-Q plot for species richness
qqnorm(meta_scaled$DO_Percent_Local, pch = 1, frame = FALSE)
qqline(meta_scaled$DO_Percent_Local, col = "steelblue", lwd = 2)

shapiro.test(meta_scaled$ORP_mV) # what is the p-value? p-value = 3.323e-12
hist(meta_scaled$ORP_mV, col="blue")
# visualize Q-Q plot for species richness
qqnorm(meta_scaled$ORP_mV, pch = 1, frame = FALSE)
qqline(meta_scaled$ORP_mV, col = "steelblue", lwd = 2)

shapiro.test(meta_scaled$Temp_DegC) # what is the p-value? p-value = 3.562e-06
hist(meta_scaled$Temp_DegC, col="blue")
# visualize Q-Q plot for species richness
qqnorm(meta_scaled$Temp_DegC, pch = 1, frame = FALSE)
qqline(meta_scaled$Temp_DegC, col = "steelblue", lwd = 2)

shapiro.test(meta_scaled$Dissolved_OrganicMatter_RFU) # what is the p-value? p-value = 1.997e-07
hist(meta_scaled$Dissolved_OrganicMatter_RFU, col="blue")
# visualize Q-Q plot for species richness
qqnorm(meta_scaled$Dissolved_OrganicMatter_RFU, pch = 1, frame = FALSE)
qqline(meta_scaled$Dissolved_OrganicMatter_RFU, col = "steelblue", lwd = 2)

shapiro.test(meta_scaled$Sulfate_milliM) # what is the p-value?
# my p-value was p-value =  0.006965
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(meta_scaled$Sulfate_milliM, col="blue")

# visualize Q-Q plot for species richness
qqnorm(meta_scaled$Sulfate_milliM, pch = 1, frame = FALSE)
qqline(meta_scaled$Sulfate_milliM, col = "steelblue", lwd = 2)

shapiro.test(meta_scaled$Sulfide_microM) # what is the p-value?
# my p-value was p-value =  5.934e-12
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(meta_scaled$Sulfide_microM, col="blue")

# visualize Q-Q plot for species richness
qqnorm(meta_scaled$Sulfide_microM, pch = 1, frame = FALSE)
qqline(meta_scaled$Sulfide_microM, col = "steelblue", lwd = 2)

shapiro.test(meta_scaled$Turbidity_FNU) # what is the p-value?  p-value = 0.0005629
hist(meta_scaled$Turbidity_FNU, col="blue")
# visualize Q-Q plot for species richness
qqnorm(meta_scaled$Turbidity_FNU, pch = 1, frame = FALSE)
qqline(meta_scaled$Turbidity_FNU, col = "steelblue", lwd = 2)

shapiro.test(meta_scaled$Chlorophyll_RFU) # what is the p-value? p-value = 1.044e-11
hist(meta_scaled$Chlorophyll_RFU, col="blue")
# visualize Q-Q plot for species richness
qqnorm(meta_scaled$Chlorophyll_RFU, pch = 1, frame = FALSE)
qqline(meta_scaled$Chlorophyll_RFU, col = "steelblue", lwd = 2)

#### PCA w/ STF Variables ####
head(SurfTypFreq)

STFs<-SurfTypFreq[,c(3:12)]
head(STFs)

# NOTE: PCA requires normally distributed data, so we are clr transforming STF data before scaling
# this is because scaled STF data alone is not normally distributed
# we can normalize the data using a clr transformation
STF.clr<-decostand(STFs,method = "clr", pseudocount = 1) # clr transformation of STF data before scaling
STF.clr[1:4,1:4]

# check rownames of clr transformed STF data & metadata
rownames(STF.clr) %in% rownames(dust.meta.surf)
dust.meta.surf=dust.meta.surf[rownames(STF.clr),] ## reorder metadata to match order of clr data

# calculate our Euclidean distance matrix using clr data
STF.euc_dist <- dist(STF.clr, method = "euclidean")

# creating our hierarcical clustering dendrogram
STF.euc_clust <- hclust(STF.euc_dist, method="ward.D2")

# let's make it a little nicer...
STF.euc_dend <- as.dendrogram(STF.euc_clust, hang=0.2)
STF.dend_cols <- as.character(dust_meta$SampDate_Color[order.dendrogram(STF.euc_dend)])
labels_colors(STF.euc_dend) <- STF.dend_cols

plot(STF.euc_dend, ylab="clr Euclidean Distance", horiz=TRUE,cex=03) + title(main = "Surface Type Frequencies Clustering Dendrogram", font.main= 1, cex.sub = 0.8, font.sub = 3)
#legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c( "#ef781c","#03045e","#32cbff","#059c3f"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# let's use our Euclidean distance matrix from before
STF.pca <- prcomp(STF.clr, center=TRUE, scale=TRUE) # pca of euclidean distance matrix = PCA of euclidean distance matrix
STF.pca$x # where sites fall on PC axes
STF.pca$rotation # variables on PC axes
summary(STF.pca)$importance
# The proportion of variances

# extract principal coordinates
STF.pca.vectors<-data.frame(STF.pca$x)
STF.pca.vectors$SampleID<-rownames(STF.pca$x)

# merge pca coordinates w/ metadata
STF.pca.meta<-merge(STF.pca.vectors, dust.meta.surf, by.x="SampleID", by.y="SampleID")
STF.pca.meta$STF_Date

head(STF.pca.meta)
summary(STF.pca)$importance # percentage of variation explained for pca below

# create pca ggplot fig
pca1a<-ggplot(STF.pca.meta, aes(x=PC1, y=PC2)) +geom_point(aes(color=factor(SampDate),shape=Site), size=5)+theme_bw()+
  labs(title="PCA: Surface Type Frequencies in Salton Sea Dust",subtitle="Using CLR Transformed Surface Type Frequencies",color="Collection Date",shape="Site")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Date",values=unique(STF.pca.meta$SampDate_Color[order(STF.pca.meta$SampDate)]),
                     labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  scale_shape_manual(values = c(7,8, 15,16,17,18)) +
  xlab("PC1 [63.78%]") + ylab("PC2 [17.03%]")

ggsave(pca1a,filename = "figures/SurfaceTypeFrequencies/SSD_CLR_STFs_PCA.png", width=15, height=10, dpi=600)

pca1b<-ggplot(STF.pca.meta, aes(x=PC1, y=PC2)) +geom_point(aes(color=factor(SampDate),shape=Site,size=CollectionYear))+theme_bw()+
  labs(title="PCA: Surface Type Frequencies in Salton Sea Dust",subtitle="Using CLR Transformed Surface Type Frequencies",color="Collection Date",shape="Site",size="Collection Year")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  scale_color_manual(name ="Collection Date",values=unique(STF.pca.meta$SampDate_Color[order(STF.pca.meta$SampDate)]),
                     labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  scale_shape_manual(values = c(7,8, 15,16,17,18)) +
  xlab("PC1 [63.78%]") + ylab("PC2 [17.03%]") + scale_size_manual(values = c("2020" = 7, "2021"=4)) +   guides(shape = guide_legend(override.aes = list(size = 5)))

ggsave(pca1b,filename = "figures/SurfaceTypeFrequencies/SSD_CLR_STFs_PCA_v2.png", width=15, height=10, dpi=600)

pca1c<-ggplot(STF.pca.meta, aes(x=PC1, y=PC2)) +geom_point(aes(color=factor(STF_Date),shape=Site,size=CollectionYear))+theme_bw()+
  labs(title="PCA: Surface Type Frequencies in Salton Sea Dust",subtitle="Using CLR Transformed Surface Type Frequencies",color="Surface Date",shape="Site",size="Collection Year")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  scale_shape_manual(values = c(7,8, 15,16,17,18)) + scale_color_manual(name ="Date Ranges",values=brewer.pal(n = 7, name = "Dark2"), labels=c("5/13/20 - 7/10/20", "7/10/20 - 8/30/20", "8/30/20 - 10/10/20","10/10/20 - 11/6/20",
                                                                                                     "6/5/21 - 8/19/21", "8/19/21 - 10/1/21", "10/1/21 - 12/8/21")) +
  xlab("PC1 [63.78%]") + ylab("PC2 [17.03%]") + scale_size_manual(values = c("2020" = 7, "2021"=4)) +   guides(shape = guide_legend(override.aes = list(size = 5)))

ggsave(pca1c,filename = "figures/SurfaceTypeFrequencies/SSD_CLR_STFs_PCA_v3.png", width=15, height=10, dpi=600)

# sample month shape, depth color
# pca2<-ggplot(STF.pca.meta, aes(x=PC1, y=PC2)) +
#   geom_point(aes(color=as.numeric(Depth_m),shape=SampleMonth), size=5)+theme_bw()+
#   labs(title="PCA: STFironmental Variables in Salton Seawater",subtitle="Using clr Transformed Data",xlab="Axis 1", ylab="Axis 2",color="Depth (m)")+
#   theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
#   scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
#   xlab("Axis 1 [63.78%]") + ylab("Axis 2 [17.03%]")
#
# ggsave(pca2,filename = "figures/STFVariablesOnly/SSD_clrSTFOnly_PCA_Depth_SampDate.png", width=12, height=10, dpi=600)

#### Compare Env Samples Variance by Sample Date ####
# Kruskal-Wallis test = nonparametric one-way ANOVA
kruskal.test(DO_Percent_Local ~ SampleMonth, data = meta_scaled)
pairwise.wilcox.test(meta_scaled$DO_Percent_Local, meta_scaled$SampleMonth, p.adjust.method = "bonf") # returns p values
#         August  December
#December 0.00029 -
#April    0.00704 0.00051

kruskal.test(ORP_mV ~ SampleMonth, data = meta_scaled)
pairwise.wilcox.test(meta_scaled$ORP_mV, meta_scaled$SampleMonth, p.adjust.method = "bonf") # returns p values
#         August  December
#December 0.00096 -
#April    0.00027 4.2e-06

kruskal.test(Temp_DegC ~ SampleMonth, data = meta_scaled)
pairwise.wilcox.test(meta_scaled$Temp_DegC, meta_scaled$SampleMonth, p.adjust.method = "bonf") # returns p values
#         August  December
#December 0.00029 -
#April  0.00029 4.3e-06


kruskal.test(Dissolved_OrganicMatter_RFU ~ SampleMonth, data = meta_scaled)
pairwise.wilcox.test(meta_scaled$Dissolved_OrganicMatter_RFU, meta_scaled$SampleMonth, p.adjust.method = "bonf") # returns p values
#         August  December
#December 0.00028 -
#April    0.00029 4.2e-06

kruskal.test(Sulfate_milliM ~ SampleMonth, data = meta_scaled)
pairwise.wilcox.test(meta_scaled$Sulfate_milliM, meta_scaled$SampleMonth, p.adjust.method = "bonf") # returns p values
#         August  December
#December 0.00029 -
#April   1.00000 4.5e-06

kruskal.test(Sulfide_microM ~ SampleMonth, data = meta_scaled)
pairwise.wilcox.test(meta_scaled$Sulfide_microM, meta_scaled$SampleMonth, p.adjust.method = "bonf") # returns p values
#         August  December
#December 0.00705 -
#April    0.44826 0.00015

#### Do Env Variables Correlate? ####
head(meta_scaled)
# check for colinearity among env variables themselves
heatmap(abs(cor(meta_scaled[,c(8,10:12,15:17)])),
        # Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(6)),
        Colv = NA, Rowv = NA)
legend("topleft",
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(6)))
dev.off()

# Visualize with a corrplot
cor_mat.env1 <- cor(meta_scaled[,c(8,10:12,15:17)], method='pearson')
cor_mat.env1

symnum(cor_mat.env1)

png('figures/EnvVariablesOnly/SSD_ScaledCentered_EnvVarOnly_CorrPlot.png')
crrplt1<-corrplot.mixed(cor_mat.env1, tl.pos='lt', tl.cex=0.7, sig.level = 0.05, number.cex=0.8,
               diag='l',lower.col = 'black', cl.ratio = 0.2, tl.srt = 45)
# env variables with a correlation of <|0.7| is a good threshold for determining if predictors correlate
dev.off()

# DO %
cor.test(meta_scaled$DO_Percent_Local, meta_scaled$ORP_mV, method="pearson") # ***
# r = 0.5619567, p = 0.000161 --> medium corr, & it's significant
cor.test(meta_scaled$DO_Percent_Local, meta_scaled$Temp_DegC, method="pearson") # ***
# r = -0.7973317, p-value = 7.386e-10 --> strong negative correlation, significant
cor.test(meta_scaled$DO_Percent_Local, meta_scaled$Dissolved_OrganicMatter_RFU, method="pearson") # ***
# r = -0.6131097, p-value = 2.6e-05 --> medium correlation & significant
cor.test(meta_scaled$DO_Percent_Local, meta_scaled$Depth.num, method="pearson") #
# r = -0.3850655 , p = 0.01414 --> not strong negative correlation,significant
cor.test(meta_scaled$DO_Percent_Local, meta_scaled$Salinity_ppt, method="pearson") # ****
# r = 0.8283734, p-value = 4.192e-11 --> strong positive correlation & significant

plot(x=meta_scaled$DO_Percent_Local, y=meta_scaled$Dissolved_OrganicMatter_RFU, col=meta_scaled$SCY_Color)
plot(x=meta_scaled$DO_Percent_Local, y=meta_scaled$ORP_mV, col=meta_scaled$SCY_Color)
plot(x=meta_scaled$DO_Percent_Local, y=meta_scaled$Temp_DegC, col=meta_scaled$SCY_Color)
plot(x=meta_scaled$DO_Percent_Local, y=meta_scaled$Depth.num, col=meta_scaled$SCY_Color)
plot(x=meta_scaled$DO_Percent_Local, y=meta_scaled$Salinity_ppt, col=meta_scaled$SCY_Color)

# ORP
cor.test(meta_scaled$ORP_mV, meta_scaled$Temp_DegC, method="pearson") # ***
# r = -0.513023, p-value = 0.0007117 --> medium negative correlation & significant
cor.test(meta_scaled$ORP_mV, meta_scaled$Dissolved_OrganicMatter_RFU, method="pearson") # ***
# r = -0.380499, p-value = 0.00833 --> lame negative correlation, it's significant
cor.test(meta_scaled$ORP_mV, meta_scaled$Depth.num, method="pearson")
# r = -0.2568054, p-value =  0.1097 --> not sig, not strong corr
cor.test(meta_scaled$ORP_mV, meta_scaled$Salinity_ppt, method="pearson") # ***
# r = 0.5235966, p-value =0.0005261 --> medium corr, significant

plot(x=meta_scaled$ORP_mV, y=meta_scaled$Dissolved_OrganicMatter_RFU, col=meta_scaled$SampDate_Color)
plot(x=meta_scaled$ORP_mV, y=meta_scaled$Temp_DegC, col=meta_scaled$SampDate_Color)
plot(x=meta_scaled$ORP_mV, y=meta_scaled$Depth.num, col=meta_scaled$SampDate_Color)
plot(x=meta_scaled$ORP_mV, y=meta_scaled$Salinity_ppt, col=meta_scaled$SampDate_Color)

# Dissolved Organic Matter
cor.test(meta_scaled$Dissolved_OrganicMatter_RFU, meta_scaled$Temp_DegC, method="pearson") # ***
# r = 0.4585448, p-value = 0.002923 # medium corr, significant
cor.test(meta_scaled$Dissolved_OrganicMatter_RFU, meta_scaled$Depth.num, method="pearson")
# r = 0.3099062 , p = 0.05165 # not sig, lame corr
cor.test(meta_scaled$Dissolved_OrganicMatter_RFU, meta_scaled$Salinity_ppt, method="pearson")
# r = -0.6411353, p-value = 8.303e-06

plot(x=meta_scaled$Dissolved_OrganicMatter_RFU, y=meta_scaled$Temp_DegC, col=meta_scaled$SampDate_Color)
plot(x=meta_scaled$Dissolved_OrganicMatter_RFU, y=meta_scaled$Depth.num, col=meta_scaled$SampDate_Color)
plot(x=meta_scaled$Dissolved_OrganicMatter_RFU, y=meta_scaled$Salinity_ppt, col=meta_scaled$SampDate_Color)

# Sulfate (milliM)
cor.test(meta_scaled$Sulfate_milliM, meta_scaled$ORP_mV, method="pearson")
# r = 0.1102182, p = 0.4984 --> not a strong correlation & it's significant
cor.test(meta_scaled$Sulfate_milliM, meta_scaled$Temp_DegC, method="pearson") # ***
# r = -0.5614189, p = 0.0001639 --> strong negative correlation, significant
cor.test(meta_scaled$Sulfate_milliM, meta_scaled$DO_Percent_Local, method="pearson") # ****
# r = 0.6452846, p = 6.944e-06 --> strong correlation & significant
cor.test(meta_scaled$Sulfate_milliM, meta_scaled$Dissolved_OrganicMatter_RFU, method="pearson")
# r = -0.03455734 , p = 0.8323 --> no corr, not sig
cor.test(meta_scaled$Sulfate_milliM, meta_scaled$Sulfide_microM, method="pearson")
# r = -0.1571489 , p = 0.3328 --> no corr, not sig
cor.test(meta_scaled$Sulfate_milliM, meta_scaled$Salinity_ppt, method="pearson")
# r = 0.515623, p-value = 0.0006613
cor.test(meta_scaled$Sulfate_milliM, meta_scaled$Depth.num, method="pearson")
# r = -0.1800819, p-value = 0.2662

plot(x=meta_scaled$Sulfate_milliM, y=meta_scaled$DO_Percent_Local, col=meta_scaled$SampDate_Color)
plot(x=meta_scaled$Sulfate_milliM, y=meta_scaled$Dissolved_OrganicMatter_RFU, col=meta_scaled$SampDate_Color)
plot(x=meta_scaled$Sulfate_milliM, y=meta_scaled$ORP_mV, col=meta_scaled$SampDate_Color)
plot(x=meta_scaled$Sulfate_milliM, y=meta_scaled$Temp_DegC, col=meta_scaled$SampDate_Color)
plot(x=meta_scaled$Sulfate_milliM, y=meta_scaled$Depth.num, col=meta_scaled$SampDate_Color)
plot(x=meta_scaled$Sulfate_milliM, y=meta_scaled$Salinity_ppt, col=meta_scaled$SampDate_Color)

# Sulfide (microM)
cor.test(meta_scaled$Sulfide_microM, meta_scaled$ORP_mV, method="pearson") # ******
# r = -0.979198, p < 2.2e-16 --> STRONG & significant negative correlation
cor.test(meta_scaled$Sulfide_microM, meta_scaled$Temp_DegC, method="pearson") # ***
# r = 0.5532523 , p = 0.0002134 --> medium correlation, significant
cor.test(meta_scaled$Sulfide_microM, meta_scaled$DO_Percent_Local, method="pearson") # ****
# r = -0.6286855, p = 1.398e-05 --> medium-strong negative correlation & significant
cor.test(meta_scaled$Sulfide_microM, meta_scaled$Dissolved_OrganicMatter_RFU, method="pearson") # ****
# r = 0.621356 , p = 1.88e-05 --> medium to strong correlation, significant
cor.test(meta_scaled$Sulfide_microM, meta_scaled$Salinity_ppt, method="pearson")
# r = -0.5595043, p-value = 0.0001745 # medium neg corr, significant
cor.test(meta_scaled$Sulfide_microM, meta_scaled$Depth.num, method="pearson")
# r = 0.2837005 , p-value = 0.07606 # not sig, no corr

plot(x=meta_scaled$Sulfide_microM, y=meta_scaled$DO_Percent_Local, col=meta_scaled$SampDate_Color)
plot(x=meta_scaled$Sulfide_microM, y=meta_scaled$Dissolved_OrganicMatter_RFU, col=meta_scaled$SampDate_Color)
plot(x=meta_scaled$Sulfide_microM, y=meta_scaled$ORP_mV, col=meta_scaled$SampDate_Color)
plot(x=meta_scaled$Sulfide_microM, y=meta_scaled$Temp_DegC, col=meta_scaled$SampDate_Color)
plot(x=meta_scaled$Sulfide_microM, y=meta_scaled$Depth.num, col=meta_scaled$SampDate_Color)
plot(x=meta_scaled$Sulfide_microM, y=meta_scaled$Salinity_ppt, col=meta_scaled$SampDate_Color)

# Chlorophyll
#cor.test(meta_scaled$Chlorophyll_RFU, meta_scaled$Dissolved_OrganicMatter_RFU, method="pearson") # ****
# r = -0.70582, p-value = 3.002e-08 --> strong correlation & significant
#cor.test(meta_scaled$Chlorophyll_RFU, meta_scaled$Temp_DegC, method="pearson")

#### Do Env Data Vary Significantly By Group?#####


salf1<-aov(Salinity_ppt ~ SampDate, data=meta_scaled)
#pairwise.adonis(bac.div.metadat$Bac_Species_Richness, bac.div.metadat$Depth_m, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(salf1)
#Df           Sum Sq Mean Sq    F value   Pr(>F)
#SampDate     2 22.887  11.444    2128 <2e-16 ***
#Residuals   21  0.113   0.005
Tuk1<-TukeyHSD(salf1)
Tuk1$SampDate
#                               diff        lwr        upr        p adj
# December.2021-August.2021  2.2667637  2.1743447  2.3591827 4.363176e-14
# April.2022-August.2021     1.7949426  1.7025236  1.8873615 4.363176e-14
# April.2022-December.2021  -0.4718212 -0.5642402 -0.3794022 5.833745e-11

salf2<-aov(Salinity_ppt ~ Depth_m, data=meta_scaled)
#pairwise.adonis(bac.div.metadat$Bac_Species_Richness, bac.div.metadat$Depth_m, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(salf2)
#Df           Sum Sq Mean Sq    F value   Pr(>F)
#Depth_m      7  0.013  0.0019   0.001      1
#Residuals   16 22.987  1.4367
Tuk2<-TukeyHSD(salf2)
Tuk2$Depth_m

domf1<-aov(Dissolved_OrganicMatter_RFU ~ SampDate, data=meta_scaled)
#pairwise.adonis(bac.div.metadat$Bac_Species_Richness, bac.div.metadat$Depth_m, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(domf1)
#Df           Sum Sq Mean Sq    F value   Pr(>F)
#SampDate     2  17.87   8.935   36.58 1.44e-07 ***
#Residuals   21   5.13   0.244
Tuk3<-TukeyHSD(domf1)
Tuk3$SampDate
#                             diff       lwr        upr        p adj
# December.2021-August.2021 -1.3249869 -1.947862 -0.7021123 7.340349e-05
# April.2022-August.2021    -2.0886802 -2.711555 -1.4658055 9.867955e-08
# April.2022-December.2021  -0.7636932 -1.386568 -0.1408186 1.468732e-02
## sampling depth

domf2<-aov(Dissolved_OrganicMatter_RFU ~ Depth_m, data=meta_scaled)
#pairwise.adonis(bac.div.metadat$Bac_Species_Richness, bac.div.metadat$Depth_m, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(domf2)
#Df           Sum Sq Mean Sq    F value   Pr(>F)
#Depth_m      7  2.875  0.4107   0.326  0.931
#Residuals   16 20.125  1.2578
Tuk4<-TukeyHSD(domf2)
Tuk4$Depth_m
## sampling depth

#### Does Env Predict Alpha Diversity? ####
## linear regression time!

## First let's make the alpha diversity data frame
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

# drop the outliers
bac.div.metadat2<-subset(bac.div.metadat, bac.div.metadat$Bac_Shannon_Diversity<=200)

# Linear Regression time!
## here the focus is comparing dust complexity to alpha diversity, species richness, & elevation
head(bac.div.metadat2)
s.div.lm.fit1<-lm(Bac_Shannon_Diversity ~ DO_Percent_Local, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.div.lm.fit1)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)        97.679      6.548  14.916   <2e-16 ***
#  DO_Percent_Local   -1.922      6.647  -0.289    0.774

s.div.lm.fit2<-lm(Bac_Shannon_Diversity ~ ORP_mV, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.div.lm.fit2)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)   97.674      3.709  26.333  < 2e-16 ***
#  ORP_mV       -15.535      3.720  -4.176 0.000138 ***
## ^^^ the two lms below show that this model is significant only for June & August 2021, not December & April

not_summer_months<-subset(bac.div.metadat2, SampDate=="December.2021" | SampDate=="April.2022" )

s.div.lm.fit2a<-lm(Bac_Shannon_Diversity ~ ORP_mV, data=not_summer_months) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.div.lm.fit2a)

summer_months<-subset(bac.div.metadat2, SampDate=="June.2021" | SampDate=="August.2021" )

s.div.lm.fit2b<-lm(Bac_Shannon_Diversity ~ ORP_mV, data=summer_months) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.div.lm.fit2b)

s.div.lm.fit3<-lm(Bac_Shannon_Diversity ~ Temp_DegC, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.div.lm.fit3)

s.div.lm.fit5<-lm(Bac_Shannon_Diversity ~ Dissolved_OrganicMatter_RFU, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.div.lm.fit5)




fit1<-aov(Bac_Shannon_Diversity ~ Depth_m, data=bac.div.metadat2)
pairwise.adonis(bac.div.metadat2$Bac_Shannon_Diversity, bac.div.metadat2$Depth_m, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(fit1)
#Df           Sum Sq Mean Sq    F value   Pr(>F)
#Elevation2   3 0.7277 0.24258   0.084 0.774
#Residuals   27 0.4444 0.01646
Tuk1<-TukeyHSD(fit1)
Tuk1$Depth_m

#plot(DustComplexity ~ Elevation, data=bac.div.metadat2)
#abline(aov(DustComplexity ~ Elevation, data=bac.div.metadat2))

# fit.0<-aov(DustComplexity ~ as.factor(Elevation), data=bac.div.metadat2)
# summary(fit.0)
# TukeyHSD(fit.0)
# Levene's test with one independent variable
## Levene's tests whether variances of 2 samples are equal
## we want variances to be the same -- want NON SIGNIFICANCE!
## t test assumes that variances are the same, so Levene's test needs to be non significant
fligner.test(Bac_Shannon_Diversity ~ Depth_m, data = bac.div.metadat2)
# Levenes Test for Homogeneity of Variance
#        Df  Chi square value  Pr(>F)
# group  3   1.0952   0.7411
# Which shows that the data do not deviate significantly from homogeneity.
elev<-bac.div.metadat2$Elevation2
compare_means(DustComplexity ~ Elevation2, data=bac.div.metadat2, method="anova",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input

p.adj.dc.elev<-compare_means(DustComplexity ~ Elevation2, data=bac.div.metadat2, method="t.test",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input# Note https://github.com/kassambara/ggpubr/issues/65

fit.test<-ggplot(bac.div.metadat2, aes(x = as.factor(Elevation), y = DustComplexity, fill=as.factor(Elevation))) +
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_manual(values=saturation(fair_cols, 0.9))+stat_compare_means(method = "anova",label.y=1.5) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", hide.ns = TRUE,label = "p.signif")

ggsave(fit.test,filename = "figures/EnvVariablesOnly/DustComp_by_Elevation_ALL_sigbars_5.24.21.pdf", width=10, height=8, dpi=600)

fit.testa<-ggplot(bac.div.metadat2, aes(x = as.factor(Elevation), y = DustComplexity, fill=as.factor(Elevation))) +
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_manual(values=saturation(fair_cols, 0.9))+stat_compare_means(method = "anova",label.y=1.5,mapping=aes(label = format.pval(..p.adj.., digits = 3))) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", hide.ns = TRUE,mapping=aes(label = format.pval(..p.adj.., digits = 3)))

ggsave(fit.testa,filename = "figures/EnvVariablesOnly/DustComp_by_Elevation_ALL_sigbars_5.24.21.pdf", width=10, height=8, dpi=600)

fit.test0<-ggplot(bac.div.metadat2, aes(x = as.factor(Elevation), y = DustComplexity, fill=as.factor(Elevation))) +
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_manual(values=saturation(fair_cols, 0.9))+stat_compare_means(method = "anova",label.y=1.5) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", hide.ns = TRUE,label = "p.signif")

ggsave(fit.test,filename = "figures/EnvVariablesOnly/DustComp_by_Elevation_ALL_sigbars_5.24.21.pdf", width=10, height=8, dpi=600)

fit.testa<-ggplot(bac.div.metadat2, aes(x = as.factor(Elevation), y = DustComplexity, fill=as.factor(Elevation))) +
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_manual(values=saturation(fair_cols, 0.9))

ggsave(fit.testa,filename = "figures/EnvVariablesOnly/DustComp_by_Elevation_ALL_no.sigbars_5.24.21.pdf", width=10, height=8, dpi=600)

fit.testb<-ggplot(bac.div.metadat2, aes(x = as.factor(Elevation), y = DustComplexity, fill=as.factor(Elevation))) +
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_grey(start=0.8, end=0.3)

ggsave(fit.testb,filename = "figures/EnvVariablesOnly/DustComp_by_Elevation_ALL_gray_5.24.21.pdf", width=10, height=8, dpi=600)

fit.testb.0<-ggplot(bac.div.metadat2, aes(x = as.factor(Elevation), y = DustComplexity, fill=as.factor(Elevation))) +
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_grey(start=0.8, end=0.3)+stat_compare_means(method = "anova",label.y=1.5) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", hide.ns = TRUE,label = "p.signif")

ggsave(fit.testb.0,filename = "figures/EnvVariablesOnly/DustComp_by_Elevation_ALL_gray_sigbars_5.24.21.pdf", width=10, height=8, dpi=600)

### Fungi comparisons first
# Dust Comp x ITS1 Shannon diversity
hist(bac.div.metadat2$ITS1_Shannon_Diversity) # NOT normally distributed
hist(bac.div.metadat2$DustComplexity) # somewhat normally distributed
chisq.test(bac.div.metadat2$ITS1_Shannon_Diversity, bac.div.metadat2$DustComplexity)

its1.fit1<-lm(DustComplexity ~ ITS1_Shannon_Diversity, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
summary(its1.fit1)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)             0.614082   0.050858  12.074  7.8e-13 ***
#  ITS1_Shannon_Diversity -0.001942   0.001416  -1.372    0.181

its1.fit1.p<-p.adjust(coef(summary(its1.fit1))[8], method="bonferroni") # pvalue

plot(DustComplexity ~ ITS1_Shannon_Diversity, data=bac.div.metadat2)
abline(its1.fit1)


#leveneTest(bac.div.metadat2$DustComplexity,
#            bac.div.metadat2$ITS1_Shannon_Diversity,
#            location = c("median"),
#            trim.alpha = 0.25)
# Levenes Test for Homogeneity of Variance
#        Df  F value  Pr(>F)
# group  3   2.3415   0.0818
# Which shows that the data do not deviate significantly from homogeneity.

fig.its1.fit1<-ggplot(its1_div_meta, aes(x = ITS1_Shannon_Diversity, y = DustComplexity)) +
  geom_point(aes(color=Elev.num), size=3) + theme_classic() + saturation(scale_colour_gradientn(colours=fair_cols,limits=c(400,2700),breaks = c(500,1250,2000,2600),labels=c("400","1100","2000","2700")), 0.9) +
  stat_smooth(method = "glm", col = "black", se=FALSE, size=1)+ labs(color="Elevation (ft)")+ylab("Dust Complexity")+xlab("ITS1 Shannon Diversity")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  stat_cor(label.y = 1, label.x=75) +
  stat_regline_equation(aes(label=paste(..adj.rr.label..)),label.y = 1.20,label.x=75)
## use summary(its1.fit1) to double check that stat_cor gives same p value as linear regression!

ggsave(fig.its1.fit1,filename = "figures/EnvVariablesOnly/DustComp_by_ITS1_ShanDiv_ALL_1.4.22.pdf", width=10, height=8, dpi=600)

fig.its1.fit1<-ggplot(its1_div_meta, aes(x = ITS1_Shannon_Diversity, y = DustComplexity)) +
  geom_point(aes(color=Elev.num),size=3) + theme_classic() + saturation(scale_colour_gradientn(colours=fair_cols,limits=c(400,2700),breaks = c(500,1250,2000,2600),labels=c("400","1100","2000","2700")), 0.9) +
  stat_smooth(method = "glm", col = "black", se=FALSE, size=1)+ labs(title="Dust Complexity x ITS1 Shannon Diversity", color="Elevation (ft)")+ylab("Dust Complexity")+xlab("ITS1 Shannon Diversity")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  stat_cor(label.y = 1, label.x=75) +
  stat_regline_equation(aes(label=paste(..adj.rr.label..)),label.y = 1.05,label.x=75)
## use summary(its1.fit1) to double check that stat_cor gives same p value as linear regression!

#fig.its1.fit2<-ggplot(bac.div.metadat2, aes(x = ITS1_Shannon_Diversity, y = DustComplexity)) +
#  geom_point(aes(color=Elevation), size=3) + theme_classic() + saturation(scale_colour_gradientn(colours=fair_cols,limits=c(400,2700),breaks = c(500,1250,2000,2600),labels=c("400","1100","2000","2700")), 0.9) +
#  stat_smooth(method = "glm", col = "black", se=FALSE, size=1)+ labs(title="Dust Complexity x ITS1 Shannon Diversity", color="Elevation (ft)")+ylab("Dust Complexity")+xlab("ITS1 Shannon Diversity")+
#  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
#  stat_cor(label.y = 1, label.x=75) +
#  stat_regline_equation(label.y = 1.05,label.x=75)

#ggsave(fig.its1.fit2,filename = "figures/EnvVariablesOnly/DustComp_by_ITS1_Shan_Div_ALL_5.19.21.pdf", width=10, height=8, dpi=600)


# DustComp x ITS1 Species Richness
hist(bac.div.metadat2$ITS1_Species_Richness)
its1.sr.fit1<-lm(DustComplexity ~ ITS1_Species_Richness, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
summary(its1.sr.fit1)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)
# (Intercept)            0.7761611  0.0757511  10.246 3.79e-11 ***
#  ITS1_Species_Richness -0.0004551  0.0001475  -3.084  0.00445 **
coef(summary(its1.sr.fit1)) # pvalue

its1.sr.fit2<-lm(DustComplexity ~ ITS1_Species_Richness, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
summary(its1.sr.fit2)
coef(summary(its1.sr.fit2))
p.adjust(coef(summary(its1.sr.fit2))[,4], method="bonferroni") # pvalue

its1.sr.fit2<-glm(DustComplexity ~ ITS1_Species_Richness, data=bac.div.metadat2, family=poisson)
its1.sr.fit3<-glm.nb(DustComplexity ~ ITS1_Species_Richness, data=bac.div.metadat2)

summary(its1.sr.fit2)
dispersiontest(its1.sr.fit2)
# null hypothesis is that equidispersion exists; alternative hypothesis is overdispersion
# if overdispersion, use negative binomial not Poisson
## Poisson distribution implies that the mean and variance are equal --> little dispersion
# negative binomial means # of observations is not fixed, whereas binomial means observations are a fixed #

# z = -16.609, p-value = 1 (cannot reject null)
# alternative hypothesis: true dispersion is greater than 1
# sample estimates:
#   dispersion
# 0.0495281 -- equidispersion exists

plot(DustComplexity ~ ITS1_Species_Richness, data=bac.div.metadat2)
abline(its1.fit2)

fig.its1.sr.fit1<-ggplot(its1_div_meta, aes(x = ITS1_Species_Richness, y = DustComplexity)) +
  geom_point(aes(color=fair_cols), size=3) + theme_classic() + saturation(scale_colour_gradientn(colours=fair_cols,limits=c(400,2700),breaks = c(500,1250,2000,2600),labels=c("400","1100","2000","2700")), 0.9) +
  stat_smooth(method = "glm", col = "black", se=FALSE, size=1)+ labs(title="Dust Complexity x ITS1 Species Richness", color="Elevation (ft)")+ylab("Dust Complexity")+xlab("ITS1 Species Richness")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  stat_cor(label.y = 1, label.x=700) +
  stat_regline_equation(aes(label=paste(..adj.rr.label..)),label.y = 1.05,label.x=700)

## use summary(its1.sr.fit1) to double check that stat_cor gives same p value as linear regression!

ggsave(fig.its1.sr.fit1,filename = "figures/EnvVariablesOnly/DustComp_by_ITS1_Spec_Richness_ALL_1.4.22.pdf", width=10, height=8, dpi=600)



#### Plots of Env Variables ####

# Compare all variables across Depths
dep.dom<-ggplot(meta_scaled, aes(x=Depth_m, y=Dissolved_OrganicMatter_RFU,color=SampDate,group=SampDate)) + geom_point(size=3) + geom_line() + theme_bw()+
  labs(title="Dissolved Organic Matter (DOM) by Depth & Sample Date",subtitle="Using Scaled DOM RFU Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(meta_scaled$SampDate_Color[order(meta_scaled$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Depth (m)") + ylab("DOM (RFU)")+coord_flip()+ scale_x_discrete(limits=rev)

ggsave(dep.dom,filename = "figures/EnvVariablesOnly/SSD_DOM_Depth_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

dep.orp<-ggplot(meta_scaled, aes(x=Depth_m, y=ORP_mV,color=SampDate,group=SampDate)) + geom_point(size=3) + geom_line() + theme_bw()+
  labs(title="Oxidative-Reduction Potential by Depth & Sample Date",subtitle="Using Scaled ORP (mV) Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(meta_scaled$SampDate_Color[order(meta_scaled$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Depth (m)") + ylab("ORP (mV)")+coord_flip()+ scale_x_discrete(limits=rev)

ggsave(dep.orp,filename = "figures/EnvVariablesOnly/SSD_ORP_Depth_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

dep.sal<-ggplot(meta_scaled, aes(x=Depth_m, y=Salinity_ppt,color=SampDate,group=SampDate)) + geom_point(size=3) + geom_line() + theme_bw()+
  labs(title="Salinity by Depth & Sample Date",subtitle="Using Scaled Salinity (PPT) Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(meta_scaled$SampDate_Color[order(meta_scaled$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Depth (m)") + ylab("Salinity (PPT)")+coord_flip()+ scale_x_discrete(limits=rev)

ggsave(dep.orp,filename = "figures/EnvVariablesOnly/SSD_Salinity_Depth_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

dep.sulf<-ggplot(meta_scaled, aes(x=Depth_m, y=Sulfate_milliM,color=SampDate,group=SampDate)) + geom_point(size=3) + geom_line() + theme_bw()+
  labs(title="Sulfate by Depth & Sample Date",subtitle="Using Scaled Sulfate (milliM) Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(meta_scaled$SampDate_Color[order(meta_scaled$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Depth (m)") + ylab("Sulfate (milliM)")+coord_flip()+ scale_x_discrete(limits=rev)

ggsave(dep.sulf,filename = "figures/EnvVariablesOnly/SSD_Sulfate_Depth_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

dep.do<-ggplot(meta_scaled, aes(x=Depth_m, y=DO_Percent_Local,color=SampDate,group=SampDate)) + geom_point(size=3) + geom_line() + theme_bw()+
  labs(title="Dissolved Oxygen by Depth & Sample Date",subtitle="Using Scaled DO (%) Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(meta_scaled$SampDate_Color[order(meta_scaled$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Depth (m)") + ylab("DO (%)")+coord_flip()+ scale_x_discrete(limits=rev)

ggsave(dep.do,filename = "figures/EnvVariablesOnly/SSD_DO_Depth_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

dep.hs<-ggplot(meta_scaled, aes(x=Depth_m, y=Sulfide_microM,color=SampDate,group=SampDate)) + geom_point(size=3) + geom_line() + theme_bw()+
  labs(title="Sulfide by Depth & Sample Date",subtitle="Using Scaled Sulfate (microM) Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(meta_scaled$SampDate_Color[order(meta_scaled$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Depth (m)") + ylab("Sulfide (microM)")+coord_flip()+ scale_x_discrete(limits=rev)

ggsave(dep.hs,filename = "figures/EnvVariablesOnly/SSD_Sulfide_Depth_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

dep.chlr<-ggplot(meta_scaled, aes(x=Depth_m, y=Chlorophyll_RFU,color=SampDate,group=SampDate)) + geom_point(size=3) + geom_line() + theme_bw()+
  labs(title="Chlorophyll by Depth & Sample Date",subtitle="Using Scaled Chlorophyll (RFU) Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(meta_scaled$SampDate_Color[order(meta_scaled$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Depth (m)") + ylab("Chlorophyll (RFU)")+coord_flip()+ scale_x_discrete(limits=rev)

ggsave(dep.chlr,filename = "figures/EnvVariablesOnly/SSD_Chlorophyll_Depth_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

dep.temp<-ggplot(meta_scaled, aes(x=Depth_m, y=Temp_DegC,color=SampDate,group=SampDate)) + geom_point(size=3) + geom_line() + theme_bw()+
  labs(title="Temperature by Depth & Sample Date",subtitle="Using Scaled Chlorophyll (RFU) Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(meta_scaled$SampDate_Color[order(meta_scaled$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Depth (m)") + ylab("Temperature (C)")+coord_flip()+ scale_x_discrete(limits=rev)

ggsave(dep.temp,filename = "figures/EnvVariablesOnly/SSD_Temp_Depth_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

# Compare variables to each other
dom.hs<-ggplot(meta_scaled, aes(x=Dissolved_OrganicMatter_RFU, y=Sulfide_microM,color=SampDate,group=SampDate)) + geom_point(size=3) + geom_line() + theme_bw()+
  labs(title="Sulfide ~ DOM",subtitle="Using Scaled Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(meta_scaled$SampDate_Color[order(meta_scaled$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("DOM (RFU)") + ylab("Sulfide (microM)")

ggsave(dom.hs,filename = "figures/EnvVariablesOnly/SSD_Sulfide_DOM_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

do.hs<-ggplot(meta_scaled, aes(x=DO_Percent_Local, y=Sulfide_microM,color=SampDate,group=SampDate)) + geom_point(size=3) + geom_line() + theme_bw()+
  labs(title="Sulfide ~ DO%",subtitle="Using Scaled Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(meta_scaled$SampDate_Color[order(meta_scaled$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("DO (%)") + ylab("Sulfide (microM)")

ggsave(do.hs,filename = "figures/EnvVariablesOnly/SSD_Sulfide_DO.Percent_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

orp.hs<-ggplot(meta_scaled, aes(x=ORP_mV, y=Sulfide_microM,color=SampDate,group=SampDate)) + geom_point(size=3) + geom_line() + theme_bw()+
  labs(title="Sulfide ~ ORP (mV)",subtitle="Using Scaled Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(meta_scaled$SampDate_Color[order(meta_scaled$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("ORP (mV") + ylab("Sulfide (microM)")

ggsave(orp.hs,filename = "figures/EnvVariablesOnly/SSD_Sulfide_ORP_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

orp.dom<-ggplot(meta_scaled, aes(x=ORP_mV, y=Dissolved_OrganicMatter_RFU,color=SampDate,group=SampDate)) + geom_point(size=3) + geom_line() + theme_bw()+
  labs(title="DOM (RFU) ~ ORP (mV)",subtitle="Using Scaled Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(meta_scaled$SampDate_Color[order(meta_scaled$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("DOM (RFU") + ylab("Sulfide (microM)")

ggsave(orp.dom,filename = "figures/EnvVariablesOnly/SSD_DOM_ORP_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

orp.do<-ggplot(meta_scaled, aes(x=ORP_mV, y=DO_Percent_Local,color=SampDate,group=SampDate)) + geom_point(size=3) + geom_line() + theme_bw()+
  labs(title="DO% ~ ORP (mV)",subtitle="Using Scaled Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(meta_scaled$SampDate_Color[order(meta_scaled$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("ORP (mV") + ylab("DO %")

ggsave(orp.do,filename = "figures/EnvVariablesOnly/SSD_DO.Percent_ORP_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)
