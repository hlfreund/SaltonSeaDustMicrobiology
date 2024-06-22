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
  library(ggcorrplot)
  library(ggpmisc)
})

#### Load Global Env to Import Count/ASV Tables ####
load("data/Amplicon/SSDust_16S.V3V4_W23_Data_Ready.Rdata") # save global env to Rdata file

head(b.dust.all)
bac.ASV_table[1:4,1:4]
bac.ASV_table[(nrow(bac.ASV_table)-4):(nrow(bac.ASV_table)),(ncol(bac.ASV_table)-4):(ncol(bac.ASV_table))] # last 4 rows & cols
head(meta.all.scaled)
head(SurfTypFreq)

#### Look at Surface Type Frequencies by Sample for Modeling Later ####

# melt down data so that we can make a stacked barplot
STF.melt<-melt(SurfTypFreq[,c(1,3:12,15:16)],by=c("Site","SampleID","STF_Date"))
colnames(STF.melt)[which(names(STF.melt) == "variable")] <- "SurfaceType"
colnames(STF.melt)[which(names(STF.melt) == "value")] <- "Frequency"
head(STF.melt)

# create palette
colorset9 = melt(c("BarrenLand"="peachpuff2","CropLand"="gold1","Developed"="gray","Forest"="darkgreen","Herbaceous"="limegreen",
                   "Mexico"="red1","OpenWater"="mediumblue","Others"="black","SaltonSea"="darkturquoise","Shrub"="saddlebrown"))
colorset9$SurfaceType<-rownames(colorset9)
colnames(colorset9)[which(names(colorset9) == "value")] <- "ST_Color"
colorset9

# merge color palette and STF data together
STF.melt<-merge(STF.melt, colorset9, by="SurfaceType")
head(STF.melt)
STF.melt$ST_Color <- as.character(STF.melt$ST_Color)
head(SurfTypFreq)

# create factors for organizing data in plot
STF.melt$Site<-factor(STF.melt$Site,levels=c("PD","BDC","DP","WI"))
unique(STF.melt$STF_Date)
STF.melt$SampleID = factor(STF.melt$SampleID, levels=unique(STF.melt$SampleID[order(STF.melt$Site,STF.melt$STF_Date)]), ordered=TRUE)

# plot time
stfplt1<-ggplot(STF.melt, aes(x=SampleID, y=Frequency, fill=SurfaceType))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Surface Type Frequencies in Salton Sea Dust", x="SampleID", y="Frequency", subtitle="",fill="Surface Type")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(name ="Surface Type",values=unique(STF.melt$ST_Color[order(STF.melt$SurfaceType)]),labels=c("Barren Land","Crop Land","Developed","Forest","Herbaceous","Mexico","Open Water","Others","Salton Sea","Shrub"))

ggsave(stfplt1,filename = "figures/SurfaceTypeFrequencies/SSD_STFs_Barplot.png", width=15, height=10, dpi=600)

stfplt2<-ggplot(STF.melt, aes(x=SampleID, y=Frequency, fill=SurfaceType))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Surface Type Frequencies in Salton Sea Dust", x="SampleID", y="Frequency", subtitle="",fill="Surface Type")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(name ="Surface Type",values=unique(STF.melt$ST_Color[order(STF.melt$SurfaceType)]),labels=c("Barren Land","Crop Land","Developed","Forest","Herbaceous","Mexico","Open Water","Others","Salton Sea","Shrub")) +
  facet_wrap(vars(Site), scales = "free")

ggsave(stfplt2,filename = "figures/SurfaceTypeFrequencies/SSD_STFs_Barplot2.png", width=15, height=15, dpi=600)

# #### Using Shapiro-Wilk test for Normality ####
#
# shapiro.test(SurfTypFreq$BarrenLand) # what is the p-value?
# # W = 0.94805, p-value = 0.1769
# # p > 0.05 states distribution of data are not significantly different from normal distribution
# # p < 0.05 means that data is significantly different from a normal distribution
# hist(SurfTypFreq$BarrenLand, col="blue")
#
# # visualize Q-Q plot for species richness
# qqnorm(meta.all.scaled$DO_Percent_Local, pch = 1, frame = FALSE)
# qqline(meta.all.scaled$DO_Percent_Local, col = "steelblue", lwd = 2)
#
# shapiro.test(meta.all.scaled$ORP_mV) # what is the p-value? p-value = 3.323e-12
# hist(meta.all.scaled$ORP_mV, col="blue")
# # visualize Q-Q plot for species richness
# qqnorm(meta.all.scaled$ORP_mV, pch = 1, frame = FALSE)
# qqline(meta.all.scaled$ORP_mV, col = "steelblue", lwd = 2)
#
# shapiro.test(meta.all.scaled$Temp_DegC) # what is the p-value? p-value = 3.562e-06
# hist(meta.all.scaled$Temp_DegC, col="blue")
# # visualize Q-Q plot for species richness
# qqnorm(meta.all.scaled$Temp_DegC, pch = 1, frame = FALSE)
# qqline(meta.all.scaled$Temp_DegC, col = "steelblue", lwd = 2)
#
# shapiro.test(meta.all.scaled$Dissolved_OrganicMatter_RFU) # what is the p-value? p-value = 1.997e-07
# hist(meta.all.scaled$Dissolved_OrganicMatter_RFU, col="blue")
# # visualize Q-Q plot for species richness
# qqnorm(meta.all.scaled$Dissolved_OrganicMatter_RFU, pch = 1, frame = FALSE)
# qqline(meta.all.scaled$Dissolved_OrganicMatter_RFU, col = "steelblue", lwd = 2)
#
# shapiro.test(meta.all.scaled$Sulfate_milliM) # what is the p-value?
# # my p-value was p-value =  0.006965
# # p > 0.05 states distribution of data are not significantly different from normal distribution
# # p < 0.05 means that data is significantly different from a normal distribution
# hist(meta.all.scaled$Sulfate_milliM, col="blue")
#
# # visualize Q-Q plot for species richness
# qqnorm(meta.all.scaled$Sulfate_milliM, pch = 1, frame = FALSE)
# qqline(meta.all.scaled$Sulfate_milliM, col = "steelblue", lwd = 2)
#
# shapiro.test(meta.all.scaled$Sulfide_microM) # what is the p-value?
# # my p-value was p-value =  5.934e-12
# # p > 0.05 states distribution of data are not significantly different from normal distribution
# # p < 0.05 means that data is significantly different from a normal distribution
# hist(meta.all.scaled$Sulfide_microM, col="blue")
#
# # visualize Q-Q plot for species richness
# qqnorm(meta.all.scaled$Sulfide_microM, pch = 1, frame = FALSE)
# qqline(meta.all.scaled$Sulfide_microM, col = "steelblue", lwd = 2)
#
# shapiro.test(meta.all.scaled$Turbidity_FNU) # what is the p-value?  p-value = 0.0005629
# hist(meta.all.scaled$Turbidity_FNU, col="blue")
# # visualize Q-Q plot for species richness
# qqnorm(meta.all.scaled$Turbidity_FNU, pch = 1, frame = FALSE)
# qqline(meta.all.scaled$Turbidity_FNU, col = "steelblue", lwd = 2)
#
# shapiro.test(meta.all.scaled$Chlorophyll_RFU) # what is the p-value? p-value = 1.044e-11
# hist(meta.all.scaled$Chlorophyll_RFU, col="blue")
# # visualize Q-Q plot for species richness
# qqnorm(meta.all.scaled$Chlorophyll_RFU, pch = 1, frame = FALSE)
# qqline(meta.all.scaled$Chlorophyll_RFU, col = "steelblue", lwd = 2)
#
#### PCA w/ STF Variables ####
head(meta.all.scaled[,c(4,7:8,10:12,38:47)])

STFs<-meta.all.scaled[,c(38:47)]
head(STFs)

# Note: log transforming STFs because this is what we do with them statistically later
STF.log<-decostand(STFs,method = "log") # log transformation of STF data before scaling
STF.log[1:4,1:4]

# check rownames of log transformed STF data & metadata
rownames(STF.log) %in% rownames(dust.meta.surf)
dust.meta.surf=dust.meta.surf[rownames(STF.log),] ## reorder metadata to match order of log data

# calculate our Euclidean distance matrix using log data
STF.euc_dist <- dist(STF.log, method = "euclidean")

# creating our hierarcical clustering dendrogram
STF.euc_clust <- hclust(STF.euc_dist, method="ward.D2")

# let's make it a little nicer...
STF.euc_dend <- as.dendrogram(STF.euc_clust, hang=0.2)
STF.dend_cols <- as.character(meta.all.scaled$SampDate_Color[order.dendrogram(STF.euc_dend)])
labels_colors(STF.euc_dend) <- STF.dend_cols

plot(STF.euc_dend, ylab="log Euclidean Distance", horiz=TRUE,cex=03) + title(main = "Surface Type Frequencies Clustering Dendrogram", font.main= 1, cex.sub = 0.8, font.sub = 3)
#legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c( "#ef781c","#03045e","#32cbff","#059c3f"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# let's use our Euclidean distance matrix from before
STF.pca <- prcomp(STF.log, center=TRUE, scale=TRUE) # pca of euclidean distance matrix = PCA of euclidean distance matrix
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
  labs(title="PCA: Surface Type Frequencies in Salton Sea Dust",subtitle="Using Log Transformed Surface Type Frequencies",color="Collection Date",shape="Site")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Date",values=unique(STF.pca.meta$SampDate_Color[order(STF.pca.meta$SampDate)]),
                     labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  scale_shape_manual(values = c(7,10, 15,16)) +
  xlab("PC1 [63.58%]") + ylab("PC2 [18.92]")

ggsave(pca1a,filename = "figures/SurfaceTypeFrequencies/SSD_log_STFs_PCA.png", width=15, height=10, dpi=600)

pca1b<-ggplot(STF.pca.meta, aes(x=PC1, y=PC2)) +geom_point(aes(color=factor(SampDate),shape=Site,size=CollectionYear))+theme_bw()+
  labs(title="PCA: Surface Type Frequencies in Salton Sea Dust",subtitle="Using log Transformed Surface Type Frequencies",color="Collection Date",shape="Site",size="Collection Year")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  scale_color_manual(name ="Collection Date",values=unique(STF.pca.meta$SampDate_Color[order(STF.pca.meta$SampDate)]),
                     labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  scale_shape_manual(values = c(7,10, 15,16,17,18)) +
  xlab("PC1 [63.58%]") + ylab("PC2 [18.92]") + scale_size_manual(values = c("2020" = 7, "2021"=4)) +   guides(shape = guide_legend(override.aes = list(size = 5)))

ggsave(pca1b,filename = "figures/SurfaceTypeFrequencies/SSD_log_STFs_PCA_v2.png", width=15, height=10, dpi=600)

pca1c<-ggplot(STF.pca.meta, aes(x=PC1, y=PC2)) +geom_point(aes(color=factor(STF_Date),shape=Site,size=CollectionYear))+theme_bw()+
  labs(title="PCA: Surface Type Frequencies in Salton Sea Dust",subtitle="Using log Transformed Surface Type Frequencies",color="Surface Date",shape="Site",size="Collection Year")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  scale_shape_manual(values = c(7,10, 15,16,17,18)) + scale_color_manual(name ="Date Ranges",values=brewer.pal(n = 7, name = "Dark2"), labels=c("5/13/20 - 7/10/20", "7/10/20 - 8/30/20", "8/30/20 - 10/10/20","10/10/20 - 11/6/20",
                                                                                                     "6/5/21 - 8/19/21", "8/19/21 - 10/1/21", "10/1/21 - 12/8/21")) +
  xlab("PC1 [63.58%]") + ylab("PC2 [18.92]") + scale_size_manual(values = c("2020" = 7, "2021"=4)) +   guides(shape = guide_legend(override.aes = list(size = 5)))

ggsave(pca1c,filename = "figures/SurfaceTypeFrequencies/SSD_log_STFs_PCA_v3.png", width=15, height=10, dpi=600)

# sample month shape, depth color
# pca2<-ggplot(STF.pca.meta, aes(x=PC1, y=PC2)) +
#   geom_point(aes(color=as.numeric(Depth_m),shape=SampleMonth), size=5)+theme_bw()+
#   labs(title="PCA: Surface Type Frequencies in Salton Seawater",subtitle="Using log Transformed Data",xlab="Axis 1", ylab="Axis 2",color="Depth (m)")+
#   theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
#   scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="sample site") +
#   xlab("Axis 1 [63.58%]") + ylab("Axis 2 [18.92]")
#
# ggsave(pca2,filename = "figures/STFVariablesOnly/SSD_logSTFOnly_PCA_Depth_SampDate.png", width=12, height=10, dpi=600)

#### PCA w/ Wind Condition Variables ####
head(meta.all.scaled[,c(4,7:8,10:12,38:47)])

clim.scaled<-meta.all.scaled[,c(4,7:8,10:12)]
head(clim.scaled)

# Note: log transforming clim.scaled because this is what we do with them statistically later
clim.log<-decostand(clim.scaled,method = "log") # log transformation of clim data before scaling
clim.log[1:4,1:4]

# check rownames of log transformed STF data & metadata
rownames(clim.log) %in% rownames(dust.meta.surf)
dust.meta.surf=dust.meta.surf[rownames(clim.log),] ## reorder metadata to match order of log data

# calculate our Euclidean distance matrix using log data
clim.euc_dist <- dist(clim.log, method = "euclidean")

# creating our hierarcical clustering dendrogram
clim.euc_clust <- hclust(clim.euc_dist, method="ward.D2")

# let's make it a little nicer...
clim.euc_dend <- as.dendrogram(clim.euc_clust, hang=0.2)
clim.dend_cols <- as.character(meta.all.scaled$SampDate_Color[order.dendrogram(clim.euc_dend)])
labels_colors(clim.euc_dend) <- clim.dend_cols

plot(clim.euc_dend, ylab="log Euclidean Distance", horiz=TRUE,cex=03) + title(main = "Surface Type Frequencies Clustering Dendrogram", font.main= 1, cex.sub = 0.8, font.sub = 3)
#legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c( "#ef781c","#03045e","#32cbff","#059c3f"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# let's use our Euclidean distance matrix from before
clim.pca <- prcomp(clim.log, center=TRUE, scale=TRUE) # pca of euclidean distance matrix = PCA of euclidean distance matrix
clim.pca$x # where sites fall on PC axes
clim.pca$rotation # variables on PC axes
summary(clim.pca)$importance
# The proportion of variances
# PC1 = 33.49%, PC2 = 19.84%

# extract principal coordinates
clim.pca.vectors<-data.frame(clim.pca$x)
clim.pca.vectors$SampleID<-rownames(clim.pca$x)

# merge pca coordinates w/ metadata
clim.pca.meta<-merge(clim.pca.vectors, dust.meta.surf, by.x="SampleID", by.y="SampleID")
clim.pca.meta$STF_Date

head(clim.pca.meta)
summary(clim.pca)$importance # percentage of variation explained for pca below

# create pca ggplot fig
pca1a<-ggplot(clim.pca.meta, aes(x=PC1, y=PC2)) +geom_point(aes(color=factor(SampDate),shape=Site), size=5)+theme_bw()+
  labs(title="PCA: Surface Type Frequencies in Salton Sea Dust",subtitle="Using Log Transformed Climate Data",color="Collection Date",shape="Site")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Date",values=unique(clim.pca.meta$SampDate_Color[order(clim.pca.meta$SampDate)]),
                     labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  scale_shape_manual(values = c(7,10, 15,16)) +
  xlab("PC1 [33.49%]") + ylab("PC2 [19.84]")

ggsave(pca1a,filename = "figures/SurfaceTypeFrequencies/SSD_log_clim.scaled_PCA.png", width=15, height=10, dpi=600)

pca1b<-ggplot(clim.pca.meta, aes(x=PC1, y=PC2)) +geom_point(aes(color=factor(SampDate),shape=Site,size=CollectionYear))+theme_bw()+
  labs(title="PCA: Surface Type Frequencies in Salton Sea Dust",subtitle="Using log Transformed Climate Data",color="Collection Date",shape="Site",size="Collection Year")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  scale_color_manual(name ="Collection Date",values=unique(clim.pca.meta$SampDate_Color[order(clim.pca.meta$SampDate)]),
                     labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  scale_shape_manual(values = c(7,10, 15,16,17,18)) +
  xlab("PC1 [33.49%]") + ylab("PC2 [19.84]") + scale_size_manual(values = c("2020" = 7, "2021"=4)) +   guides(shape = guide_legend(override.aes = list(size = 5)))

ggsave(pca1b,filename = "figures/SurfaceTypeFrequencies/SSD_log_clim.scaled_PCA_v2.png", width=15, height=10, dpi=600)

pca1c<-ggplot(clim.pca.meta, aes(x=PC1, y=PC2)) +geom_point(aes(color=factor(STF_Date),shape=Site,size=CollectionYear))+theme_bw()+
  labs(title="PCA: Surface Type Frequencies in Salton Sea Dust",subtitle="Using log Transformed Climate Data",color="Surface Date",shape="Site",size="Collection Year")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  scale_shape_manual(values = c(7,10, 15,16,17,18)) + scale_color_manual(name ="Date Ranges",values=brewer.pal(n = 7, name = "Dark2"), labels=c("5/13/20 - 7/10/20", "7/10/20 - 8/30/20", "8/30/20 - 10/10/20","10/10/20 - 11/6/20",
                                                                                                                                                "6/5/21 - 8/19/21", "8/19/21 - 10/1/21", "10/1/21 - 12/8/21")) +
  xlab("PC1 [33.49%]") + ylab("PC2 [19.84]") + scale_size_manual(values = c("2020" = 7, "2021"=4)) +   guides(shape = guide_legend(override.aes = list(size = 5)))

ggsave(pca1c,filename = "figures/SurfaceTypeFrequencies/SSD_log_clim.scaled_PCA_v3.png", width=15, height=10, dpi=600)

# sample month shape, depth color
# pca2<-ggplot(clim.pca.meta, aes(x=PC1, y=PC2)) +
#   geom_point(aes(color=as.numeric(Depth_m),shape=SampleMonth), size=5)+theme_bw()+
#   labs(title="PCA: Surface Type Frequencies in Salton Seawater",subtitle="Using log Transformed Data",xlab="Axis 1", ylab="Axis 2",color="Depth (m)")+
#   theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
#   scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="sample site") +
#   xlab("PC1 [33.49%]") + ylab("PC2 [19.84]")
#
# ggsave(pca2,filename = "figures/STFVariablesOnly/SSD_logSTFOnly_PCA_Depth_SampDate.png", width=12, height=10, dpi=600)

#### Do Env Variables Correlate? ####
head(meta.all.scaled)
# check for colinearity among env variables themselves
heatmap(abs(cor(meta.all.scaled[,c(4,7:8,10:12,38:47)])),
        # Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(6)),
        Colv = NA, Rowv = NA)
legend("topleft",
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       bty = "n",
       fill = rev(heat.colors(6)))
dev.off()

# Visualize with a corrplot - ALL climate variables together
cor_mat.env1 <- cor(meta.all.scaled[,c(4,7:8,10:12,38:47)], method='pearson')
cor_mat.env1

symnum(cor_mat.env1)

p.env <- cor_pmat(meta.all.scaled[,c(4,7:8,10:12,38:47)])

symnum(cor_mat.env1)

ggcorrplot(cor_mat.env1,method="square",lab = T,p.mat=p.env,
           hc.order=TRUE,outline.color="white",type = "lower")

png('figures/EnvVariablesOnly/SSD_ScaledCentered_EnvVarOnly_CorrPlot.png',width = 1100, height = 1100)
ggcorrplot(cor_mat.env1,p.mat=p.env,hc.order=TRUE,
           method="square",outline.color="white",type = "lower")

# env variables with a correlation of <|0.7| is a good threshold for determining if predictors correlate
dev.off()

png('figures/EnvVariablesOnly/SSD_ScaledCentered_EnvVarOnly_CorrPlot_labeled.png',width = 1100, height = 1100)
ggcorrplot(cor_mat.env1,method="square",lab = T,p.mat=p.env,sig.level = 0.05,insig="blank",
           hc.order=TRUE,outline.color="white",type = "lower")

# env variables with a correlation of <|0.7| is a good threshold for determining if predictors correlate
dev.off()

# save corrplots into 1 plot
env.corrplot1<-ggcorrplot(cor_mat.env1,p.mat=p.env,hc.order=TRUE,
                          method="square",outline.color="white",type = "lower")

env.corrplot2<-ggcorrplot(cor_mat.env1,method="square",lab = T,p.mat=p.env,sig.level = 0.05,insig="blank",
                          hc.order=TRUE,outline.color="white",type = "lower")

env.corrplot.together<-ggarrange(env.corrplot1,env.corrplot2, ncol = 1, nrow = 2,common.legend=TRUE)

ggsave(env.corrplot.together,filename = "figures/EnvVariablesOnly/SSD_ScaledCentered_EnvVarOnly_CorrPlot_Combined.png", width=25, height=35, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).


# subset corr and p value results
cor.envs.df<-as.data.frame(cor_mat.env1)
p.env.df<-as.data.frame(p.env)

# recA.corr<-data.frame(SampleID=rownames(cor.envs.df),envCorr.recA=cor.envs.df$`recA; recombination protein RecA`)
# recA.pval<-data.frame(SampleID=rownames(p.env.df),envCorrPval.recA=p.env.df$`recA; recombination protein RecA`)
#
# recA.corr.df<-merge(recA.corr,recA.pval,by="SampleID")
# recA.corr.df<-recA.corr.df[order(recA.corr.df$envCorrPval.recA),]


#### Homogeneity of Variance & PERMANOVA tests - STFs Only by Groups ####
## betadisper to look at homogeneity of group dispersions (aka variance) when considering multiple variables
# multivariate analogue to Levene's test of homogeneity of variances
# program finds spatial median or centroid of the group, & compare distances of group to centroid/spatial median via ANOVA

#While PERMANOVA tests differences in group means (analogous to MANOVA),
## a related test called PERMDISP can be used to evaluate homogeneity of group dispersion
#(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
## * need a distance matrix!

head(meta.all.scaled[,c(4,7:8,10:12,38:47)])

# log transform the scaled STFs
STF.log<-decostand(meta.all.scaled[,c(38:47)],method = "log")

# create Euclidean distance matrix from log transformed, scaled STF data (STF data, all scaled)
STF.dist<-dist(STF.log,method="euclidean")
# first by compare dispersions by sampling date
STF.disper1<-betadisper(STF.dist, meta.all.scaled$Site)
STF.disper1

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:

permutest(STF.disper1, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:

anova(STF.disper1) # p = 0.7037 --> accept the Null H, spatial medians (a measure of dispersion) are NOT significantly difference across sample sites
# ANOVA adjusted p-value
aov.beta.p1<-anova(STF.disper1)[["Pr(>F)"]] # get p values from ANOVA
p.adjust(aov.beta.p1,method="bonferroni",n=length(aov.beta.p1))

TukeyHSD(STF.disper1) # tells us which sample sites/category's dispersion MEANS are significantly different than each other

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

pnova1<-adonis2(STF.dist ~ Site,data=meta.all.scaled,by="terms",permutations=10000)
pnova1 # p-value = 0.000999
p.adjust(pnova1$`Pr(>F)`,method="bonferroni",n=length(pnova1$`Pr(>F)`)) # adjusted pval
# 0.002997003

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

pair.mod1<-pairwise.adonis(STF.dist,meta.all.scaled$Site, p.adjust.m='bonferroni',perm=9999) # shows us variation for each sample to see which ones are different
pair.mod1
#       pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
# 1 BDC vs DP  1 138017552 10.5827274 0.46862043   0.001      0.006   *
# 2 BDC vs PD  1   7719598  0.7262013 0.05706348   0.438      1.000
# 3 BDC vs WI  1 253642722 14.7717374 0.55176611   0.001      0.006   *
# 4  DP vs PD  1 117686302  8.6941739 0.42012665   0.001      0.006   *
# 5  DP vs WI  1  39521423  1.9684993 0.14092418   0.166      0.996
# 6  PD vs WI  1 223488508 12.6512985 0.51321023   0.001      0.006   *

# Visualize dispersions
png('figures/SurfaceTypeFrequencies/SSD_STFsOnly_pcoa_betadispersion_site.png',width = 700, height = 600, res=100)
plot(STF.disper1,main = "Centroids and Dispersion based on Eucldiean Distance of Scaled Surface Type Frequencies", col=colorset6$Site_Color)
dev.off()

png('figures/SurfaceTypeFrequencies/SSD_STFsOnly__boxplot_centroid_distance_site.png',width = 700, height = 600, res=100)
boxplot(STF.disper1,xlab="By Site", main = "Distance to Centroid by Category", sub="Based on Weighted Unifrac Distance", col=colorset6$Site_Color)
dev.off()

### now compare dispersions by site + year
STF.disper2<-betadisper(STF.dist, interaction(meta.all.scaled$Site,meta.all.scaled$CollectionYear,sep="."))
STF.disper2
# NOTE: SB and RHB have less samples than other sites since they are supposed to represent similar locations
# maybe we should remove these sites and then rerun the PERMANOVA?

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:

permutest(STF.disper2, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:

anova(STF.disper2) # p = 0.8886 --> accept the Null H, spatial medians (a measure of dispersion) are NOT significantly difference across sample sites
# ANOVA adjusted p-value
aov.beta.p3<-anova(STF.disper2)[["Pr(>F)"]] # get p values from ANOVA
p.adjust(aov.beta.p3,method="bonferroni",n=length(aov.beta.p3))
# 1

TukeyHSD(STF.disper2) # tells us which sample sites/category's dispersion MEANS are significantly different than each other
#                     diff       lwr       upr     p adj

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

pnova3<-adonis2(STF.dist ~ Site*CollectionYear,data=meta.all.scaled,by="terms",permutations=10000)
pnova3
#                       Df SumOfSqs      R2      F  Pr(>F)
# Site                 3  138.244 0.51201 7.8074 0.000999 ***
# CollectionYear       1    8.186 0.03032 1.3869 0.240759
# Site:CollectionYear  3    5.525 0.02046 0.3120 0.973027
# Residual            20  118.045 0.43720
# Total               27  270.000 1.00000

p.adjust(pnova3$`Pr(>F)`,method="bonferroni",n=length(pnova3$`Pr(>F)`)) # adjusted pval
#0.004995005 1.000000000 1.000000000  NA

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

pair.mod3<-pairwise.adonis(STF.dist,interaction(meta.all.scaled$Site,meta.all.scaled$CollectionYear,sep="."), p.adjust.m='bonferroni',perm=9999) # shows us variation for each sample to see which ones are different
pair.mod3

# Visualize dispersions
png('figures/SurfaceTypeFrequencies/SSD_STFsOnly__pcoa_betadispersion_site_by_year.png',width = 700, height = 600, res=100)
plot(STF.disper2,main = "Centroids and Dispersion based on Weighted Unifrac Distance", col=colorset6$Site_Color)
dev.off()

png('figures/SurfaceTypeFrequencies/SSD_STFsOnly__boxplot_centroid_distance_site_by_year.png',width = 900, height = 600, res=100)
boxplot(STF.disper2,xlab="By Site x Collection Year", main = "Distance to Centroid by Category", sub="Based on Weighted Unifrac Distance", col=colorset6$Site_Color)
dev.off()

#### Homogeneity of Variance & PERMANOVA tests - Climate Data Only by Groups ####
## betadisper to look at homogeneity of group dispersions (aka variance) when considering multiple variables
# multivariate analogue to Levene's test of homogeneity of variances
# program finds spatial median or centroid of the group, & compare distances of group to centroid/spatial median via ANOVA

#While PERMANOVA tests differences in group means (analogous to MANOVA),
## a related test called PERMDISP can be used to evaluate homogeneity of group dispersion
#(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
## * need a distance matrix!

head(meta.all.scaled[,c(4,7:8,10:12,38:47)])

# log transform the scaled STFs
clim.log<-decostand(meta.all.scaled[,c(4,7:8,10:12)],method = "log")

# create Euclidean distance matrix from log transformed, scaled climate data (climate data, all scaled)
clim.dist<-dist(clim.log,method="euclidean")
# first by compare dispersions by sampling date
clim.disper1<-betadisper(clim.dist, meta.all.scaled$Site)
clim.disper1

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:

permutest(clim.disper1, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#       PD     BDC      DP    WI

anova(clim.disper1) # p = 0.8567 --> accept the Null H, spatial medians (a measure of dispersion) are NOT significantly difference across sample sites
# ANOVA adjusted p-value
aov.beta.p1<-anova(clim.disper1)[["Pr(>F)"]] # get p values from ANOVA
p.adjust(aov.beta.p1,method="bonferroni",n=length(aov.beta.p1))

TukeyHSD(clim.disper1) # tells us which sample sites/category's dispersion MEANS are significantly different than each other
#         diff       lwr       upr     p adj

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

pnova2<-adonis2(clim.dist ~ Site,data=meta.all.scaled,by="terms",permutations=10000)
pnova2 # p-value = 0.0000999
p.adjust(pnova2$`Pr(>F)`,method="bonferroni",n=length(pnova2$`Pr(>F)`)) # adjusted pval
# 0.0002997003

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

pair.mod2<-pairwise.adonis(clim.dist,meta.all.scaled$Site, p.adjust.m='bonferroni',perm=9999) # shows us variation for each sample to see which ones are different
pair.mod2
#       pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
# 1 BDC vs DP  1  6293.282  8.503451 0.4147327  0.0009     0.0054   *
# 2 BDC vs PD  1  8699.037 10.694396 0.4712351  0.0010     0.0060   *
# 3 BDC vs WI  1  4955.641  5.799380 0.3258192  0.0036     0.0216   .
# 4  DP vs PD  1  6929.736  8.067000 0.4020033  0.0004     0.0024   *
# 5  DP vs WI  1  6049.102  6.720368 0.3589870  0.0006     0.0036   *
# 6  PD vs WI  1  6656.352  6.837904 0.3629864  0.0030     0.0180   .

# Visualize dispersions
png('figures/SurfaceTypeFrequencies/SSD_ClimDataOnly_pcoa_betadispersion_site.png',width = 700, height = 600, res=100)
plot(clim.disper1,main = "Centroids and Dispersion based on Eucldiean Distance of Scaled climironmental Variables", col=colorset6$Site_Color)
dev.off()

png('figures/SurfaceTypeFrequencies/SSD_ClimDataOnly_boxplot_centroid_distance_site.png',width = 700, height = 600, res=100)
boxplot(clim.disper1,xlab="By Site", main = "Distance to Centroid by Category", sub="Based on Weighted Unifrac Distance", col=colorset6$Site_Color)
dev.off()

### now compare dispersions by site + year
clim.disper2<-betadisper(clim.dist, interaction(meta.all.scaled$Site,meta.all.scaled$CollectionYear,sep="."))
clim.disper2
# NOTE: SB and RHB have less samples than other sites since they are supposed to represent similar locations
# maybe we should remove these sites and then rerun the PERMANOVA?

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:

permutest(clim.disper2, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:

anova(clim.disper2) # p = 0.9512 --> accept the Null H, spatial medians (a measure of dispersion) are NOT significantly difference across sample sites
# ANOVA adjusted p-value
aov.beta.p3<-anova(clim.disper2)[["Pr(>F)"]] # get p values from ANOVA
p.adjust(aov.beta.p3,method="bonferroni",n=length(aov.beta.p3))
# 1

TukeyHSD(clim.disper2) # tells us which sample sites/category's dispersion MEANS are significantly different than each other
#                     diff       lwr       upr     p adj

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

pnova3<-adonis2(clim.dist ~ Site*CollectionYear,data=meta.all.scaled,by="terms",permutations=10000)
pnova3
#                       Df SumOfSqs      R2      F  Pr(>F)
# Site                 3   197.19 0.52168 8.0956 0.000999 ***
# CollectionYear       1    11.63 0.03076 1.4319 0.223776
# Site:CollectionYear  3     6.79 0.01797 0.2788 1.000000
# Residual            20   162.39 0.42960
# Total               27   378.00 1.00000

p.adjust(pnova3$`Pr(>F)`,method="bonferroni",n=length(pnova3$`Pr(>F)`)) # adjusted pval
# 0.004995005 0.954045954 1.000000000  NA

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

pair.mod3<-pairwise.adonis(clim.dist,interaction(meta.all.scaled$Site,meta.all.scaled$CollectionYear,sep="."), p.adjust.m='bonferroni',perm=9999) # shows us variation for each sample to see which ones are different
pair.mod3

# Visualize dispersions
png('figures/SurfaceTypeFrequencies/SSD_ClimDataOnly_pcoa_betadispersion_site_by_year.png',width = 700, height = 600, res=100)
plot(clim.disper2,main = "Centroids and Dispersion based on Weighted Unifrac Distance", col=colorset6$Site_Color)
dev.off()

png('figures/SurfaceTypeFrequencies/SSD_ClimDataOnly_boxplot_centroid_distance_site_by_year.png',width = 900, height = 600, res=100)
boxplot(clim.disper2,xlab="By Site x Collection Year", main = "Distance to Centroid by Category", sub="Based on Weighted Unifrac Distance", col=colorset6$Site_Color)
dev.off()

#### Homogeneity of Variance & PERMANOVA tests - ALL Env Vars by Groups ####
## betadisper to look at homogeneity of group dispersions (aka variance) when considering multiple variables
# multivariate analogue to Levene's test of homogeneity of variances
# program finds spatial median or centroid of the group, & compare distances of group to centroid/spatial median via ANOVA

#While PERMANOVA tests differences in group means (analogous to MANOVA),
## a related test called PERMDISP can be used to evaluate homogeneity of group dispersion
#(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
## * need a distance matrix!

head(meta.all.scaled[,c(4,7:8,10:12,38:47)])

# create log distance matrix of scaeld env data (for comparisons)
env.log<-decostand(meta.all.scaled[,c(4,7:8,10:12,38:47)],method = "log", pseudocount = 2) # log transformation of STF data before scaling
env.log[1:4,1:4]
env.log.dist<-dist(env.log,method="euclidean")

# create Euclidean distance matrix from scaled env data (incldues climate + STF data, all scaled)
env.dist<-dist(meta.all.scaled[,c(4,7:8,10:12,38:47)],method="euclidean")
# first by compare dispersions by sampling date
env.disper1<-betadisper(env.dist, meta.all.scaled$Site)
env.disper1

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:

permutest(env.disper1, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#       PD     BDC      DP    WI
# PD          0.86200 0.99600 0.121
# BDC 0.86457         0.88400 0.106
# DP  0.99446 0.87099         0.140
# WI  0.12595 0.11096 0.14121

anova(env.disper1) # p = 0.1533 --> accept the Null H, spatial medians (a measure of dispersion) are NOT significantly difference across sample sites
# ANOVA adjusted p-value
aov.beta.p1<-anova(env.disper1)[["Pr(>F)"]] # get p values from ANOVA
p.adjust(aov.beta.p1,method="bonferroni",n=length(aov.beta.p1))

TukeyHSD(env.disper1) # tells us which sample sites/category's dispersion MEANS are significantly different than each other
#         diff       lwr       upr     p adj

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

pnova1<-adonis2(env.dist ~ Site,data=meta.all.scaled,by="terms",permutations=10000)
pnova1 # p-value = 0.000999
p.adjust(pnova1$`Pr(>F)`,method="bonferroni",n=length(pnova1$`Pr(>F)`)) # adjusted pval
# 0.002997003

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

pair.mod1<-pairwise.adonis(env.dist,meta.all.scaled$Site, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
pair.mod1
#       pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
# 1 BDC vs DP  1  59.30312 10.200572 0.4594734   0.001      0.006   *
# 2 BDC vs PD  1  48.28358  8.444903 0.4130566   0.001      0.006   *
# 3 BDC vs WI  1 104.15556 11.202059 0.4828045   0.001      0.006   *
# 4  DP vs PD  1  52.47065  9.094897 0.4311420   0.002      0.012   .
# 5  DP vs WI  1  36.81390  3.937462 0.2470570   0.034      0.204
# 6  PD vs WI  1  93.36193 10.089440 0.4567540   0.002      0.012   .

# Visualize dispersions
png('figures/SurfaceTypeFrequencies/SSD_pcoa_betadispersion_site.png',width = 700, height = 600, res=100)
plot(env.disper1,main = "Centroids and Dispersion based on Eucldiean Distance of Scaled Environmental Variables", col=colorset6$Site_Color)
dev.off()

png('figures/SurfaceTypeFrequencies/SSD_boxplot_centroid_distance_site.png',width = 700, height = 600, res=100)
boxplot(env.disper1,xlab="By Site", main = "Distance to Centroid by Category", sub="Based on Weighted Unifrac Distance", col=colorset6$Site_Color)
dev.off()

### now compare dispersions by site + year
env.disper2<-betadisper(env.dist, interaction(meta.all.scaled$Site,meta.all.scaled$CollectionYear,sep="."))
env.disper2
# NOTE: SB and RHB have less samples than other sites since they are supposed to represent similar locations
# maybe we should remove these sites and then rerun the PERMANOVA?

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:

permutest(env.disper2, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:

anova(env.disper2) # p = 0.9512 --> accept the Null H, spatial medians (a measure of dispersion) are NOT significantly difference across sample sites
# ANOVA adjusted p-value
aov.beta.p3<-anova(env.disper2)[["Pr(>F)"]] # get p values from ANOVA
p.adjust(aov.beta.p3,method="bonferroni",n=length(aov.beta.p3))
# 1

TukeyHSD(env.disper2) # tells us which sample sites/category's dispersion MEANS are significantly different than each other
#                     diff       lwr       upr     p adj

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

pnova3<-adonis2(env.dist ~ Site*CollectionYear,data=meta.all.scaled,by="terms",permutations=10000)
pnova3
#                       Df SumOfSqs      R2      F  Pr(>F)
# Site                 3   197.19 0.52168 8.0956 0.000999 ***
# CollectionYear       1    11.63 0.03076 1.4319 0.223776
# Site:CollectionYear  3     6.79 0.01797 0.2788 1.000000
# Residual            20   162.39 0.42960
# Total               27   378.00 1.00000

p.adjust(pnova3$`Pr(>F)`,method="bonferroni",n=length(pnova3$`Pr(>F)`)) # adjusted pval
# 0.004995005 0.954045954 1.000000000  NA

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

pair.mod3<-pairwise.adonis(env.dist,interaction(meta.all.scaled$Site,meta.all.scaled$CollectionYear,sep="."), p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
pair.mod3

# Visualize dispersions
png('figures/SurfaceTypeFrequencies/SSD_pcoa_betadispersion_site_by_year.png',width = 700, height = 600, res=100)
plot(env.disper2,main = "Centroids and Dispersion based on Weighted Unifrac Distance", col=colorset6$Site_Color)
dev.off()

png('figures/SurfaceTypeFrequencies/SSD_boxplot_centroid_distance_site_by_year.png',width = 900, height = 600, res=100)
boxplot(env.disper2,xlab="By Site x Collection Year", main = "Distance to Centroid by Category", sub="Based on Weighted Unifrac Distance", col=colorset6$Site_Color)
dev.off()


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
bac.div.metadat <- merge(d.r_16s,meta.all.scaled, by.x="SampleID", by.y="SampleID")
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
dep.dom<-ggplot(meta.all.scaled, aes(x=Depth_m, y=Dissolved_OrganicMatter_RFU,color=SampDate,group=SampDate)) + geom_point(size=3) + geom_line() + theme_bw()+
  labs(title="Dissolved Organic Matter (DOM) by Depth & sample site",subtitle="Using Scaled DOM RFU Data",color="sample site")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="sample site",values=unique(meta.all.scaled$SampDate_Color[order(meta.all.scaled$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Depth (m)") + ylab("DOM (RFU)")+coord_flip()+ scale_x_discrete(limits=rev)

ggsave(dep.dom,filename = "figures/EnvVariablesOnly/SSD_DOM_Depth_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

dep.orp<-ggplot(meta.all.scaled, aes(x=Depth_m, y=ORP_mV,color=SampDate,group=SampDate)) + geom_point(size=3) + geom_line() + theme_bw()+
  labs(title="Oxidative-Reduction Potential by Depth & sample site",subtitle="Using Scaled ORP (mV) Data",color="sample site")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="sample site",values=unique(meta.all.scaled$SampDate_Color[order(meta.all.scaled$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Depth (m)") + ylab("ORP (mV)")+coord_flip()+ scale_x_discrete(limits=rev)

ggsave(dep.orp,filename = "figures/EnvVariablesOnly/SSD_ORP_Depth_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

dep.sal<-ggplot(meta.all.scaled, aes(x=Depth_m, y=Salinity_ppt,color=SampDate,group=SampDate)) + geom_point(size=3) + geom_line() + theme_bw()+
  labs(title="Salinity by Depth & sample site",subtitle="Using Scaled Salinity (PPT) Data",color="sample site")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="sample site",values=unique(meta.all.scaled$SampDate_Color[order(meta.all.scaled$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Depth (m)") + ylab("Salinity (PPT)")+coord_flip()+ scale_x_discrete(limits=rev)

ggsave(dep.orp,filename = "figures/EnvVariablesOnly/SSD_Salinity_Depth_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

dep.sulf<-ggplot(meta.all.scaled, aes(x=Depth_m, y=Sulfate_milliM,color=SampDate,group=SampDate)) + geom_point(size=3) + geom_line() + theme_bw()+
  labs(title="Sulfate by Depth & sample site",subtitle="Using Scaled Sulfate (milliM) Data",color="sample site")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="sample site",values=unique(meta.all.scaled$SampDate_Color[order(meta.all.scaled$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Depth (m)") + ylab("Sulfate (milliM)")+coord_flip()+ scale_x_discrete(limits=rev)

ggsave(dep.sulf,filename = "figures/EnvVariablesOnly/SSD_Sulfate_Depth_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

dep.do<-ggplot(meta.all.scaled, aes(x=Depth_m, y=DO_Percent_Local,color=SampDate,group=SampDate)) + geom_point(size=3) + geom_line() + theme_bw()+
  labs(title="Dissolved Oxygen by Depth & sample site",subtitle="Using Scaled DO (%) Data",color="sample site")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="sample site",values=unique(meta.all.scaled$SampDate_Color[order(meta.all.scaled$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Depth (m)") + ylab("DO (%)")+coord_flip()+ scale_x_discrete(limits=rev)

ggsave(dep.do,filename = "figures/EnvVariablesOnly/SSD_DO_Depth_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

dep.hs<-ggplot(meta.all.scaled, aes(x=Depth_m, y=Sulfide_microM,color=SampDate,group=SampDate)) + geom_point(size=3) + geom_line() + theme_bw()+
  labs(title="Sulfide by Depth & sample site",subtitle="Using Scaled Sulfate (microM) Data",color="sample site")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="sample site",values=unique(meta.all.scaled$SampDate_Color[order(meta.all.scaled$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Depth (m)") + ylab("Sulfide (microM)")+coord_flip()+ scale_x_discrete(limits=rev)

ggsave(dep.hs,filename = "figures/EnvVariablesOnly/SSD_Sulfide_Depth_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

dep.chlr<-ggplot(meta.all.scaled, aes(x=Depth_m, y=Chlorophyll_RFU,color=SampDate,group=SampDate)) + geom_point(size=3) + geom_line() + theme_bw()+
  labs(title="Chlorophyll by Depth & sample site",subtitle="Using Scaled Chlorophyll (RFU) Data",color="sample site")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="sample site",values=unique(meta.all.scaled$SampDate_Color[order(meta.all.scaled$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Depth (m)") + ylab("Chlorophyll (RFU)")+coord_flip()+ scale_x_discrete(limits=rev)

ggsave(dep.chlr,filename = "figures/EnvVariablesOnly/SSD_Chlorophyll_Depth_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

dep.temp<-ggplot(meta.all.scaled, aes(x=Depth_m, y=Temp_DegC,color=SampDate,group=SampDate)) + geom_point(size=3) + geom_line() + theme_bw()+
  labs(title="Temperature by Depth & sample site",subtitle="Using Scaled Chlorophyll (RFU) Data",color="sample site")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="sample site",values=unique(meta.all.scaled$SampDate_Color[order(meta.all.scaled$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Depth (m)") + ylab("Temperature (C)")+coord_flip()+ scale_x_discrete(limits=rev)

ggsave(dep.temp,filename = "figures/EnvVariablesOnly/SSD_Temp_Depth_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

# Compare variables to each other
dom.hs<-ggplot(meta.all.scaled, aes(x=Dissolved_OrganicMatter_RFU, y=Sulfide_microM,color=SampDate,group=SampDate)) + geom_point(size=3) + geom_line() + theme_bw()+
  labs(title="Sulfide ~ DOM",subtitle="Using Scaled Data",color="sample site")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="sample site",values=unique(meta.all.scaled$SampDate_Color[order(meta.all.scaled$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("DOM (RFU)") + ylab("Sulfide (microM)")

ggsave(dom.hs,filename = "figures/EnvVariablesOnly/SSD_Sulfide_DOM_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

do.hs<-ggplot(meta.all.scaled, aes(x=DO_Percent_Local, y=Sulfide_microM,color=SampDate,group=SampDate)) + geom_point(size=3) + geom_line() + theme_bw()+
  labs(title="Sulfide ~ DO%",subtitle="Using Scaled Data",color="sample site")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="sample site",values=unique(meta.all.scaled$SampDate_Color[order(meta.all.scaled$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("DO (%)") + ylab("Sulfide (microM)")

ggsave(do.hs,filename = "figures/EnvVariablesOnly/SSD_Sulfide_DO.Percent_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

orp.hs<-ggplot(meta.all.scaled, aes(x=ORP_mV, y=Sulfide_microM,color=SampDate,group=SampDate)) + geom_point(size=3) + geom_line() + theme_bw()+
  labs(title="Sulfide ~ ORP (mV)",subtitle="Using Scaled Data",color="sample site")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="sample site",values=unique(meta.all.scaled$SampDate_Color[order(meta.all.scaled$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("ORP (mV") + ylab("Sulfide (microM)")

ggsave(orp.hs,filename = "figures/EnvVariablesOnly/SSD_Sulfide_ORP_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

orp.dom<-ggplot(meta.all.scaled, aes(x=ORP_mV, y=Dissolved_OrganicMatter_RFU,color=SampDate,group=SampDate)) + geom_point(size=3) + geom_line() + theme_bw()+
  labs(title="DOM (RFU) ~ ORP (mV)",subtitle="Using Scaled Data",color="sample site")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="sample site",values=unique(meta.all.scaled$SampDate_Color[order(meta.all.scaled$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("DOM (RFU") + ylab("Sulfide (microM)")

ggsave(orp.dom,filename = "figures/EnvVariablesOnly/SSD_DOM_ORP_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

orp.do<-ggplot(meta.all.scaled, aes(x=ORP_mV, y=DO_Percent_Local,color=SampDate,group=SampDate)) + geom_point(size=3) + geom_line() + theme_bw()+
  labs(title="DO% ~ ORP (mV)",subtitle="Using Scaled Data",color="sample site")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="sample site",values=unique(meta.all.scaled$SampDate_Color[order(meta.all.scaled$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("ORP (mV") + ylab("DO %")

ggsave(orp.do,filename = "figures/EnvVariablesOnly/SSD_DO.Percent_ORP_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)
