#### Set WD & Load Libraries ####
getwd() # use setwd("path/to/files") if you are not in the right directory
#setwd("/Volumes/HLF_SSD/Aronson_Lab_Data/Salton_Sea/SaltonSeaDust")
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
  library(fst)
  library(plotly)
  library(htmlwidgets)
})

#### Load Global Env to Import Count/ASV Tables ####
load("data/Amplicon/SSDust_16S.V3V4_W23_Data_Ready.Rdata") # save global env to Rdata file
#load("data/Amplicon/SSD_16S_CLR_EucDist_Ready.Rdata")

head(b.dust.all)
bac.ASV_table[1:4,1:4]
bac.ASV_table[(nrow(bac.ASV_table)-4):(nrow(bac.ASV_table)),(ncol(bac.ASV_table)-4):(ncol(bac.ASV_table))] # last 4 rows & cols
head(dust_meta)

# create metadata dfs by categories for downstream comparisons
summer_meta<-dust_meta[dust_meta$Season_General=="Summer",]
meta.2020<-dust_meta[dust_meta$CollectionYear=="2020",]
meta.2021<-dust_meta[dust_meta$CollectionYear=="2021",]

#### PCoA on Raw Data - Sanity Check ####

# check rownames of CLR transformed ASV data & metadata
rownames(bac.ASV_table) %in% rownames(dust_meta)
dust_meta=dust_meta[rownames(bac.ASV_table[,-1]),] ## reorder metadata to match order of CLR data

# calculate our Euclidean distance matrix using CLR data
raw.b.euc_dist <- dist(bac.ASV_table[,-1], method = "euclidean")

# creating our hierarcical clustering dendrogram
raw.b.euc_clust <- hclust(raw.b.euc_dist, method="ward.D2")

# let's make it a little nicer...
raw.b.euc_dend <- as.dendrogram(raw.b.euc_clust, hang=0.2)
raw.b.dend_cols <- as.character(dust_meta$SampDate_Color[order.dendrogram(raw.b.euc_dend)])
labels_colors(raw.b.euc_dend) <- raw.b.dend_cols

## DO NOT RUN THIS LINE, THIS IS YOUR COLOR REFERENCE!!!!
#(August.2021="#ef781c",December.2021="#03045e",April.2022="#059c3f")

plot(raw.b.euc_dend, ylab="CLR Euclidean Distance",cex = 0.5) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 2)
#legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c("#ef781c","#03045e","#059c3f"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# PCOA w/ Euclidean distance matrix (of CLR data)
raw.b.pcoa <- pcoa(raw.b.euc_dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
#save.image("data/Amplicon/SSD_16S_CLR_EucDist_Ready.Rdata")

# The proportion of variances explained is in its element values$Relative_eig
raw.b.pcoa$values

# extract principal coordinates
raw.b.pcoa.vectors<-data.frame(raw.b.pcoa$vectors)
raw.b.pcoa.vectors$SampleID<-rownames(raw.b.pcoa$vectors)

# merge pcoa coordinates w/ metadata
raw.b.pcoa.meta<-merge(raw.b.pcoa.vectors, dust_meta, by.x="SampleID", by.y="SampleID")
raw.b.pcoa.meta$SampleMonth
raw.b.pcoa.meta$SampDate

head(raw.b.pcoa.meta)

head(raw.b.pcoa$values) # pull out Relative (Relative_eig) variation % to add to axes labels
# PC1 = 50.78%, PC2 = 18.29%

# drop outliars
#outliars<-c("BDC.D.7.27.21","PD.D.7.27.21","WI.D.9.18.21","WI.D.7.10.20")
#raw.b.pcoa.meta2<-raw.b.pcoa.meta[!(raw.b.pcoa.meta$SampleID %in% outliars),]

#### Visualize RAW PCoAs####

# create PCoA ggplot fig
ggplot(raw.b.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(SampleMonth),shape=CollectionYear), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Date",values=unique(raw.b.pcoa.meta$SampMonth_Color[order(raw.b.pcoa.meta$SampleMonth)]),labels=c(unique(raw.b.pcoa.meta$SampleMonth[order(raw.b.pcoa.meta$SampleMonth)]))) +
  xlab("PC1 [48.24%]") + ylab("PC2 [17.3%]")

ggplot(raw.b.pcoa.meta2, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site),shape=CollectionYear), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Date",values=unique(raw.b.pcoa.meta$Site_Color[order(raw.b.pcoa.meta$Site)]),labels=c(unique(raw.b.pcoa.meta$Site[order(raw.b.pcoa.meta$Site)]))) +
  xlab("PC1 [48.24%]") + ylab("PC2 [10.31%]")

# by collection period & site
ggplot(raw.b.pcoa.meta2, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year),shape=Site), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(raw.b.pcoa.meta$SCY_Color[order(raw.b.pcoa.meta$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
  xlab("PC1 [48.24%]") + ylab("PC2 [10.31%]")


#### CLR Transform All Comp Data ####
rownames(bac.ASV_table)
bac.ASV_table[1:4,1:4]

# CLR transformation of ASV table
# df must have rownames are SampleIDs, columns are ASV IDs for vegan functions below
b.clr<-decostand(bac.ASV_table[,-1],method = "clr", pseudocount = 1) #CLR transformation
b.clr[1:4,1:4]

#### Beta Diversity - All Data ####

# check rownames of CLR transformed ASV data & metadata
rownames(b.clr) %in% rownames(dust_meta)
dust_meta=dust_meta[rownames(b.clr),] ## reorder metadata to match order of CLR data

# calculate our Euclidean distance matrix using CLR data
b.euc_dist <- dist(b.clr, method = "euclidean")

# creating our hierarcical clustering dendrogram
b.euc_clust <- hclust(b.euc_dist, method="ward.D2")

# let's make it a little nicer...
b.euc_dend <- as.dendrogram(b.euc_clust, hang=0.02)
b.dend_cols <- as.character(dust_meta$SampDate_Color[order.dendrogram(b.euc_dend)])
labels_colors(b.euc_dend) <- b.dend_cols

colorset8 # color dendrogram by collection date
png(filename="figures/BetaDiversity/SSD_16S_CLR_EucDist_Dendrogram1.png",width = 7, height = 7, units = "in",res = 800)
plot(b.euc_dend, ylab="CLR Euclidean Distance", cex = 0.1,horiz=TRUE) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 2)
#legend("topleft",legend = colorset8$SampDate,cex=.8,col = colorset8$SampDate_Color,pch = 15, bty = "n")
dev.off()

# b.euc_dend1 <- as.dendrogram(b.euc_clust, hang=0.06)
# b.dend_cols1 <- as.character(dust_meta$SampMonth_Color[order.dendrogram(b.euc_dend1)])
# labels_colors(b.euc_dend1) <- b.dend_cols1
#
# colorset2 # color dendrogram by month of collection
# png(filename="figures/BetaDiversity/SSD_16S_CLR_EucDist_Dendrogram2.png",width = 7, height = 7, units = "in",res = 800)
# plot(b.euc_dend1, ylab="CLR Euclidean Distance", cex = 0.6) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 2)
# legend("topright",legend = colorset2$SampleMonth,cex=.8,col = colorset2$SampMonth_Color,pch = 15, bty = "n")
# dev.off()

# PCOA w/ Euclidean distance matrix (of CLR data)
b.pcoa <- pcoa(b.euc_dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
#save.image("data/Amplicon/SSD_16S_CLR_EucDist_Ready.Rdata")

# The proportion of variances explained is in its element values$Relative_eig
b.pcoa$values

# extract principal coordinates
b.pcoa.vectors<-data.frame(b.pcoa$vectors)
b.pcoa.vectors$SampleID<-rownames(b.pcoa$vectors)

# merge pcoa coordinates w/ metadata
b.pcoa.meta<-merge(b.pcoa.vectors, dust_meta, by.x="SampleID", by.y="SampleID")
b.pcoa.meta$SampleMonth
b.pcoa.meta$SampDate

head(b.pcoa.meta)

head(b.pcoa$values) # pull out Relative (Relative_eig) variation % to add to axes labels
# PC1 = 22.18%, PC2 = 10.31%

save.image("data/Amplicon/SSD_16S_CLR_EucDist_Ready.Rdata")

#### Visualize PCoAs - All Data ####

## *** all figures that end in _PCOA1 come from the same single PCoA
#data is just subsetted to understand how points are related to each other w/in & across timepoints

# create PCoA ggplot fig
pcoa1<-ggplot(b.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=SampDate,shape=Site,size=CollectionYear))+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Sample Date",size="Collection Year")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  scale_size_manual(values = c("2020" = 7, "2021"=4),labels=c("2020","2021")) +
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Date",values=unique(b.pcoa.meta$SampDate_Color[order(b.pcoa.meta$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  scale_shape_manual(values = c(7,10, 15,16)) +
  xlab("PC1 [22.18%]") + ylab("PC2 [10.31%]")

ggsave(pcoa1,filename = "figures/BetaDiversity/SSD_16S_CLR_SampDate_Site_Year_PCOA1.png", width=14, height=10, dpi=600)

# pcoa1a<-ggplot(b.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site),shape=CollectionYear), size=5)+theme_bw()+
#   labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Collection Date",values=unique(b.pcoa.meta$Site_Color[order(b.pcoa.meta$Site)]),labels=c(unique(b.pcoa.meta$Site[order(b.pcoa.meta$Site)]))) +
#   xlab("PC1 [22.18%]") + ylab("PC2 [10.31%]")
#
# ggsave(pcoa1a,filename = "figures/BetaDiversity/SSD_16S_CLR_Site_Year_PCOA1.png", width=12, height=10, dpi=600)
#
# # specific season
# pcoa3<-ggplot(b.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Season_Specific),shape=CollectionYear), size=5)+theme_bw()+
#   labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Season")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Season",values=unique(b.pcoa.meta$SeasonSpec_Color[order(b.pcoa.meta$Season_Specific)]),
#                      labels=c("Early.Summer"="Early Summer","Late.Summer"="Late Summer","Early.Fall"="Early Fall","Late.Fall"="Late Fall","Fall.Winter"="Fall-Winter")) +
#   xlab("PC1 [22.18%]") + ylab("PC2 [10.31%]")
#
# ggsave(pcoa3,filename = "figures/BetaDiversity/SSD_16S_CLR_season_specific_PCOA1.png", width=12, height=10, dpi=600)
#
# # by collection year & site
# pcoa4a<-ggplot(b.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(CollectionYear),shape=Site), size=5)+theme_bw()+
#   labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Collection Year")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Collection Year",values=unique(b.pcoa.meta$Year_Color[order(b.pcoa.meta$CollectionYear)])) +
#   xlab("PC1 [22.18%]") + ylab("PC2 [10.31%]")
#
# ggsave(pcoa4a,filename = "figures/BetaDiversity/SSD_16S_CLR_CollectionYear_Site_PCOA1.png", width=12, height=10, dpi=600)
#
# # by collection period & site
# pcoa5<-ggplot(b.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year),shape=Site), size=5)+theme_bw()+
#   labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Collection Period",values=unique(b.pcoa.meta$SCY_Color[order(b.pcoa.meta$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   xlab("PC1 [22.18%]") + ylab("PC2 [10.31%]")
#
# ggsave(pcoa5,filename = "figures/BetaDiversity/SSD_16S_CLR_SeasCollYr_Site_PCOA1.png", width=12, height=10, dpi=600)

# 3D PCoA

pltly.all.a<-plot_ly(b.pcoa.meta, x=~Axis.1,y=~Axis.2,z=~Axis.3, color = ~SampDate, colors = c(unique(b.pcoa.meta$SampDate_Color[order(b.pcoa.meta$SampDate)])),
        symbol=~Site,symbols = c("square-open", "circle","triangle-up-dot","diamond")) %>%
  layout(scene = list(xaxis = list(title = 'PC1 22.18%'),
                      yaxis = list(title = 'PC2 10.31%'),
                      zaxis = list(title = 'PC3 8.25%')))

saveWidget(widget = pltly.all.a, #the plotly object,
           file = "figures/BetaDiversity/SSD_16S_CLR_SampDate_CollYr_Site_3D_PCOA1.html", #the path & file name
           selfcontained = TRUE #creates a single html file
)

## visualizing only specific sites

# WI - by collection period
wi.pc1<-ggplot(b.pcoa.meta[b.pcoa.meta$Site=="WI",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year)), size=5)+theme_bw()+
  labs(title="PCoA: Wister Microbial Composition",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.pcoa.meta$SCY_Color[order(b.pcoa.meta$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
  xlab("PC1 [22.18%]") + ylab("PC2 [10.31%]")

ggsave(wi.pc1,filename = "figures/BetaDiversity/SSD_16S_CLR_SeasCollYr_Wister_PCOA1.png", width=12, height=10, dpi=600)

dp.pc1<-ggplot(b.pcoa.meta[b.pcoa.meta$Site=="DP",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year)), size=5)+theme_bw()+
  labs(title="PCoA: Dos Palmas Microbial Composition",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.pcoa.meta$SCY_Color[order(b.pcoa.meta$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
  xlab("PC1 [22.18%]") + ylab("PC2 [10.31%]")

ggsave(dp.pc1,filename = "figures/BetaDiversity/SSD_16S_CLR_SeasCollYr_DosPalmas_PCOA1.png", width=12, height=10, dpi=600)

s.r.pc1<-ggplot(b.pcoa.meta[b.pcoa.meta$Site=="SB" | b.pcoa.meta$Site=="RHB",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year),shape=Site), size=5)+theme_bw()+
  labs(title="PCoA: Sonny Bono & Red Hill Bay Microbial Composition",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.pcoa.meta$SCY_Color[order(b.pcoa.meta$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
  xlab("PC1 [22.18%]") + ylab("PC2 [10.31%]")

ggsave(s.r.pc1,filename = "figures/BetaDiversity/SSD_16S_CLR_SeasCollYr_SB_and_RHB_PCOA1.png", width=12, height=10, dpi=600)

pd.pc1<-ggplot(b.pcoa.meta[b.pcoa.meta$Site=="PD",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year)), size=5)+theme_bw()+
  labs(title="PCoA: Palm Desert Microbial Composition",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.pcoa.meta$SCY_Color[order(b.pcoa.meta$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
  xlab("PC1 [22.18%]") + ylab("PC2 [10.31%]")

ggsave(pd.pc1,filename = "figures/BetaDiversity/SSD_16S_CLR_SeasCollYr_PD_PCOA1.png", width=12, height=10, dpi=600)

bdc.pc1<-ggplot(b.pcoa.meta[b.pcoa.meta$Site=="BDC",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year)), size=5)+theme_bw()+
  labs(title="PCoA: Boyd Deep Canyon Microbial Composition",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.pcoa.meta$SCY_Color[order(b.pcoa.meta$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
  xlab("PC1 [22.18%]") + ylab("PC2 [10.31%]")

ggsave(bdc.pc1,filename = "figures/BetaDiversity/SSD_16S_CLR_SeasCollYr_BDC_PCOA1.png", width=12, height=10, dpi=600)

# visualizing pcoa by month per year
unique(b.pcoa.meta$SampDate)

j2020.p1<-ggplot(b.pcoa.meta[b.pcoa.meta$SampDate=="July.2020",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: July 2020 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.pcoa.meta$Site_Color[order(b.pcoa.meta$Site)])) +
  xlab("PC1 [22.18%]") + ylab("PC2 [10.31%]")

ggsave(j2020.p1,filename = "figures/BetaDiversity/SSD_16S_CLR_AllSites_July2020_PCOA1.png", width=12, height=10, dpi=600)

a2020.p1<-ggplot(b.pcoa.meta[b.pcoa.meta$SampDate=="August.2020",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: August 2020 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.pcoa.meta$Site_Color[order(b.pcoa.meta$Site)])) +
  xlab("PC1 [22.18%]") + ylab("PC2 [10.31%]")

ggsave(a2020.p1,filename = "figures/BetaDiversity/SSD_16S_CLR_AllSites_Aug2020_PCOA1.png", width=12, height=10, dpi=600)

o2020.p1<-ggplot(b.pcoa.meta[b.pcoa.meta$SampDate=="October.2020",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: October 2020 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.pcoa.meta$Site_Color[order(b.pcoa.meta$Site)])) +
  xlab("PC1 [22.18%]") + ylab("PC2 [10.31%]")

ggsave(o2020.p1,filename = "figures/BetaDiversity/SSD_16S_CLR_AllSites_Oct2020_PCOA1.png", width=12, height=10, dpi=600)

n2020.p1<-ggplot(b.pcoa.meta[b.pcoa.meta$SampDate=="November.2020",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: November 2020 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.pcoa.meta$Site_Color[order(b.pcoa.meta$Site)])) +
  xlab("PC1 [22.18%]") + ylab("PC2 [10.31%]")

ggsave(n2020.p1,filename = "figures/BetaDiversity/SSD_16S_CLR_AllSites_Nov2020_PCOA1.png", width=12, height=10, dpi=600)

ja2021.p1<-ggplot(b.pcoa.meta[b.pcoa.meta$SampDate=="July.2021" | b.pcoa.meta$SampDate=="August.2021",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: July & August 2021 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.pcoa.meta$Site_Color[order(b.pcoa.meta$Site)])) +
  xlab("PC1 [22.18%]") + ylab("PC2 [10.31%]")

ggsave(ja2021.p1,filename = "figures/BetaDiversity/SSD_16S_CLR_AllSites_July_Aug_2021_PCOA1.png", width=12, height=10, dpi=600)

s2021.p1<-ggplot(b.pcoa.meta[b.pcoa.meta$SampDate=="September.2021",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: September 2021 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.pcoa.meta$Site_Color[order(b.pcoa.meta$Site)])) +
  xlab("PC1 [22.18%]") + ylab("PC2 [10.31%]")

ggsave(s2021.p1,filename = "figures/BetaDiversity/SSD_16S_CLR_AllSites_Sept2021_PCOA1.png", width=12, height=10, dpi=600)

d2021.p1<-ggplot(b.pcoa.meta[b.pcoa.meta$SampDate=="December.2021",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: December 2021 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.pcoa.meta$Site_Color[order(b.pcoa.meta$Site)])) +
  xlab("PC1 [22.18%]") + ylab("PC2 [10.31%]")

ggsave(d2021.p1,filename = "figures/BetaDiversity/SSD_16S_CLR_AllSites_Dec2021_PCOA1.png", width=12, height=10, dpi=600)

## visualizing by year
b.pcoa1.20<-b.pcoa.meta[b.pcoa.meta$CollectionYear=="2020",]
b.pcoa1.21<-b.pcoa.meta[b.pcoa.meta$CollectionYear=="2021",]

# 2020
twntytwnty.pc1<-ggplot(b.pcoa1.20, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year),shape=Site), size=5)+theme_bw()+
  labs(title="PCoA: Microbial Composition in 2020",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.pcoa1.20$SCY_Color[order(b.pcoa1.20$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
  xlab("PC1 [22.18%]") + ylab("PC2 [10.31%]")

ggsave(twntytwnty.pc1,filename = "figures/BetaDiversity/SSD_16S_CLR_SeasCollYr_2020_PCOA1.png", width=12, height=10, dpi=600)

# 2021
twntytwnty1.pc1<-ggplot(b.pcoa1.21, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year),shape=Site), size=5)+theme_bw()+
  labs(title="PCoA: Microbial Composition in 2021",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.pcoa1.21$SCY_Color[order(b.pcoa1.21$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
  xlab("PC1 [22.18%]") + ylab("PC2 [10.31%]")

ggsave(twntytwnty1.pc1,filename = "figures/BetaDiversity/SSD_16S_CLR_SeasCollYr_2021_PCOA1.png", width=12, height=10, dpi=600)

#### Visualize Location & Date with PC Axes ####
# this idea was suggested by Dr. Will Porter as a way to see if we can identify compositional similiarites or dissimilarites by site and sample collection date
# using PC axes as input data, where x axis is site and y axis is the sample date in the heat map

dust.time.site<-subset(dust_meta, select=c(SampleID, Site, SampDate))
b.pcoa.dts<-merge(b.pcoa.vectors, dust.time.site, by.x="SampleID", by.y="SampleID")

head(b.pcoa.dts)

ggplot(b.pcoa.dts, aes(Site, SampDate)) +
  geom_tile(aes(fill = Axis.1)) +
  geom_text(aes(label = round(Axis.1, 1))) +
  scale_fill_gradient(low = "white", high = "red")

# For heatmap color gradient, PC1
max(b.pcoa.dts$Axis.1, na.rm=TRUE)
max(b.pcoa.dts$Axis.1, na.rm=TRUE)/2
min(b.pcoa.dts$Axis.1, na.rm=TRUE)

pcoa.axis1.hm<-ggplot(b.pcoa.dts, aes(Site, SampDate, Axis.1, fill=Axis.1)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="blue", high="red",labels=c("25","-36.5","-95"),breaks=c(25,-36,-95)) + labs(title="Salton Sea Dust Microbial PCoA Axis 1 by Sample Date & Site",subtitle="Using CLR-Transformed 16S Data",fill="PC1 Values") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0)) + geom_text(aes(label = round(Axis.1, 2)),size=9)

ggsave(pcoa.axis1.hm,filename = "figures/BetaDiversity/SSD_16S_CLR_PCoA_PC1_Site_by_SampDate.png", width=18, height=13, dpi=600)

# For heatmap color gradient, PC2
max(b.pcoa.dts$Axis.2, na.rm=TRUE)
max(b.pcoa.dts$Axis.2, na.rm=TRUE)/2
min(b.pcoa.dts$Axis.2, na.rm=TRUE)

pcoa.axis2.hm<-ggplot(b.pcoa.dts, aes(Site, SampDate, Axis.2, fill=Axis.2)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="blue", high="red",labels=c("45","-36","-115"),breaks=c(45,-36,-115)) + labs(title="Salton Sea Dust Microbial PCoA Axis 2 by Sample Date & Site",subtitle="Using CLR-Transformed 16S Data",fill="PC2 Values") +
  geom_text(aes(label = round(Axis.2, 2)),size=9) +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(pcoa.axis2.hm,filename = "figures/BetaDiversity/SSD_16S_CLR_PCoA_PC2_Site_by_SampDate.png", width=18, height=13, dpi=600)


s.t.pcoa<-melt(b.pcoa.dts[,-1], by=c("Site","SampDate"))
colnames(s.t.pcoa)[which(names(s.t.pcoa) == "variable")] <- "PCoA_Axis"
colnames(s.t.pcoa)[which(names(s.t.pcoa) == "value")] <- "PC_Axis_Value"

ggplot(s.t.pcoa, aes(Site, SampDate, PCoA_Axis, fill=PC_Axis_Value)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="blue", high="red",labels=c(),breaks=c()) + labs(title="PCoA by Sample Date & Site",subtitle="Using CLR-Transformed 16S Data",fill="PC Axis Value") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(sulf.hm1a,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_SampID_by_Function_heatmap.png", width=18, height=13, dpi=600)

#### Homogeneity of Variance & PERMANOVA tests - Composition by Groups ####
## betadisper to look at homogeneity of group dispersions (aka variance) when considering multiple variables
# multivariate analogue to Levene's test of homogeneity of variances
# program finds spatial median or centroid of the group, & compare distances of group to centroid/spatial median via ANOVA

#While PERMANOVA tests differences in group means (analogous to MANOVA),
## a related test called PERMDISP can be used to evaluate homogeneity of group dispersion
#(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
## * need a distance matrix!

head(dust_meta)
head(b.clr)
rownames(dust_meta) %in% rownames(b.clr) #b.clr was used to make the distance matrix b.euc_dist

# first by compare dispersions by sampling date
b.disper1<-betadisper((vegdist(b.clr,method="euclidean")), dust_meta$Site)
b.disper1
# NOTE: SB and RHB have less samples than other sites since they are supposed to represent similar locations
# maybe we should remove these sites and then rerun the PERMANOVA?

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:

permutest(b.disper1, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#         PD     BDC      DP    WI
# PD          0.16400 0.70100 0.295
# BDC 0.17952         0.31300 0.868
# DP  0.71671 0.33040         0.511
# WI  0.31744 0.86631 0.49913

anova(b.disper1) # p = 0.5283 --> accept the Null H, spatial medians (a measure of dispersion) are NOT significantly difference across sample dates
# ANOVA adjusted p-value
aov.beta.p1<-anova(b.disper1)[["Pr(>F)"]] # get p values from ANOVA
p.adjust(aov.beta.p1,method="bonferroni",n=length(aov.beta.p1))

TukeyHSD(b.disper1) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other
#         diff       lwr       upr     p adj
# BDC-PD -25.942841 -81.04497 29.15928 0.5725376
# DP-PD   -7.476742 -62.57887 47.62538 0.9816864
# WI-PD  -22.533183 -77.63531 32.56894 0.6761934
# DP-BDC  18.466099 -36.63603 73.56822 0.7920997
# WI-BDC   3.409658 -51.69247 58.51178 0.9981808
# WI-DP  -15.056441 -70.15857 40.04568 0.8741546

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

pnova1<-adonis2(b.clr ~ Site,data=dust_meta,method = "euclidean",by="terms",permutations=1000)
pnova1 # p-value = 0.05694 .
p.adjust(pnova1$`Pr(>F)`,method="bonferroni",n=length(pnova1$`Pr(>F)`)) # adjusted pval
# 0.1708292

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

b.clr.dist = (vegdist(b.clr, "euclidean", na.rm = TRUE)) #distance matrix using Bray's dissimilarity index for trait distribution (traits of interest only)
pair.mod1<-pairwise.adonis(b.clr.dist,dust_meta$Site, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
pair.mod1
#     pairs Df SumsOfSqs  F.Model        R2     p.value     p.adjusted sig
# 1 BDC vs DP  1 11164.940 1.287096 0.09686810   0.095      0.570
# 2 BDC vs PD  1 12777.753 1.365205 0.10214622   0.137      0.822
# 3 BDC vs WI  1  7719.261 1.050797 0.08051591   0.292      1.000
# 4  DP vs PD  1 15806.183 1.401118 0.10455235   0.119      0.714
# 5  DP vs WI  1  9748.730 1.051911 0.08059441   0.299      1.000
# 6  PD vs WI  1 15866.910 1.594232 0.11727266   0.080      0.480

# Visualize dispersions
png('figures/BetaDiversity/SSD_pcoa_betadispersion_site.png',width = 700, height = 600, res=100)
plot(b.disper1,main = "Centroids and Dispersion based on Aitchison Distance", col=colorset6$Site_Color)
dev.off()

png('figures/BetaDiversity/SSD_boxplot_centroid_distance_site.png',width = 700, height = 600, res=100)
boxplot(b.disper1,xlab="By Site", main = "Distance to Centroid by Category", sub="Based on Aitchison Distance", col=colorset6$Site_Color)
dev.off()

# Next compare dispersions by sampdate
b.disper2<-betadisper((vegdist(b.clr,method="euclidean")), dust_meta$SampDate)
b.disper2

permutest(b.disper2, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons

anova(b.disper2) # p = 0.3231 --> accept the Null H, spatial medians are NOT significantly difference across sample dates

# ANOVA adjusted p-value
aov.beta.p2<-anova(b.disper2)["Pr(>F)"] # get p values from ANOVA
p.adjust(aov.beta.p2,method="bonferroni",n=length(aov.beta.p1))

TukeyHSD(b.disper2) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other
# no sig results
#                                   diff        lwr       upr     p adj
# August.2020-July.2020         -27.872187 -127.70519  71.96081 0.9777273
# October.2020-July.2020        -41.096631 -140.92963  58.73637 0.8523812
# November.2020-July.2020        19.080104  -80.75290 118.91310 0.9976177
# July.2021-July.2020           -14.885255 -122.71722  92.94671 0.9997059
# August.2021-July.2020         -84.565119 -242.41495  73.28471 0.6234200
# September.2021-July.2020       -2.336578 -102.16958  97.49642 1.0000000
# December.2021-July.2020         2.605036  -97.22796 102.43804 1.0000000
# October.2020-August.2020      -13.224445 -113.05744  86.60856 0.9997760
# November.2020-August.2020      46.952291  -52.88071 146.78529 0.7531556
# July.2021-August.2020          12.986931  -94.84503 120.81890 0.9998811
# August.2021-August.2020       -56.692932 -214.54276 101.15690 0.9193239
# September.2021-August.2020     25.535609  -74.29739 125.36861 0.9863526
# December.2021-August.2020      30.477223  -69.35578 130.31022 0.9640410
# November.2020-October.2020     60.176735  -39.65626 160.00973 0.4866636
# July.2021-October.2020         26.211376  -81.62059 134.04334 0.9898224
# August.2021-October.2020      -43.468487 -201.31832 114.38134 0.9793488
# September.2021-October.2020    38.760054  -61.07295 138.59305 0.8850569
# December.2021-October.2020     43.701667  -56.13133 143.53467 0.8110073
# July.2021-November.2020       -33.965359 -141.79732  73.86660 0.9577071
# August.2021-November.2020    -103.645223 -261.49505  54.20461 0.3847021
# September.2021-November.2020  -21.416682 -121.24968  78.41632 0.9951533
# December.2021-November.2020   -16.475068 -116.30807  83.35793 0.9990597
# August.2021-July.2021         -69.679863 -232.70647  93.34674 0.8279773
# September.2021-July.2021       12.548678  -95.28329 120.38064 0.9999055
# December.2021-July.2021        17.490291  -90.34167 125.32226 0.9991586
# September.2021-August.2021     82.228541  -75.62129 240.07837 0.6537122
# December.2021-August.2021      87.170155  -70.67968 245.01999 0.5894398
# December.2021-September.2021    4.941614  -94.89139 104.77461 0.9999997

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

pnova2<-adonis2(b.clr ~ SampDate,data=dust_meta,method = "euclidean",by="terms",permutations=1000)
pnova2 # p-value = 0.1798
p.adjust(pnova2$`Pr(>F)`,method="bonferroni",n=length(pnova2$`Pr(>F)`)) # adjusted pval

#b.clr.dist = (vegdist(b.clr, "euclidean", na.rm = TRUE)) #distance matrix using Bray's dissimilarity index for trait distribution (traits of interest only)
pair.mod2<-pairwise.adonis(b.clr.dist,dust_meta$SampDate, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
pair.mod2
# none are significantly different
#                               pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
# 1    October.2020 vs November.2020  1 15381.677 1.6867130 0.21943228   0.052          1
# 2    October.2020 vs December.2021  1  7729.997 1.0831064 0.15291404   0.215          1
# 3        October.2020 vs July.2021  1  6411.434 1.1850124 0.19159419   0.192          1
# 4        October.2020 vs July.2020  1  6036.391 0.8133577 0.11937693   0.611          1
# 5      October.2020 vs August.2020  1  3288.238 0.6644848 0.09970535   0.934          1
# 6   October.2020 vs September.2021  1  8709.898 1.2741996 0.17516699   0.198          1
# 7      October.2020 vs August.2021  1 15472.762 4.2634927 0.58697556   0.200          1
# 8   November.2020 vs December.2021  1 13812.571 1.0938857 0.15420120   0.355          1
# 9       November.2020 vs July.2021  1 15196.800 1.2665406 0.20211161   0.209          1
# 10      November.2020 vs July.2020  1 14110.899 1.0928717 0.15408028   0.282          1
# 11    November.2020 vs August.2020  1 14048.043 1.3457597 0.18320225   0.177          1
# 12 November.2020 vs September.2021  1 13539.700 1.0984866 0.15474941   0.285          1
# 13    November.2020 vs August.2021  1 17233.485 1.1796070 0.28222917   0.400          1
# 14      December.2021 vs July.2021  1  9086.798 0.9445996 0.15890046   0.427          1
# 15      December.2021 vs July.2020  1  9408.330 0.8608341 0.12547076   0.691          1
# 16    December.2021 vs August.2020  1  7948.537 0.9399542 0.13544098   0.394          1
# 17 December.2021 vs September.2021  1 10131.283 0.9794989 0.14033943   0.468          1
# 18    December.2021 vs August.2021  1 15053.373 1.4141758 0.32037143   0.200          1
# 19          July.2021 vs July.2020  1  8460.775 0.8493588 0.14520546   0.627          1
# 20        July.2021 vs August.2020  1  6698.761 0.9578216 0.16076708   0.472          1
# 21     July.2021 vs September.2021  1  9396.839 1.0149767 0.16874158   0.399          1
# 22        July.2021 vs August.2021  1 15255.929 1.8875497 0.48553712   0.250          1
# 23        July.2020 vs August.2020  1  5269.560 0.6028559 0.09130229   0.795          1
# 24     July.2020 vs September.2021  1  9213.081 0.8668667 0.12623905   0.510          1
# 25        July.2020 vs August.2021  1 16005.070 1.4272385 0.32237670   0.400          1
# 26   August.2020 vs September.2021  1  7724.031 0.9471518 0.13633671   0.445          1
# 27      August.2020 vs August.2021  1 15520.164 2.4761037 0.45216523   0.400          1
# 28   September.2021 vs August.2021  1 15213.271 1.5149584 0.33554206   0.400          1

# Visualize dispersions
png('figures/BetaDiversity/SSD_pcoa_betadispersion_SampDate.png',width = 700, height = 600, res=100)
plot(b.disper2,main = "Centroids and Dispersion based on Aitchison Distance", col=colorset8$SampDate_Color)
dev.off()

png('figures/BetaDiversity/SSD_boxplot_centroid_distance_SampDate.png',width = 1200, height = 600, res=100)
boxplot(b.disper2,xlab="Sample Collection Date", main = "Distance to Centroid by Category", sub="Based on Aitchison Distance", col=colorset8$SampDate_Color)
dev.off()

### now compare dispersions by site + year
b.disper3<-betadisper((vegdist(b.clr,method="euclidean")), interaction(dust_meta$Site,dust_meta$CollectionYear,sep="."))
b.disper3
# NOTE: SB and RHB have less samples than other sites since they are supposed to represent similar locations
# maybe we should remove these sites and then rerun the PERMANOVA?

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:

permutest(b.disper3, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#             PD.2020  BDC.2020   DP.2020   WI.2020   PD.2021  BDC.2021   DP.2021 WI.2021
# PD.2020            0.1930000 0.5050000 0.1500000 0.6210000 0.0060000 0.5360000   0.313
# BDC.2020 0.1768387           0.6470000 0.9790000 0.6730000 0.6930000 0.3580000   0.467
# DP.2020  0.5043294 0.6649208           0.6630000 0.9610000 0.8630000 0.6940000   0.868
# WI.2020  0.1412316 0.9707172 0.6328157           0.6430000 0.6530000 0.3230000   0.402
# PD.2021  0.6391761 0.6765114 0.9541869 0.6490555           0.8280000 0.8250000   0.930
# BDC.2021 0.0055398 0.7135129 0.8675635 0.6615068 0.8420818           0.1340000   0.337
# DP.2021  0.5120869 0.3497882 0.7127328 0.3034816 0.8164188 0.1170939             0.704
# WI.2021  0.2884601 0.4605955 0.8486044 0.4114531 0.9284479 0.3315710 0.6936403

anova(b.disper3) # p = 0.8844 --> accept the Null H, spatial medians (a measure of dispersion) are NOT significantly difference across sample dates
# ANOVA adjusted p-value
aov.beta.p3<-anova(b.disper3)[["Pr(>F)"]] # get p values from ANOVA
p.adjust(aov.beta.p3,method="bonferroni",n=length(aov.beta.p3))

TukeyHSD(b.disper3) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other
#                     diff       lwr       upr     p adj
# BDC.2020-PD.2020  -35.219153 -128.55932  58.12102 0.8987716
# DP.2020-PD.2020   -19.275390 -112.61556  74.06478 0.9961562
# WI.2020-PD.2020   -36.399334 -129.73950  56.94084 0.8828134
# PD.2021-PD.2020   -16.525579 -117.34449  84.29333 0.9990997
# BDC.2021-PD.2020  -24.836004 -125.65491  75.98290 0.9890033
# DP.2021-PD.2020    -6.571516 -107.39042  94.24739 0.9999982
# WI.2021-PD.2020   -12.642545 -113.46145  88.17636 0.9998443
# DP.2020-BDC.2020   15.943763  -77.39641 109.28393 0.9988276
# WI.2020-BDC.2020   -1.180181  -94.52035  92.15999 1.0000000
# PD.2021-BDC.2020   18.693573  -82.12533 119.51248 0.9980291
# BDC.2021-BDC.2020  10.383149  -90.43576 111.20206 0.9999585
# DP.2021-BDC.2020   28.647636  -72.17127 129.46654 0.9754786
# WI.2021-BDC.2020   22.576608  -78.24230 123.39551 0.9937256
# WI.2020-DP.2020   -17.123944 -110.46411  76.21623 0.9981564
# PD.2021-DP.2020     2.749810  -98.06910 103.56872 1.0000000
# BDC.2021-DP.2020   -5.560614 -106.37952  95.25829 0.9999994
# DP.2021-DP.2020    12.703873  -88.11503 113.52278 0.9998393
# WI.2021-DP.2020     6.632845  -94.18606 107.45175 0.9999981
# PD.2021-WI.2020    19.873755  -80.94515 120.69266 0.9971130
# BDC.2021-WI.2020   11.563331  -89.25558 112.38224 0.9999143
# DP.2021-WI.2020    29.827818  -70.99109 130.64672 0.9695284
# WI.2021-WI.2020    23.756789  -77.06212 124.57570 0.9915131
# BDC.2021-PD.2021   -8.310424 -116.09037  99.46952 0.9999942
# DP.2021-PD.2021     9.954063  -97.82588 117.73401 0.9999802
# WI.2021-PD.2021     3.883035 -103.89691 111.66298 1.0000000
# DP.2021-BDC.2021   18.264487  -89.51546 126.04443 0.9988856
# WI.2021-BDC.2021   12.193459  -95.58649 119.97340 0.9999218
# WI.2021-DP.2021    -6.071028 -113.85097 101.70892 0.9999993

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

pnova3<-adonis2(b.clr ~ Site*CollectionYear,data=dust_meta,method = "euclidean",by="terms",permutations=1000)
pnova3
#                       Df SumOfSqs      R2      F  Pr(>F)
# Site                 3    36542 0.14051 1.3025 0.06793 .
# CollectionYear       1    10382 0.03992 1.1101 0.23776
# Site:CollectionYear  3    26105 0.10038 0.9304 0.57343
# Residual            20   187040 0.71920
# Total               27   260069 1.00000

p.adjust(pnova3$`Pr(>F)`,method="bonferroni",n=length(pnova3$`Pr(>F)`)) # adjusted pval
# 0.3396603 1.0000000 1.0000000        NA        NA

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

b.clr.dist = (vegdist(b.clr, "euclidean", na.rm = TRUE)) #distance matrix using Bray's dissimilarity index for trait distribution (traits of interest only)
pair.mod3<-pairwise.adonis(b.clr.dist,interaction(dust_meta$Site,dust_meta$CollectionYear,sep="."), p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
pair.mod3
#                   pairs Df SumsOfSqs  F.Model        R2     p.value     p.adjusted sig
# 1  BDC.2020 vs BDC.2021  1  6157.411 0.8960010 0.1519676   0.528          1
# 2   BDC.2020 vs DP.2020  1  6893.020 0.8708630 0.1267473   0.740          1
# 3   BDC.2020 vs DP.2021  1 12916.728 1.4862864 0.2291429   0.085          1
# 4   BDC.2020 vs PD.2020  1 14642.101 1.5596363 0.2063110   0.058          1
# 5   BDC.2020 vs PD.2021  1  7781.978 0.9108059 0.1540917   0.484          1
# 6   BDC.2020 vs WI.2020  1  5893.646 0.9744444 0.1397164   0.477          1
# 7   BDC.2020 vs WI.2021  1  9244.302 1.1439425 0.1861903   0.277          1
# 8   BDC.2021 vs DP.2020  1  8597.794 0.9803923 0.1639344   0.469          1
# 9   BDC.2021 vs DP.2021  1 11351.797 1.1405742 0.2218768   0.300          1
# 10  BDC.2021 vs PD.2020  1 12777.920 1.2126325 0.1951882   0.259          1
# 11  BDC.2021 vs PD.2021  1  6031.447 0.6173741 0.1337067   0.900          1
# 12  BDC.2021 vs WI.2020  1  7071.317 1.0829943 0.1780364   0.276          1
# 13  BDC.2021 vs WI.2021  1  8736.029 0.9505179 0.1920037   0.600          1
# 14   DP.2020 vs DP.2021  1 10635.080 1.0044226 0.1672805   0.409          1
# 15   DP.2020 vs PD.2020  1 15559.884 1.4184650 0.1912073   0.193          1
# 16   DP.2020 vs PD.2021  1  9177.831 0.8789591 0.1495093   0.648          1
# 17   DP.2020 vs WI.2020  1  7371.709 0.9662006 0.1386984   0.437          1
# 18   DP.2020 vs WI.2021  1  9115.108 0.9134531 0.1544703   0.600          1
# 19   DP.2021 vs PD.2020  1 16531.591 1.3379575 0.2111023   0.201          1
# 20   DP.2021 vs PD.2021  1 11090.241 0.9209140 0.1871429   0.600          1
# 21   DP.2021 vs WI.2020  1 12357.660 1.4803289 0.2284342   0.053          1
# 22   DP.2021 vs WI.2021  1 10707.981 0.9340578 0.1893082   0.700          1
# 23   PD.2020 vs PD.2021  1 10750.622 0.8805278 0.1497362   0.400          1
# 24   PD.2020 vs WI.2020  1 17885.412 1.9648737 0.2466924   0.113          1
# 25   PD.2020 vs WI.2021  1 14517.712 1.2359363 0.1981958   0.198          1
# 26   PD.2021 vs WI.2020  1  8643.953 1.0539649 0.1740950   0.339          1
# 27   PD.2021 vs WI.2021  1  9129.881 0.8093326 0.1682838   0.900          1
# 28   WI.2020 vs WI.2021  1  8943.307 1.1557049 0.1877453   0.170          1

# Visualize dispersions
png('figures/BetaDiversity/SSD_pcoa_betadispersion_site_by_year.png',width = 700, height = 600, res=100)
plot(b.disper3,main = "Centroids and Dispersion based on Aitchison Distance", col=colorset6$Site_Color)
dev.off()

png('figures/BetaDiversity/SSD_boxplot_centroid_distance_site_by_year.png',width = 900, height = 600, res=100)
boxplot(b.disper3,xlab="By Site x Collection Year", main = "Distance to Centroid by Category", sub="Based on Aitchison Distance", col=colorset6$Site_Color)
dev.off()

### now compare dispersions within summer only, site + year

# pulling out summer samples only
b.clr[which(rownames(b.clr) %in% rownames(summer_meta)),][1:10,1:5] # does our indexing idea work? yes!
summ.clr<-b.clr[which(rownames(b.clr) %in% rownames(summer_meta)),]

b.disper4<-betadisper((vegdist(summ.clr,method="euclidean")), interaction(summer_meta$Site,summer_meta$CollectionYear,sep="."))
b.disper4
# NOTE: SB and RHB have less samples than other sites since they are supposed to represent similar locations
# maybe we should remove these sites and then rerun the PERMANOVA?

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:

permutest(b.disper4, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#             PD.2020  BDC.2020   DP.2020   WI.2020   PD.2021  BDC.2021   DP.2021 WI.2021
# PD.2020            0.0260000 0.1620000 0.1920000 0.9130000 0.1250000 0.2750000   0.161
# BDC.2020 0.0082168           0.5290000 0.7150000 0.0100000 0.0640000 0.0210000   0.035
# DP.2020  0.1684466 0.5682766           0.9150000 0.2470000 0.6970000 0.4480000   0.572
# WI.2020  0.1990247 0.7277465 0.9183499           0.2850000 0.6940000 0.4750000   0.567
# PD.2021  0.9238252 0.0050571 0.2397645 0.2846356           0.0010000 0.0010000   0.001
# BDC.2021 0.1004492 0.0425528 0.7163219 0.6935400 0.0000000           0.0010000   0.002
# DP.2021  0.2610152 0.0155220 0.4729980 0.4946793 0.0000000 0.0000000             0.001
# WI.2021  0.1579286 0.0249566 0.5873339 0.5895183 0.0000000 0.0000000 0.0000000

anova(b.disper4) # p = 0.2164 --> accept the Null H, spatial medians (a measure of dispersion) are NOT significantly difference across sample dates
# ANOVA adjusted p-value
aov.beta.p4<-anova(b.disper4)[["Pr(>F)"]] # get p values from ANOVA
p.adjust(aov.beta.p4,method="bonferroni",n=length(aov.beta.p4))

TukeyHSD(b.disper4) # tells us which summer sites + years /category's dispersion MEANS are significantly different than each other
#                     diff       lwr       upr     p adj
# BDC.2020-PD.2020  -53.821951 -129.49858  21.85468 0.2529067
# DP.2020-PD.2020   -39.840723 -115.51735  35.83590 0.5714260
# WI.2020-PD.2020   -43.604545 -119.28117  32.07208 0.4720854
# PD.2021-PD.2020     1.265656  -83.34339  85.87470 1.0000000
# BDC.2021-PD.2020  -28.614849 -113.22389  55.99419 0.9094776
# DP.2021-PD.2020   -16.834014 -101.44306  67.77503 0.9944790
# WI.2021-PD.2020   -22.810660 -107.41970  61.79838 0.9697823
# DP.2020-BDC.2020   13.981227  -61.69540  89.65785 0.9964611
# WI.2020-BDC.2020   10.217406  -65.45922  85.89403 0.9995079
# PD.2021-BDC.2020   55.087607  -29.52143 139.69665 0.3394887
# BDC.2021-BDC.2020  25.207101  -59.40194 109.81614 0.9499478
# DP.2021-BDC.2020   36.987936  -47.62111 121.59698 0.7521783
# WI.2021-BDC.2020   31.011290  -53.59775 115.62033 0.8721998
# WI.2020-DP.2020    -3.763822  -79.44045  71.91281 0.9999994
# PD.2021-DP.2020    41.106380  -43.50266 125.71542 0.6550537
# BDC.2021-DP.2020   11.225874  -73.38317  95.83492 0.9995605
# DP.2021-DP.2020    23.006709  -61.60233 107.61575 0.9684095
# WI.2021-DP.2020    17.030063  -67.57898 101.63910 0.9940872
# PD.2021-WI.2020    44.870201  -39.73884 129.47924 0.5635184
# BDC.2021-WI.2020   14.989696  -69.61935  99.59874 0.9972598
# DP.2021-WI.2020    26.770531  -57.83851 111.37957 0.9332230
# WI.2021-WI.2020    20.793885  -63.81516 105.40293 0.9815588
# BDC.2021-PD.2021  -29.880506 -122.56507  62.80406 0.9271368
# DP.2021-PD.2021   -18.099671 -110.78423  74.58489 0.9950595
# WI.2021-PD.2021   -24.076317 -116.76088  68.60825 0.9751454
# DP.2021-BDC.2021   11.780835  -80.90373 104.46540 0.9996677
# WI.2021-BDC.2021    5.804189  -86.88037  98.48875 0.9999972
# WI.2021-DP.2021    -5.976646  -98.66121  86.70792 0.9999966

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

pnova4<-adonis2(summ.clr ~ Site*CollectionYear,data=summer_meta,method = "euclidean",by="terms",permutations=1000)
pnova4
#                       Df SumOfSqs      R2      F  Pr(>F)
# Site                 3    31186 0.20099 1.3628 0.07193 .
# CollectionYear       1    10451 0.06736 1.3701 0.13387
# Site:CollectionYear  3    21992 0.14173 0.9610 0.49650
# Residual            12    91534 0.58992
# Total               19   155163 1.00000

p.adjust(pnova4$`Pr(>F)`,method="bonferroni",n=length(pnova4$`Pr(>F)`)) # adjusted pval

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

summ.clr.dist = (vegdist(summ.clr, "euclidean", na.rm = TRUE)) #distance matrix using Bray's dissimilarity index for trait distribution (traits of interest only)
pair.mod4<-pairwise.adonis(summ.clr.dist,interaction(summer_meta$Site,summer_meta$CollectionYear,sep="."), p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
pair.mod4
#                   pairs Df SumsOfSqs  F.Model        R2     p.value     p.adjusted sig
# 1  BDC.2020 vs BDC.2021  1  6617.165 1.6784236 0.3587584 0.1000000          1
# 2   BDC.2020 vs DP.2020  1  3784.424 1.1387840 0.2216057 0.1000000          1
# 3   BDC.2020 vs DP.2021  1 10539.027 2.1043550 0.4122666 0.1000000          1
# 4   BDC.2020 vs PD.2020  1 14581.317 2.0093638 0.3343721 0.1000000          1
# 5   BDC.2020 vs PD.2021  1 13063.954 1.8646618 0.3833076 0.1000000          1
# 6   BDC.2020 vs WI.2020  1  4157.002 1.2875871 0.2435113 0.1000000          1
# 7   BDC.2020 vs WI.2021  1  7141.092 1.6067557 0.3487825 0.1000000          1
# 8   BDC.2021 vs DP.2020  1  5607.339 1.0001559 0.2500292 0.4000000          1
# 9   BDC.2021 vs DP.2021  1  8701.261 0.9380843 0.3192843 0.6666667          1
# 10  BDC.2021 vs PD.2020  1 10482.593 0.9660406 0.2435781 0.4000000          1
# 11  BDC.2021 vs PD.2021  1  7777.384 0.6337310 0.2406210 0.6666667          1
# 12  BDC.2021 vs WI.2020  1  6222.593 1.1354671 0.2745680 0.2000000          1
# 13  BDC.2021 vs WI.2021  1  6792.955 0.8058175 0.2871953 1.0000000          1
# 14   DP.2020 vs DP.2021  1  8758.931 1.3127548 0.3043889 0.2000000          1
# 15   DP.2020 vs PD.2020  1 11476.255 1.3494069 0.2522536 0.3000000          1
# 16   DP.2020 vs PD.2021  1 10905.372 1.2578213 0.2954143 0.2000000          1
# 17   DP.2020 vs WI.2020  1  4100.205 0.9159392 0.1863203 0.5000000          1
# 18   DP.2020 vs WI.2021  1  5977.509 0.9785730 0.2459608 0.3000000          1
# 19   DP.2021 vs PD.2020  1 13853.945 1.1625558 0.2792889 0.3000000          1
# 20   DP.2021 vs PD.2021  1 12072.352 0.8703341 0.3032170 0.6666667          1
# 21   DP.2021 vs WI.2020  1  9668.344 1.4770039 0.3299090 0.2000000          1
# 22   DP.2021 vs WI.2021  1  8869.441 0.8844273 0.3066214 0.6666667          1
# 23   PD.2020 vs PD.2021  1 10322.990 0.7418781 0.1982636 0.5000000          1
# 24   PD.2020 vs WI.2020  1 14160.930 1.6838261 0.2962487 0.3000000          1
# 25   PD.2020 vs WI.2021  1 10095.950 0.8892746 0.2286479 0.5000000          1
# 26   PD.2021 vs WI.2020  1 12527.185 1.4662329 0.3282930 0.2000000          1
# 27   PD.2021 vs WI.2021  1  9882.088 0.7586860 0.2750172 0.6666667          1
# 28   WI.2020 vs WI.2021  1  6743.746 1.1273143 0.2731351 0.2000000          1

# Visualize dispersions
png('figures/BetaDiversity/SSD_pcoa_betadispersion_site_by_year_summer_only.png',width = 700, height = 600, res=100)
plot(b.disper4,main = "Centroids and Dispersion based on Aitchison Distance", col=colorset6$Site_Color)
dev.off()

png('figures/BetaDiversity/SSD_boxplot_centroid_distance_site_by_year_summer_only.png',width = 700, height = 600, res=100)
boxplot(b.disper4,xlab="By Site x Collection Year (Summer Only)", main = "Distance to Centroid by Category", sub="Based on Aitchison Distance", col=colorset6$Site_Color)
dev.off()


### compare only sites in 2020

# pulling out 2020 samples only
b.clr[which(rownames(b.clr) %in% rownames(meta.2020)),][1:10,1:5] # does our indexing idea work? yes!
clr.2020<-b.clr[which(rownames(b.clr) %in% rownames(meta.2020)),]

b.disper5<-betadisper((vegdist(clr.2020,method="euclidean")), meta.2020$Site)
b.disper5
# NOTE: SB and RHB have less samples than other sites since they are supposed to represent similar locations
# maybe we should remove these sites and then rerun the PERMANOVA?

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:

permutest(b.disper5, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#         PD     BDC      DP    WI
# PD          0.19200 0.49700 0.149
# BDC 0.17683         0.65300 0.973
# DP  0.50433 0.66492         0.623
# WI  0.14123 0.97072 0.63281

anova(b.disper5) # p = 0.5764 --> accept the Null H, spatial medians (a measure of dispersion) are NOT significantly difference across sample dates
# ANOVA adjusted p-value
aov.beta.p5<-anova(b.disper5)[["Pr(>F)"]] # get p values from ANOVA
p.adjust(aov.beta.p5,method="bonferroni",n=length(aov.beta.p5))

TukeyHSD(b.disper5) # tells us which summer sites + years /category's dispersion MEANS are significantly different than each other
#            diff       lwr       upr     p adj
# BDC-PD -35.219152 -121.47415  51.03585 0.6311958
# DP-PD  -19.275391 -105.53039  66.97961 0.9087550
# WI-PD  -36.399334 -122.65433  49.85567 0.6074864
# DP-BDC  15.943761  -70.31124 102.19876 0.9451076
# WI-BDC  -1.180182  -87.43518  85.07482 0.9999745
# WI-DP  -17.123943 -103.37894  69.13106 0.9333724

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

pnova5<-adonis2(clr.2020 ~ Site,data=meta.2020,method = "euclidean",by="terms",permutations=1000)
pnova5
#           Df SumOfSqs      R2      F  Pr(>F)
# Site      3    34123 0.25048 1.3368 0.08292 .
# Residual 12   102106 0.74952
# Total    15   136229 1.00000

p.adjust(pnova5$`Pr(>F)`,method="bonferroni",n=length(pnova5$`Pr(>F)`)) # adjusted pval
# [1] 0.2487512        NA        NA

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

clr.2020.dist = (vegdist(clr.2020, "euclidean", na.rm = TRUE)) #distance matrix using Bray's dissimilarity index for trait distribution (traits of interest only)
pair.mod5<-pairwise.adonis(clr.2020.dist,meta.2020$Site, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
pair.mod5
#       pairs Df SumsOfSqs  F.Model        R2     p.value     p.adjusted sig
# 1 BDC vs DP  1  6893.020 0.8708630 0.1267473   0.745      1.000
# 2 BDC vs PD  1 14642.101 1.5596363 0.2063110   0.065      0.390
# 3 BDC vs WI  1  5893.646 0.9744444 0.1397164   0.440      1.000
# 4  DP vs PD  1 15559.884 1.4184650 0.1912073   0.173      1.000
# 5  DP vs WI  1  7371.709 0.9662006 0.1386984   0.440      1.000
# 6  PD vs WI  1 17885.412 1.9648737 0.2466924   0.113      0.678

# Visualize dispersions
png('figures/BetaDiversity/SSD_pcoa_betadispersion_site_2020_only.png',width = 700, height = 600, res=100)
plot(b.disper5,main = "Centroids and Dispersion based on Aitchison Distance (2020 Only)", col=colorset6$Site_Color)
dev.off()

png('figures/BetaDiversity/SSD_boxplot_centroid_distance_site_2020_only.png',width = 700, height = 600, res=100)
boxplot(b.disper5,xlab="By Site (2020 Only)", main = "Distance to Centroid by Category", sub="Based on Aitchison Distance", col=colorset6$Site_Color)
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
dust_meta$Depth.num<-as.numeric(as.character(dust_meta$Depth_m))

# now make sure your data frames you're comparing are in the same exact order!!
rownames(b.clr) %in% rownames(dust_meta)
dust_meta=dust_meta[rownames(b.clr),] ## reorder metadata to match order of CLR data
perm <- with(dust_meta, how(nperm = 1000)) # using SampDate as block because there is a significant difference between sample dates, trying to remove this effect when looking at permanovas

pnova1<-adonis2(b.clr ~ DO_Percent_Local*ORP_mV*Temp_DegC*Dissolved_OrganicMatter_RFU*Depth.num*Sulfate_milliM*Sulfide_microM,data=dust_meta,method = "euclidean",by="terms",permutations=perm)
pnova1
# nothing

adonis2(b.clr ~ DO_Percent_Local*ORP_mV*Temp_DegC*Dissolved_OrganicMatter_RFU*Depth.num*Sulfate_milliM*Sulfide_microM,data=dust_meta,method = "euclidean",by=NULL,permutations=perm)
#         Df SumOfSqs     R2    F Pr(>F)
#Model    23    25343  1
#Residual  0        0  0
#Total    23    25343  1

pnova2<-adonis2(b.clr ~ DO_Percent_Local*ORP_mV*Temp_DegC*Dissolved_OrganicMatter_RFU*Sulfate_milliM*Sulfide_microM,data=dust_meta,method = "euclidean",by="terms",permutations=perm)
pnova2
# nothing significant

adonis2(b.clr ~ DO_Percent_Local*ORP_mV*Temp_DegC*Dissolved_OrganicMatter_RFU*Sulfate_milliM*Sulfide_microM,data=dust_meta,method = "euclidean",by=NULL,permutations=perm)
#         Df SumOfSqs      R2      F   Pr(>F)
#Model    23    34412 0.73114 1.8918 0.4615
#Residual 16    12654 0.26886
#Total    39    47066 1.00000

pnova3<-adonis2(b.clr ~ DO_Percent_Local*ORP_mV*Temp_DegC*Dissolved_OrganicMatter_RFU*Sulfate_milliM*Sulfide_microM,data=dust_meta,method = "euclidean",by="terms",permutations=perm)
pnova3

adonis2(b.clr ~ ORP_mV*DO_Percent_Local*Dissolved_OrganicMatter_RFU*Sulfate_milliM*Sulfide_microM,data=dust_meta,method = "euclidean",by=NULL,permutations=perm)

pnova4<-adonis2(b.clr ~ ORP_mV*DO_Percent_Local*Dissolved_OrganicMatter_RFU*Sulfate_milliM,data=dust_meta,method = "euclidean",by="terms",permutations=perm)
pnova4
#                                         Df SumOfSqs      R2       F   Pr(>F)
# ORP_mV                                                              1   3435.1 0.13555  8.0320 0.000999 ***
# DO_Percent_Local                                                    1   4009.5 0.15821  9.3749 0.000999 ***
# Dissolved_OrganicMatter_RFU                                         1   5843.0 0.23056 13.6620 0.000999 ***
# Sulfate_milliM                                                      1    961.4 0.03794  2.2480 0.042957 *
# ORP_mV:DO_Percent_Local                                             1   1322.9 0.05220  3.0932 0.002997 **
# DO_Percent_Local:Dissolved_OrganicMatter_RFU                        1   1413.3 0.05577  3.3045 0.005994 **
# DO_Percent_Local:Sulfate_milliM                                     1    794.0 0.03133  1.8566 0.072927 .
p.adjust(pnova4$`Pr(>F)`,method="bonferroni",n=length(pnova4$`Pr(>F)`)) # adjusted pval

adonis2(b.clr ~ ORP_mV*DO_Percent_Local*Dissolved_OrganicMatter_RFU*Sulfate_milliM,data=dust_meta,method = "euclidean",by=NULL,permutations=perm)
#         Df SumOfSqs      R2      F   Pr(>F)
#Model    15  21921.1 0.86499 3.417 0.000999 ***
#Residual  8   3421.5 0.13501
#Total    23  25342.6 1.00000

pnova4b<-adonis2(b.clr ~ DO_Percent_Local*Dissolved_OrganicMatter_RFU*ORP_mV,data=dust_meta,method = "euclidean",by="terms",permutations=perm)
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

p.adjust(pnova4b$`Pr(>F)`,method="bonferroni",n=length(pnova4b$`Pr(>F)`)) # adjusted pval

pnova5<-adonis2(b.clr ~ DO_Percent_Local*Dissolved_OrganicMatter_RFU,data=dust_meta,method = "euclidean",by="terms",permutations=perm)
pnova5
#                                               Df SumOfSqs      R2       F   Pr(>F)
# DO_Percent_Local                              1   5830.3 0.23006  9.5297 0.000999 ***
#   Dissolved_OrganicMatter_RFU                   1   6344.6 0.25036 10.3704 0.000999 ***
#   DO_Percent_Local:Dissolved_OrganicMatter_RFU  1    931.5 0.03676  1.5226 0.134865
# Residual                                     20  12236.1 0.48283
# Total                                        23  25342.6 1.00000

p.adjust(pnova5$`Pr(>F)`,method="bonferroni",n=length(pnova5$`Pr(>F)`)) # adjusted pval

adonis2(b.clr ~ DO_Percent_Local*Dissolved_OrganicMatter_RFU,data=dust_meta,method = "euclidean",by=NULL,permutations=perm)
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
save.image("data/SSea Dust_BetaDiv_Data.Rdata")
