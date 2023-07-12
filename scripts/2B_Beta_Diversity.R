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
load("data/SSDust_16S.V3V4_W23_Data_Ready.Rdata") # save global env to Rdata file
#load("data/SSD_16S_CLR_EucDist_Ready.Rdata")

head(b.dust.all)
bac.ASV_table[1:4,1:4]
bac.ASV_table[(nrow(bac.ASV_table)-4):(nrow(bac.ASV_table)),(ncol(bac.ASV_table)-4):(ncol(bac.ASV_table))] # last 4 rows & cols
head(dust_meta)

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
raw.b.dend_cols <- as.character(dust_meta$SampMonth_Color[order.dendrogram(raw.b.euc_dend)])
labels_colors(raw.b.euc_dend) <- raw.b.dend_cols

## DO NOT RUN THIS LINE, THIS IS YOUR COLOR REFERENCE!!!!
#(August.2021="#ef781c",December.2021="#03045e",April.2022="#059c3f")

plot(raw.b.euc_dend, ylab="CLR Euclidean Distance",cex = 0.5) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 2)
#legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c("#ef781c","#03045e","#059c3f"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# PCOA w/ Euclidean distance matrix (of CLR data)
raw.b.pcoa <- pcoa(raw.b.euc_dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
#save.image("data/SSD_16S_CLR_EucDist_Ready.Rdata")

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
# PC1 = 48.24%, PC2 = 17.3%

# drop outliars
outliars<-c("BDC.D.7.27.21","PD.D.7.27.21","WI.D.9.18.21","WI.D.7.10.20")
raw.b.pcoa.meta2<-raw.b.pcoa.meta[!(raw.b.pcoa.meta$SampleID %in% outliars),]

#### Visualize RAW PCoAs####

# create PCoA ggplot fig
ggplot(raw.b.pcoa.meta2, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(SampleMonth),shape=CollectionYear), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Date",values=unique(raw.b.pcoa.meta$SampMonth_Color[order(raw.b.pcoa.meta$SampleMonth)]),labels=c(unique(raw.b.pcoa.meta$SampleMonth[order(raw.b.pcoa.meta$SampleMonth)]))) +
  xlab("PC1 [48.24%]") + ylab("PC2 [17.3%]")

ggplot(raw.b.pcoa.meta2, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site),shape=CollectionYear), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Date",values=unique(raw.b.pcoa.meta$Site_Color[order(raw.b.pcoa.meta$Site)]),labels=c(unique(raw.b.pcoa.meta$Site[order(raw.b.pcoa.meta$Site)]))) +
  xlab("PC1 [48.24%]") + ylab("PC2 [9.73%]")

# by collection period & site
ggplot(raw.b.pcoa.meta2, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year),shape=Site), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(raw.b.pcoa.meta$SCY_Color[order(raw.b.pcoa.meta$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
  xlab("PC1 [48.24%]") + ylab("PC2 [9.73%]")


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
b.euc_dend <- as.dendrogram(b.euc_clust, hang=0.2)
b.dend_cols <- as.character(dust_meta$SampMonth_Color[order.dendrogram(b.euc_dend)])
labels_colors(b.euc_dend) <- b.dend_cols

## DO NOT RUN THIS LINE, THIS IS YOUR COLOR REFERENCE!!!!
#(August.2021="#ef781c",December.2021="#03045e",April.2022="#059c3f")

plot(b.euc_dend, ylab="CLR Euclidean Distance",cex = 0.5) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 2)
#legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c("#ef781c","#03045e","#059c3f"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# PCOA w/ Euclidean distance matrix (of CLR data)
b.pcoa <- pcoa(b.euc_dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
#save.image("data/SSD_16S_CLR_EucDist_Ready.Rdata")

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
# PC1 = 23.63%, PC2 = 9.73%

save.image("data/SSD_16S_CLR_EucDist_Ready.Rdata")

#### Visualize PCoAs - All Data ####

## *** all figures that end in _PCOA1 come from the same single PCoA
#data is just subsetted to understand how points are related to each other w/in & across timepoints

# create PCoA ggplot fig
pcoa1<-ggplot(b.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(SampleMonth),shape=CollectionYear), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Date",values=unique(b.pcoa.meta$SampMonth_Color[order(b.pcoa.meta$SampleMonth)]),labels=c(unique(b.pcoa.meta$SampleMonth[order(b.pcoa.meta$SampleMonth)]))) +
  xlab("PC1 [23.63%]") + ylab("PC2 [9.73%]")

ggsave(pcoa1,filename = "figures/BetaDiversity/SSD_16S_CLR_SampMonth_Year_PCOA1.png", width=12, height=10, dpi=600)

pcoa1a<-ggplot(b.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site),shape=CollectionYear), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Date",values=unique(b.pcoa.meta$Site_Color[order(b.pcoa.meta$Site)]),labels=c(unique(b.pcoa.meta$Site[order(b.pcoa.meta$Site)]))) +
  xlab("PC1 [23.63%]") + ylab("PC2 [9.73%]")

ggsave(pcoa1a,filename = "figures/BetaDiversity/SSD_16S_CLR_Site_Year_PCOA1.png", width=12, height=10, dpi=600)

# specific season
pcoa3<-ggplot(b.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Season_Specific),shape=CollectionYear), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Season")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Season",values=unique(b.pcoa.meta$SeasonSpec_Color[order(b.pcoa.meta$Season_Specific)]),
                     labels=c("Early.Summer"="Early Summer","Late.Summer"="Late Summer","Early.Fall"="Early Fall","Late.Fall"="Late Fall","Fall.Winter"="Fall-Winter")) +
  xlab("PC1 [23.63%]") + ylab("PC2 [9.73%]")

ggsave(pcoa3,filename = "figures/BetaDiversity/SSD_16S_CLR_season_specific_PCOA1.png", width=12, height=10, dpi=600)

# by collection year & site
pcoa4a<-ggplot(b.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(CollectionYear),shape=Site), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Collection Year")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Year",values=unique(b.pcoa.meta$Year_Color[order(b.pcoa.meta$CollectionYear)])) +
  xlab("PC1 [23.63%]") + ylab("PC2 [9.73%]")

ggsave(pcoa4a,filename = "figures/BetaDiversity/SSD_16S_CLR_CollectionYear_Site_PCOA1.png", width=12, height=10, dpi=600)

# by collection period & site
pcoa5<-ggplot(b.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year),shape=Site), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.pcoa.meta$SCY_Color[order(b.pcoa.meta$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
  xlab("PC1 [23.63%]") + ylab("PC2 [9.73%]")

ggsave(pcoa5,filename = "figures/BetaDiversity/SSD_16S_CLR_SeasCollYr_Site_PCOA1.png", width=12, height=10, dpi=600)

# 3D PCoA

plot_ly(b.pcoa.meta, x=~Axis.1,y=~Axis.2,z=~Axis.3, color = ~Season_Specific, colors = c(unique(b.pcoa.meta$SeasonSpec_Color[order(b.pcoa.meta$Season_Specific)])),
        symbol=~CollectionYear,symbols = c("square", "circle")) %>%
  layout(scene = list(xaxis = list(title = 'PC1 23.63%'),
                      yaxis = list(title = 'PC2 9.73%'),
                      zaxis = list(title = 'PC3 6.97%')))

pltly.all.a<-plot_ly(b.pcoa.meta, x=~Axis.1,y=~Axis.2,z=~Axis.3, color = ~factor(Seas_Coll_Year), colors = c(unique(b.pcoa.meta$SCY_Color[order(b.pcoa.meta$Seas_Coll_Year)])),
                    symbol=~Site,symbols = c("square", "circle","diamond-open","x","circle-open")) %>%
  layout(scene = list(xaxis = list(title = 'PC1 23.63%'),
                      yaxis = list(title = 'PC2 9.73%'),
                      zaxis = list(title = 'PC3 6.97%')))
saveWidget(widget = pltly.all.a, #the plotly object,
           file = "figures/BetaDiversity/SSD_16S_CLR_SeasCollYr_Site_3D_PCOA1.html", #the path & file name
           selfcontained = TRUE #creates a single html file
)

## visualizing only specific sites

# WI - by collection period
wi.pc1<-ggplot(b.pcoa.meta[b.pcoa.meta$Site=="WI",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year)), size=5)+theme_bw()+
  labs(title="PCoA: Wister Microbial Composition",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.pcoa.meta$SCY_Color[order(b.pcoa.meta$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
  xlab("PC1 [23.63%]") + ylab("PC2 [9.73%]")

ggsave(wi.pc1,filename = "figures/BetaDiversity/SSD_16S_CLR_SeasCollYr_Wister_PCOA1.png", width=12, height=10, dpi=600)

dp.pc1<-ggplot(b.pcoa.meta[b.pcoa.meta$Site=="DP",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year)), size=5)+theme_bw()+
  labs(title="PCoA: Dos Palmas Microbial Composition",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.pcoa.meta$SCY_Color[order(b.pcoa.meta$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
  xlab("PC1 [23.63%]") + ylab("PC2 [9.73%]")

ggsave(dp.pc1,filename = "figures/BetaDiversity/SSD_16S_CLR_SeasCollYr_DosPalmas_PCOA1.png", width=12, height=10, dpi=600)

s.r.pc1<-ggplot(b.pcoa.meta[b.pcoa.meta$Site=="SB" | b.pcoa.meta$Site=="RHB",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year),shape=Site), size=5)+theme_bw()+
  labs(title="PCoA: Sonny Bono & Red Hill Bay Microbial Composition",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.pcoa.meta$SCY_Color[order(b.pcoa.meta$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
  xlab("PC1 [23.63%]") + ylab("PC2 [9.73%]")

ggsave(s.r.pc1,filename = "figures/BetaDiversity/SSD_16S_CLR_SeasCollYr_SB_and_RHB_PCOA1.png", width=12, height=10, dpi=600)

pd.pc1<-ggplot(b.pcoa.meta[b.pcoa.meta$Site=="PD",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year)), size=5)+theme_bw()+
  labs(title="PCoA: Palm Desert Microbial Composition",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.pcoa.meta$SCY_Color[order(b.pcoa.meta$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
  xlab("PC1 [23.63%]") + ylab("PC2 [9.73%]")

ggsave(pd.pc1,filename = "figures/BetaDiversity/SSD_16S_CLR_SeasCollYr_PD_PCOA1.png", width=12, height=10, dpi=600)

bdc.pc1<-ggplot(b.pcoa.meta[b.pcoa.meta$Site=="BDC",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year)), size=5)+theme_bw()+
  labs(title="PCoA: Boyd Deep Canyon Microbial Composition",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.pcoa.meta$SCY_Color[order(b.pcoa.meta$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
  xlab("PC1 [23.63%]") + ylab("PC2 [9.73%]")

ggsave(bdc.pc1,filename = "figures/BetaDiversity/SSD_16S_CLR_SeasCollYr_BDC_PCOA1.png", width=12, height=10, dpi=600)

# visualizing pcoa by month per year
unique(b.pcoa.meta$SampDate)

j2020.p1<-ggplot(b.pcoa.meta[b.pcoa.meta$SampDate=="July.2020",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: July 2020 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.pcoa.meta$Site_Color[order(b.pcoa.meta$Site)])) +
  xlab("PC1 [23.63%]") + ylab("PC2 [9.73%]")

ggsave(j2020.p1,filename = "figures/BetaDiversity/SSD_16S_CLR_AllSites_July2020_PCOA1.png", width=12, height=10, dpi=600)

a2020.p1<-ggplot(b.pcoa.meta[b.pcoa.meta$SampDate=="August.2020",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: August 2020 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.pcoa.meta$Site_Color[order(b.pcoa.meta$Site)])) +
  xlab("PC1 [23.63%]") + ylab("PC2 [9.73%]")

ggsave(a2020.p1,filename = "figures/BetaDiversity/SSD_16S_CLR_AllSites_Aug2020_PCOA1.png", width=12, height=10, dpi=600)

o2020.p1<-ggplot(b.pcoa.meta[b.pcoa.meta$SampDate=="October.2020",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: October 2020 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.pcoa.meta$Site_Color[order(b.pcoa.meta$Site)])) +
  xlab("PC1 [23.63%]") + ylab("PC2 [9.73%]")

ggsave(o2020.p1,filename = "figures/BetaDiversity/SSD_16S_CLR_AllSites_Oct2020_PCOA1.png", width=12, height=10, dpi=600)

n2020.p1<-ggplot(b.pcoa.meta[b.pcoa.meta$SampDate=="November.2020",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: November 2020 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.pcoa.meta$Site_Color[order(b.pcoa.meta$Site)])) +
  xlab("PC1 [23.63%]") + ylab("PC2 [9.73%]")

ggsave(n2020.p1,filename = "figures/BetaDiversity/SSD_16S_CLR_AllSites_Nov2020_PCOA1.png", width=12, height=10, dpi=600)

ja2021.p1<-ggplot(b.pcoa.meta[b.pcoa.meta$SampDate=="July.2021" | b.pcoa.meta$SampDate=="August.2021",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: July & August 2021 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.pcoa.meta$Site_Color[order(b.pcoa.meta$Site)])) +
  xlab("PC1 [23.63%]") + ylab("PC2 [9.73%]")

ggsave(ja2021.p1,filename = "figures/BetaDiversity/SSD_16S_CLR_AllSites_July_Aug_2021_PCOA1.png", width=12, height=10, dpi=600)

s2021.p1<-ggplot(b.pcoa.meta[b.pcoa.meta$SampDate=="September.2021",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: September 2021 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.pcoa.meta$Site_Color[order(b.pcoa.meta$Site)])) +
  xlab("PC1 [23.63%]") + ylab("PC2 [9.73%]")

ggsave(s2021.p1,filename = "figures/BetaDiversity/SSD_16S_CLR_AllSites_Sept2021_PCOA1.png", width=12, height=10, dpi=600)

d2021.p1<-ggplot(b.pcoa.meta[b.pcoa.meta$SampDate=="December.2021",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: December 2021 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.pcoa.meta$Site_Color[order(b.pcoa.meta$Site)])) +
  xlab("PC1 [23.63%]") + ylab("PC2 [9.73%]")

ggsave(d2021.p1,filename = "figures/BetaDiversity/SSD_16S_CLR_AllSites_Dec2021_PCOA1.png", width=12, height=10, dpi=600)

#### CLR Transform 2020 Comp Data ####
d.meta_20<-dust_meta[which(dust_meta$CollectionYear=="2020"),]

asv.table.20<-bac.ASV_table[which(bac.ASV_table$SampleID %in% d.meta_20$SampleID),]
asv.table.20[1:4,1:4]
rownames(asv.table.20)

# CLR transformation of ASV table
# df must have rownames are SampleIDs, columns are ASV IDs for vegan functions below
b.clr.20<-decostand(asv.table.20[,-1],method = "clr", pseudocount = 1) #CLR transformation
b.clr.20[1:4,1:4]

#### Beta Diversity - 2020 Data ####

# check rownames of CLR transformed ASV data & metadata
rownames(b.clr.20) %in% rownames(d.meta_20)
d.meta_20=d.meta_20[rownames(b.clr.20),] ## reorder metadata to match order of CLR data

# calculate our Euclidean distance matrix using CLR data
b.20.euc_dist <- dist(b.clr.20, method = "euclidean")

# creating our hierarcical clustering dendrogram
b.20.euc_clust <- hclust(b.20.euc_dist, method="ward.D2")

# let's make it a little nicer...
b.20.euc_dend <- as.dendrogram(b.20.euc_clust, hang=0.2)
b.20.dend_cols <- as.character(d.meta_20$SampMonth_Color[order.dendrogram(b.20.euc_dend)])
labels_colors(b.20.euc_dend) <- b.20.dend_cols

## DO NOT RUN THIS LINE, THIS IS YOUR COLOR REFERENCE!!!!
#(August.2021="#ef781c",December.2021="#03045e",April.2022="#059c3f")

plot(b.20.euc_dend, ylab="CLR Euclidean Distance",cex = 0.5) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 2)
#legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c("#ef781c","#03045e","#059c3f"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# PCOA w/ Euclidean distance matrix (of CLR data)
b.20.pcoa <- pcoa(b.20.euc_dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
#save.image("data/SSD_16S_CLR_EucDist_Ready.Rdata")

# The proportion of variances explained is in its element values$Relative_eig
b.20.pcoa$values

# extract principal coordinates
b.20.pcoa.vectors<-data.frame(b.20.pcoa$vectors)
b.20.pcoa.vectors$SampleID<-rownames(b.20.pcoa$vectors)

# merge pcoa coordinates w/ metadata
b.20.pcoa.meta<-merge(b.20.pcoa.vectors, d.meta_20, by.x="SampleID", by.y="SampleID")
b.20.pcoa.meta$SampleMonth
b.20.pcoa.meta$SampDate

head(b.20.pcoa.meta)

head(b.20.pcoa$values) # pull out Relative (Relative_eig) variation % to add to axes labels
# PC1 = 23.63%, PC2 = 9.73%

save.image("data/SSD_16S_2020_CLR_EucDist_Ready.Rdata")

#### Visualize PCoAs - 2020 Data ####

# create PCoA ggplot fig
pcoa1<-ggplot(b.20.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Season_Specific),shape=Site), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust, 2020",subtitle="Using Centered-Log Ratio Data",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Date",values=unique(b.20.pcoa.meta$SeasonSpec_Color[order(b.20.pcoa.meta$Season_Specific)]),labels=c(unique(b.20.pcoa.meta$Season_Specific[order(b.20.pcoa.meta$Season_Specific)]))) +
  xlab("PC1 [31.85%]") + ylab("PC2 [15.26%]")

ggsave(pcoa1,filename = "figures/BetaDiversity/2020/SSD_16S_2020_CLR_SampMonth_Site_PCOA1.png", width=12, height=10, dpi=600)

# collection period
pcoa1a<-ggplot(b.20.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=Seas_Coll_Year,shape=Site), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust, 2020",subtitle="Using Centered-Log Ratio Data",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.20.pcoa.meta$SCY_Color[order(b.20.pcoa.meta$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1","S.2.2020"="Summer #2","S.3.2020"="Summer #3","F.1.2020"="Fall #1")) +
  xlab("PC1 [31.85%]") + ylab("PC2 [15.26%]")

ggsave(pcoa1a,filename = "figures/BetaDiversity/2020/SSD_16S_2020_CLR_Site_SeasCollYr_PCOA1.png", width=12, height=10, dpi=600)

# specific season
pcoa3<-ggplot(b.20.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Season_Specific),shape=Site), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust, 2020",subtitle="Using Centered-Log Ratio Data",color="Season")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Season",values=unique(b.20.pcoa.meta$SeasonSpec_Color[order(b.20.pcoa.meta$Season_Specific)]),
                     labels=c("Early.Summer"="Early Summer","Late.Summer"="Late Summer","Early.Fall"="Early Fall","Late.Fall"="Late Fall","Fall.Winter"="Fall-Winter")) +
  xlab("PC1 [31.85%]") + ylab("PC2 [15.26%]")

ggsave(pcoa3,filename = "figures/BetaDiversity/2020/SSD_16S_2020_CLR_season_specific_PCOA1.png", width=12, height=10, dpi=600)

# 3D PCoA
pltly.20.a<-plot_ly(b.20.pcoa.meta, x=~Axis.1,y=~Axis.2,z=~Axis.3, color = ~Season_Specific, colors = c(unique(b.20.pcoa.meta$SeasonSpec_Color[order(b.20.pcoa.meta$Season_Specific)])),
                    symbol=~Site,symbols = c("square", "circle","diamond-open","x","circle-open")) %>%
  layout(scene = list(xaxis = list(title = 'PC1 31.85%'),
                      yaxis = list(title = 'PC2 15.26%'),
                      zaxis = list(title = 'PC3 9.33%')))
saveWidget(widget = pltly.20.a, #the plotly object,
           file = "figures/BetaDiversity/2020/SSD_16S_2020_CLR_SeasonSpec_Site_3D_PCOA.html", #the path & file name
           selfcontained = TRUE #creates a single html file
)
pltly.20.b<-plot_ly(b.20.pcoa.meta, x=~Axis.1,y=~Axis.2,z=~Axis.3, color = ~factor(Seas_Coll_Year), colors = c(unique(b.20.pcoa.meta$SCY_Color[order(b.20.pcoa.meta$Seas_Coll_Year)])),
                    symbol=~Site,symbols = c("square", "circle","diamond-open","x","circle-open")) %>%
  layout(scene = list(xaxis = list(title = 'PC1 31.85%'),
                      yaxis = list(title = 'PC2 15.26%'),
                      zaxis = list(title = 'PC3 9.33%')))
saveWidget(widget = pltly.20.b, #the plotly object,
           file = "figures/BetaDiversity/2020/SSD_16S_2020_CLR_SeasCollYr_Site_3D_PCOA.html", #the path & file name
           selfcontained = TRUE #creates a single html file
)
#### CLR Transform 2021 Comp Data ####
d.meta_21<-dust_meta[which(dust_meta$CollectionYear=="2021"),]

asv.table.21<-bac.ASV_table[which(bac.ASV_table$SampleID %in% d.meta_21$SampleID),]
asv.table.21[1:4,1:4]
rownames(asv.table.21)

# CLR transformation of ASV table
# df must have rownames are SampleIDs, columns are ASV IDs for vegan functions below
b.clr.21<-decostand(asv.table.21[,-1],method = "clr", pseudocount = 1) #CLR transformation
b.clr.21[1:4,1:4]

#### Beta Diversity - 2021 Data ####

# check rownames of CLR transformed ASV data & metadata
rownames(b.clr.21) %in% rownames(d.meta_21)
d.meta_21=d.meta_21[rownames(b.clr.21),] ## reorder metadata to match order of CLR data

# calculate our Euclidean distance matrix using CLR data
b.21.euc_dist <- dist(b.clr.21, method = "euclidean")

# creating our hierarcical clustering dendrogram
b.21.euc_clust <- hclust(b.21.euc_dist, method="ward.D2")

# let's make it a little nicer...
b.21.euc_dend <- as.dendrogram(b.21.euc_clust, hang=0.2)
b.21.dend_cols <- as.character(d.meta_21$SampMonth_Color[order.dendrogram(b.21.euc_dend)])
labels_colors(b.21.euc_dend) <- b.21.dend_cols

## DO NOT RUN THIS LINE, THIS IS YOUR COLOR REFERENCE!!!!
#(August.2121="#ef781c",December.2121="#03045e",April.2122="#059c3f")

plot(b.21.euc_dend, ylab="CLR Euclidean Distance",cex = 0.5) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 2)
#legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c("#ef781c","#03045e","#059c3f"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# PCOA w/ Euclidean distance matrix (of CLR data)
b.21.pcoa <- pcoa(b.21.euc_dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
#save.image("data/SSD_16S_CLR_EucDist_Ready.Rdata")

# The proportion of variances explained is in its element values$Relative_eig
b.21.pcoa$values

# extract principal coordinates
b.21.pcoa.vectors<-data.frame(b.21.pcoa$vectors)
b.21.pcoa.vectors$SampleID<-rownames(b.21.pcoa$vectors)

# merge pcoa coordinates w/ metadata
b.21.pcoa.meta<-merge(b.21.pcoa.vectors, d.meta_21, by.x="SampleID", by.y="SampleID")
b.21.pcoa.meta$SampleMonth
b.21.pcoa.meta$SampDate

head(b.21.pcoa.meta)

head(b.21.pcoa$values) # pull out Relative (Relative_eig) variation % to add to axes labels
# PC1 = 24.72%, PC2 = 15.13%

save.image("data/SSD_16S_2021_CLR_EucDist_Ready.Rdata")

#### Visualize PCoAs - 2021 Data ####

# create PCoA ggplot fig
pcoa1<-ggplot(b.21.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Season_Specific),shape=Site), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust, 2021",subtitle="Using Centered-Log Ratio Data",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Date",values=unique(b.21.pcoa.meta$SeasonSpec_Color[order(b.21.pcoa.meta$Season_Specific)]),labels=c(unique(b.21.pcoa.meta$Season_Specific[order(b.21.pcoa.meta$Season_Specific)]))) +
  xlab("PC1 [24.72%]") + ylab("PC2 [15.13%]")

ggsave(pcoa1,filename = "figures/BetaDiversity/2021/SSD_16S_2021_CLR_SampMonth_Site_PCOA1.png", width=12, height=10, dpi=600)

pcoa1a<-ggplot(b.21.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=Seas_Coll_Year,shape=Site), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust, 2021",subtitle="Using Centered-Log Ratio Data",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.21.pcoa.meta$SCY_Color[order(b.21.pcoa.meta$Seas_Coll_Year)]),labels=c("S.1.2021"="Summer #1","S.2.2021"="Summer #2","S.3.2021"="Summer #3","F.1.2021"="Fall #1")) +
  xlab("PC1 [24.72%]") + ylab("PC2 [15.13%]")

ggsave(pcoa1a,filename = "figures/BetaDiversity/2021/SSD_16S_2021_CLR_Site_SeasCollYr_PCOA1.png", width=12, height=10, dpi=600)

# specific season
pcoa3<-ggplot(b.21.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Season_Specific),shape=Site), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust, 2021",subtitle="Using Centered-Log Ratio Data",color="Season")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Season",values=unique(b.21.pcoa.meta$SeasonSpec_Color[order(b.21.pcoa.meta$Season_Specific)]),
                     labels=c("Early.Summer"="Early Summer","Late.Summer"="Late Summer","Early.Fall"="Early Fall","Late.Fall"="Late Fall","Fall.Winter"="Fall-Winter")) +
  xlab("PC1 [24.72%]") + ylab("PC2 [15.13%]")

ggsave(pcoa3,filename = "figures/BetaDiversity/2021/SSD_16S_2021_CLR_season_specific_PCOA1.png", width=12, height=10, dpi=600)

# 3D PCoAs
pltly.21.a<-plot_ly(b.21.pcoa.meta, x=~Axis.1,y=~Axis.2,z=~Axis.3, color = ~factor(Season_Specific), colors = c(unique(b.21.pcoa.meta$SeasonSpec_Color[order(b.21.pcoa.meta$Season_Specific)])),
                    symbol=~Site,symbols = c("square", "circle","diamond-open","x","circle-open")) %>%
  layout(scene = list(xaxis = list(title = 'PC1 24.72%'),
                      yaxis = list(title = 'PC2 15.13%'),
                      zaxis = list(title = 'PC3 12.79%')))
saveWidget(widget = pltly.21.a, #the plotly object,
           file = "figures/BetaDiversity/2021/SSD_16S_2021_CLR_SeasonSpec_Site_3D_PCOA.html", #the path & file name
           selfcontained = TRUE #creates a single html file
)
pltly.21.b<-plot_ly(b.21.pcoa.meta, x=~Axis.1,y=~Axis.2,z=~Axis.3, color = ~factor(Seas_Coll_Year), colors = c(unique(b.21.pcoa.meta$SCY_Color[order(b.21.pcoa.meta$Seas_Coll_Year)])),
                    symbol=~Site,symbols = c("square", "circle","diamond-open","x","circle-open")) %>%
  layout(scene = list(xaxis = list(title = 'PC1 24.72%'),
                      yaxis = list(title = 'PC2 15.13%'),
                      zaxis = list(title = 'PC3 12.79%')))
saveWidget(widget = pltly.21.b, #the plotly object,
           file = "figures/BetaDiversity/2021/SSD_16S_2021_CLR_SeasCollYr_Site_3D_PCOA.html", #the path & file name
           selfcontained = TRUE #creates a single html file
)


#### Separate All Data by Timepoints ####

head(dust_meta)

date_list<-unique(dust_meta$SampDate) #define an array of string values
# go through metadata & create a list of data frames
## when metadata$Variable == element in date_list (aka x in this case), subset metadata by said element into elements of a list
date_list
# here the function(x) is using date_list aka x to subset metadata, when $Variable column == date_list
# Run the function so it's stored in Global Env
date_subsets<-lapply(date_list, function(x) {subset(dust_meta, SampDate==x)})
# date_subsets is a list; each element is the sample date

date_subsets # sanity check1 (should see all elements in list)
date_subsets[[1]] # sanity check2 (see 1st element in list)
date_subsets[[2]]
date_subsets[[3]]
#rename the list elements

# name each element in list
names(date_subsets)<-date_list # * only do this if the order of names in date_list match order of the elements in date_subsets!
date_subsets$July.2021 # sanity check3 - should be able to pull dataframes by names rather than index now

# example of subsetting
date_subsets[[2]][1:6]
date_subsets$November.2020 # should produce same ouptut as line above

date_subsets[[2]][1:2,1:2] # another example of indexing w/ list

# ^ subsetting list to [[second dataframe]], [[row #, column #]]
date_subsets[[2]][[1,2]] # [[second dataframe]], [[row 1, column 2]]

# set up the function and run this to store fxn in our Global environment
df_specific.subset<-function(var_vec,var_subsets){
  # var_vec = vector of variable elements from specific categorical variable;
  ## e.g. vector of names from date categorical variable (metadata dates)
  # var_subsets = list of dataframes subsetted by column$element from original dataframe;
  ## e.g. list of dataframes (each df = element of list) subsetted from metadata using vector of metadata$date names
  for(i in seq_along(var_vec)){
    # print(var_vec[i]) -- var_vec[i] = each element in var_vec
    # print(var_subsets[[i]]) -- var_subsets[[i]] = each sub
    df<-paste(var_vec[i])
    #print(df)
    assign(df, var_subsets[[i]], envir = .GlobalEnv)
    print(paste("Dataframe", var_vec[i] ,"done"))

  }

}

# run the function
df_specific.subset(date_list, date_subsets) # used scaled metadata quantitative values

head(August.2021) # sanity check
August.2021[1:5,] # double check that our new Variable (here SampDate) data frames still have scaled chemical data
rownames(August.2021)

# matching data with user defined function -- here is the function, must run to store function in Global env
match_dat<-function(compdata, subset_metadata){
  subset_comp_data = pullrow<-(is.element(row.names(compdata), row.names(subset_metadata)))
  ### * comp data and metadata need to have (same) row names - rownames should be Sample IDs
  subset_comp_data=compdata[pullrow,]
  return(subset_comp_data)
}

# double check that our data frames are ready for this function, aka that they both have the same rownames
## row #s do not have to be the same, but their row names should be in the same format and be able to match up
rownames(b.clr)
rownames(August.2021)

# run the match_dat function in loop
# using elements in date_list to look for object w/ same name in global env, then run match_dat function on that object to create new object
for (x in date_list) {
  #print(x)
  if(x %in% names(which(unlist(eapply(.GlobalEnv,is.data.frame))))){
    print(x)
    df<-paste0("b.clr_",x)
    df2<-match_dat(b.clr,get(x))
    assign(df, df2, envir = .GlobalEnv)
  }
}

# did the function work the way we wanted it to?
b.clr_July.2020[,1:4]
b.clr_July.2021[,1:4]

# function version of the match dat + loop in nested function
loop_match_dat<-function(obj_list, compdata){
  # must run match_dat function first so it's in Global Env
  match_dat<-function(compdata, subset_metadata){
    subset_comp_data = pullrow<-(is.element(row.names(compdata), row.names(subset_metadata)))
    ### * comp data and metadata need to have (same) row names - rownames should be Sample IDs
    subset_comp_data=compdata[pullrow,]
    return(subset_comp_data)
  }
  for (x in obj_list) {
    #print(x)
    if(x %in% names(which(unlist(eapply(.GlobalEnv,is.data.frame))))){
      print(x)
      df<-paste0("b.clr_",x) # create name of df
      df2<-match_dat(compdata,get(x)) # create df2 w/ custom match_dat function
      assign(df, df2, envir = .GlobalEnv) # assign new df the values of df2 from match_dat function
    }
  }
  return(df)
}

#### Separate All Data by Collections ####

head(dust_meta)

coll_list<-unique(dust_meta$Seas_Coll_Year) #define an array of string values
# go through metadata & create a list of data frames
## when metadata$Variable == element in coll_list (aka x in this case), subset metadata by said element into elements of a list
coll_list
# here the function(x) is using coll_list aka x to subset metadata, when $Variable column == coll_list
# Run the function so it's stored in Global Env
coll_subsets<-lapply(coll_list, function(x) {subset(dust_meta, Seas_Coll_Year==x)})
# coll_subsets is a list; each element is the sample date

coll_subsets # sanity check1 (should see all elements in list)
coll_subsets[[1]] # sanity check2 (see 1st element in list)
coll_subsets[[2]]
coll_subsets[[3]]
#rename the list elements

# name each element in list
names(coll_subsets)<-coll_list # * only do this if the order of names in coll_list match order of the elements in coll_subsets!
coll_subsets$July.2021 # sanity check3 - should be able to pull dataframes by names rather than index now

# example of subsetting
coll_subsets[[2]][1:6]
coll_subsets$November.2020 # should produce same ouptut as line above

coll_subsets[[2]][1:2,1:2] # another example of indexing w/ list

# ^ subsetting list to [[second dataframe]], [[row #, column #]]
coll_subsets[[2]][[1,2]] # [[second dataframe]], [[row 1, column 2]]

# set up the function and run this to store fxn in our Global environment
# df_specific.subset<-function(var_vec,var_subsets){
#   # var_vec = vector of variable elements from specific categorical variable;
#   ## e.g. vector of names from coll categorical variable (metadata colls)
#   # var_subsets = list of dataframes subsetted by column$element from original dataframe;
#   ## e.g. list of dataframes (each df = element of list) subsetted from metadata using vector of metadata$coll names
#   for(i in seq_along(var_vec)){
#     # print(var_vec[i]) -- var_vec[i] = each element in var_vec
#     # print(var_subsets[[i]]) -- var_subsets[[i]] = each sub
#     df<-paste(var_vec[i])
#     #print(df)
#     assign(df, var_subsets[[i]], envir = .GlobalEnv)
#     print(paste("Dataframe", var_vec[i] ,"done"))
#
#   }
#
# }

# run the function
df_specific.subset(coll_list, coll_subsets) # used scaled metadata quantitative values

head(S.1.2020) # sanity check
S.1.2020[1:5,] # double check that our new Variable (here SampDate) data frames still have scaled chemical data
rownames(S.1.2020)

# matching data with user defined function -- here is the function, must run to store function in Global env
# match_dat<-function(compdata, subset_metadata){
#   subset_comp_data = pullrow<-(is.element(row.names(compdata), row.names(subset_metadata)))
#   ### * comp data and metadata need to have (same) row names - rownames should be Sample IDs
#   subset_comp_data=compdata[pullrow,]
#   return(subset_comp_data)
# }

# double check that our data frames are ready for this function, aka that they both have the same rownames
## row #s do not have to be the same, but their row names should be in the same format and be able to match up
rownames(b.clr)
rownames(S.1.2020)

# run the match_dat function in loop
# using elements in coll_list to look for object w/ same name in global env, then run match_dat function on that object to create new object
for (x in coll_list) {
  #print(x)
  if(x %in% names(which(unlist(eapply(.GlobalEnv,is.data.frame))))){
    print(x)
    df<-paste0("b.clr_",x)
    df2<-match_dat(b.clr,get(x))
    assign(df, df2, envir = .GlobalEnv)
  }
}

# did the function work the way we wanted it to?
b.clr_S.1.2020[,1:4]
b.clr_S.1.2021[,1:4]

#### Beta Diversity - Summer Early 2020 ####

# check rownames of CLR transformed ASV data & metadata
rownames(b.clr_S.1.2020) %in% rownames(dust_meta)
meta_s120=dust_meta[rownames(b.clr_S.1.2020),] ## reorder metadata to match order of CLR data

# calculate our Euclidean distance matrix using CLR data
b.s120.euc_dist <- dist(b.clr_S.1.2020, method = "euclidean")

# creating our hierarcical clustering dendrogram
b.s120.euc_clust <- hclust(b.s120.euc_dist, method="ward.D2")

# let's make it a little nicer...
b.s120.euc_dend <- as.dendrogram(b.s120.euc_clust, hang=0.2)
b.s120.dend_cols <- as.character(meta_s120$Site_Color[order.dendrogram(b.s120.euc_dend)])
labels_colors(b.s120.euc_dend) <- b.s120.dend_cols

## DO NOT RUN THIS LINE, THIS IS YOUR COLOR REFERENCE!!!!
#(August.2021="#ef781c",December.2021="#03045e",April.2022="#059c3f")

plot(b.s120.euc_dend, ylab="CLR Euclidean Distance",cex = 0.5) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 2)
#legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c("#ef781c","#03045e","#059c3f"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# PCOA w/ Euclidean distance matrix (of CLR data)
b.s120.pcoa <- pcoa(b.s120.euc_dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
#save.image("data/SSD_16S_CLR_EucDist_Ready.Rdata")

# The proportion of variances explained is in its element values$Relative_eig
b.s120.pcoa$values

# extract principal coordinates
b.s120.pcoa.vectors<-data.frame(b.s120.pcoa$vectors)
b.s120.pcoa.vectors$SampleID<-rownames(b.s120.pcoa$vectors)

# merge pcoa coordinates w/ metadata
b.s120.pcoa.meta<-merge(b.s120.pcoa.vectors, dust_meta, by.x="SampleID", by.y="SampleID")
b.s120.pcoa.meta$Site

head(b.s120.pcoa.meta)

head(b.s120.pcoa$values) # pull out Relative (Relative_eig) variation % to add to axes labels
# PC1 = 62.22%, PC2 = 25.75%

#### Visualize PCoAs - S.1.2020 Data ####
# create PCoA ggplot fig

pcoa1a<-ggplot(b.s120.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Date",values=unique(b.s120.pcoa.meta$Site_Color[order(b.s120.pcoa.meta$Site)]),labels=c(unique(b.s120.pcoa.meta$Site[order(b.s120.pcoa.meta$Site)]))) +
  xlab("PC1 [23.63%]") + ylab("PC2 [9.73%]")

ggsave(pcoa1a,filename = "figures/BetaDiversity/SSD_16S_CLR_Site_Year_PCOA1.png", width=12, height=10, dpi=600)


#### Beta Diversity - Exclude Control Sites ####
rownames(bac.ASV_table)
bac.ASV_table[1:4,1:4]

# remove control sites from ASV table
bac.ASV_table1<-bac.ASV_table[!grepl(c("BDC"), bac.ASV_table$SampleID),] #remove BDC
bac.ASV_tab.noctrl<-bac.ASV_table1[!grepl(c("PD"), bac.ASV_table1$SampleID),] #remove PD
unique(bac.ASV_tab.noctrl$SampleID) # sanity check

# CLR transformation of ASV table
# df must have rownames are SampleIDs, columns are ASV IDs for vegan functions below
b.near.clr<-decostand(bac.ASV_tab.noctrl[,-1],method = "clr", pseudocount = 1) #CLR transformation
b.near.clr[1:4,1:4]
dim(b.near.clr)

# check rownames of CLR transformed ASV data & metadata
rownames(b.near.clr) %in% rownames(dust_meta)
dust_meta_near=dust_meta[rownames(b.near.clr),] ## reorder metadata to match order of CLR data
dim(dust_meta_near)

# calculate our Euclidean distance matrix using CLR data
b.near.euc_dist <- dist(b.near.clr, method = "euclidean")

# creating our hierarcical clustering dendrogram
b.near.euc_clust <- hclust(b.near.euc_dist, method="ward.D2")

# let's make it a little nicer...
b.near.euc_dend <- as.dendrogram(b.near.euc_clust, hang=0.2)
b.near.dend_cols <- as.character(dust_meta_near$SeasonSpec_Color[order.dendrogram(b.near.euc_dend)])
labels_colors(b.near.euc_dend) <- b.near.dend_cols

plot(b.near.euc_dend, ylab="CLR Euclidean Distance",cex = 0.5) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 2)
dev.off()

# PCOA w/ Euclidean distance matrix (of CLR data)
b.near.pcoa <- pcoa(b.near.euc_dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix

# The proportion of variances explained is in its element values$Relative_eig
b.near.pcoa$values

# extract principal coordinates
b.near.pcoa.vectors<-data.frame(b.near.pcoa$vectors)
b.near.pcoa.vectors$SampleID<-rownames(b.near.pcoa$vectors)

# merge pcoa coordinates w/ metadata
b.near.pcoa.meta<-merge(b.near.pcoa.vectors, dust_meta_near, by.x="SampleID", by.y="SampleID")
b.near.pcoa.meta$SampleMonth
b.near.pcoa.meta$SampDate

head(b.near.pcoa.meta)

head(b.near.pcoa$values) # pull out Relative (Relative_eig) variation % to add to axes labels
# PC1 = 24.03%, PC2 = 14.26%

save.image("data/SSD_NearSitesOnly_16S_CLR_EucDist_Ready.Rdata")

#### Visualize PCoAs - Exclude Control Sites Data ####
plot_ly(b.near.pcoa.meta, x=~Axis.1,y=~Axis.2,z=~Axis.3, color = ~Season_Specific, colors = c(unique(b.near.pcoa.meta$SeasonSpec_Color[order(b.near.pcoa.meta$Season_Specific)])),
        symbol=~CollectionYear,symbols = c("square", "circle")) %>%
  layout(scene = list(xaxis = list(title = 'PC1 24.03%'),
                      yaxis = list(title = 'PC2 14.26%'),
                      zaxis = list(title = 'PC3 11.16%')))

# create PCoA ggplot fig
pcoa.n.1<-ggplot(b.near.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(SampleMonth),shape=CollectionYear), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Date",values=unique(b.near.pcoa.meta$SampMonth_Color[order(b.near.pcoa.meta$SampleMonth)]),labels=c(unique(b.near.pcoa.meta$SampleMonth[order(b.near.pcoa.meta$SampleMonth)]))) +
  xlab("PC1 [24.03%]") + ylab("PC2 [14.26%]")

ggsave(pcoa.n.1,filename = "figures/BetaDiversity/SSD_16S_NearSitesOnly_CLR_SampMonth_year_PCOA1.png", width=12, height=10, dpi=600)

pcoa.n.1a<-ggplot(b.near.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site),shape=CollectionYear), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Date",values=unique(b.near.pcoa.meta$Site_Color[order(b.near.pcoa.meta$Site)]),labels=c(unique(b.near.pcoa.meta$Site[order(b.near.pcoa.meta$Site)]))) +
  xlab("PC1 [24.03%]") + ylab("PC2 [14.26%]")

ggsave(pcoa.n.1a,filename = "figures/BetaDiversity/SSD_16S_NearSitesOnly_CLR_Site_Year_PCOA1.png", width=12, height=10, dpi=600)

# general season
pcoa.n.2<-ggplot(b.near.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Season_General),shape=CollectionYear), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Season")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Season",values=unique(b.near.pcoa.meta$SeasonGen_Color[order(b.near.pcoa.meta$Season_General)]),labels=c(unique(b.near.pcoa.meta$Season_General[order(b.near.pcoa.meta$Season_General)]))) +
  xlab("PC1 [24.03%]") + ylab("PC2 [14.26%]")

ggsave(pcoa.n.2,filename = "figures/BetaDiversity/SSD_16S_NearSitesOnly_CLR_season_general_PCOA1.png", width=12, height=10, dpi=600)

# specific season
pcoa.n.3<-ggplot(b.near.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Season_Specific),shape=CollectionYear), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Season")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Season",values=unique(b.near.pcoa.meta$SeasonSpec_Color[order(b.near.pcoa.meta$Season_Specific)]),
                     labels=c("Early.Summer"="Early Summer","Late.Summer"="Late Summer","Early.Fall"="Early Fall","Late.Fall"="Late Fall","Fall.Winter"="Fall-Winter")) +
  xlab("PC1 [24.03%]") + ylab("PC2 [14.26%]")

ggsave(pcoa.n.3,filename = "figures/BetaDiversity/SSD_16S_NearSitesOnly_CLR_season_specific_PCOA1.png", width=12, height=10, dpi=600)

# by collection year
pcoa.n.4<-ggplot(b.near.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(CollectionYear)), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Collection Year")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Year",values=unique(b.near.pcoa.meta$Year_Color[order(b.near.pcoa.meta$CollectionYear)])) +
  xlab("PC1 [24.03%]") + ylab("PC2 [14.26%]")

ggsave(pcoa.n.4,filename = "figures/BetaDiversity/SSD_16S_NearSitesOnly_CLR_CollectionYear_PCOA1.png", width=12, height=10, dpi=600)

# by collection year & site
pcoa.n.4a<-ggplot(b.near.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(CollectionYear),shape=Site), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Collection Year")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Year",values=unique(b.near.pcoa.meta$Year_Color[order(b.near.pcoa.meta$CollectionYear)])) +
  xlab("PC1 [24.03%]") + ylab("PC2 [14.26%]")

ggsave(pcoa.n.4a,filename = "figures/BetaDiversity/SSD_16S_NearSitesOnly_CLR_CollectionYear_Site_PCOA1.png", width=12, height=10, dpi=600)

### PCOAs by Collection Year

# 2020
pcoa.n.2020.1<-ggplot(b.pcoa.meta[b.pcoa.meta$CollectionYear=="2020",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(SampleMonth),shape=Site), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Date",values=unique(b.pcoa.meta$SampMonth_Color[order(b.pcoa.meta$SampleMonth)]),labels=c(unique(b.pcoa.meta$SampleMonth[order(b.pcoa.meta$SampleMonth)]))) +
  xlab("PC1 [23.63%]") + ylab("PC2 [9.73%]")

ggsave(pcoa.n.2020.1,filename = "figures/BetaDiversity/SSD_16S_NearSitesOnly_2020_CLR_sampdate_PCOA1.png", width=12, height=10, dpi=600)

# general season
pcoa.n.2020.2<-ggplot(b.pcoa.meta[b.pcoa.meta$CollectionYear=="2020",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Season_General),shape=CollectionYear), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Season")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Season",values=unique(b.pcoa.meta$SeasonGen_Color[order(b.pcoa.meta$Season_General)]),labels=c(unique(b.pcoa.meta$Season_General[order(b.pcoa.meta$Season_General)]))) +
  xlab("PC1 [23.63%]") + ylab("PC2 [9.73%]")

ggsave(pcoa.n.2020.2,filename = "figures/BetaDiversity/SSD_16S_NearSitesOnly_2020_CLR_season_general_PCOA1.png", width=12, height=10, dpi=600)

# specific season
pcoa.n.2020.3<-ggplot(b.pcoa.meta[b.pcoa.meta$CollectionYear=="2020",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Season_Specific),shape=CollectionYear), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Season")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Season",values=unique(b.pcoa.meta$SeasonSpec_Color[order(b.pcoa.meta$Season_Specific)]),
                     labels=c("Early.Summer"="Early Summer","Late.Summer"="Late Summer","Early.Fall"="Early Fall","Late.Fall"="Late Fall","Fall.Winter"="Fall-Winter")) +
  xlab("PC1 [23.63%]") + ylab("PC2 [9.73%]")

ggsave(pcoa.n.2020.3,filename = "figures/BetaDiversity/SSD_16S_NearSitesOnly_2020_CLR_season_specific_PCOA1.png", width=12, height=10, dpi=600)

# 2021
pcoa.n.2021.1<-ggplot(b.pcoa.meta[b.pcoa.meta$CollectionYear=="2021",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(SampleMonth),shape=Site), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Date",values=unique(b.pcoa.meta$SampMonth_Color[order(b.pcoa.meta$SampleMonth)]),labels=c(unique(b.pcoa.meta$SampleMonth[order(b.pcoa.meta$SampleMonth)]))) +
  xlab("PC1 [23.63%]") + ylab("PC2 [9.73%]")

ggsave(pcoa.n.2021.1,filename = "figures/BetaDiversity/SSD_16S_NearSitesOnly_2021_CLR_sampdate_PCOA1.png", width=12, height=10, dpi=600)

# general season
pcoa.n.2021.2<-ggplot(b.pcoa.meta[b.pcoa.meta$CollectionYear=="2021",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Season_General),shape=CollectionYear), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Season")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Season",values=unique(b.pcoa.meta$SeasonGen_Color[order(b.pcoa.meta$Season_General)]),labels=c(unique(b.pcoa.meta$Season_General[order(b.pcoa.meta$Season_General)]))) +
  xlab("PC1 [23.63%]") + ylab("PC2 [9.73%]")

ggsave(pcoa.n.2021.2,filename = "figures/BetaDiversity/SSD_16S_NearSitesOnly_2021_CLR_season_general_PCOA1.png", width=12, height=10, dpi=600)

# specific season
pcoa.n.2021.3<-ggplot(b.pcoa.meta[b.pcoa.meta$CollectionYear=="2021",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Season_Specific),shape=CollectionYear), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Season")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Season",values=unique(b.pcoa.meta$SeasonSpec_Color[order(b.pcoa.meta$Season_Specific)]),
                     labels=c("Early.Summer"="Early Summer","Late.Summer"="Late Summer","Early.Fall"="Early Fall","Late.Fall"="Late Fall","Fall.Winter"="Fall-Winter")) +
  xlab("PC1 [23.63%]") + ylab("PC2 [9.73%]")

ggsave(pcoa.n.2021.3,filename = "figures/BetaDiversity/SSD_16S_NearSitesOnly_2021_CLR_season_specific_PCOA1.png", width=12, height=10, dpi=600)

#### Homogeneity of Variance & PERMANOVA tests - Composition by Groups ####
## betadisper to look at homogeneity of group dispersions (aka variance) when considering multiple variables
# multivariate analogue to Levene's test of homogeneity of variances
# program finds spatial median or centroid of the group, & compare distances of group to centroid/spatial median via ANOVA

#While PERMANOVA tests differences in group means (analogous to MANOVA),
## a related test called PERMDISP can be used to evaluate homogeneity of group dispersion
#(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
## * need a distance matrix!

rownames(dust_meta) %in% rownames(b.clr) #b.clr was used to make the distance matrix b.euc_dist

# first by compare dispersions by sampling date
b.disper1<-betadisper((vegdist(b.clr,method="euclidean")), dust_meta$SampDate)
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
# ANOVA adjusted p-value
aov.beta.p1<-summary(anova(b.disper1))[[1]][["Pr(>F)"]] # get p values from ANOVA
p.adjust(aov.beta.p1,method="bonferroni",n=length(aov.beta.p1))

TukeyHSD(b.disper1) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other
#                               diff       lwr       upr     p adj
# December.2021-August.2021 -4.556897 -7.631983 -1.481811 0.0033548
# April.2022-August.2021    -5.605532 -8.680617 -2.530446 0.0004428
# April.2022-December.2021  -1.048634 -4.123720  2.026451 0.6710007

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

pnova1<-adonis2(b.clr ~ SampDate,data=dust_meta,method = "euclidean",by="terms",permutations=1000)
pnova1 # p-value = 0.000999
p.adjust(pnova1$`Pr(>F)`,method="bonferroni",n=length(pnova1$`Pr(>F)`)) # adjusted pval

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

b.clr.dist = (vegdist(b.clr, "euclidean", na.rm = TRUE)) #distance matrix using Bray's dissimilarity index for trait distribution (traits of interest only)
pair.mod1<-pairwise.adonis(b.clr.dist,dust_meta$SampDate, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
pair.mod1
#                           pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
# 1  December.2021 vs April.2022  1  7006.287 17.13323 0.5503196   0.001      0.003   *
# 2 December.2021 vs August.2021  1  7416.531 13.48028 0.4905437   0.001      0.003   *
# 3    April.2022 vs August.2021  1  7987.429 15.15788 0.5198554   0.001      0.003   *

# Visualize dispersions
png('figures/BetaDiversity/pcoa_betadispersion_sampledate.png',width = 700, height = 600, res=100)
plot(b.disper1,main = "Centroids and Dispersion based on Aitchison Distance", col=colorset1$SampMonth_Color)
dev.off()

png('boxplot_centroid_distance_sampledate.png',width = 700, height = 600, res=100)
boxplot(b.disper1,xlab="Sample Collection Date", main = "Distance to Centroid by Category", sub="Based on Aitchison Distance", col=colorset1$SampMonth_Color)
dev.off()

# Next compare dispersions by depth
b.disper2<-betadisper((vegdist(b.clr,method="euclidean")), dust_meta$Depth_m)
b.disper2

permutest(b.disper2, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons

anova(b.disper2) # p = 0.6277 --> accept the Null H, spatial medians are NOT significantly difference across sample dates

# ANOVA adjusted p-value
aov.beta.p2<-summary(anova(b.disper2))[[1]][["Pr(>F)"]] # get p values from ANOVA
p.adjust(aov.beta.p2,method="bonferroni",n=length(aov.beta.p1))

TukeyHSD(b.disper2) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other
# no sig results

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

pnova2<-adonis2(b.clr ~ Depth_m,data=dust_meta,method = "euclidean",by="terms",permutations=1000)
pnova2 # p-value = 1
p.adjust(pnova2$`Pr(>F)`,method="bonferroni",n=length(pnova2$`Pr(>F)`)) # adjusted pval

#b.clr.dist = (vegdist(b.clr, "euclidean", na.rm = TRUE)) #distance matrix using Bray's dissimilarity index for trait distribution (traits of interest only)
pair.mod2<-pairwise.adonis(b.clr.dist,dust_meta$Depth_m, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
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
