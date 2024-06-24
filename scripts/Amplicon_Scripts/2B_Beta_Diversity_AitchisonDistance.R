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
  library("factoextra")
  library(reticulate)
})

#### Load Global Env to Import Count/ASV Tables ####
load("data/Amplicon/SSDust_16S.V3V4_W23_Data_Ready.Rdata") # save global env to Rdata file
#load("data/Amplicon/SSD_16S_CLR_EucDist_Ready.Rdata")

# NOTE about ASV table:
# ASV counts were divided by the number of deployment days per sample collector, then multiplied by 100 & rounded
# this is so that all ASV counts are scaled by deployment days
# bac.ASV_round.table is the ASV table with scaled ASV counts that have been * 100 and rounded

head(b.dust.all)
bac.ASV_round.table[1:4,1:4]
bac.ASV_round.table[(nrow(bac.ASV_round.table)-4):(nrow(bac.ASV_round.table)),(ncol(bac.ASV_round.table)-4):(ncol(bac.ASV_round.table))] # last 4 rows & cols
head(dust_meta)

# create metadata dfs by categories for downstream comparisons
summer_meta<-dust_meta[dust_meta$Season_General=="Summer",]
meta.2020<-dust_meta[dust_meta$CollectionYear=="2020",]
meta.2021<-dust_meta[dust_meta$CollectionYear=="2021",]

# merge precip data and dust_meta to create rain vs no rain category (for PERMANOVAs later)

precip.meta<-merge(dust_meta,precip.data,by=c("Deploy_dth","Collect_dth","Precip_STID"))
precip.meta$RainCat<-ifelse(precip.meta$precip_24hr_accum>0,"Rain","No Rain")

# #### PCoA on Raw Data - Sanity Check ####
#
# # check rownames of CLR transformed ASV data & metadata
# rownames(bac.ASV_round.table) %in% rownames(dust_meta)
# dust_meta=dust_meta[rownames(bac.ASV_round.table[,-1]),] ## reorder metadata to match order of CLR data
#
# # calculate our Euclidean distance matrix using CLR data
# raw.b.euc_dist <- dist(bac.ASV_round.table[,-1], method = "euclidean")
#
# # creating our hierarcical clustering dendrogram
# raw.b.euc_clust <- hclust(raw.b.euc_dist, method="ward.D2")
#
# # let's make it a little nicer...
# raw.b.euc_dend <- as.dendrogram(raw.b.euc_clust, hang=0.2)
# raw.b.dend_cols <- as.character(dust_meta$SampDate_Color[order.dendrogram(raw.b.euc_dend)])
# labels_colors(raw.b.euc_dend) <- raw.b.dend_cols
#
# ## DO NOT RUN THIS LINE, THIS IS YOUR COLOR REFERENCE!!!!
# #(August.2021="#ef781c",December.2021="#03045e",April.2022="#059c3f")
#
# plot(raw.b.euc_dend, ylab="CLR Euclidean Distance", cex = 0.1,horiz=TRUE) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 2)
# #legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c("#ef781c","#03045e","#059c3f"),pch = 15, bty = "n")
# # Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
# dev.off()
#
# # PCOA w/ Euclidean distance matrix (of CLR data)
# raw.b.pcoa <- pcoa(raw.b.euc_dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
# #save.image("data/Amplicon/SSD_16S_CLR_EucDist_Ready.Rdata")
#
# # The proportion of variances explained is in its element values$Relative_eig
# raw.b.pcoa$values
#
# # extract principal coordinates
# raw.b.pcoa.vectors<-data.frame(raw.b.pcoa$vectors)
# raw.b.pcoa.vectors$SampleID<-rownames(raw.b.pcoa$vectors)
#
# # merge pcoa coordinates w/ metadata
# raw.b.pcoa.meta<-merge(raw.b.pcoa.vectors, dust_meta, by.x="SampleID", by.y="SampleID")
# raw.b.pcoa.meta$SampleMonth
# raw.b.pcoa.meta$SampDate
#
# head(raw.b.pcoa.meta)
#
# head(raw.b.pcoa$values) # pull out Relative (Relative_eig) variation % to add to axes labels
# # PC1 = 50.78%, PC2 = 18.29%
#
# # drop outliars
# #outliars<-c("BDC.D.7.27.21","PD.D.7.27.21","WI.D.9.18.21","WI.D.7.10.20")
# #raw.b.pcoa.meta2<-raw.b.pcoa.meta[!(raw.b.pcoa.meta$SampleID %in% outliars),]
#
# #### Visualize RAW PCoAs####
#
# # create PCoA ggplot fig
# ggplot(raw.b.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(SampleMonth),shape=CollectionYear), size=5)+theme_bw()+
#   labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Collection Date",values=unique(raw.b.pcoa.meta$SampMonth_Color[order(raw.b.pcoa.meta$SampleMonth)]),labels=c(unique(raw.b.pcoa.meta$SampleMonth[order(raw.b.pcoa.meta$SampleMonth)]))) +
#   xlab("PC1 [48.24%]") + ylab("PC2 [17.3%]")
#
# ggplot(raw.b.pcoa.meta2, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site),shape=CollectionYear), size=5)+theme_bw()+
#   labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Collection Date",values=unique(raw.b.pcoa.meta$Site_Color[order(raw.b.pcoa.meta$Site)]),labels=c(unique(raw.b.pcoa.meta$Site[order(raw.b.pcoa.meta$Site)]))) +
#   xlab("PC1 [48.24%]") + ylab("PC2 [8.31%]")
#
# # by collection period & site
# ggplot(raw.b.pcoa.meta2, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year),shape=Site), size=5)+theme_bw()+
#   labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Collection Period",values=unique(raw.b.pcoa.meta$SCY_Color[order(raw.b.pcoa.meta$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   xlab("PC1 [48.24%]") + ylab("PC2 [8.31%]")
#

#### CLR Transform All Comp Data ####
rownames(bac.ASV_round.table)
bac.ASV_round.table[1:4,1:4]

# CLR transformation of ASV table
# df must have rownames are SampleIDs, columns are ASV IDs for vegan functions below
b.clr<-decostand(bac.ASV_round.table[,-1],method = "clr", pseudocount = 1) #CLR transformation
b.clr[1:4,1:4]

#### Beta Diversity - All Data ####

# check rownames of CLR transformed ASV data & metadata
rownames(b.clr) %in% rownames(dust_meta)
dust_meta=dust_meta[rownames(b.clr),] ## reorder metadata to match order of CLR data

# calculate our Euclidean distance matrix using CLR data
b.euc_dist <- dist(b.clr, method = "euclidean")
# b.euc_dist2 <- vegdist(b.clr, method = "euclidean") # did this to compare the two matrices to see ifthey're identical
fviz_dist(b.euc_dist, gradient = list(low = "blue", mid = "white", high = "red"))

# creating our hierarcical clustering dendrogram
b.euc_clust <- hclust(b.euc_dist, method="ward.D2")

# let's make it a little nicer...
b.euc_dend <- as.dendrogram(b.euc_clust, hang=0.02)
b.dend_cols <- as.character(dust_meta$SampDate_Color[order.dendrogram(b.euc_dend)])
labels_colors(b.euc_dend) <- b.dend_cols

colorset8 # color dendrogram by collection date
png(filename="figures/BetaDiversity/Aitchison/SSD_16S_CLR_EucDist_Dendrogram1.png",width = 7, height = 7, units = "in",res = 800)
plot(b.euc_dend, ylab="CLR Euclidean Distance", cex = 0.1,horiz=TRUE) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 2)
#legend("topleft",legend = colorset8$SampDate,cex=.8,col = colorset8$SampDate_Color,pch = 15, bty = "n")
dev.off()

# let's try another hierarchical cluster dendrogram with eclust()
b.clr.hc<-eclust(b.clr, "hclust",nboot = 500) # use eclust() for "enhanced" hierarchical clustering; more here: https://www.datanovia.com/en/blog/cluster-analysis-in-r-simplified-and-enhanced/
fviz_dend(b.clr.hc, rect = TRUE) # plot the dendrogam

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
rownames(b.pcoa.meta)<-b.pcoa.meta$SampleID

head(b.pcoa$values) # pull out Relative (Relative_eig) variation % to add to axes labels
# PC1 = 24.47%, PC2 = 8.31%

# K-means Clustering - will show us the groups that also appear in our PCoA
## more on K-means clustering: http://www.sthda.com/english/wiki/factoextra-r-package-easy-multivariate-data-analyses-and-elegant-visualization
## more on K-means clustering: https://uc-r.github.io/kmeans_clustering#:~:text=The%20gap%20statistic%20compares%20the,simulations%20of%20the%20sampling%20process.

## first we are going to play around with different values for k (aka # of clusters) and plot them to see what seems best
k2 <- kmeans(b.clr, centers = 2, nstart = 25)
k3 <- kmeans(b.clr, centers = 3, nstart = 25)
k4 <- kmeans(b.clr, centers = 4, nstart = 25)
k5 <- kmeans(b.clr, centers = 5, nstart = 25)

# plots to compare different values for k (aka # clusters)
p1 <- fviz_cluster(k2, geom = "point",  data = b.clr) + ggtitle("k = 2")
p2 <- fviz_cluster(k3, geom = "point",  data = b.clr) + ggtitle("k = 3")
p3 <- fviz_cluster(k4, geom = "point",  data = b.clr) + ggtitle("k = 4")
p4 <- fviz_cluster(k5, geom = "point",  data = b.clr) + ggtitle("k = 5")

# next we can calculate the K-means clustering and get back the ideal # of clusters
b.clr.km <- eclust(b.clr, "kmeans", hc_metric="euclid",nboot = 1000)
b.clr.km

# Visualize kmeans clustering with output form eclust()
# use repel = TRUE to avoid overplotting
fviz_cluster(b.clr.km, b.clr, ellipse = TRUE,
             ellipse.level = 0.95, ellipse.alpha = 0.2,ellipse.type = "convex",outlier.color = "black")

kmeans.plot1<-fviz_cluster(b.clr.km, b.clr, ellipse = TRUE,
                          ellipse.level = 0.95, ellipse.alpha = 0.2,ellipse.type = "euclid",outlier.color = "black") +
  ggtitle("K-means = 3")
ggsave(kmeans.plot1,filename = "figures/BetaDiversity/Aitchison/SSD_16S_CLR_KMeans_Cluster_Euclid_3clusters_PCoA1.png", width=14, height=10, dpi=600,create.dir = TRUE)

kmeans.plot2<-fviz_cluster(b.clr.km, b.clr, ellipse = TRUE,
                           ellipse.level = 0.95, ellipse.alpha = 0.2,ellipse.type = "convex",outlier.color = "black") +
  ggtitle("K-means = 3")
ggsave(kmeans.plot2,filename = "figures/BetaDiversity/Aitchison/SSD_16S_CLR_KMeans_Cluster_Convex_3clusters_PCoA2.png", width=14, height=10, dpi=600,create.dir = TRUE)


# Gap statistic plot to determine the ideal # of clusters k
# The gap statistic compares the total intracluster variation for different values of k with their expected values under null reference distribution of the data (i.e. a distribution with no obvious clustering)
fviz_gap_stat(b.clr.km$gap_stat) # shows 3 clusters
fviz_nbclust(b.clr, kmeans, method = "gap_stat") # also shows 3 clusters
fviz_nbclust(b.clr, kmeans, method = "silhouette") # shows 2 clusters

#fviz_silhouette(b.clr.km)

# Optimal number of clusters using gap statistics
b.clr.km$nbclust

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
  scale_shape_manual(values = c(0,1,16,15)) +
  xlab("PC1 [24.47%]") + ylab("PC2 [8.31%]")

ggsave(pcoa1,filename = "figures/BetaDiversity/Aitchison/SSD_16S_CLR_SampDate_Site_Year_PCOA1.png", width=14, height=10, dpi=600,create.dir = TRUE)

# pcoa1a<-ggplot(b.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site),shape=CollectionYear), size=5)+theme_bw()+
#   labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Collection Date",values=unique(b.pcoa.meta$Site_Color[order(b.pcoa.meta$Site)]),labels=c(unique(b.pcoa.meta$Site[order(b.pcoa.meta$Site)]))) +
#   xlab("PC1 [24.47%]") + ylab("PC2 [8.31%]")
#
# ggsave(pcoa1a,filename = "figures/BetaDiversity/Aitchison/SSD_16S_CLR_Site_Year_PCOA1.png", width=12, height=10, dpi=600,create.dir = TRUE)
#
# # specific season
# pcoa3<-ggplot(b.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Season_Specific),shape=CollectionYear), size=5)+theme_bw()+
#   labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Season")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Season",values=unique(b.pcoa.meta$SeasonSpec_Color[order(b.pcoa.meta$Season_Specific)]),
#                      labels=c("Early.Summer"="Early Summer","Late.Summer"="Late Summer","Early.Fall"="Early Fall","Late.Fall"="Late Fall","Fall.Winter"="Fall-Winter")) +
#   xlab("PC1 [24.47%]") + ylab("PC2 [8.31%]")
#
# ggsave(pcoa3,filename = "figures/BetaDiversity/Aitchison/SSD_16S_CLR_season_specific_PCOA1.png", width=12, height=10, dpi=600,create.dir = TRUE)
#
# # by collection year & site
# pcoa4a<-ggplot(b.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(CollectionYear),shape=Site), size=5)+theme_bw()+
#   labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Collection Year")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Collection Year",values=unique(b.pcoa.meta$Year_Color[order(b.pcoa.meta$CollectionYear)])) +
#   xlab("PC1 [24.47%]") + ylab("PC2 [8.31%]")
#
# ggsave(pcoa4a,filename = "figures/BetaDiversity/Aitchison/SSD_16S_CLR_CollectionYear_Site_PCOA1.png", width=12, height=10, dpi=600,create.dir = TRUE)
#
# # by collection period & site
# pcoa5<-ggplot(b.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year),shape=Site), size=5)+theme_bw()+
#   labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Collection Period",values=unique(b.pcoa.meta$SCY_Color[order(b.pcoa.meta$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   xlab("PC1 [24.47%]") + ylab("PC2 [8.31%]")
#
# ggsave(pcoa5,filename = "figures/BetaDiversity/Aitchison/SSD_16S_CLR_SeasCollYr_Site_PCOA1.png", width=12, height=10, dpi=600,create.dir = TRUE)

# 3D PCoA

pltly.all.a<-plot_ly(b.pcoa.meta, x=~Axis.1,y=~Axis.2,z=~Axis.3, color = ~SampDate, colors = c(unique(b.pcoa.meta$SampDate_Color[order(b.pcoa.meta$SampDate)])),
        symbol=~Site,symbols = c("square-open", "circle-open","circle","square")) %>%
  layout(scene = list(xaxis = list(title = 'PC1 24.47%'),
                      yaxis = list(title = 'PC2 8.31%'),
                      zaxis = list(title = 'PC3 7.31%')))

# before you can run save_image(), run the following lines; follow instructions: https://search.r-project.org/CRAN/refmans/plotly/html/save_image.html
#install.packages('reticulate')
# reticulate::install_miniconda()
# reticulate::conda_install('r-reticulate', 'python-kaleido')
# reticulate::conda_install('r-reticulate', 'plotly', channel = 'plotly')
reticulate::use_miniconda('r-reticulate')

save_image(pltly.all.a, "figures/BetaDiversity/Aitchison/SSD_16S_CLR_SampDate_CollYr_Site_3D_Aitchison_PCOA1.png",width=1200,height=1000)
save_image(pltly.all.a, "figures/BetaDiversity/Aitchison/SSD_16S_CLR_SampDate_CollYr_Site_3D_Aitchison_PCOA2.png",width=1400,height=1100)


# save 3D plot as an HTml
saveWidget(widget = pltly.all.a, #the plotly object,
           file = "figures/BetaDiversity/Aitchison/SSD_16S_CLR_SampDate_CollYr_Site_3D_Aitchison_PCOA1.html", #the path & file name
           selfcontained = TRUE #creates a single html file
)

## visualizing only specific sites

# WI - by collection period
wi.pc1<-ggplot(b.pcoa.meta[b.pcoa.meta$Site=="WI",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year)), size=5)+theme_bw()+
  labs(title="PCoA: Wister Microbial Composition",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.pcoa.meta$SCY_Color[order(b.pcoa.meta$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
  xlab("PC1 [24.47%]") + ylab("PC2 [8.31%]")

ggsave(wi.pc1,filename = "figures/BetaDiversity/Aitchison/SSD_16S_CLR_SeasCollYr_Wister_PCOA1.png", width=12, height=10, dpi=600,create.dir = TRUE)

dp.pc1<-ggplot(b.pcoa.meta[b.pcoa.meta$Site=="DP",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year)), size=5)+theme_bw()+
  labs(title="PCoA: Dos Palmas Microbial Composition",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.pcoa.meta$SCY_Color[order(b.pcoa.meta$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
  xlab("PC1 [24.47%]") + ylab("PC2 [8.31%]")

ggsave(dp.pc1,filename = "figures/BetaDiversity/Aitchison/SSD_16S_CLR_SeasCollYr_DosPalmas_PCOA1.png", width=12, height=10, dpi=600,create.dir = TRUE)

s.r.pc1<-ggplot(b.pcoa.meta[b.pcoa.meta$Site=="SB" | b.pcoa.meta$Site=="RHB",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year),shape=Site), size=5)+theme_bw()+
  labs(title="PCoA: Sonny Bono & Red Hill Bay Microbial Composition",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.pcoa.meta$SCY_Color[order(b.pcoa.meta$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
  xlab("PC1 [24.47%]") + ylab("PC2 [8.31%]")

ggsave(s.r.pc1,filename = "figures/BetaDiversity/Aitchison/SSD_16S_CLR_SeasCollYr_SB_and_RHB_PCOA1.png", width=12, height=10, dpi=600,create.dir = TRUE)

pd.pc1<-ggplot(b.pcoa.meta[b.pcoa.meta$Site=="PD",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year)), size=5)+theme_bw()+
  labs(title="PCoA: Palm Desert Microbial Composition",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.pcoa.meta$SCY_Color[order(b.pcoa.meta$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
  xlab("PC1 [24.47%]") + ylab("PC2 [8.31%]")

ggsave(pd.pc1,filename = "figures/BetaDiversity/Aitchison/SSD_16S_CLR_SeasCollYr_PD_PCOA1.png", width=12, height=10, dpi=600,create.dir = TRUE)

bdc.pc1<-ggplot(b.pcoa.meta[b.pcoa.meta$Site=="BDC",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year)), size=5)+theme_bw()+
  labs(title="PCoA: Boyd Deep Canyon Microbial Composition",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.pcoa.meta$SCY_Color[order(b.pcoa.meta$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
  xlab("PC1 [24.47%]") + ylab("PC2 [8.31%]")

ggsave(bdc.pc1,filename = "figures/BetaDiversity/Aitchison/SSD_16S_CLR_SeasCollYr_BDC_PCOA1.png", width=12, height=10, dpi=600,create.dir = TRUE)

# visualizing pcoa by month per year
unique(b.pcoa.meta$SampDate)

j2020.p1<-ggplot(b.pcoa.meta[b.pcoa.meta$SampDate=="July.2020",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: July 2020 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.pcoa.meta$Site_Color[order(b.pcoa.meta$Site)])) +
  xlab("PC1 [24.47%]") + ylab("PC2 [8.31%]")

ggsave(j2020.p1,filename = "figures/BetaDiversity/Aitchison/SSD_16S_CLR_AllSites_July2020_PCOA1.png", width=12, height=10, dpi=600,create.dir = TRUE)

a2020.p1<-ggplot(b.pcoa.meta[b.pcoa.meta$SampDate=="August.2020",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: August 2020 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.pcoa.meta$Site_Color[order(b.pcoa.meta$Site)])) +
  xlab("PC1 [24.47%]") + ylab("PC2 [8.31%]")

ggsave(a2020.p1,filename = "figures/BetaDiversity/Aitchison/SSD_16S_CLR_AllSites_Aug2020_PCOA1.png", width=12, height=10, dpi=600,create.dir = TRUE)

o2020.p1<-ggplot(b.pcoa.meta[b.pcoa.meta$SampDate=="October.2020",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: October 2020 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.pcoa.meta$Site_Color[order(b.pcoa.meta$Site)])) +
  xlab("PC1 [24.47%]") + ylab("PC2 [8.31%]")

ggsave(o2020.p1,filename = "figures/BetaDiversity/Aitchison/SSD_16S_CLR_AllSites_Oct2020_PCOA1.png", width=12, height=10, dpi=600,create.dir = TRUE)

n2020.p1<-ggplot(b.pcoa.meta[b.pcoa.meta$SampDate=="November.2020",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: November 2020 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.pcoa.meta$Site_Color[order(b.pcoa.meta$Site)])) +
  xlab("PC1 [24.47%]") + ylab("PC2 [8.31%]")

ggsave(n2020.p1,filename = "figures/BetaDiversity/Aitchison/SSD_16S_CLR_AllSites_Nov2020_PCOA1.png", width=12, height=10, dpi=600,create.dir = TRUE)

ja2021.p1<-ggplot(b.pcoa.meta[b.pcoa.meta$SampDate=="July.2021" | b.pcoa.meta$SampDate=="August.2021",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: July & August 2021 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.pcoa.meta$Site_Color[order(b.pcoa.meta$Site)])) +
  xlab("PC1 [24.47%]") + ylab("PC2 [8.31%]")

ggsave(ja2021.p1,filename = "figures/BetaDiversity/Aitchison/SSD_16S_CLR_AllSites_July_Aug_2021_PCOA1.png", width=12, height=10, dpi=600,create.dir = TRUE)

s2021.p1<-ggplot(b.pcoa.meta[b.pcoa.meta$SampDate=="September.2021",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: September 2021 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.pcoa.meta$Site_Color[order(b.pcoa.meta$Site)])) +
  xlab("PC1 [24.47%]") + ylab("PC2 [8.31%]")

ggsave(s2021.p1,filename = "figures/BetaDiversity/Aitchison/SSD_16S_CLR_AllSites_Sept2021_PCOA1.png", width=12, height=10, dpi=600,create.dir = TRUE)

d2021.p1<-ggplot(b.pcoa.meta[b.pcoa.meta$SampDate=="December.2021",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: December 2021 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.pcoa.meta$Site_Color[order(b.pcoa.meta$Site)])) +
  xlab("PC1 [24.47%]") + ylab("PC2 [8.31%]")

ggsave(d2021.p1,filename = "figures/BetaDiversity/Aitchison/SSD_16S_CLR_AllSites_Dec2021_PCOA1.png", width=12, height=10, dpi=600,create.dir = TRUE)

## visualizing by year
b.pcoa1.20<-b.pcoa.meta[b.pcoa.meta$CollectionYear=="2020",]
b.pcoa1.21<-b.pcoa.meta[b.pcoa.meta$CollectionYear=="2021",]

# 2020
twntytwnty.pc1<-ggplot(b.pcoa1.20, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year),shape=Site), size=5)+theme_bw()+
  labs(title="PCoA: Microbial Composition in 2020",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.pcoa1.20$SCY_Color[order(b.pcoa1.20$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
  xlab("PC1 [24.47%]") + ylab("PC2 [8.31%]")

ggsave(twntytwnty.pc1,filename = "figures/BetaDiversity/Aitchison/SSD_16S_CLR_SeasCollYr_2020_PCOA1.png", width=12, height=10, dpi=600,create.dir = TRUE)

# 2021
twntytwnty1.pc1<-ggplot(b.pcoa1.21, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year),shape=Site), size=5)+theme_bw()+
  labs(title="PCoA: Microbial Composition in 2021",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.pcoa1.21$SCY_Color[order(b.pcoa1.21$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
  xlab("PC1 [24.47%]") + ylab("PC2 [8.31%]")

ggsave(twntytwnty1.pc1,filename = "figures/BetaDiversity/Aitchison/SSD_16S_CLR_SeasCollYr_2021_PCOA1.png", width=12, height=10, dpi=600,create.dir = TRUE)

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

ggsave(pcoa.axis1.hm,filename = "figures/BetaDiversity/Aitchison/SSD_16S_CLR_PCoA_PC1_Site_by_SampDate.png", width=18, height=13, dpi=600,create.dir = TRUE)

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

ggsave(pcoa.axis2.hm,filename = "figures/BetaDiversity/Aitchison/SSD_16S_CLR_PCoA_PC2_Site_by_SampDate.png", width=18, height=13, dpi=600,create.dir = TRUE)

## Loop to Generate Heat Map for Each PC Axis

pc.plot.list<-list() # create empty list for each plot to be stored in
pc.axes<-names(b.pcoa.dts)[grepl("Axis",names(b.pcoa.dts))] # pull out names of columns in df that contain "Axis" in name

# heatmap function that you will use to generate each heatmap
hm.fxn<-function(df, x_var, y_var, f_var) {
  a <- ggplot(df, aes(x = x_var, y = y_var, fill = f_var)) +
    geom_tile(color = "black") +
    scale_fill_gradient(low = "blue", high = "red") +
    geom_text(aes(label = round(f_var,2)), color = "white", size = 4) +
    coord_fixed() +
    theme_classic() +
    theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
          axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
          axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),text=element_text(size=14)) +
    labs(x = "",
         y = "",
         fill = "R",
         title = as.character(f_var))


  return(a)
}

# loop through variable containing col names in df of interest
# create heatmap, then adjust title and legend title, add plot to a list of plots, then save plot to file
for (i in pc.axes) {
  pc.heatmap=hm.fxn(b.pcoa.dts, b.pcoa.dts$Site, b.pcoa.dts$SampDate, b.pcoa.dts[,i])
  hm_titled = pc.heatmap + ggtitle(as.character(i)) + guides(fill=guide_legend(title="PC Values"))
  pc.plot.list[[i]]=hm_titled
  ggsave(hm_titled,filename = paste("figures/BetaDiversity/PCoA_Axes_Heatmaps/SSD_16S_PCoA_Site_SampDate_heatmap_",i,".png",sep=""), width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

}


# call each plot in list by index to view plots
pc.plot.list[[1]] # PC 1
pc.plot.list[[2]] # PC 2
pc.plot.list[[9]] # PC 9
pc.plot.list[[27]] # PC 27

#
# loop.hm<-function(df, x_var, y_var, f_var) {
#   a <- ggplot(df, aes(x = x_var, y = y_var, fill = f_var)) +
#     geom_tile(color = "black") +
#     scale_fill_gradient(low = "blue", high = "red") +
#     geom_text(aes(label = round(f_var,2)), color = "white", size = 4) +
#     coord_fixed() +
#     theme_minimal() +
#     labs(x = "",
#          y = "",
#          fill = "R", # Want the legend title to be each of the column names that are looped
#          title = as.character(f_var))
#
#   #ggsave(a, file = paste0("figures/BetaDiversity/Aitchison/SSD_16S.PCoA_heatmap_", f_var,".png"), device = png, width = 15, height = 15, units = "cm")
#
#   return(a)
# }
# # loop with heatmap function to create heatmap


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

#ggsave(sulf.hm1a,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_SampID_by_Function_heatmap.png", width=18, height=13, dpi=600,create.dir = TRUE)

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
# PD          0.11500 0.25000 0.224
# BDC 0.11197         0.51500 0.773
# DP  0.25518 0.52841         0.801
# WI  0.21506 0.78889 0.77683

anova(b.disper1) # p = 0.3424 --> accept the Null H, spatial medians (a measure of dispersion) are NOT significantly difference across sample dates
# ANOVA adjusted p-value
aov.beta.p1<-anova(b.disper1)[["Pr(>F)"]] # get p values from ANOVA
p.adjust(aov.beta.p1,method="bonferroni",n=length(aov.beta.p1))

TukeyHSD(b.disper1) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other
#         diff       lwr       upr     p adj
# BDC-PD -47.052166 -121.70958 27.60524 0.3266531
# DP-PD  -31.777129 -106.43454 42.88028 0.6484299
# WI-PD  -39.515430 -114.17284 35.14198 0.4760740
# DP-BDC  15.275037  -59.38237 89.93245 0.9416395
# WI-BDC   7.536736  -67.12067 82.19415 0.9922594
# WI-DP   -7.738301  -82.39571 66.91911 0.9916363

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

pnova1<-adonis2(b.clr ~ Site,data=dust_meta,method = "euclidean",by="terms",permutations= 10000)
pnova1 # p-value = 0.0498
p.adjust(pnova1$`Pr(>F)`,method="bonferroni",n=length(pnova1$`Pr(>F)`)) # adjusted pval
# 0.1493851

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

b.clr.dist = (vegdist(b.clr, "euclidean", na.rm = TRUE)) #distance matrix using Bray's dissimilarity index for trait distribution (traits of interest only)
pair.mod1<-pairwise.adonis(b.clr.dist,dust_meta$Site, p.adjust.m='bonferroni',perm=9999) # shows us variation for each sample to see which ones are different
pair.mod1
#     pairs Df SumsOfSqs  F.Model        R2     p.value     p.adjusted sig
# 1 BDC vs DP  1  15216.07 1.225157 0.09263834  0.1349     0.8094
# 2 BDC vs PD  1  24238.00 1.420882 0.10587098  0.1490     0.8940
# 3 BDC vs WI  1  12253.51 1.025153 0.07870566  0.3068     1.0000
# 4  DP vs PD  1  28382.01 1.523036 0.11262533  0.0944     0.5664
# 5  DP vs WI  1  13648.12 1.008762 0.07754479  0.3512     1.0000
# 6  PD vs WI  1  30062.35 1.654658 0.12117904  0.0732     0.4392

# Visualize dispersions
png('figures/BetaDiversity/Aitchison/SSD_pcoa_betadispersion_site.png',width = 700, height = 600, res=100)
plot(b.disper1,main = "Centroids and Dispersion based on Aitchison Distance", col=colorset6$Site_Color)
dev.off()

png('figures/BetaDiversity/Aitchison/SSD_boxplot_centroid_distance_site.png',width = 700, height = 600, res=100)
boxplot(b.disper1,xlab="By Site", main = "Distance to Centroid by Category", sub="Based on Aitchison Distance", col=colorset6$Site_Color)
dev.off()


### now compare dispersions by site + year
b.disper2<-betadisper((vegdist(b.clr,method="euclidean")), interaction(dust_meta$Site,dust_meta$CollectionYear,sep="."))
b.disper2

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:

permutest(b.disper2, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:

anova(b.disper2) # p = 0.7639 --> accept the Null H, spatial medians (a measure of dispersion) are NOT significantly difference across sample dates
# ANOVA adjusted p-value
aov.beta.p2<-anova(b.disper2)[["Pr(>F)"]] # get p values from ANOVA
p.adjust(aov.beta.p2,method="bonferroni",n=length(aov.beta.p2))

TukeyHSD(b.disper2) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other
#                     diff       lwr       upr     p adj

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

pnova2<-adonis2(b.clr ~ Site*CollectionYear,data=dust_meta,method = "euclidean",by="terms",permutations= 10000)
pnova2
#                       Df SumOfSqs      R2      F  Pr(>F)
# Site                 2    26542 0.14051 1.2025 0.06792 .
# CollectionYear       1    10282 0.02992 1.1101 0.22776
# Site:CollectionYear  2    26105 0.10028 0.9204 0.57242
# Residual            20   187040 0.71920
# Total               27   260069 1.00000

p.adjust(pnova2$`Pr(>F)`,method="bonferroni",n=length(pnova2$`Pr(>F)`)) # adjusted pval
# 0.2296602 1.0000000 1.0000000        NA        NA

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

b.clr.dist = (vegdist(b.clr, "euclidean", na.rm = TRUE)) #distance matrix using Bray's dissimilarity index for trait distribution (traits of interest only)
pair.mod2<-pairwise.adonis(b.clr.dist,interaction(dust_meta$Site,dust_meta$CollectionYear,sep="."), p.adjust.m='bonferroni',perm=999) # shows us variation for each sample to see which ones are different
pair.mod2

# Visualize dispersions
png('figures/BetaDiversity/Aitchison/SSD_pcoa_betadispersion_site_by_year.png',width = 700, height = 600, res=100)
plot(b.disper2,main = "Centroids and Dispersion based on Aitchison Distance", col=colorset6$Site_Color)
dev.off()

png('figures/BetaDiversity/Aitchison/SSD_boxplot_centroid_distance_site_by_year.png',width = 900, height = 600, res=100)
boxplot(b.disper2,xlab="By Site x Collection Year", main = "Distance to Centroid by Category", sub="Based on Aitchison Distance", col=colorset6$Site_Color)
dev.off()

### now compare dispersions by rain vs no rain
b.disper3<-betadisper((vegdist(b.clr,method="euclidean")), precip.meta$RainCat)
b.disper3

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:

permutest(b.disper3, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:

anova(b.disper3) # p = 0.9402 --> accept the Null H, spatial medians (a measure of dispersion) are NOT significantly difference across sample dates
# ANOVA adjusted p-value
aov.beta.p3<-anova(b.disper3)[["Pr(>F)"]] # get p values from ANOVA
p.adjust(aov.beta.p3,method="bonferroni",n=length(aov.beta.p3))

TukeyHSD(b.disper3) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other
#                     diff       lwr       upr     p adj

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

pnova3<-adonis2(b.clr ~ RainCat,data=precip.meta,method = "euclidean",by="terms",permutations= 10000)
pnova3
# not sig
p.adjust(pnova3$`Pr(>F)`,method="bonferroni",n=length(pnova3$`Pr(>F)`)) # adjusted pval

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

b.clr.dist = (vegdist(b.clr, "euclidean", na.rm = TRUE)) #distance matrix using Bray's dissimilarity index for trait distribution (traits of interest only)
pair.mod3<-pairwise.adonis(b.clr.dist,precip.meta$RainCat, p.adjust.m='bonferroni',perm=9999) # shows us variation for each sample to see which ones are different
pair.mod3

# Visualize dispersions
png('figures/BetaDiversity/Aitchison/SSD_pcoa_betadispersion_site_Rain_vs_NoRain.png',width = 700, height = 600, res=100)
plot(b.disper3,main = "Centroids and Dispersion based on Aitchison Distance")
dev.off()

png('figures/BetaDiversity/Aitchison/SSD_boxplot_centroid_distance_site_Rain_vs_NoRain.png',width = 900, height = 600, res=100)
boxplot(b.disper3,xlab="By Site x Collection Year", main = "Distance to Centroid by Category", sub="Based on Aitchison Distance")
dev.off()

### now compare dispersions by site and rain vs no rain
b.disper4<-betadisper((vegdist(b.clr,method="euclidean")), interaction(precip.meta$Site,precip.meta$RainCat,sep="."))
b.disper4

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:

permutest(b.disper4, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:

anova(b.disper4) # p =  --> accept the Null H, spatial medians (a measure of dispersion) are NOT significantly difference across sample dates
# ANOVA adjusted p-value
aov.beta.p4<-anova(b.disper4)[["Pr(>F)"]] # get p values from ANOVA
p.adjust(aov.beta.p4,method="bonferroni",n=length(aov.beta.p4))

TukeyHSD(b.disper4) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other
#                     diff       lwr       upr     p adj

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

pnova4<-adonis2(b.clr ~ Site*RainCat,data=precip.meta,method = "euclidean",by="terms",permutations= 10000)
pnova4
# not sig
p.adjust(pnova4$`Pr(>F)`,method="bonferroni",n=length(pnova4$`Pr(>F)`)) # adjusted pval

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

b.clr.dist = (vegdist(b.clr, "euclidean", na.rm = TRUE)) #distance matrix using Bray's dissimilarity index for trait distribution (traits of interest only)
pair.mod4<-pairwise.adonis(b.clr.dist,interaction(precip.meta$Site,precip.meta$RainCat,sep="."), p.adjust.m='bonferroni',perm=9999) # shows us variation for each sample to see which ones are different
pair.mod4

# Visualize dispersions
png('figures/BetaDiversity/Aitchison/SSD_pcoa_betadispersion_site_Rain_vs_NoRain.png',width = 700, height = 600, res=100)
plot(b.disper4,main = "Centroids and Dispersion based on Aitchison Distance", col=colorset6$Site_Color)
dev.off()

png('figures/BetaDiversity/Aitchison/SSD_boxplot_centroid_distance_site_Rain_vs_NoRain.png',width = 900, height = 600, res=100)
boxplot(b.disper4,xlab="By Site x Collection Year", main = "Distance to Centroid by Category", sub="Based on Aitchison Distance", col=colorset6$Site_Color)
dev.off()

#### Rank Distance Comparison with ANOSIM ####
# NOTES:
# The differences observed between ANOSIM and PERMANOVA demonstrate the importance of understanding the nuances of different statistical methods and their underlying assumptions.
# While PERMANOVA assesses the multivariate dispersion between groups as a whole, ANOSIM focuses on rank distances between pairs of samples.
# The pairwise approach of ANOSIM allows for a more detailed examination of specific comparisons.

# Null hypothesis: There is no difference between the means of two or more groups of (ranked) dissimilarities.
# The R-statistic in ANOSIM is a ratio between within-group and between-group dissimilarities
# An R value close to "1.0" suggests dissimilarity between groups while an R value close to "0" suggests an even distribution of high and low ranks within and between groups.
# the higher the R value, the more dissimilar your groups are in terms of microbial community composition
# The P-value is the proportion of permutations that resulted in a value of R as large or larger than that calculated using the actual grouping factor

# check if rownames of metadata are in same order as the distance matrix
rownames(meta.all.scaled) %in% rownames(as.matrix(b.euc_dist))

anosim1<-anosim(x = b.euc_dist, grouping = meta.all.scaled$Site, permutations = 9999)
anosim1
# ANOSIM statistic R: 0.0847; Significance: 0.03741
# ^ sites are slightly dissimilar, result is significant...cannot accept null hypothesis
plot(anosim1)

anosim2<-anosim(x = b.euc_dist, grouping = interaction(meta.all.scaled$Site,meta.all.scaled$CollectionYear), permutations = 9999)
# Null hypothesis: There is no difference between the means of two or more groups of (ranked) dissimilarities.
# An R value close to "1.0" suggests dissimilarity between groups while an R value close to "0" suggests an even distribution of high and low ranks within and between groups.
# The P-value is the proportion of permutations that resulted in a value of R as large or larger than that calculated using the actual grouping factor

plot(anosim2)

#### Species Contributions to Dissimilarity with SIMPER ####
# NOTES: (some from here https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/simper/)
# When there are multiple sample units, there are also multiple dissimilarities (i.e., pairwise combinations) to consider – each species can contribute to the dissimilarity between each pair of sample units.
## If our purpose is to quantify the contribution of each species to the differences between two groups, we have to consider all of these dissimilarities.
# The average dissimilarity for each pairwise combination can be calculated directly via the vegan::meandist() function

meandist(dist = b.euc_dist, grouping = meta.all.scaled$Site)
# The distances shown here are a symmetric square matrix where each value is the average distance between two sample units.
# Values on the diagonal are the average distances between observations in the same group;
# all other values are the average distances between an observation in one group and an observation in another group.

# now for SIMPER analysis...
# Similarity percentage (SIMPER) partitions the dissimilarity matrix for every pair of sample units, and then calculates the average contribution of each species to the difference between the sample units.
# These contributions are relativized so that the average contributions of all species sum to 1.
# Statistical significance of these contributions is assessed by permuting the group identities.
# Note that the title is misleading: a high value for SIMPER means that a species has a high contribution to the difference between the two groups, NOT that it has high similarity between them!

simper1<-simper(b.clr,
       meta.all.scaled$Site,
       permutations = 999
)


simper1.summary<-summary(simper1, ordered = TRUE,
        digits = max(3,getOption("digits") - 3))
# create vector of comparison names from SIMPER
simper1.names<-c(names(simper1.summary))

for (i in seq_along(simper1.names)){
  print(simper1.summary[i])
}



for (i in seq_along(simper1.names)){
  print(as.data.frame(simper1.summary[i]))
}

subset.simper<-function(simper_object){
  # comp_names = list of the comparison names from SIMPER output (vegan::simper())
  ## e.g. $BDC-WI, $PD-WI, etc
  # var_subsets = list of dataframes subsetted by column$element from original dataframe;
  ## e.g. list of dataframes (each df = element of list) subsetted from metadata using vector of metadata$Site names
  simper.sum<-summary(simper_object, ordered = TRUE,
                           digits = max(3,getOption("digits") - 3)) # create table of SIMPER outputs
  comp_names<-c(names(simper.sum)) # create vector of comparison names in simper object
  for(i in seq_along(comp_names)){
    # print(comp_names[i]) -- comp_names[i] = each element in comp_names
    # print(simper_object[i]) -- simper_object[i] = each comparison in simper object
    df<-as.data.frame(simper.sum[i])
    df$ASV_ID<-rownames(df)
    df2<-merge(df,bac.ASV_tax,by="ASV_ID")
    #df2.order<-df2[order(df2[,grepl(".p",names(df2))],decreasing=FALSE)]
    #print(df)
    assign(paste0(comp_names[i],"_SIMPER.results"), df, envir = .GlobalEnv)
    assign(paste0(comp_names[i],"_SIMPER_taxaIDs"), df2, envir = .GlobalEnv)

    print(paste("Dataframe", comp_names[i] ,"done"))

  }

}

subset.simper(simper1)

# description of output of SIMPER...
# species are in descending order from highest to lowest contribution of differences between sample units
# columns are in this order: average, sd, ratio, ava/avb, cumsum, p
# ave = average contribution of species to average dissimilarity between observations from two groups
# sd = standard deviation of contribution of this species (based on its contribution to all dissimilarities betwen obsevations from two groups)
# ava/avb = average abundance of this species in each of the two groups; only included if there was a group variable used
# cusum = cumulative contribution fo this and all previous species in list; based on average but expressed as a proportion of the average dissimilarity
# p = permutation based p value aka probability of getting a larger or equal ave contribution for each species if the grouping factor was randomly permuted

head(simper1$BDC_DP$cusum)
head(sort(simper1$BDC_DP$p,decreasing=FALSE))

head(simper1$BDC_PD$cusum)
head(simper1$BDC_WI$cusum)

head(simper1$DP_PD$cusum)
head(simper1$DP_WI$cusum)

head(simper1$PD_WI$cusum)


#### Using Shapiro-Wilk test for Checking Normality of PCoA Axes ####
shapiro.test(b.pcoa.vectors$Axis.1) # what is the p-value?
# p-value = 0.5816
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(b.pcoa.vectors$Axis.1, col="blue") # with outliars

# visualize Q-Q plot for alpha div
# The Q-Q plot, or quantile-quantile plot, is a graphical tool to help us assess if a set of data plausibly came from some theoretical distribution such as a normal or exponential.
# For example, if we run a statistical analysis that assumes our residuals are normally distributed, we can use a normal Q-Q plot to check that assumption
# more on Q-Q plots here: https://data.library.virginia.edu/understanding-q-q-plots/
# more here too: https://grodri.github.io/glms/notes/c2s9#:~:text=8%20The%20Q%2DQ%20Plot,versus%20quantiles%20of%20a%20distribution.
qqnorm(b.pcoa.vectors$Axis.1, pch = 1, frame = FALSE)
qqline(b.pcoa.vectors$Axis.1, col = "red", lwd = 2)

#### Linear Models with Surface Type Frequencies & PC Axes ####
head(b.pcoa.meta)

SurfTypFreq[,3:12]

# create dfs of only surface type freq data and only the pcoa axes of interest
STFs_only<-SurfTypFreq[,3:12]
head(STFs_only)
pcoa.axes<-b.pcoa.vectors[,-(ncol(b.pcoa.vectors))]
head(pcoa.axes)

dim(STFs_only) # confirming that both data frames have the same # of rows
dim(pcoa.axes)

rownames(STFs_only) # check rownames to see if they are in the same order in both data frames
rownames(pcoa.axes)

# reorder data frames so they are in the same order by row (SampleID)
STFs_only=STFs_only[rownames(pcoa.axes),] ## reorder metadata to match order of CLR data

rownames(STFs_only) # check rownames to see if they are in the same order in both data frames after reordering
rownames(pcoa.axes)

# glm_<- vector('list', ncol(pcoa.axes) * ncol(STFs_only)) # create empty list where the GLM output is stored
# results_<- vector('list', ncol(pcoa.axes) * ncol(STFs_only)) # create an empty list where the GLM summaries are stored
# sig.results<-vector('list', ncol(pcoa.axes) * ncol(STFs_only))
# mdlnum <- 1 # counting our model numbers for indexes purposes in the loop
#
# # use a loop to run a bunch of GLMs
# ## pcoa.axes[i] is dependent variable (y), STFs_only[j] is independent variable (x) in GLM
# for (i in 1:ncol(pcoa.axes)){ # for each column in pcoa.axes
#   for (j in 1:ncol(STFs_only)){ # for each column in STFs_only
#     glm_[[mdlnum]] <-glm(pcoa.axes[,i]~STFs_only[,j], family=gaussian) # run the GLM with the gaussian distribution, where df1[i] is your dependent variable and df2[j] is your independent variable
#     results_[[mdlnum]] <-summary(glm_[[mdlnum]]) # save results of glm into list called results
#     names(results_)[mdlnum]<-paste(names(pcoa.axes)[i],"~",names(STFs_only)[j]) # rename list element to contain the name of the columns used in the model
#     mdlnum <- mdlnum + 1 # add 1 to modelnumber so we keep track of # of models (for indexing purposes in list)
#
#   }
# }
#
# ## pcoa.axes[i] is dependent variable (y), STFs_only[j] is independent variable (x) in GLM
# for (i in 1:ncol(pcoa.axes)){ # for each column in pcoa.axes
#   for (j in 1:ncol(STFs_only)){ # for each column in STFs_only
#     glm_[[mdlnum]] <-glm(pcoa.axes[,i]~STFs_only[,j], family=gaussian) # run the GLM with the gaussian distribution, where df1[i] is your dependent variable and df2[j] is your independent variable
#     results_[[mdlnum]] <-summary(glm_[[mdlnum]]) # save results of glm into list called results
#     names(results_)[mdlnum]<-paste(names(pcoa.axes)[i],"~",names(STFs_only)[j]) # rename list element to contain the name of the columns used in the model
#
#     ifelse(coef(results_[[mdlnum]])[,4] < 0.05, sig.results[[mdlnum]]<-results_[[mdlnum]], "Not Sig")
#     names(sig.results)[mdlnum]<-paste(names(pcoa.axes)[i],"~",names(STFs_only)[j])
#     mdlnum <- mdlnum + 1 # add 1 to modelnumber so we keep track of # of models (for indexing purposes in list)
#
#   }
# }
# sig.results[sapply(sig.results, is.null)] <- NULL

multi.univar.glm.fxn<-function(dep.var.df,indep.var.df,distro){
  # create empty lists to store stuff & model number (mdlnum) to keep track of models each iteration of loop in fxn
  glm_<- vector('list', ncol(dep.var.df) * ncol(indep.var.df)) # create empty list where the GLM output is stored
  results_<- vector('list', ncol(dep.var.df) * ncol(indep.var.df)) # create an empty list where the GLM summaries are stored
  sig.results<-vector('list', ncol(dep.var.df) * ncol(indep.var.df))
  mdlnum <- 1 # counting our model numbers for indexes purposes in the loop

  # run the nested loop that generates GLMs from each data frame
  ## dep.var.df[i] is dependent variable (y), indep.var.df[j] is independent variable (x) in GLM
  for (i in 1:ncol(dep.var.df)){ # for each column in dep.var.df
    for (j in 1:ncol(indep.var.df)){ # for each column in indep.var.df
      glm_[[mdlnum]] <-glm(dep.var.df[,i]~indep.var.df[,j], family=distro) # run the GLM with the gaussian distribution, where df1[i] is your dependent variable and df2[j] is your independent variable
      results_[[mdlnum]] <-summary(glm_[[mdlnum]]) # save results of glm into list called results
      names(results_)[mdlnum]<-paste(names(dep.var.df)[i],"~",names(indep.var.df)[j]) # rename list element to contain the name of the columns used in the model

      # save only significant GLMs to another list called sig.results
      ## if p-value < 0.05, save to sig.results list
      ifelse(coef(results_[[mdlnum]])[,4] < 0.05, sig.results[[mdlnum]]<-results_[[mdlnum]], "Not Sig")
      names(sig.results)[mdlnum]<-paste(names(dep.var.df)[i],"~",names(indep.var.df)[j])
      mdlnum <- mdlnum + 1 # add 1 to modelnumber so we keep track of # of models (for indexing purposes in list)

    }
  }

  # drop all NULL elements from sig.results list so it only includes significant GLMs
  sig.results[sapply(sig.results, is.null)] <- NULL

  # assign lists to global env so they are saved there are function ends
  assign("results.glms", results_,envir = .GlobalEnv)
  assign("sig.results.glms", sig.results,envir = .GlobalEnv)

}

multi.univar.glm.fxn(pcoa.axes,STFs_only,gaussian) # test the function!

#### Save Everything ####
save.image("data/SSea Dust_BetaDiv_Data.Rdata")
