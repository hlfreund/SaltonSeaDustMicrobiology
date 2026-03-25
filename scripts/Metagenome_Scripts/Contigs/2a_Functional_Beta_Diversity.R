#### Set WD & Load Libraries ####
getwd() # use setwd("path/to/files") if you are not in the right directory
#setwd("/bigdata/Metagenomesaronsonlab/shared/SaltonSea/Metagenomes/SeqCenter_3.30.2023/MGM_Analyses/Contigs")
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
  library(rgl)
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
})

#### Load Data & See Info About Data ####
#load("data/Metagenomes/SSD_mgm_analysis_data.Rdata") # load Rdata to global env
load("data/Metagenomes/SSD_MGM_Fxn_BetaDiv.Rdata")

#load("data/Metagenomes/Analysis/SSD_MGM_FxnBetaDiv.Rdata")
# save.image("data/MetagenomesSSD_mgm_FxnBetaDiv_data_Ready.Rdata)
head(dust_mgm_meta)
arsen.fxns[1:4,]

# ABOUT THE DATA:
# Before transformations (i.e., VST, CLR, etc) were done, the following was performed:
# contigs were co-assembled, so genes are found in same assembly (same set of contigs)
# non-normalized reads were mapped to genes in contigs (reads mapped by metagenome aka sample)
# featureCounts counted reads that mapped to genes on contigs
# Reads mapped to genes were divided by gene length for all genes across all samples because same KO can be assigned to multiple genes
# Gene coverage was then divided by the # of days deployed for each sample
# Scaled coverages were multiplied by 100 to not skew normalization & stats
# Scaled gene coverage was then added together for each KO ID, since multiple genes were assigned the same KO ID
# Summed coverage per KO was then transformed and/or normalized via median-ratio normalization, vst, and clr

## For pathway analyses -- after gene coverage was calculated and added together per KO ID, they were added together for each pathway
## summed coverages per KO ID, then per pathway were transformed by CLR

# Notes about which objects are which..
mgm_genes.cov.clean # has all read counts, coverages, scaled coverages, and scaled up coverages by gene
mgm_gene.cov_table[,1:4] # uses scaled up, relative coverage per genes (which were also scaled by deployment duration before being scaled up)
mgm_fxn.cov_table[,1:4] # summed scaled up relative coverages per KO
mgm.mr[1:4,1:4] # median of ratio normalized summed KO coverages

# NOTE about CLR transformation:
## uses a pseudocount of 1 to replace 0s, which is why not all 0s are treated equally
## need to look into robustCLR, which uses CLR transformation without 0s. Need more info on this methodology...

#### Calculate Mean, Scaled Relative Coverage by Gene & KO ####
colMeans(mgm_genes.cov.clean[,c(6,8:10)])

#colMeans(mgm_gene.cov_table[,-1]) # mean scaled relative coverage by gene
mean(colMeans(mgm_gene.cov_table[,-1])) # mean of all mean scaled relative gene coverages

#colMeans(mgm_fxn.cov_table[,-1]) # mean summed scaled relative coverage by KO
mean(colMeans(mgm_fxn.cov_table[,-1])) # mean of all mean summed scaled relative KO coverages

grep("SampleID", names(mgm.mr)) # find index for column we want aka "SampleID" column
#colMeans(mgm.mr[,-8778]) # mean summed scaled relative coverage by KO
mean(colMeans(mgm.mr[,-8778])) # mean of all mean summed scaled relative KO coverages

#### Clustering by Features Across Samples ####
#
# # using CLR data first
# mgm.clr[1:4,1:4]
#
# # calculate our Euclidean distance matrix using CLR data
# mgm.euc_dist1 <- dist(mgm.clr, method = "euclidean")
#
# # creating our hierarcical clustering dendrogram
# mgm.euc_clust1 <- hclust(mgm.euc_dist1, method="ward.D2")
#
# # let's make it a little nicer...
# mgm.euc_dend <- as.dendrogram(mgm.euc_clust1, hang=0.2)
# mgm.dend_cols <- as.character(meta.all.scaled$SampDate_Color[order.dendrogram(mgm.euc_dend)])
# labels_colors(mgm.euc_dend) <- mgm.dend_cols
#
# plot(mgm.euc_dend, ylab="CLR Euclidean Distance",cex = 0.5) + title(main = "Functional Trait Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 3)
# legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c( "#ef781c","#03045e","#059c3f"),pch = 15, bty = "n")
#
# dev.off()

#### Functional Beta Diversity - MRN data ####
# MR = median-ratio normalization
mgm.mr[1:4,1:4] # sample IDs are rows, genes (KOs) are columns
mgm_fxn.cov_table[1:4,1:4] # sanity check --> KOs with low coverage are still in this df
### NOTE ^ the table mgm_fxn.cov_table.no_lows will be used for transformations & normalizations that are NOT done using DESeq2
### this is because the DESeq2 object prep step allows us to drop the KOs with the lowest summed gene coverages

grep("SampleID", names(mgm.mr)) # find index for column we want
grep("SampleID", names(mgm.mr[,-8778])) # sanity check that we have the right column

# check rownames of MR & MR normalized feature count data & metadata
rownames(mgm.mr) %in% rownames(meta.all.scaled)

## PCOA with median of ratio normalized data first
# calculate our Euclidean distance matrix using VST data
mgm.euc_dist.mr <- dist(mgm.mr[,-8778], method = "euclidean")

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
##save.image("data/Metagenomesssd_mr.euc.dist1_3.7.23.Rdata")

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

#### Visualize PCoA with All Functions ####
# create PCoA ggplot fig
pcoa1<-ggplot(mgm.pcoa.mr.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(SampDate),shape=Site), size=4)+theme_bw()+
  labs(title="PCoA: Bacterial Functions in Salton Sea Dust",subtitle="Using Median-Ratio Normalized Feature Data",xlab="PC1 [50.41%]", ylab="PC2 [54.27%]",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=27),axis.title.y = element_text(size=27),legend.title.align=0.5, legend.title = element_text(size=27),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(mgm.pcoa.mr.meta$SampDate_Color[order(mgm.pcoa.mr.meta$SampDate)]),labels=c("July 2020", "August 2020", "November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  xlab("PC1 [50.41%]") + ylab("PC2 [54.27%]") + scale_shape_manual(values = c(7,10, 15,35))

ggsave(pcoa1,filename = "figures/MGM_Figs/Contigs/FxnDiv/PCoAs/MedianRatioNormalization/SSD_MGM_pcoa_MR_sampdate.png", width=12, height=10, dpi=600,create.dir = TRUE)

pcoa1a<-ggplot(mgm.pcoa.mr.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(SampDate),shape=Site,size=CollectionYear))+theme_bw()+
  labs(title="PCoA: Bacterial Functions in Salton Sea Dust",subtitle="Using Median-Ratio Normalized Feature Data",xlab="PC1 [50.41%]", ylab="PC2 [54.27%]",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=27),axis.title.y = element_text(size=27),legend.title.align=0.5, legend.title = element_text(size=27),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_size_manual(values = c("2020" = 7, "2021"=4),labels=c("2020","2021")) +
  scale_color_manual(name ="Sample Date",values=unique(mgm.pcoa.mr.meta$SampDate_Color[order(mgm.pcoa.mr.meta$SampDate)]),labels=c("July 2020", "August 2020", "November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  xlab("PC1 [50.41%]") + ylab("PC2 [54.27%]") + scale_shape_manual(values = c(7,10, 15,35))

ggsave(pcoa1a,filename = "figures/MGM_Figs/Contigs/FxnDiv/PCoAs/MedianRatioNormalization/SSD_MGM_pcoa_MR_sampdate_v2.png", width=12, height=10, dpi=600,create.dir = TRUE)

pcoa2<-ggplot(mgm.pcoa.mr.meta, aes(x=Axis.2, y=Axis.3)) +geom_point(aes(color=factor(SampDate),shape=Site,size=CollectionYear))+theme_bw()+
  labs(title="PCoA: Bacterial Functions in Salton Sea Dust",subtitle="Using Median-Ratio Normalized Feature Data",xlab="PC1 [50.41%]", ylab="PC2 [54.27%]",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=27),axis.title.y = element_text(size=27),legend.title.align=0.5, legend.title = element_text(size=27),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_size_manual(values = c("2020" = 7, "2021"=4),labels=c("2020","2021")) +
  scale_color_manual(name ="Sample Date",values=unique(mgm.pcoa.mr.meta$SampDate_Color[order(mgm.pcoa.mr.meta$SampDate)]),labels=c("July 2020", "August 2020", "November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  xlab("PC2 [54.27%]") + ylab("PC3 [6.42%]") + scale_shape_manual(values = c(7,10, 15,35))

ggsave(pcoa2,filename = "figures/MGM_Figs/Contigs/FxnDiv/PCoAs/MedianRatioNormalization/SSD_MGM_pcoa2_MR_sampdate.png", width=12, height=10, dpi=600,create.dir = TRUE)

pltly.mgm.mr<-plot_ly(mgm.pcoa.mr.meta, x=~Axis.1,y=~Axis.2,z=~Axis.3, color = ~Site, colors = c("seagreen","orange","red","purple","darkgreen","darkorange3","red4","purple4"),
                     symbol=~Site,symbols = c("square-open", "circle-open","circle","diamond")) %>%
  layout(scene = list(xaxis = list(title = 'PC1 50.41%'),
                      yaxis = list(title = 'PC2 54.27%'),
                      zaxis = list(title = 'PC3 6.42%')))

saveWidget(widget = pltly.mgm.mr, #the plotly object,
           file = "figures/MGM_Figs/Contigs/FxnDiv/PCoAs/MedianRatioNormalization/SSD_MGM_MR_pcoa_3d.html", #the path & file name
           selfcontained = TRUE #creates a single html file
)

#saveopts <- options(rgl.useNULL = TRUE)
#plot3d(x=mgm.pcoa.mr.meta$Axis.1,y=mgm.pcoa.mr.meta$Axis.2,z=mgm.pcoa.mr.meta$Axis.3,
#       xlab="PC1 50.41%",ylab="PC2 54.27%",zlab="PC3 6.42%",
#       col=mgm.pcoa.mr.meta$Site_Color, type="s",size=1)
#options(saveopts)
#movie3d(spin3d(axis=c(0,1,0),rpm=4),duration=15,dir="figures/MGM_Figs/Contigs/FxnDiv/PCoAs/MedianRatioNormalization/")

pltly.mgm.mr.noWI<-plot_ly(mgm.pcoa.mr.meta[mgm.pcoa.mr.meta$Site!="WI",], x=~Axis.1,y=~Axis.2,z=~Axis.3, color = ~Site, colors = c("seagreen","orange","red","purple","darkgreen","darkorange3","red4","purple4"),
        symbol=~Site,symbols = c("square-open", "circle-open","circle","diamond")) %>%
  layout(scene = list(xaxis = list(title = 'PC1 50.41%'),
                      yaxis = list(title = 'PC2 54.27%'),
                      zaxis = list(title = 'PC3 6.42%')))

saveWidget(widget = pltly.mgm.mr.noWI, #the plotly object,
           file = "figures/MGM_Figs/Contigs/FxnDiv/PCoAs/MedianRatioNormalization/SSD_MGM_MR_NoWI_pcoa_3d.html", #the path & file name
           selfcontained = TRUE #creates a single html file
)

#### Visualize Location & Date with PC Axes ####
# this idea was suggested by Dr. Will Porter as a way to see if we can identify compositional similiarites or dissimilarites by site and sample collection date
# using PC axes as input data, where x axis is site and y axis is the sample date in the heat map

dust.time.site<-subset(meta.all.scaled, select=c(SampleID, Site, SampDate))
mgm.pcoa.mr.dts<-merge(mgm.pcoa.mr.vectors, dust.time.site, by.x="SampleID", by.y="SampleID")

head(mgm.pcoa.mr.dts)

ggplot(mgm.pcoa.mr.dts, aes(Site, SampDate)) +
  geom_tile(aes(fill = Axis.1)) +
  geom_text(aes(label = round(Axis.1, 1))) +
  scale_fill_gradient(low = "white", high = "red")

# For heatmap color gradient, PC1
max(mgm.pcoa.mr.dts$Axis.1, na.rm=TRUE)
max(mgm.pcoa.mr.dts$Axis.1, na.rm=TRUE)/2
min(mgm.pcoa.mr.dts$Axis.1, na.rm=TRUE)

pcoa.axis1.hm<-ggplot(mgm.pcoa.mr.dts, aes(Site, SampDate, Axis.1, fill=Axis.1)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="blue", high="red") + labs(title="Salton Sea Dust Microbial PCoA Axis 1 by Sample Date & Site",subtitle="Using Median-Ratio Normalized KO Coverages ",fill="PC1 Values") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0)) + geom_text(aes(label = round(Axis.1, 2)),size=9)

ggsave(pcoa.axis1.hm,filename = "figures/MGM_Figs/Contigs/FxnDiv/PCoAs/MedianRatioNormalization/SSD_MGM_Fxn_MR_PCoA_PC1_Site_by_SampDate.png", width=18, height=27, dpi=600,create.dir = TRUE)

# For heatmap color gradient, PC2
max(mgm.pcoa.mr.dts$Axis.2, na.rm=TRUE)
max(mgm.pcoa.mr.dts$Axis.2, na.rm=TRUE)/2
min(mgm.pcoa.mr.dts$Axis.2, na.rm=TRUE)

pcoa.axis2.hm<-ggplot(mgm.pcoa.mr.dts, aes(Site, SampDate, Axis.2, fill=Axis.2)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="blue", high="red") + labs(title="Salton Sea Dust Microbial PCoA Axis 2 by Sample Date & Site",subtitle="Using Median-Ratio Normalized KO Coverages ",fill="PC2 Values") +
  geom_text(aes(label = round(Axis.2, 2)),size=9) +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(pcoa.axis2.hm,filename = "figures/MGM_Figs/Contigs/FxnDiv/PCoAs/MedianRatioNormalization/SSD_MGM_Fxn_MR_PCoA_PC2_Site_by_SampDate.png", width=18, height=27, dpi=600,create.dir = TRUE)

## Loop to Generate Heat Map for Each PC Axis

pc.plot.list<-list() # create empty list for each plot to be stored in
pc.axes<-names(mgm.pcoa.mr.dts)[grepl("Axis",names(mgm.pcoa.mr.dts))] # pull out names of columns in df that contain "Axis" in name

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
  pc.heatmap=hm.fxn(mgm.pcoa.mr.dts, mgm.pcoa.mr.dts$Site, mgm.pcoa.mr.dts$SampDate, mgm.pcoa.mr.dts[,i])
  hm_titled = pc.heatmap + ggtitle(as.character(i)) + guides(fill=guide_legend(title="PC Values"))
  pc.plot.list[[i]]=hm_titled
  ggsave(hm_titled,filename = paste("figures/MGM_Figs/Contigs/FxnDiv/PCoAs/MedianRatioNormalization/PCoA_Axes_Heatmaps/SSD_MGM_Fxn_MR_PCoA_Site_SampDate_heatmap_",i,".png",sep=""), width=18, height=27, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

}


# call each plot in list by index to view plots
pc.plot.list[[1]] # PC 1
pc.plot.list[[2]] # PC 2
pc.plot.list[[9]] # PC 9
pc.plot.list[[23]] # PC 23

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
#   #ggsave(a, file = paste0("figures/BetaDiversity/WeightedUnifrac/SSD_16S.PCoA_heatmap_", f_var,".png"), device = png, width = 15, height = 15, units = "cm")
#
#   return(a)
# }
# # loop with heatmap function to create heatmap


#### Homogeneity of Variance & PERMANOVA tests - Functions by Groups ####

## betadisper to look at homogeneity of group dispersions (aka variance) when considering multiple variables
# multivariate analogue to Levene's test of homogeneity of variances
# program finds spatial median or centroid of the group, & compare distances of group to centroid/spatial median via ANOVA

#While PERMANOVA tests differences in group means (analogous to MANOVA),
## a related test called PERMDISP can be used to evaluate homogeneity of group dispersion
#(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
## * need a distance matrix!

# check rownames of MR & MR normalized feature count data & metadata
# rownames(mgm.mr) %in% rownames(meta.all.scaled)
# mgm.euc_dist.mr came from mgm.mr (pre-Euclidean distance calculation)

# first by Site
mgm.disper1<-betadisper(mgm.euc_dist.mr, meta.all.scaled$Site)
mgm.disper1

permutest(mgm.disper1, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)

anova(mgm.disper1) # p = 0.132 --> accept the Null H, spatial medians are NOT significantly difference across sample dates

TukeyHSD(mgm.disper1) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other

#                             diff        lwr       upr     p adj

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

pnova1<-adonis2(mgm.mr[,-8778] ~ Site,data=meta.all.scaled,method = "euclidean",by="terms",permutations= 10000)
pnova1 # p-value = 0.117
p.adjust(pnova1$`Pr(>F)`,method="bonferroni",n=length(pnova1$`Pr(>F)`)) # adjusted pval
# 0.3509649

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

b.mr.dist = (vegdist(mgm.mr[,-8778], "euclidean", na.rm = TRUE)) #distance matrix using Euclidean distance matrix for trait distribution (traits of interest only)
pair.mod1<-pairwise.adonis(b.mr.dist,meta.all.scaled$Site, p.adjust.m='bonferroni',perm=9999) # shows us variation for each sample to see which ones are different
pair.mod1

# Visualize dispersions
png('figures/MGM_Figs/Contigs/FxnDiv/PCoAs/MedianRatioNormalization/SSD_MGM_pcoa_MR_betadispersion_sampledate.png',width = 700, height = 600, res=100)
plot(mgm.disper1,label=FALSE,main = "Centroids and Dispersion (Median-Ratio Data)", col=colorset6$Site_Color)
legend("bottomleft",legend = c("BDC","DP","PD","WI"),cex=.8,col = c( "#390099","#ffbd00","#eb5e28","#008000"),pch = 15, bty = "n")

dev.off()

png('figures/MGM_Figs/Contigs/FxnDiv/PCoAs/MedianRatioNormalization/SSD_MGM_boxplot_MR_centroid_distance_sampledate.png',width = 700, height = 600, res=100)
boxplot(mgm.disper1,xlab="Site", main = "Distance to Centroid by Category (Median-Ratio Data)", sub="Euclidean Distance of Median-Ratio Normalized Data", col=colorset6$Site_Color)
dev.off()

# next by site & collection year
mgm.disper2<-betadisper(mgm.euc_dist.mr, interaction(meta.all.scaled$Site,meta.all.scaled$CollectionYear))
mgm.disper2

permutest(mgm.disper2, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:

anova(mgm.disper2) # p = 0.4337 --> accept the Null H, spatial medians are NOT significantly difference across sample dates

TukeyHSD(mgm.disper2) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other

pnova2<-adonis2(mgm.mr[,-8778] ~ interaction(meta.all.scaled$Site,meta.all.scaled$CollectionYear),data=meta.all.scaled,method = "euclidean",by="terms",permutations= 10000)
pnova2 # p-value = 0.6475
#p.adjust(pnova2$`Pr(>F)`,method="bonferroni",n=length(pnova2$`Pr(>F)`)) # adjusted pval

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

b.mr.dist = (vegdist(mgm.mr[,-8778], "euclidean", na.rm = TRUE)) #distance matrix using Euclidean distance matrix for trait distribution (traits of interest only)
pair.mod2<-pairwise.adonis(b.mr.dist,interaction(meta.all.scaled$Site,meta.all.scaled$CollectionYear), p.adjust.m='bonferroni',perm=9999) # shows us variation for each sample to see which ones are different
pair.mod2

# Visualize dispersions
png('figures/MGM_Figs/Contigs/FxnDiv/PCoAs/MedianRatioNormalization/SSD_MGM_pcoa_MR_betadispersion_Site_Year.png',width = 700, height = 600, res=100)
plot(mgm.disper2,label=FALSE,main = "Centroids and Dispersion (Median-Ratio Data)", col=colorset8$SampDate_Color)
dev.off()

png('figures/MGM_Figs/Contigs/FxnDiv/PCoAs/MedianRatioNormalization/SSD_MGM_boxplot_MR_centroid_distance_Site_Year.png',width = 700, height = 600, res=100)
boxplot(mgm.disper2,xlab="Collection Date + Site", main = "Distance to Centroid by Category (Median-Ratio Data)", sub="Euclidean Distance of Median-Ratio Normalized Data", col=colorset6$Site_Color)
dev.off()

#### Pull Out LPS Metabolic Fxns from MR data ####
## heatmaps of traits of interest
mgm.mr.lps<-merge(lps.fxns, mgm_mr_melt, by=c("KO_ID"),allow.cartesian = TRUE)
head(mgm.mr.lps)

mgm.mr[1:4,1:4]

# pull out LPS functions from MR normalized, summed coverages (summed gene coverage per KO)
lps.ko<-mgm.mr[,which(colnames(mgm.mr) %in% lps.fxns$KO_ID)] # merge MR data w/ S fxns found in contigs from KOFamScan
lps.ko$SampleID<-rownames(lps.ko)
lps.ko.melt<-melt(lps.ko, by="SampleID")
colnames(lps.ko.melt)[which(names(lps.ko.melt) == "variable")] <- "KO_ID"
colnames(lps.ko.melt)[which(names(lps.ko.melt) == "value")] <- "MR_SumCovPerKO"
head(lps.ko.melt) #sanity check

mr.lps.ko<-merge(lps.ko.melt,lps.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(mr.lps.ko)
colnames(mr.lps.ko)[which(names(mr.lps.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
mr.cov.sum.lps.ko<-as.data.frame(dcast(mr.lps.ko, SampleID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(mr.cov.sum.lps.ko)<-mr.cov.sum.lps.ko$SampleID
mr.cov.sum.lps.ko[1:4,]

# sanity check
mr.cov.sum.lps.ko$`lpxL, htrB; Kdo2-lipid IVA lauroyltransferase/acyltransferase [EC:2.3.1.241 2.3.1.-]`[1:4]
head(mr.lps.ko)

#### LPS Heat Maps ####
# see max & mean of summed
max(mr.cov.sum.lps.ko[,-1])
min(mr.cov.sum.lps.ko[,-1])
max(mr.cov.sum.lps.ko[,-1])/2

# first heat map of sulfur KOs
#heatmap(as.matrix(mr.cov.sum.lps.ko[,-1]), scale = "none")

colSums(mr.cov.sum.lps.ko[,-1])
#mr.cov.sum.lps.ko2 <- mr.cov.sum.lps.ko[,which(colSums(mr.cov.sum.lps.ko[,-1])>10)]

#heatmap(as.matrix(mr.cov.sum.lps.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap -- using MR data that includes NAs so they are blocked out on heatmaps
mr.lps.ko[1:4,]
mr.lps.all<-merge(mr.lps.ko,meta.all.scaled,by="SampleID")
head(mr.lps.all)
mr.lps.all$PlotID<-mr.lps.all$SampleID
mr.lps.all$PlotID = factor(mr.lps.all$PlotID, levels=unique(mr.lps.all$PlotID[order(mr.lps.all$Site,mr.lps.all$SampDate)]), ordered=TRUE)

unique(mr.lps.all$LPS_Structure)
# mr.lps.all$LPS_Structure<-mr.lps.all$LPS_Structure
# mr.lps.all$LPS_Structure[(mr.lps.all$LPS_Structure) == "KDO2-lipid A biosynthesis, Raetz LPS_Structure"] <- "Lipid A"
# mr.lps.all$LPS_Structure[(mr.lps.all$LPS_Structure) == "CMP-KDO biosynthesis"] <- "CMP-KDO"
# mr.lps.all$LPS_Structure[(mr.lps.all$LPS_Structure) == "KDO2-lipid A biosynthesis, Raetz LPS_Structure, non-LpxL-LpxM type"] <- "Lipid A, non-LpxL-LpxM"
# mr.lps.all$LPS_Structure[(mr.lps.all$LPS_Structure) == "ADP-L-glycero-D-manno-heptose biosynthesis"] <- "ALgDmh"
# mr.lps.all$LPS_Structure[(mr.lps.all$LPS_Structure) == "KDO2-lipid A modification LPS_Structure"] <- "Lipid A Mod"
#
# mr.lps.all$LPS_Structure<-factor(mr.lps.all$LPS_Structure,levels=c("KDO2-lipid A biosynthesis, Raetz LPS_Structure","CMP-KDO biosynthesis","KDO2-lipid A biosynthesis, Raetz LPS_Structure, non-LpxL-LpxM type","ADP-L-glycero-D-manno-heptose biosynthesis","KDO2-lipid A modification LPS_Structure"))
# mr.lps.all$LPS_Structure<-factor(mr.lps.all$LPS_Structure,levels=c("Lipid A","CMP-KDO","Lipid A, non-LpxL-LpxM","ALgDmh","Lipid A Mod"))

#mr.lps.all$LPS_StructureSpecific<-factor(mr.lps.all$LPS_StructureSpecific,levels=c("Assimilatory Sulfate Reduction","Dissimilatory Sulfate Redox","Multiple LPS_Structures","SOX","LPS Disproportionation","Sulfide Oxidation","Sulfite Oxidation","Thiosulfate Oxidation"))
#mr.lps.all$PathSpecShort<-factor(mr.lps.all$PathSpecShort,levels=c("A.SO4 Red","D.SO4 RedOx","Multi","Dispro","H2S Ox","SO3 Ox","S2O3 Ox"))

#mr.lps.all$KO_Function.KEGG = factor(mr.lps.all$KO_Function.KEGG, levels=unique(mr.lps.all$KO_Function.KEGG[order(mr.lps.all$LPS_Structure)]), ordered=TRUE)

head(mr.lps.all)

# convert all 0s to NAs so they appear gray on the heatmap
mr.lps.all$MR_SumCovPerKO_NA<-ifelse(mr.lps.all$MR_SumCovPerKO==0,NA,mr.lps.all$MR_SumCovPerKO)
head(mr.lps.all)

# For heatmap color gradient
max(mr.lps.all$MR_SumCovPerKO, na.rm=TRUE)
max(mr.lps.all$MR_SumCovPerKO, na.rm=TRUE)/2
min(mr.lps.all$MR_SumCovPerKO, na.rm=TRUE)

# Find mean coverage of fxns
LPS.mgm.means<-data.frame(MeanScaledCov=colMeans(mr.cov.sum.lps.ko[,-1]))
LPS.mgm.means$KO_ID<-rownames(LPS.mgm.means)
LPS.mgm.means$Category<-"LPS"
LPS.mgm.means

# Figures below
# by SampleID

lps.hm1a<-ggplot(mr.lps.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="LPS Biosynthesis in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(lps.hm1a,filename = "figures/MGM_Figs/Contigs/FxnDiv/LPS_Biosynthesis/LPS_KOFxns_MGMs_SampID_by_Function_heatmap.png", width=35, height=20, dpi=600,create.dir = TRUE)

lps.hm1a2<-ggplot(mr.lps.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="LPS Biosynthesis in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(lps.hm1a2,filename = "figures/MGM_Figs/Contigs/FxnDiv/LPS_Biosynthesis/LPS_KOFxns_MGMs_SampID_by_Function_Site_heatmap.png", width=35, height=20, dpi=600,create.dir = TRUE)

lps.hm1a3<-ggplot(mr.lps.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="LPS Biosynthesis in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(LPS_Structure~Site,scales="free", space = "free")

ggsave(lps.hm1a3,filename = "figures/MGM_Figs/Contigs/FxnDiv/LPS_Biosynthesis/LPS_KOFxns_MGMs_SampleID_by_Function_LPS_Structure_Site_heatmap.png", width=35, height=20, dpi=600,create.dir = TRUE)

lps.hm1a3b<-ggplot(mr.lps.all[mr.lps.all$MR_SumCovPerKO>30,], aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) + geom_text(label=round(mr.lps.all$MR_SumCovPerKO[mr.lps.all$MR_SumCovPerKO>30],2)) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="LPS Biosynthesis in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(LPS_Structure~Site,scales="free", space = "free")

ggsave(lps.hm1a3b,filename = "figures/MGM_Figs/Contigs/FxnDiv/LPS_Biosynthesis/LPS_KOFxns_MGMs_SampleID_by_Function_LPS_Structure_Site_HighCov_labeled_heatmap.png", width=35, height=20, dpi=600,create.dir = TRUE)

# lps.hm1a4<-ggplot(mr.lps.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.15) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="LPS Biosynthesis in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(LPS_Structure~.,scales="free_y", space = "free")
#
# ggsave(lps.hm1a4,filename = "figures/MGM_Figs/Contigs/FxnDiv/LPS_Biosynthesis/LPS_KOFxns_MGMs_SampID_by_Function_LPS_Structure_heatmap2.png", width=17, height=15, dpi=600,create.dir = TRUE)

lps.hm1a5<-ggplot(mr.lps.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="LPS Biosynthesis in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(LPS_Structure~SampDate, scales="free", space = "free")

ggsave(lps.hm1a5,filename = "figures/MGM_Figs/Contigs/FxnDiv/LPS_Biosynthesis/LPS_KOFxns_MGMs_SampleID_by_Function_LPS_Structure_SampDate_heatmap.png", width=35, height=20, dpi=600,create.dir = TRUE)

lps.hm1a5a<-ggplot(mr.lps.all[mr.lps.all$MR_SumCovPerKO>30,], aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) + geom_text(label=round(mr.lps.all$MR_SumCovPerKO[mr.lps.all$MR_SumCovPerKO>30],2)) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="LPS Biosynthesis in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Summed Coverages >= 30)",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(LPS_Structure~SampDate, scales="free", space = "free")

ggsave(lps.hm1a5a,filename = "figures/MGM_Figs/Contigs/FxnDiv/LPS_Biosynthesis/LPS_KOFxns_MGMs_SampleID_by_Function_LPS_Structure_SampDate_HighCov_labeled_heatmap.png", width=38, height=20, dpi=600,create.dir = TRUE)


# by LPS structure

lps.hm2a<-ggplot(mr.lps.all[mr.lps.all$LPS_Structure=="O-antigen Repeat Unit",], aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="O-Antigen (LPS) Biosynthesis Genes in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(lps.hm2a,filename = "figures/MGM_Figs/Contigs/FxnDiv/LPS_Biosynthesis/LPS_KOFxns_MGMs_SampID_by_Function_OAntigen_Only_Site_heatmap.png", width=35, height=15, dpi=600,create.dir = TRUE)

lps.hm2a1<-ggplot(mr.lps.all[mr.lps.all$LPS_Structure=="Lipid A",], aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="O-Antigen (LPS) Biosynthesis Genes in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(lps.hm2a1,filename = "figures/MGM_Figs/Contigs/FxnDiv/LPS_Biosynthesis/LPS_KOFxns_MGMs_SampID_by_Function_LipidA_Only_Site_heatmap.png", width=35, height=15, dpi=600,create.dir = TRUE)

lps.hm2a2<-ggplot(mr.lps.all[mr.lps.all$LPS_Structure=="Core Region",], aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="LPS Core Region Biosynthesis Genes in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(lps.hm2a2,filename = "figures/MGM_Figs/Contigs/FxnDiv/LPS_Biosynthesis/LPS_KOFxns_MGMs_SampID_by_Function_CoreRegion_Only_Site_heatmap.png", width=35, height=15, dpi=600,create.dir = TRUE)

#### Pull Out Quorum Sensing Metabolic Fxns from MR data ####
## heatmaps of traits of interest
mgm.mr.quorsens<-merge(QuorSens.fxns, mgm_mr_melt, by=c("KO_ID"),allow.cartesian = TRUE)
head(mgm.mr.quorsens)

mgm.mr[1:4,1:4]

# pull out QuorumSensing functions from MR normalized, summed coverages (summed gene coverage per KO)
quorsens.ko<-mgm.mr[,which(colnames(mgm.mr) %in% QuorSens.fxns$KO_ID)] # merge MR data w/ S fxns found in contigs from KOFamScan
quorsens.ko$SampleID<-rownames(quorsens.ko)
quorsens.ko.melt<-melt(quorsens.ko, by="SampleID")
colnames(quorsens.ko.melt)[which(names(quorsens.ko.melt) == "variable")] <- "KO_ID"
colnames(quorsens.ko.melt)[which(names(quorsens.ko.melt) == "value")] <- "MR_SumCovPerKO"
head(quorsens.ko.melt) #sanity check

mr.quorsens.ko<-merge(quorsens.ko.melt,QuorSens.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(mr.quorsens.ko)
colnames(mr.quorsens.ko)[which(names(mr.quorsens.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
mr.cov.sum.quorsens.ko<-as.data.frame(dcast(mr.quorsens.ko, SampleID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(mr.cov.sum.quorsens.ko)<-mr.cov.sum.quorsens.ko$SampleID
mr.cov.sum.quorsens.ko[1:4,]

# sanity check
mr.cov.sum.quorsens.ko$`secA; preprotein translocase subunit SecA [EC:7.4.2.8]`[1:4]
head(mr.quorsens.ko)

#### Quorum Sensing Heat Maps ####
# see max & mean of summed
max(mr.cov.sum.quorsens.ko[,-1])
min(mr.cov.sum.quorsens.ko[,-1])
max(mr.cov.sum.quorsens.ko[,-1])/2

# first heat map of sulfur KOs
#heatmap(as.matrix(mr.cov.sum.quorsens.ko[,-1]), scale = "none")

colSums(mr.cov.sum.quorsens.ko[,-1])
#mr.cov.sum.quorsens.ko2 <- mr.cov.sum.quorsens.ko[,which(colSums(mr.cov.sum.quorsens.ko[,-1])>10)]

#heatmap(as.matrix(mr.cov.sum.quorsens.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap -- using MR data that includes NAs so they are blocked out on heatmaps
mr.quorsens.ko[1:4,]
mr.quorsens.all<-merge(mr.quorsens.ko,meta.all.scaled,by="SampleID")
head(mr.quorsens.all)
mr.quorsens.all$PlotID<-mr.quorsens.all$SampleID
mr.quorsens.all$PlotID = factor(mr.quorsens.all$PlotID, levels=unique(mr.quorsens.all$PlotID[order(mr.quorsens.all$Site,mr.quorsens.all$SampDate)]), ordered=TRUE)

#mr.quorsens.all$KO_Function.KEGG = factor(mr.quorsens.all$KO_Function.KEGG, levels=unique(mr.quorsens.all$KO_Function.KEGG[order(mr.quorsens.all$Pathway)]), ordered=TRUE)

head(mr.quorsens.all)

# convert all 0s to NAs so they appear gray on the heatmap
mr.quorsens.all$MR_SumCovPerKO_NA<-ifelse(mr.quorsens.all$MR_SumCovPerKO==0,NA,mr.quorsens.all$MR_SumCovPerKO)
head(mr.quorsens.all)

# For heatmap color gradient
max(mr.quorsens.all$MR_SumCovPerKO, na.rm=TRUE)
max(mr.quorsens.all$MR_SumCovPerKO, na.rm=TRUE)/2
min(mr.quorsens.all$MR_SumCovPerKO, na.rm=TRUE)

# Find mean coverage of fxns
QS.mgm.means<-data.frame(MeanScaledCov=colMeans(mr.cov.sum.quorsens.ko[,-1]))
QS.mgm.means$KO_ID<-rownames(QS.mgm.means)
QS.mgm.means$Category<-"Quorum Sensing"
QS.mgm.means

# Figures below
# by SampleID

quorsens.hm1a<-ggplot(mr.quorsens.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Quorum Sensing in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(quorsens.hm1a,filename = "figures/MGM_Figs/Contigs/FxnDiv/QuorumSensing/QuorumSensing_KOFxns_MGMs_SampID_by_Function_heatmap.png", width=28, height=27, dpi=600,create.dir = TRUE)

quorsens.hm1a2<-ggplot(mr.quorsens.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Quorum Sensing in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(quorsens.hm1a2,filename = "figures/MGM_Figs/Contigs/FxnDiv/QuorumSensing/QuorumSensing_KOFxns_MGMs_SampID_by_Function_Site_heatmap.png", width=28, height=15, dpi=600,create.dir = TRUE)

quorsens.hm1a3<-ggplot(mr.quorsens.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Quorum Sensing in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Cell.Type~Site,scales="free", space = "free")

ggsave(quorsens.hm1a3,filename = "figures/MGM_Figs/Contigs/FxnDiv/QuorumSensing/QuorumSensing_KOFxns_MGMs_SampID_by_Function_CellType_Site_heatmap.png", width=28, height=15, dpi=600,create.dir = TRUE)

quorsens.hm1a3a<-ggplot(mr.quorsens.all[mr.quorsens.all$MR_SumCovPerKO>30,], aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) + geom_text(label=round(mr.quorsens.all$MR_SumCovPerKO[mr.quorsens.all$MR_SumCovPerKO>30],2)) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Quorum Sensing in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Summed Coverages >= 30)",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Cell.Type~Site, scales="free", space = "free")

ggsave(quorsens.hm1a3a,filename = "figures/MGM_Figs/Contigs/FxnDiv/QuorumSensing/QuorSens_KOFxns_MGMs_SampleID_by_Function_CellType_Site_HighCov_labeled_heatmap.png", width=28, height=20, dpi=600,create.dir = TRUE)


quorsens.hm1a4<-ggplot(mr.quorsens.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) + geom_text(label=round(mr.quorsens.all$MR_SumCovPerKO,2)) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Quorum Sensing in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Cell.Type~SampDate,scales="free", space = "free")

ggsave(quorsens.hm1a4,filename = "figures/MGM_Figs/Contigs/FxnDiv/QuorumSensing/QuorumSensing_KOFxns_MGMs_SampID_by_Function_CellType_SampDate_heatmap.png", width=28, height=15, dpi=600,create.dir = TRUE)

quorsens.hm1a4a<-ggplot(mr.quorsens.all[mr.quorsens.all$MR_SumCovPerKO>30,], aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) + geom_text(label=round(mr.quorsens.all$MR_SumCovPerKO[mr.quorsens.all$MR_SumCovPerKO>30],2)) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Quorum Sensing in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Summed Coverages >= 30)",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Cell.Type~SampDate, scales="free", space = "free")

ggsave(quorsens.hm1a4a,filename = "figures/MGM_Figs/Contigs/FxnDiv/QuorumSensing/QuorSens_KOFxns_MGMs_SampleID_by_Function_CellType_SampDate_HighCov_labeled_heatmap.png", width=28, height=20, dpi=600,create.dir = TRUE)

#### Pull Out Sporulation Fxns from MR data ####
## heatmaps of traits of interest
# mgm.mr.spor<-merge(spor.fxns, mgm_mr_melt, by=c("KO_ID"),allow.cartesian = TRUE)
# head(mgm.mr.spor)

mgm.mr[1:4,1:4]

# pull out Sporulation functions from MR normalized, summed coverages (summed gene coverage per KO)
spor.ko<-mgm.mr[,which(colnames(mgm.mr) %in% spor.fxns$KO_ID)] # merge MR data w/ S fxns found in contigs from KOFamScan
spor.ko$SampleID<-rownames(spor.ko)
spor.ko.melt<-melt(spor.ko, by="SampleID")
colnames(spor.ko.melt)[which(names(spor.ko.melt) == "variable")] <- "KO_ID"
colnames(spor.ko.melt)[which(names(spor.ko.melt) == "value")] <- "MR_SumCovPerKO"
head(spor.ko.melt) #sanity check

mr.spor.ko<-merge(spor.ko.melt,spor.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(mr.spor.ko)
colnames(mr.spor.ko)[which(names(mr.spor.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
mr.cov.sum.spor.ko<-as.data.frame(dcast(mr.spor.ko, SampleID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(mr.cov.sum.spor.ko)<-mr.cov.sum.spor.ko$SampleID
mr.cov.sum.spor.ko[1:4,]

# sanity check
mr.cov.sum.spor.ko$`spoVS; stage V sporulation protein S`[1:4]
head(mr.spor.ko)

#### Sporulation Heat Maps ####
# see max & mean of summed
max(mr.cov.sum.spor.ko[,-1])
min(mr.cov.sum.spor.ko[,-1])
max(mr.cov.sum.spor.ko[,-1])/2

# first heat map of sulfur KOs
#heatmap(as.matrix(mr.cov.sum.spor.ko[,-1]), scale = "none")

colSums(mr.cov.sum.spor.ko[,-1])
#mr.cov.sum.spor.ko2 <- mr.cov.sum.spor.ko[,which(colSums(mr.cov.sum.spor.ko[,-1])>10)]

#heatmap(as.matrix(mr.cov.sum.spor.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap -- using MR data that includes NAs so they are blocked out on heatmaps
mr.spor.ko[1:4,]
mr.spor.all<-merge(mr.spor.ko,meta.all.scaled,by="SampleID")
head(mr.spor.all)
mr.spor.all$PlotID<-mr.spor.all$SampleID
mr.spor.all$PlotID = factor(mr.spor.all$PlotID, levels=unique(mr.spor.all$PlotID[order(mr.spor.all$Site,mr.spor.all$SampDate)]), ordered=TRUE)

head(mr.spor.all)

# convert all 0s to NAs so they appear gray on the heatmap
mr.spor.all$MR_SumCovPerKO_NA<-ifelse(mr.spor.all$MR_SumCovPerKO==0,NA,mr.spor.all$MR_SumCovPerKO)
head(mr.spor.all)

# For heatmap color gradient
max(mr.spor.all$MR_SumCovPerKO, na.rm=TRUE)
max(mr.spor.all$MR_SumCovPerKO, na.rm=TRUE)/2
min(mr.spor.all$MR_SumCovPerKO, na.rm=TRUE)

# Find mean coverage of fxns
Spor.mgm.means<-data.frame(MeanScaledCov=colMeans(mr.cov.sum.spor.ko[,-1]))
Spor.mgm.means$KO_ID<-rownames(Spor.mgm.means)
Spor.mgm.means$Category<-"Sporulation"
Spor.mgm.means

# Figures below
# by SampleID

spor.hm1a<-ggplot(mr.spor.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sporulation Functions in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(spor.hm1a,filename = "figures/MGM_Figs/Contigs/FxnDiv/Sporulation/Sporulation_KOFxns_MGMs_SampID_by_Function_heatmap.png", width=18, height=27, dpi=600,create.dir = TRUE)

spor.hm1a2<-ggplot(mr.spor.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sporulation Functions in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(spor.hm1a2,filename = "figures/MGM_Figs/Contigs/FxnDiv/Sporulation/Sporulation_KOFxns_MGMs_SampID_by_Function_Site_heatmap.png", width=18, height=15, dpi=600,create.dir = TRUE)

spor.hm1a2a<-ggplot(mr.spor.all[mr.spor.all$MR_SumCovPerKO>10,], aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) + geom_text(label=round(mr.spor.all$MR_SumCovPerKO[mr.spor.all$MR_SumCovPerKO>10],2)) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sporulation Functions in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Summed Coverages >= 10)",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site, scales="free", space = "free")

ggsave(spor.hm1a2a,filename = "figures/MGM_Figs/Contigs/FxnDiv/Sporulation/Sporulation_KOFxns_MGMs_SampleID_by_Function_Site_HighCov_labeled_heatmap.png", width=28, height=20, dpi=600,create.dir = TRUE)

#### Pull Out Temp Shock Metabolic Fxns from MRN data ####
## heatmaps of traits of interest

mgm.mr[1:4,1:4]

# pull out heat shock functions from MR normalized, summed coverages (summed gene coverage per KO)
tempshock.ko<-mgm.mr[,which(colnames(mgm.mr) %in% tempshock.fxns$KO_ID)] # merge MRN data w/ S fxns found in contigs from KOFamScan
tempshock.ko$SampleID<-rownames(tempshock.ko)
tempshock.ko.melt<-melt(tempshock.ko, by="SampleID")
colnames(tempshock.ko.melt)[which(names(tempshock.ko.melt) == "variable")] <- "KO_ID"
colnames(tempshock.ko.melt)[which(names(tempshock.ko.melt) == "value")] <- "MR_SumCovPerKO"
head(tempshock.ko.melt) #sanity check

mr.tempshock.ko<-merge(tempshock.ko.melt,tempshock.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(mr.tempshock.ko)
colnames(mr.tempshock.ko)[which(names(mr.tempshock.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
mr.cov.sum.tempshock.ko<-as.data.frame(dcast(mr.tempshock.ko, SampleID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(mr.cov.sum.tempshock.ko)<-mr.cov.sum.tempshock.ko$SampleID
mr.cov.sum.tempshock.ko[1:4,]

# sanity check
mr.cov.sum.tempshock.ko$`hslR; ribosome-associated heat shock protein Hsp15`[1:4]
head(mr.tempshock.ko)


#### Temp Shock KO Heat Maps ####
# see max & mean of summed
max(mr.cov.sum.tempshock.ko[,-1])
min(mr.cov.sum.tempshock.ko[,-1])
max(mr.cov.sum.tempshock.ko[,-1])/2

# first heat map of tempshock KOs
#heatmap(as.matrix(mr.cov.sum.tempshock.ko[,-1]), scale = "none")

colSums(mr.cov.sum.tempshock.ko[,-1])
#mr.cov.sum.tempshock.ko2 <- mr.cov.sum.tempshock.ko[,which(colSums(mr.cov.sum.tempshock.ko[,-1])>10)]

##heatmap(as.matrix(mr.cov.sum.tempshock.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap -- using MR data that includes NAs so they are blocked out on heatmaps
mr.tempshock.ko[1:4,]
mr.tempshock.all<-merge(mr.tempshock.ko,meta.all.scaled,by="SampleID")
head(mr.tempshock.all)
mr.tempshock.all$PlotID<-mr.tempshock.all$SampleID
mr.tempshock.all$PlotID = factor(mr.tempshock.all$PlotID, levels=unique(mr.tempshock.all$PlotID[order(mr.tempshock.all$Site,mr.tempshock.all$SampDate)]), ordered=TRUE)

head(mr.tempshock.all)

# convert all 0s to NAs so they appear gray on the heatmap
mr.tempshock.all$MR_SumCovPerKO_NA<-ifelse(mr.tempshock.all$MR_SumCovPerKO==0,NA,mr.tempshock.all$MR_SumCovPerKO)
head(mr.tempshock.all)

# For heatmap color gradient
max(mr.tempshock.all$MR_SumCovPerKO, na.rm=TRUE)
max(mr.tempshock.all$MR_SumCovPerKO, na.rm=TRUE)/2
min(mr.tempshock.all$MR_SumCovPerKO, na.rm=TRUE)

# Find mean coverage of fxns
TempShock.mgm.means<-data.frame(MeanScaledCov=colMeans(mr.cov.sum.tempshock.ko[,-1]))
TempShock.mgm.means$KO_ID<-rownames(TempShock.mgm.means)
TempShock.mgm.means$Category<-"Temperature Shock"
TempShock.mgm.means

# Figures below
# by SampleID

tempshock.hm1a<-ggplot(mr.tempshock.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Temperature Shock Proteins in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(tempshock.hm1a,filename = "figures/MGM_Figs/Contigs/FxnDiv/TempShock/TempShock_KOFxns_MGMs_SampID_by_Function_heatmap.png", width=25, height=27, dpi=600,create.dir = TRUE)

tempshock.hm1a2<-ggplot(mr.tempshock.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Temperature Shock Proteins in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(tempshock.hm1a2,filename = "figures/MGM_Figs/Contigs/FxnDiv/TempShock/TempShock_KOFxns_MGMs_SampID_by_Function_Site_heatmap.png", width=25, height=27, dpi=600,create.dir = TRUE)

tempshock.hm1a2a<-ggplot(mr.tempshock.all[mr.tempshock.all$MR_SumCovPerKO>20,], aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) + geom_text(label=round(mr.tempshock.all$MR_SumCovPerKO[mr.tempshock.all$MR_SumCovPerKO>20],2)) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Temperature Shock Proteins in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Summed Coverages >= 20)",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site, scales="free", space = "free")

ggsave(tempshock.hm1a2a,filename = "figures/MGM_Figs/Contigs/FxnDiv/TempShock/TempShock_KOFxns_MGMs_SampleID_by_Function_Site_HighCov_labeled_heatmap.png", width=28, height=20, dpi=600,create.dir = TRUE)

tempshock.hm1a3<-ggplot(mr.tempshock.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Temperature Shock Proteins in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate,scales="free_x", space = "free")

ggsave(tempshock.hm1a3,filename = "figures/MGM_Figs/Contigs/FxnDiv/TempShock/TempShock_KOFxns_MGMs_SampID_by_Function_SampDate_heatmap.png", width=28, height=27, dpi=600,create.dir = TRUE)

#### Pull Out Aerotaxis Fxns from MRN data ####
## heatmaps of traits of interest

mgm.mr[1:4,1:4]

# pull out heat shock functions from MR normalized, summed coverages (summed gene coverage per KO)
aerotaxis.ko<-data.frame(K03776=mgm.mr[,which(colnames(mgm.mr) %in% "K03776")],SampleID=rownames(mgm.mr)) # merge MRN data w/ S fxns found in contigs from KOFamScan
aerotaxis.ko.melt<-melt(aerotaxis.ko, by="SampleID")
colnames(aerotaxis.ko.melt)[which(names(aerotaxis.ko.melt) == "variable")] <- "KO_ID"
colnames(aerotaxis.ko.melt)[which(names(aerotaxis.ko.melt) == "value")] <- "MR_SumCovPerKO"
head(aerotaxis.ko.melt) #sanity check

#mr.aerotaxis.ko<-merge(aerotaxis.ko.melt,heatshock.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
#colnames(mr.aerotaxis.ko)[which(names(mr.aerotaxis.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
aerotaxis.ko.melt$KO_Function.KEGG<-"aer; aerotaxis receptor"
mr.cov.sum.aerotaxis.ko<-as.data.frame(dcast(aerotaxis.ko.melt, SampleID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(mr.cov.sum.aerotaxis.ko)<-mr.cov.sum.aerotaxis.ko$SampleID
mr.cov.sum.aerotaxis.ko[1:4,]

# sanity check
mr.cov.sum.aerotaxis.ko$`aer; aerotaxis receptor`[1:4]
head(mr.cov.sum.aerotaxis.ko)

#### Aerotaxis KO Heat Maps ####
# see max & mean of summed
max(mr.cov.sum.aerotaxis.ko[,-1])
min(mr.cov.sum.aerotaxis.ko[,-1])
max(mr.cov.sum.aerotaxis.ko[,-1])/2

# first heat map of sulfur KOs
# #heatmap(as.matrix(mr.cov.sum.aerotaxis.ko[,-1]), scale = "none")
#
# colSums(mr.cov.sum.aerotaxis.ko[,-1])
# #mr.cov.sum.aerotaxis.ko2 <- mr.cov.sum.aerotaxis.ko[,which(colSums(mr.cov.sum.aerotaxis.ko[,-1])>10)]

#heatmap(as.matrix(mr.cov.sum.aerotaxis.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap -- using MR data that includes NAs so they are blocked out on heatmaps
aerotaxis.ko.melt[1:4,]
mr.aerotaxis.all<-merge(aerotaxis.ko.melt,meta.all.scaled,by="SampleID")
head(mr.aerotaxis.all)
mr.aerotaxis.all$PlotID<-mr.aerotaxis.all$SampleID
mr.aerotaxis.all$PlotID = factor(mr.aerotaxis.all$PlotID, levels=unique(mr.aerotaxis.all$PlotID[order(mr.aerotaxis.all$Site,mr.aerotaxis.all$SampDate)]), ordered=TRUE)

head(mr.aerotaxis.all)

# convert all 0s to NAs so they appear gray on the heatmap
mr.aerotaxis.all$MR_SumCovPerKO_NA<-ifelse(mr.aerotaxis.all$MR_SumCovPerKO==0,NA,mr.aerotaxis.all$MR_SumCovPerKO)

# For heatmap color gradient
max(aerotaxis.ko.melt$MR_SumCovPerKO, na.rm=TRUE)
max(aerotaxis.ko.melt$MR_SumCovPerKO, na.rm=TRUE)/2
min(aerotaxis.ko.melt$MR_SumCovPerKO, na.rm=TRUE)

# Find mean coverage of fxns
mean(mr.cov.sum.aerotaxis.ko$`aer; aerotaxis receptor`)

# Figures below
# by SampleID

aerotaxis.hm1a<-ggplot(mr.aerotaxis.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Aerotaxis Receptor in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(aerotaxis.hm1a,filename = "figures/MGM_Figs/Contigs/FxnDiv/Aerotaxis/Aerotaxis_KOFxns_MGMs_SampID_by_Function_heatmap.png", width=18, height=27, dpi=600,create.dir = TRUE)

aerotaxis.hm1a2<-ggplot(mr.aerotaxis.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Aerotaxis Receptors in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(aerotaxis.hm1a2,filename = "figures/MGM_Figs/Contigs/FxnDiv/Aerotaxis/Aerotaxis_KOFxns_MGMs_SampID_by_Function_Site_heatmap.png", width=17, height=15, dpi=600,create.dir = TRUE)

#### Pull Out UV-Damaged DNA Repair Fxns from MRN data ####
## heatmaps of traits of interest

mgm.mr[1:4,1:4]

# pull out heat shock functions from MR normalized, summed coverages (summed gene coverage per KO)
UVrep.ko<-mgm.mr[,which(colnames(mgm.mr) %in% UVrep.fxns$KO_ID)] # merge MRN data w/ S fxns found in contigs from KOFamScan
UVrep.ko$SampleID<-rownames(UVrep.ko)
UVrep.ko.melt<-melt(UVrep.ko, by="SampleID")
colnames(UVrep.ko.melt)[which(names(UVrep.ko.melt) == "variable")] <- "KO_ID"
colnames(UVrep.ko.melt)[which(names(UVrep.ko.melt) == "value")] <- "MR_SumCovPerKO"
head(UVrep.ko.melt) #sanity check

mr.UVrep.ko<-merge(UVrep.ko.melt,UV.DNA.rep.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(mr.UVrep.ko)
colnames(mr.UVrep.ko)[which(names(mr.UVrep.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
mr.cov.sum.UVrep.ko<-as.data.frame(dcast(mr.UVrep.ko, SampleID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(mr.cov.sum.UVrep.ko)<-mr.cov.sum.UVrep.ko$SampleID
mr.cov.sum.UVrep.ko[1:4,]

# sanity check
mr.cov.sum.UVrep.ko$`uvrC; excinuclease ABC subunit C`[1:4]
head(mr.UVrep.ko)

#### UV-Damage DNA Repair KO Heat Maps ####
# see max & mean of summed
max(mr.cov.sum.UVrep.ko[,-1])
min(mr.cov.sum.UVrep.ko[,-1])
max(mr.cov.sum.UVrep.ko[,-1])/2

# first heat map of UVrep KOs
##heatmap(as.matrix(mr.cov.sum.UVrep.ko[,-1]), scale = "none")

colSums(mr.cov.sum.UVrep.ko[,-1])
#mr.cov.sum.UVrep.ko2 <- mr.cov.sum.UVrep.ko[,which(colSums(mr.cov.sum.UVrep.ko[,-1])>10)]

##heatmap(as.matrix(mr.cov.sum.UVrep.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap -- using MR data that includes NAs so they are blocked out on heatmaps
mr.UVrep.ko[1:4,]
mr.UVrep.all<-merge(mr.UVrep.ko,meta.all.scaled,by="SampleID")
head(mr.UVrep.all)
#mr.UVrep.all$PlotID<-mr.UVrep.all$SampleID
mr.UVrep.all$PlotID = factor(mr.UVrep.all$SampleID, levels=unique(mr.UVrep.all$SampleID[order(mr.UVrep.all$Site,mr.UVrep.all$SampDate)]), ordered=TRUE)

#mr.UVrep.all$PathShort<-mr.UVrep.all$Pathway
# mr.UVrep.all$PathShort[(mr.UVrep.all$PathShort) == "Dissimilatory Sulfate Redox"] <- "D.SO4 RedOx"

head(mr.UVrep.all)

# convert all 0s to NAs so they appear gray in heat map
mr.UVrep.all$MR_SumCovPerKO_NA<-ifelse(mr.UVrep.all$MR_SumCovPerKO==0,NA,mr.UVrep.all$MR_SumCovPerKO)
head(mr.UVrep.all)

# For heatmap color gradient
max(mr.UVrep.all$MR_SumCovPerKO, na.rm=TRUE)
max(mr.UVrep.all$MR_SumCovPerKO, na.rm=TRUE)/2
min(mr.UVrep.all$MR_SumCovPerKO, na.rm=TRUE)

# Find mean coverage of fxns
UVrep.mgm.means<-data.frame(MeanScaledCov=colMeans(mr.cov.sum.UVrep.ko[,-1]))
UVrep.mgm.means$KO_ID<-rownames(UVrep.mgm.means)
UVrep.mgm.means$Category<-"UV-Damaged DNA Repair"
UVrep.mgm.means

unique(mr.UVrep.all$SampDate)
mr.UVrep.all$SampDatePlot<-as.character(mr.UVrep.all$SampDate)
mr.UVrep.all$SampDatePlot[which(mr.UVrep.all$SampDatePlot == "August.2020")] <- "Aug 2020"
mr.UVrep.all$SampDatePlot[which(mr.UVrep.all$SampDatePlot == "November.2020")] <- "Nov 2020"
mr.UVrep.all$SampDatePlot[which(mr.UVrep.all$SampDatePlot == "August.2021")] <- "Aug 2021"
mr.UVrep.all$SampDatePlot[which(mr.UVrep.all$SampDatePlot == "September.2021")] <- "Sep 2021"
mr.UVrep.all$SampDatePlot[which(mr.UVrep.all$SampDatePlot == "December.2021")] <- "Dec 2021"
mr.UVrep.all$SampDatePlot<-gsub("\\."," ",mr.UVrep.all$SampDatePlot)

unique(mr.UVrep.all$SampDatePlot)
mr.UVrep.all$SampDatePlot = factor(mr.UVrep.all$SampDatePlot, levels=unique(mr.UVrep.all$SampDatePlot[order(mr.UVrep.all$SampDate)]), ordered=TRUE)
unique(mr.UVrep.all$SampDatePlot)


# Figures below
# by SampleID

UVrep.hm1a<-ggplot(mr.UVrep.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="UV-Damage DNA Repair Functions in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(UVrep.hm1a,filename = "figures/MGM_Figs/Contigs/FxnDiv/UV_DNA_Repair/UV.DNA.Repair_KOFxns_MGMs_SampID_by_Function_heatmap.png", width=25, height=27, dpi=600,create.dir = TRUE)

UVrep.hm1a2<-ggplot(mr.UVrep.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="UV-Damage DNA Repair Functions in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(UVrep.hm1a2,filename = "figures/MGM_Figs/Contigs/FxnDiv/UV_DNA_Repair/UV.DNA.Repair_KOFxns_MGMs_SampID_by_Function_Site_heatmap.png", width=25, height=15, dpi=600,create.dir = TRUE)

UVrep.hm1a2a<-ggplot(mr.UVrep.all[mr.UVrep.all$MR_SumCovPerKO>30,], aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) + geom_text(label=round(mr.UVrep.all$MR_SumCovPerKO[mr.UVrep.all$MR_SumCovPerKO>30],2)) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="UV-Damage DNA Repair Functions in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Summed Coverages >= 30)",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site, scales="free", space = "free")

ggsave(UVrep.hm1a2a,filename = "figures/MGM_Figs/Contigs/FxnDiv/UV_DNA_Repair/UV.DNA.Repair_KOFxns_MGMs_SampleID_by_Function_Site_HighCov_labeled_heatmap.png", width=28, height=20, dpi=600,create.dir = TRUE)

UVrep.hm1a3<-ggplot(mr.UVrep.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=25),axis.title.y = element_text(size=25),legend.title.align=0.5, legend.title = element_text(size=23),
        axis.text= element_text(size=20),axis.text.x = element_text(angle=45, hjust=1,size=20),legend.text = element_text(size=18),plot.title = element_blank(),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_blank(),strip.text.x = element_text(size = 16,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+
  facet_grid(.~SampDatePlot,scales="free_x", space = "free")

ggsave(UVrep.hm1a3,filename = "figures/MGM_Figs/Contigs/FxnDiv/UV_DNA_Repair/UV.DNA.Repair_KOFxns_MGMs_SampID_by_Function_SampDate_heatmap.png", width=35, height=15, dpi=300,create.dir = TRUE)

#### Pull Out Osmoprotectant Fxns from MRN data ####
## heatmaps of traits of interest

mgm.mr[1:4,1:4]

# pull out heat shock functions from MR normalized, summed coverages (summed gene coverage per KO)
osmo.ko<-mgm.mr[,which(colnames(mgm.mr) %in% osmo.fxns$KO_ID)] # merge MRN data w/ S fxns found in contigs from KOFamScan
osmo.ko$SampleID<-rownames(osmo.ko)
osmo.ko.melt<-melt(osmo.ko, by="SampleID")
colnames(osmo.ko.melt)[which(names(osmo.ko.melt) == "variable")] <- "KO_ID"
colnames(osmo.ko.melt)[which(names(osmo.ko.melt) == "value")] <- "MR_SumCovPerKO"
head(osmo.ko.melt) #sanity check

mr.osmo.ko<-merge(osmo.ko.melt,osmo.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(mr.osmo.ko)
colnames(mr.osmo.ko)[which(names(mr.osmo.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
mr.cov.sum.osmo.ko<-as.data.frame(dcast(mr.osmo.ko, SampleID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(mr.cov.sum.osmo.ko)<-mr.cov.sum.osmo.ko$SampleID
mr.cov.sum.osmo.ko[1:4,]

# sanity check
mr.cov.sum.osmo.ko$`osmY; hyperosmotically inducible periplasmic protein`[1:4]
head(mr.osmo.ko)

#### Osmoprotectant KO Heat Maps ####
# see max & mean of summed
max(mr.cov.sum.osmo.ko[,-1])
min(mr.cov.sum.osmo.ko[,-1])
max(mr.cov.sum.osmo.ko[,-1])/2

# first heat map of metal KOs
##heatmap(as.matrix(mr.cov.sum.osmo.ko[,-1]), scale = "none")

colSums(mr.cov.sum.osmo.ko[,-1])
#mr.cov.sum.osmo.ko2 <- mr.cov.sum.osmo.ko[,which(colSums(mr.cov.sum.osmo.ko[,-1])>10)]

##heatmap(as.matrix(mr.cov.sum.osmo.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap -- using MR data that includes NAs so they are blocked out on heatmaps
mr.osmo.ko[1:4,]
mr.osmo.all<-merge(mr.osmo.ko,meta.all.scaled,by="SampleID")
head(mr.osmo.all)
#mr.osmo.all$PlotID<-mr.osmo.all$SampleID
mr.osmo.all$PlotID = factor(mr.osmo.all$SampleID, levels=unique(mr.osmo.all$SampleID[order(mr.osmo.all$Site,mr.osmo.all$SampDate)]), ordered=TRUE)

head(mr.osmo.all)

# convert all 0s to NAs so they appear gray in heat map
mr.osmo.all$MR_SumCovPerKO_NA<-ifelse(mr.osmo.all$MR_SumCovPerKO==0,NA,mr.osmo.all$MR_SumCovPerKO)
head(mr.osmo.all)

# For heatmap color gradient
max(mr.osmo.all$MR_SumCovPerKO, na.rm=TRUE)
max(mr.osmo.all$MR_SumCovPerKO, na.rm=TRUE)/2
min(mr.osmo.all$MR_SumCovPerKO, na.rm=TRUE)

# Find mean coverage of fxns
Osmo.mgm.means<-data.frame(MeanScaledCov=colMeans(mr.cov.sum.osmo.ko[,-1]))
Osmo.mgm.means$KO_ID<-rownames(Osmo.mgm.means)
Osmo.mgm.means$Category<-"Osmoprotectant"
Osmo.mgm.means

# Figures below
# by SampleID

osmo.hm1a<-ggplot(mr.osmo.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Osmoprotectant Transport in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(osmo.hm1a,filename = "figures/MGM_Figs/Contigs/FxnDiv/Osmoprotectant/Osmoprotectant_KOFxns_MGMs_SampID_by_Function_heatmap.png", width=25, height=27, dpi=600,create.dir = TRUE)

osmo.hm1a2<-ggplot(mr.osmo.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Osmoprotectant Transport in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(osmo.hm1a2,filename = "figures/MGM_Figs/Contigs/FxnDiv/Osmoprotectant/Osmoprotectant_KOFxns_MGMs_SampID_by_Function_Site_heatmap.png", width=25, height=15, dpi=600,create.dir = TRUE)

osmo.hm1a2a<-ggplot(mr.osmo.all[mr.osmo.all$MR_SumCovPerKO>20,], aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) + geom_text(label=round(mr.osmo.all$MR_SumCovPerKO[mr.osmo.all$MR_SumCovPerKO>20],2)) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Osmoprotectant Transport in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Summed Coverages >= 20)",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site, scales="free", space = "free")

ggsave(osmo.hm1a2a,filename = "figures/MGM_Figs/Contigs/FxnDiv/Osmoprotectant/Osmoprotectant_KOFxns_MGMs_SampleID_by_Function_Site_HighCov_labeled_heatmap.png", width=28, height=20, dpi=600,create.dir = TRUE)

osmo.hm1a3<-ggplot(mr.osmo.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Osmoprotectant Transport in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate,scales="free_x", space = "free")

ggsave(osmo.hm1a3,filename = "figures/MGM_Figs/Contigs/FxnDiv/Osmoprotectant/Osmoprotectant_KOFxns_MGMs_SampID_by_Function_SampDate_heatmap.png", width=35, height=15, dpi=600,create.dir = TRUE)

#### Pull Out ALL Genes of Interest per Pathway/Cycle from MR data ####
## heatmaps of traits of interest

mgm.mr[1:4,1:4]

# pull out sulfur functions from MR normalized, summed coverages (summed coverage per KO)
All_GOI.ko<-mgm.mr[,which(colnames(mgm.mr) %in% All_GOI.fxns$KO_ID)] # merge MR data w/ All_GOI-related fxns found in contigs from KOFamScan
All_GOI.ko$SampleID<-rownames(All_GOI.ko)
All_GOI.ko.melt<-melt(All_GOI.ko, by="SampleID")
colnames(All_GOI.ko.melt)[which(names(All_GOI.ko.melt) == "variable")] <- "KO_ID"
colnames(All_GOI.ko.melt)[which(names(All_GOI.ko.melt) == "value")] <- "MR_SumCovPerKO"
head(All_GOI.ko.melt) #sanity check

mr.All_GOI.ko<-merge(All_GOI.ko.melt,all_goi.kegg,by.x=c("KO_ID"),by.y=c("KO_ID"))
head(mr.All_GOI.ko)
colnames(mr.All_GOI.ko)[which(names(mr.All_GOI.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
mr.cov.sum.All_GOI.ko<-as.data.frame(dcast(mr.All_GOI.ko, SampleID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###
rownames(mr.cov.sum.All_GOI.ko)<-mr.cov.sum.All_GOI.ko$SampleID
mr.cov.sum.All_GOI.ko[1:4,1:4]

#### All KOs of Interest - KO Heat Maps ####
# see max & mean of summed
max(mr.cov.sum.All_GOI.ko[,-1])
min(mr.cov.sum.All_GOI.ko[,-1])
max(mr.cov.sum.All_GOI.ko[,-1])/2

# first heat map of All_GOI KOs
##heatmap(as.matrix(mr.cov.sum.All_GOI.ko[,-1]), scale = "none")

colSums(mr.cov.sum.All_GOI.ko[,-1])
#mr.cov.sum.All_GOI.ko2 <- mr.cov.sum.All_GOI.ko[,which(colSums(mr.cov.sum.All_GOI.ko[,-1])>10)]

##heatmap(as.matrix(mr.cov.sum.All_GOI.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap -- using MR data that includes NAs so they are blocked out on heatmaps
mr.All_GOI.ko[1:4,]
mr.All_GOI.all<-merge(mr.All_GOI.ko,meta.all.scaled,by="SampleID")
head(mr.All_GOI.all)
#mr.All_GOI.all$PlotID<-mr.All_GOI.all$SampleID
mr.All_GOI.all$PlotID = factor(mr.All_GOI.all$SampleID, levels=unique(mr.All_GOI.all$SampleID[order(mr.All_GOI.all$Site,mr.All_GOI.all$SampDate)]), ordered=TRUE)

unique(mr.All_GOI.all$Category)
mr.All_GOI.all$CatShort<-mr.All_GOI.all$Category
mr.All_GOI.all$CatShort[(mr.All_GOI.all$CatShort) == "LPS Biosynthesis"] <- "LPS"
mr.All_GOI.all$CatShort[(mr.All_GOI.all$CatShort) == "Sporulation"] <- "Spore"
mr.All_GOI.all$CatShort[(mr.All_GOI.all$CatShort) == "Temperature Shock Resistance"] <- "Temp."
mr.All_GOI.all$CatShort[(mr.All_GOI.all$CatShort) == "Osmoprotectant Transport/Accumulation"] <- "Osmo"
mr.All_GOI.all$CatShort[(mr.All_GOI.all$CatShort) == "Quorum Sensing"] <- "QS"
mr.All_GOI.all$CatShort[(mr.All_GOI.all$CatShort) == "UV-Damange DNA Repair"] <- "UV"
mr.All_GOI.all$CatShort[(mr.All_GOI.all$CatShort) == "House Keeping Gene"] <- "House"

head(mr.All_GOI.all)

mr.All_GOI.all$KOFxn_Short<-gsub("; .*","",mr.All_GOI.all$KO_Function.KEGG)

# convert all 0s to NAs so they appear gray in heat map
mr.All_GOI.all$MR_SumCovPerKO_NA<-ifelse(mr.All_GOI.all$MR_SumCovPerKO==0,NA,mr.All_GOI.all$MR_SumCovPerKO)
head(mr.All_GOI.all)

# For heatmap color gradient
max(mr.All_GOI.all$MR_SumCovPerKO, na.rm=TRUE)
max(mr.All_GOI.all$MR_SumCovPerKO, na.rm=TRUE)/2
min(mr.All_GOI.all$MR_SumCovPerKO, na.rm=TRUE)

# Find mean coverage of fxns
All_GOI.mgm.means<-data.frame(MeanScaledCov=colMeans(mr.cov.sum.All_GOI.ko[,-1]))
All_GOI.mgm.means$KO_ID<-rownames(All_GOI.mgm.means)
#All_GOI.mgm.means$Category<-"All Genes of Interest"
All_GOI.mgm.means

# save it all!
# save.image("data/MetagenomesSSD_mgm_FxnBetaDiv_data_Ready.Rdata")

# Figures below
# by SampleID

All_GOI.hm1a<-ggplot(mr.All_GOI.all, aes(PlotID, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Dust Microbial Survival in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(All_GOI.hm1a,filename = "figures/MGM_Figs/Contigs/FxnDiv/All_GOI/All_GOI_KOFxns_MGMs_SampID_by_Function_heatmap.png", width=30, height=30, dpi=600,create.dir = TRUE)

All_GOI.hm1a2<-ggplot(mr.All_GOI.all, aes(PlotID, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Dust Microbial Survival in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(All_GOI.hm1a2,filename = "figures/MGM_Figs/Contigs/FxnDiv/All_GOI/All_GOI_KOFxns_MGMs_SampID_by_Function_Site_heatmap.png", width=30, height=30, dpi=600,create.dir = TRUE)

All_GOI.hm1a2b<-ggplot(mr.All_GOI.all[mr.All_GOI.all$MR_SumCovPerKO>30,], aes(PlotID, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Dust Microbial Survival in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Summed Coverages >= 30)",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(All_GOI.hm1a2b,filename = "figures/MGM_Figs/Contigs/FxnDiv/All_GOI/All_GOI_KOFxns_MGMs_SampID_by_Function_Site_HigherCov_heatmap.png", width=30, height=30, dpi=600,create.dir = TRUE)

All_GOI.hm1a2c<-ggplot(mr.All_GOI.all[mr.All_GOI.all$MR_SumCovPerKO>30,], aes(PlotID, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Dust Microbial Survival in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Summed Coverages >= 30)",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate,scales="free_x", space = "free")

ggsave(All_GOI.hm1a2c,filename = "figures/MGM_Figs/Contigs/FxnDiv/All_GOI/All_GOI_KOFxns_MGMs_SampID_by_Function_SampDate_HigherCov_heatmap.png", width=30, height=30, dpi=600,create.dir = TRUE)

All_GOI.hm1a2d<-ggplot(mr.All_GOI.all[mr.All_GOI.all$MR_SumCovPerKO>30,], aes(PlotID, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) + geom_text(label=round(mr.All_GOI.all$MR_SumCovPerKO[mr.All_GOI.all$MR_SumCovPerKO>30],2)) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Dust Microbial Survival in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Summed Coverages >= 30)",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate,scales="free_x", space = "free")

ggsave(All_GOI.hm1a2d,filename = "figures/MGM_Figs/Contigs/FxnDiv/All_GOI/All_GOI_KOFxns_MGMs_SampID_by_Function_SampDate_HigherCov_heatmap_labeled.png", width=35, height=30, dpi=600,create.dir = TRUE)

All_GOI.hm1a3<-ggplot(mr.All_GOI.all, aes(PlotID, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Dust Microbial Survival in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(CatShort~Site, scales="free", space = "free")

ggsave(All_GOI.hm1a3,filename = "figures/MGM_Figs/Contigs/FxnDiv/All_GOI/ALL_GOI_KOFxns_MGMs_SampleID_by_Function_KO_Category_Site_heatmap.png", width=35, height=30, dpi=600,create.dir = TRUE)

All_GOI.hm1a3a<-ggplot(mr.All_GOI.all[mr.All_GOI.all$MR_SumCovPerKO>30,], aes(PlotID, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Dust Microbial Survival in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Summed Coverages >= 30)",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(CatShort~Site, scales="free", space = "free")

ggsave(All_GOI.hm1a3a,filename = "figures/MGM_Figs/Contigs/FxnDiv/All_GOI/All_GOI_KOFxns_MGMs_SampleID_by_Function_KO_Category_Site_HigherCov_heatmap.png", width=35, height=30, dpi=600,create.dir = TRUE)

All_GOI.hm1a3b<-ggplot(mr.All_GOI.all[mr.All_GOI.all$MR_SumCovPerKO>30,], aes(PlotID, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) + geom_text(label=round(mr.All_GOI.all$MR_SumCovPerKO[mr.All_GOI.all$MR_SumCovPerKO>30],2)) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Dust Microbial Survival in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Summed Coverages >= 30)",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(CatShort~Site, scales="free", space = "free")

ggsave(All_GOI.hm1a3b,filename = "figures/MGM_Figs/Contigs/FxnDiv/All_GOI/All_GOI_KOFxns_MGMs_SampleID_by_Function_KO_Category_Site_HigherCov_heatmap_labeled.png", width=35, height=30, dpi=600,create.dir = TRUE)

All_GOI.hm1a4<-ggplot(mr.All_GOI.all, aes(PlotID, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(fill="MR Coverage Per KO") +
  theme(legend.title.align=0.5, legend.title = element_text(size=22),
        axis.text = element_text(size=24),axis.text.x = element_text(angle=45, hjust=1,size=30),legend.text = element_text(size=20),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=11),strip.text = element_text(size = 25,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(CatShort~SampDate, scales="free", space = "free")

ggsave(All_GOI.hm1a4,filename = "figures/MGM_Figs/Contigs/FxnDiv/All_GOI/All_GOI_KOFxns_MGMs_SampleID_by_Function_KO_Category_SampDate_heatmap.png", width=35, height=47, dpi=600,limitsize=FALSE,create.dir = TRUE)

All_GOI.hm1a4a<-ggplot(mr.All_GOI.all[mr.All_GOI.all$MR_SumCovPerKO>30,], aes(PlotID, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Dust Microbial Survival in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Summed Coverages >= 30)",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(CatShort~SampDate, scales="free", space = "free")

ggsave(All_GOI.hm1a4a,filename = "figures/MGM_Figs/Contigs/FxnDiv/All_GOI/All_GOI_KOFxns_MGMs_SampleID_by_Function_KO_Category_SampDate_HigherCov_heatmap.png", width=35, height=30, dpi=600,create.dir = TRUE)

All_GOI.hm1a4b<-ggplot(mr.All_GOI.all[mr.All_GOI.all$MR_SumCovPerKO>30,], aes(PlotID, KOFxn_Short, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) + geom_text(label=round(mr.All_GOI.all$MR_SumCovPerKO[mr.All_GOI.all$MR_SumCovPerKO>30],2)) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Dust Microbial Survival in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Summed Coverages >= 30)",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(CatShort~SampDate, scales="free", space = "free")

ggsave(All_GOI.hm1a4b,filename = "figures/MGM_Figs/Contigs/FxnDiv/All_GOI/All_GOI_KOFxns_MGMs_SampleID_by_Function_KO_Category_SampDate_HigherCov_heatmap_labeled.png", width=35, height=30, dpi=600,create.dir = TRUE)

#### Pull Out Antibiotic Resistance Gene (ARGs) Fxns from MRN data ####
## heatmaps of traits of interest

mgm.mr[1:4,1:4]

# pull out ARG functions from MR normalized, summed coverages (summed gene coverage per KO)
ARG.ko<-mgm.mr[,which(colnames(mgm.mr) %in% ARG.fxns$KO_ID)] # merge MRN data w/ S fxns found in contigs from KOFamScan
ARG.ko$SampleID<-rownames(ARG.ko)
ARG.ko.melt<-melt(ARG.ko, by="SampleID")
colnames(ARG.ko.melt)[which(names(ARG.ko.melt) == "variable")] <- "KO_ID"
colnames(ARG.ko.melt)[which(names(ARG.ko.melt) == "value")] <- "MR_SumCovPerKO"
head(ARG.ko.melt) #sanity check

mr.ARG.ko<-merge(ARG.ko.melt,ARG.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(mr.ARG.ko)
colnames(mr.ARG.ko)[which(names(mr.ARG.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
mr.cov.sum.ARG.ko<-as.data.frame(dcast(mr.ARG.ko, SampleID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(mr.cov.sum.ARG.ko)<-mr.cov.sum.ARG.ko$SampleID
mr.cov.sum.ARG.ko[1:4,]

# sanity check
mr.cov.sum.ARG.ko$`aminoglycoside 6'-N-acetyltransferase I`[1:4]
head(mr.ARG.ko)

#### ARG KO Heat Maps ####
# see max & mean of summed
max(mr.cov.sum.ARG.ko[,-1])
min(mr.cov.sum.ARG.ko[,-1])
max(mr.cov.sum.ARG.ko[,-1])/2

# first heat map of metal KOs
##heatmap(as.matrix(mr.cov.sum.ARG.ko[,-1]), scale = "none")

colSums(mr.cov.sum.ARG.ko[,-1])
#mr.cov.sum.ARG.ko2 <- mr.cov.sum.ARG.ko[,which(colSums(mr.cov.sum.ARG.ko[,-1])>10)]

##heatmap(as.matrix(mr.cov.sum.ARG.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap -- using MR data that includes NAs so they are blocked out on heatmaps
mr.ARG.ko[1:4,]
mr.ARG.all<-merge(mr.ARG.ko,meta.all.scaled,by="SampleID")
head(mr.ARG.all)
#mr.ARG.all$PlotID<-mr.ARG.all$SampleID
mr.ARG.all$PlotID = factor(mr.ARG.all$SampleID, levels=unique(mr.ARG.all$SampleID[order(mr.ARG.all$Site,mr.ARG.all$SampDate)]), ordered=TRUE)

head(mr.ARG.all)

# convert all 0s to NAs so they appear gray in heat map
mr.ARG.all$MR_SumCovPerKO_NA<-ifelse(mr.ARG.all$MR_SumCovPerKO==0,NA,mr.ARG.all$MR_SumCovPerKO)
head(mr.ARG.all)

# For heatmap color gradient
max(mr.ARG.all$MR_SumCovPerKO, na.rm=TRUE)
max(mr.ARG.all$MR_SumCovPerKO, na.rm=TRUE)/2
min(mr.ARG.all$MR_SumCovPerKO, na.rm=TRUE)

# Find mean coverage of fxns
ARG.mgm.means<-data.frame(MeanScaledCov=colMeans(mr.cov.sum.ARG.ko[,-1]))
ARG.mgm.means$KO_ID<-rownames(ARG.mgm.means)
ARG.mgm.means$Category<-"ARG"
ARG.mgm.means

unique(mr.ARG.all$SampDate)
mr.ARG.all$SampDatePlot<-as.character(mr.ARG.all$SampDate)
mr.ARG.all$SampDatePlot[which(mr.ARG.all$SampDatePlot == "August.2020")] <- "Aug 2020"
mr.ARG.all$SampDatePlot[which(mr.ARG.all$SampDatePlot == "November.2020")] <- "Nov 2020"
mr.ARG.all$SampDatePlot[which(mr.ARG.all$SampDatePlot == "August.2021")] <- "Aug 2021"
mr.ARG.all$SampDatePlot[which(mr.ARG.all$SampDatePlot == "September.2021")] <- "Sep 2021"
mr.ARG.all$SampDatePlot[which(mr.ARG.all$SampDatePlot == "December.2021")] <- "Dec 2021"
mr.ARG.all$SampDatePlot<-gsub("\\."," ",mr.ARG.all$SampDatePlot)

unique(mr.ARG.all$SampDatePlot)
mr.ARG.all$SampDatePlot = factor(mr.ARG.all$SampDatePlot, levels=unique(mr.ARG.all$SampDatePlot[order(mr.ARG.all$SampDate)]), ordered=TRUE)
unique(mr.ARG.all$SampDatePlot)

# Figures below
# by SampleID

ARG.hm1a<-ggplot(mr.ARG.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Antibiotic Resistance Genes in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(ARG.hm1a,filename = "figures/MGM_Figs/Contigs/FxnDiv/AntibioticResistanceGenes/ARGs_KOFxns_MGMs_SampID_by_Function_heatmap.png", width=25, height=27, dpi=600,create.dir = TRUE)

ARG.hm1a2<-ggplot(mr.ARG.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="ARGs Transport in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")

ggsave(ARG.hm1a2,filename = "figures/MGM_Figs/Contigs/FxnDiv/AntibioticResistanceGenes/ARGs_KOFxns_MGMs_SampID_by_Function_Site_heatmap.png", width=25, height=15, dpi=600,create.dir = TRUE)

ARG.hm1a2a<-ggplot(mr.ARG.all[mr.ARG.all$MR_SumCovPerKO>20,], aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.25) + geom_text(label=round(mr.ARG.all$MR_SumCovPerKO[mr.ARG.all$MR_SumCovPerKO>20],2)) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="ARGs Transport in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Summed Coverages >= 20)",fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site, scales="free", space = "free")

ggsave(ARG.hm1a2a,filename = "figures/MGM_Figs/Contigs/FxnDiv/AntibioticResistanceGenes/ARGs_KOFxns_MGMs_SampleID_by_Function_Site_HighCov_labeled_heatmap.png", width=28, height=20, dpi=600,create.dir = TRUE)

ARG.hm1a3<-ggplot(mr.ARG.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(fill="MR Coverage Per KO") +
  theme(axis.title.x = element_text(size=25),axis.title.y = element_text(size=25),legend.title.align=0.5, legend.title = element_text(size=23),
        axis.text = element_text(size=20),axis.text.x = element_text(angle=45, hjust=1,size=17),legend.text = element_text(size=17),plot.title = element_blank(),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_blank(),strip.text.x = element_text(size = 17,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDatePlot,scales="free_x", space = "free")

ggsave(ARG.hm1a3,filename = "figures/MGM_Figs/Contigs/FxnDiv/AntibioticResistanceGenes/ARGs_KOFxns_MGMs_SampID_by_Function_SampDate_heatmap.png", width=35, height=15, dpi=300,create.dir = TRUE)

#### Pull out House Keeping Gene Coverage ####
# will use this later!
mgm.mr$K03553 # recA is a house keeping gene

# pull out HK functions from MR normalized, summed coverages (summed gene coverage per KO)
HK.ko<-data.frame(KO_ID="K03553",MR_SumCovPerKO=mgm.mr$K03553,SampleID=rownames(mgm.mr)) # merge MR data w/ S fxns found in contigs from KOFamScan
head(HK.ko) #sanity check

mr.HK.ko<-merge(HK.ko,Housekeep.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(mr.HK.ko)
colnames(mr.HK.ko)[which(names(mr.HK.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website

head(mr.HK.ko)

## create table of Sample x KO_ID
mr.HK.ko_table<-as.data.frame(dcast(mr.HK.ko, SampleID~KO_Function.KEGG, value.var="MR_SumCovPerKO"))
mr.HK.ko_table[1:5,]

## merge HK functions with scaled metadata
mr.HK.ko.all<-merge(meta.all.scaled,mr.HK.ko_table,by="SampleID")
mr.HK.ko.all[1:5,]

# Find mean coverage of recA
HK.mgm.means<-data.frame(MeanScaledCov=mean(mr.HK.ko.all$`recA; recombination protein RecA`),KO_ID=names(mr.HK.ko.all)[ncol(mr.HK.ko.all)])
rownames(HK.mgm.means)<-HK.mgm.means$KO_ID
HK.mgm.means$Category<-"House Keeping"
HK.mgm.means

# #### Pull Out Metal Resistance Fxns from MRN data ####
# ## heatmaps of traits of interest
#
# mgm.mr[1:4,1:4]
#
# # pull out heat shock functions from MR normalized, summed coverages (summed gene coverage per KO)
# metal.ko<-mgm.mr[,which(colnames(mgm.mr) %in% metal.fxns$KO_ID)] # merge MRN data w/ S fxns found in contigs from KOFamScan
# metal.ko$SampleID<-rownames(metal.ko)
# metal.ko.melt<-melt(metal.ko, by="SampleID")
# colnames(metal.ko.melt)[which(names(metal.ko.melt) == "variable")] <- "KO_ID"
# colnames(metal.ko.melt)[which(names(metal.ko.melt) == "value")] <- "MR_SumCovPerKO"
# head(metal.ko.melt) #sanity check
#
# mr.metal.ko<-merge(metal.ko.melt,metal.re.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
# head(mr.metal.ko)
# colnames(mr.metal.ko)[which(names(mr.metal.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
# mr.cov.sum.metal.ko<-as.data.frame(dcast(mr.metal.ko, SampleID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###just dcast, nothing is being added here!
# rownames(mr.cov.sum.metal.ko)<-mr.cov.sum.metal.ko$SampleID
# mr.cov.sum.metal.ko[1:4,]
#
# # sanity check
# mr.cov.sum.metal.ko$`copC, pcoC; copper resistance protein C`[1:4]
# head(mr.metal.ko)
#
# #### Metal Resistance KO Heat Maps ####
# # see max & mean of summed
# max(mr.cov.sum.metal.ko[,-1])
# min(mr.cov.sum.metal.ko[,-1])
# max(mr.cov.sum.metal.ko[,-1])/2
#
# # first heat map of metal KOs
# ##heatmap(as.matrix(mr.cov.sum.metal.ko[,-1]), scale = "none")
#
# colSums(mr.cov.sum.metal.ko[,-1])
# #mr.cov.sum.metal.ko2 <- mr.cov.sum.metal.ko[,which(colSums(mr.cov.sum.metal.ko[,-1])>10)]
#
# ##heatmap(as.matrix(mr.cov.sum.metal.ko[,-1]), scale = "none")
#
# # prep for ggplot2 heatmap -- using MR data that includes NAs so they are blocked out on heatmaps
# mr.metal.ko[1:4,]
# mr.metal.all<-merge(mr.metal.ko,meta.all.scaled,by="SampleID")
# head(mr.metal.all)
# #mr.metal.all$PlotID<-mr.metal.all$SampleID
# mr.metal.all$PlotID = factor(mr.metal.all$SampleID, levels=unique(mr.metal.all$SampleID[order(mr.metal.all$Site,mr.metal.all$SampDate)]), ordered=TRUE)
#
# #mr.metal.all$PathShort<-mr.metal.all$Pathway
# # mr.metal.all$PathShort[(mr.metal.all$PathShort) == "Dissimilatory Sulfate Redox"] <- "D.SO4 RedOx"
#
# head(mr.metal.all)
#
# # convert all 0s to NAs so they appear gray in heat map
# mr.metal.all$MR_SumCovPerKO_NA<-ifelse(mr.metal.all$MR_SumCovPerKO==0,NA,mr.metal.all$MR_SumCovPerKO)
# head(mr.metal.all)
#
# # For heatmap color gradient
# max(mr.metal.all$MR_SumCovPerKO, na.rm=TRUE)
# max(mr.metal.all$MR_SumCovPerKO, na.rm=TRUE)/2
# min(mr.metal.all$MR_SumCovPerKO, na.rm=TRUE)
#
# # Figures below
# # by SampleID
#
# metal.hm1a<-ggplot(mr.metal.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Metal Resistance in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))
#
# ggsave(metal.hm1a,filename = "figures/MGM_Figs/Contigs/FxnDiv/MetalResistance/MetalResistance_KOFxns_MGMs_SampID_by_Function_heatmap.png", width=25, height=27, dpi=600,create.dir = TRUE)
#
# metal.hm1a2<-ggplot(mr.metal.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.15) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Metal Resistance in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")
#
# ggsave(metal.hm1a2,filename = "figures/MGM_Figs/Contigs/FxnDiv/MetalResistance/MetalResistance_KOFxns_MGMs_SampID_by_Function_Site_heatmap.png", width=25, height=15, dpi=600,create.dir = TRUE)
#
# metal.hm1a2a<-ggplot(mr.metal.all[mr.metal.all$MR_SumCovPerKO>20,], aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.25) + geom_text(label=round(mr.metal.all$MR_SumCovPerKO[mr.metal.all$MR_SumCovPerKO>20],2)) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Metal Resistance in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO\n(Only Summed Coverages >= 20)",fill="MR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site, scales="free", space = "free")
#
# ggsave(metal.hm1a2a,filename = "figures/MGM_Figs/Contigs/FxnDiv/MetalResistance/MetalResistance_KOFxns_MGMs_SampleID_by_Function_Site_HighCov_labeled_heatmap.png", width=28, height=20, dpi=600,create.dir = TRUE)

# #### Pull Out Arsenic Metabolic Fxns from MRN data ####
# ## heatmaps of traits of interest
#
# mgm.mr[1:4,1:4]
#
# # pull out heat shock functions from MR normalized, summed coverages (summed gene coverage per KO)
# arsen.ko<-mgm.mr[,which(colnames(mgm.mr) %in% arsen.fxns$KO_ID)] # merge MRN data w/ S fxns found in contigs from KOFamScan
# arsen.ko$SampleID<-rownames(arsen.ko)
# arsen.ko.melt<-melt(arsen.ko, by="SampleID")
# colnames(arsen.ko.melt)[which(names(arsen.ko.melt) == "variable")] <- "KO_ID"
# colnames(arsen.ko.melt)[which(names(arsen.ko.melt) == "value")] <- "MR_SumCovPerKO"
# head(arsen.ko.melt) #sanity check
#
# mr.arsen.ko<-merge(arsen.ko.melt,arsen.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
# head(mr.arsen.ko)
# colnames(mr.arsen.ko)[which(names(mr.arsen.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
# mr.cov.sum.arsen.ko<-as.data.frame(dcast(mr.arsen.ko, SampleID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###just dcast, nothing is being added here!
# rownames(mr.cov.sum.arsen.ko)<-mr.cov.sum.arsen.ko$SampleID
# mr.cov.sum.arsen.ko[1:4,]
#
# # sanity check
# mr.cov.sum.arsen.ko$`arsH; arsenical resistance protein ArsH`[1:4]
# head(mr.arsen.ko)
#
# #### Arsenic KO Heat Maps ####
# # see max & mean of summed
# max(mr.cov.sum.arsen.ko[,-1])
# min(mr.cov.sum.arsen.ko[,-1])
# max(mr.cov.sum.arsen.ko[,-1])/2
#
# # first heat map of metal KOs
# ##heatmap(as.matrix(mr.cov.sum.arsen.ko[,-1]), scale = "none")
#
# colSums(mr.cov.sum.arsen.ko[,-1])
# #mr.cov.sum.arsen.ko2 <- mr.cov.sum.arsen.ko[,which(colSums(mr.cov.sum.arsen.ko[,-1])>10)]
#
# ##heatmap(as.matrix(mr.cov.sum.arsen.ko[,-1]), scale = "none")
#
# # prep for ggplot2 heatmap -- using MR data that includes NAs so they are blocked out on heatmaps
# mr.arsen.ko[1:4,]
# mr.arsen.all<-merge(mr.arsen.ko,meta.all.scaled,by="SampleID")
# head(mr.arsen.all)
# #mr.arsen.all$PlotID<-mr.arsen.all$SampleID
# mr.arsen.all$PlotID = factor(mr.arsen.all$SampleID, levels=unique(mr.arsen.all$SampleID[order(mr.arsen.all$Site,mr.arsen.all$SampDate)]), ordered=TRUE)
#
# #mr.arsen.all$PathShort<-mr.arsen.all$Pathway
# # mr.arsen.all$PathShort[(mr.arsen.all$PathShort) == "Dissimilatory Sulfate Redox"] <- "D.SO4 RedOx"
#
# head(mr.arsen.all)
#
# # convert all 0s to NAs so they appear gray in heat map
# mr.arsen.all$MR_SumCovPerKO_NA<-ifelse(mr.arsen.all$MR_SumCovPerKO==0,NA,mr.arsen.all$MR_SumCovPerKO)
# head(mr.arsen.all)
#
# # For heatmap color gradient
# max(mr.arsen.all$MR_SumCovPerKO, na.rm=TRUE)
# max(mr.arsen.all$MR_SumCovPerKO, na.rm=TRUE)/2
# min(mr.arsen.all$MR_SumCovPerKO, na.rm=TRUE)
#
# # Figures below
# # by SampleID
#
# arsen.hm1a<-ggplot(mr.arsen.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Arsenic Metabolism in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))
#
# ggsave(arsen.hm1a,filename = "figures/MGM_Figs/Contigs/FxnDiv/Arsenic/Arsenic_KOFxns_MGMs_SampID_by_Function_heatmap.png", width=25, height=27, dpi=600,create.dir = TRUE)
#
# arsen.hm1a2<-ggplot(mr.arsen.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.15) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Arsenic Metabolism in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")
#
# ggsave(arsen.hm1a2,filename = "figures/MGM_Figs/Contigs/FxnDiv/Arsenic/Arsenic_KOFxns_MGMs_SampID_by_Function_Site_heatmap.png", width=25, height=15, dpi=600,create.dir = TRUE)
#
#
# #### Pull Out Selenium Metabolic Fxns from MRN data ####
# ## heatmaps of traits of interest
#
# mgm.mr[1:4,1:4]
#
# # pull out heat shock functions from MR normalized, summed coverages (summed gene coverage per KO)
# selen.ko<-mgm.mr[,which(colnames(mgm.mr) %in% selen.fxns$KO_ID)] # merge MRN data w/ S fxns found in contigs from KOFamScan
# selen.ko$SampleID<-rownames(selen.ko)
# selen.ko.melt<-melt(selen.ko, by="SampleID")
# colnames(selen.ko.melt)[which(names(selen.ko.melt) == "variable")] <- "KO_ID"
# colnames(selen.ko.melt)[which(names(selen.ko.melt) == "value")] <- "MR_SumCovPerKO"
# head(selen.ko.melt) #sanity check
#
# mr.selen.ko<-merge(selen.ko.melt,selen.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
# head(mr.selen.ko)
# colnames(mr.selen.ko)[which(names(mr.selen.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
# mr.cov.sum.selen.ko<-as.data.frame(dcast(mr.selen.ko, SampleID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###just dcast, nothing is being added here!
# rownames(mr.cov.sum.selen.ko)<-mr.cov.sum.selen.ko$SampleID
# mr.cov.sum.selen.ko[1:4,]
#
# # sanity check
# mr.cov.sum.selen.ko$`trxB, TRR; thioredoxin reductase (NADPH) [EC:1.8.1.9]`[1:4]
# head(mr.selen.ko)
#
# #### Selenium KO Heat Maps ####
# # see max & mean of summed
# max(mr.cov.sum.selen.ko[,-1])
# min(mr.cov.sum.selen.ko[,-1])
# max(mr.cov.sum.selen.ko[,-1])/2
#
# # first heat map of metal KOs
# ##heatmap(as.matrix(mr.cov.sum.selen.ko[,-1]), scale = "none")
#
# colSums(mr.cov.sum.selen.ko[,-1])
# #mr.cov.sum.selen.ko2 <- mr.cov.sum.selen.ko[,which(colSums(mr.cov.sum.selen.ko[,-1])>10)]
#
# ##heatmap(as.matrix(mr.cov.sum.selen.ko[,-1]), scale = "none")
#
# # prep for ggplot2 heatmap -- using MR data that includes NAs so they are blocked out on heatmaps
# mr.selen.ko[1:4,]
# mr.selen.all<-merge(mr.selen.ko,meta.all.scaled,by="SampleID")
# head(mr.selen.all)
# #mr.selen.all$PlotID<-mr.selen.all$SampleID
# mr.selen.all$PlotID = factor(mr.selen.all$SampleID, levels=unique(mr.selen.all$SampleID[order(mr.selen.all$Site,mr.selen.all$SampDate)]), ordered=TRUE)
#
# #mr.selen.all$PathShort<-mr.selen.all$Pathway
# # mr.selen.all$PathShort[(mr.selen.all$PathShort) == "Dissimilatory Sulfate Redox"] <- "D.SO4 RedOx"
#
# head(mr.selen.all)
#
# # convert all 0s to NAs so they appear gray in heat map
# mr.selen.all$MR_SumCovPerKO_NA<-ifelse(mr.selen.all$MR_SumCovPerKO==0,NA,mr.selen.all$MR_SumCovPerKO)
# head(mr.selen.all)
#
# # For heatmap color gradient
# max(mr.selen.all$MR_SumCovPerKO, na.rm=TRUE)
# max(mr.selen.all$MR_SumCovPerKO, na.rm=TRUE)/2
# min(mr.selen.all$MR_SumCovPerKO, na.rm=TRUE)
#
# # Figures below
# # by SampleID
#
# selen.hm1a<-ggplot(mr.selen.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Selenium Metabolism in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))
#
# ggsave(selen.hm1a,filename = "figures/MGM_Figs/Contigs/FxnDiv/Selenium/Selenium_KOFxns_MGMs_SampID_by_Function_heatmap.png", width=25, height=27, dpi=600,create.dir = TRUE)
#
# selen.hm1a2<-ggplot(mr.selen.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.15) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Selenium Metabolism in Salton Sea Dust Metagenomes",subtitle="Using Median-Ratio Normalized Gene Coverage Summed by KO",fill="MR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~Site,scales="free_x", space = "free")
#
# ggsave(selen.hm1a2,filename = "figures/MGM_Figs/Contigs/FxnDiv/Selenium/Selenium_KOFxns_MGMs_SampID_by_Function_Site_heatmap.png", width=25, height=15, dpi=600,create.dir = TRUE)
#

# #### Pull Out Phototrophy Fxns from MR data ####
# ## heatmaps of traits of interest
#
# mgm.mr[1:4,1:4]
#
# # pull out Phototrophy functions from MR normalized, summed coverages (summed coverage per KO)
# photo.ko<-mgm.mr[,which(colnames(mgm.mr) %in% photo.fxns$KO_ID)] # merge MR data w/ photoon-related fxns found in contigs from KOFamScan
# photo.ko$SampleID<-rownames(photo.ko)
# photo.ko.melt<-melt(photo.ko, by="SampleID")
# colnames(photo.ko.melt)[which(names(photo.ko.melt) == "variable")] <- "KO_ID"
# colnames(photo.ko.melt)[which(names(photo.ko.melt) == "value")] <- "MR_SumCovPerKO"
# head(photo.ko.melt) #sanity check
#
# mr.photo.ko<-merge(photo.ko.melt,photo.kegg,by.x=c("KO_ID"),by.y=c("KO_ID"))
# head(mr.photo.ko)
# colnames(mr.photo.ko)[which(names(mr.photo.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
# mr.cov.sum.photo.ko<-as.data.frame(dcast(mr.photo.ko, SampleID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###
# rownames(mr.cov.sum.photo.ko)<-mr.cov.sum.photo.ko$SampleID
# mr.cov.sum.photo.ko[1:4,1:4]
#
# ### Phototrophy Heat Maps ####
# # see max & mean of summed
# max(mr.cov.sum.photo.ko[,-1])
# mean(as.matrix(mr.cov.sum.photo.ko[,-1]))
#
# # first heat map of sulfur KOs
# #heatmap(as.matrix(mr.cov.sum.photo.ko[,-1]), scale = "none")
#
# colSums(mr.cov.sum.photo.ko[,-1])
#
# #heatmap(as.matrix(mr.cov.sum.photo.ko[,-1]), scale = "none")
#
# # prep for ggplot2 heatmap
# mr.photo.ko[1:4,]
# mr.photo.all<-merge(mr.photo.ko,meta.all.scaled,by="SampleID")
# head(mr.photo.all)
# mr.photo.all$PlotID = factor(mr.photo.all$SampleID, levels=unique(mr.photo.all$SampleID[order(mr.photo.all$Site,mr.photo.all$SampDate)]), ordered=TRUE)
# mr.photo.all$SampDate<-gsub("\\."," ",mr.photo.all$SampDate)
# mr.photo.all$SampDate<-factor(mr.photo.all$SampDate, levels=c("August 2021","December 2021","April 2022"))
#
# unique(mr.photo.all$Pathway)
# mr.photo.all$PathShort<-mr.photo.all$Pathway
# mr.photo.all$PathShort[(mr.photo.all$PathShort) == "Proteorhodopsin"] <- "PR"
# #mr.photo.all$PathShort[(mr.photo.all$PathShort) == "Bacteriorhodopsin"] <- "BR"
# mr.photo.all$PathShort[(mr.photo.all$PathShort) == "Sensory Rhodopsin"] <- "SR"
# #mr.photo.all$PathShort[(mr.photo.all$PathShort) == "Halorhodopsin"] <- "HR"
# mr.photo.all$PathShort[(mr.photo.all$PathShort) == "Photosystem II"] <- "PS II"
# mr.photo.all$PathShort[(mr.photo.all$PathShort) == "Photosystem I"] <- "PS I"
# mr.photo.all$PathShort[(mr.photo.all$PathShort) == "Anoxygenic Photosystem II"] <- "AnOx PS"
#
# mr.photo.all$Pathway<-factor(mr.photo.all$Pathway,levels=c("Proteorhodopsin","Sensory Rhodopsin","Photosystem II","Photosystem I","Anoxygenic Photosystem II"))
# mr.photo.all$PathShort<-factor(mr.photo.all$PathShort,levels=c("PR","SR","PS II","PS I","AnOx PS"))
#
# mr.photo.all$MethShort<-mr.photo.all$Method
# mr.photo.all$MethShort[(mr.photo.all$MethShort) == "Bacterial Rhodopsin"] <- "Bac Rhod"
# mr.photo.all$MethShort[(mr.photo.all$MethShort) == "Anoxygenic Photosynthesis"] <- "AnOx PS"
# mr.photo.all$MethShort[(mr.photo.all$MethShort) == "Oxygenic Photosynthesis"] <- "Ox PS"
#
# mr.photo.all$Method<-factor(mr.photo.all$Method,levels=c("Bacterial Rhodopsin","Oxygenic Photosynthesis","Anoxygenic Photosynthesis"))
# mr.photo.all$MethShort<-factor(mr.photo.all$MethShort,levels=c("Bac Rhod","Ox PS","AnOx PS"))
#
# unique(mr.photo.all$Phototrophy)
# mr.photo.all$Phototrophy<-factor(mr.photo.all$Phototrophy,levels=c("Hetero","Auto"))
#
# mr.photo.all$KO_Function.KEGG = factor(mr.photo.all$KO_Function.KEGG, levels=unique(mr.photo.all$KO_Function.KEGG[order(mr.photo.all$Phototrophy)]), ordered=TRUE)
#
# head(mr.photo.all)
#
# # convert all 0s to NAs so they appear gray in heat map
# mr.photo.all$MR_SumCovPerKO_NA<-ifelse(mr.photo.all$MR_SumCovPerKO==0,NA,mr.photo.all$MR_SumCovPerKO)
# head(mr.photo.all)
#
# # For heatmap color gradient
# max(mr.photo.all$MR_SumCovPerKO, na.rm=TRUE)
# median(mr.photo.all$MR_SumCovPerKO, na.rm=TRUE)
# min(mr.photo.all$MR_SumCovPerKO, na.rm=TRUE)
#
# # Figures below
# # by SampleID
#
# photo.hm1a<-ggplot(mr.photo.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Phototrophy Metabolism in Salton Sea Dust Metagenomes",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))
#
# ggsave(photo.hm1a,filename = "figures/MGM_Figs/Contigs/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_SampID_by_Function_heatmap.png", width=18, height=27, dpi=600,create.dir = TRUE)
#
# photo.hm1b<-ggplot(mr.photo.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.15) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Phototrophy Metabolism in Salton Sea Dust Metagenomes",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Phototrophy~.,scales="free_y", space = "free")
#
# ggsave(photo.hm1b,filename = "figures/MGM_Figs/Contigs/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_SampID_by_Function_Phototrophy_heatmap.png", width=17, height=15, dpi=600,create.dir = TRUE)
#
# photo.hm1c<-ggplot(mr.photo.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Phototrophy Metabolism in Salton Sea Dust Metagenomes",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~Site,scales="free_x", space = "free")
#
# ggsave(photo.hm1c,filename = "figures/MGM_Figs/Contigs/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_SampleID_by_Function_Site_best_heatmap.png", width=25, height=27, dpi=600,create.dir = TRUE)
#
# photo.hm1d<-ggplot(mr.photo.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.15) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Phototrophy Metabolism in Salton Sea Dust Metagenomes",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")
#
# ggsave(photo.hm1d,filename = "figures/MGM_Figs/Contigs/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_SampID_by_Function_Phototrophy_System_heatmap2.png", width=17, height=15, dpi=600,create.dir = TRUE)
#
# photo.hm1e<-ggplot(mr.photo.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Phototrophy Metabolism in Salton Sea Dust Metagenomes",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Phototrophy~Site, scales="free", space = "free")
#
# ggsave(photo.hm1e,filename = "figures/MGM_Figs/Contigs/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_SampleID_by_Function_Phototrophy_Site_best_heatmap.png", width=25, height=15, dpi=600,create.dir = TRUE)
#
# photo.hm1f<-ggplot(mr.photo.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Phototrophy Metabolism in Salton Sea Dust Metagenomes",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~Site, scales="free", space = "free")
#
# ggsave(photo.hm1f,filename = "figures/MGM_Figs/Contigs/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_SampleID_by_Function_Phototrophy_System_Site_best_heatmap2.png", width=25, height=15, dpi=600,create.dir = TRUE)
#
# photo.hm1g<-ggplot(mr.photo.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Phototrophy Metabolism in Salton Sea Dust Metagenomes",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(MethShort~Site, scales="free", space = "free")
#
# ggsave(photo.hm1g,filename = "figures/MGM_Figs/Contigs/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_SampleID_by_Function_SampDate_Phototrophy_Method_best_heatmap.png", width=25, height=15, dpi=600,create.dir = TRUE)
#
# #### Pull Out Sulfur Metabolic Fxns from MR data ####
# ## heatmaps of traits of interest
#
# mgm.mr[1:4,1:4]
#
# # pull out sulfur functions from MR normalized, summed coverages (summed gene coverage per KO)
# sulf.ko<-mgm.mr[,which(colnames(mgm.mr) %in% sulfur.fxns$KO_ID)] # merge MR data w/ S fxns found in contigs from KOFamScan
# sulf.ko$SampleID<-rownames(sulf.ko)
# sulf.ko.melt<-melt(sulf.ko, by="SampleID")
# colnames(sulf.ko.melt)[which(names(sulf.ko.melt) == "variable")] <- "KO_ID"
# colnames(sulf.ko.melt)[which(names(sulf.ko.melt) == "value")] <- "MR_SumCovPerKO"
# head(sulf.ko.melt) #sanity check
#
# mr.sulf.ko<-merge(sulf.ko.melt,sulf.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
# head(mr.sulf.ko)
# colnames(mr.sulf.ko)[which(names(mr.sulf.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
# mr.cov.sum.sulf.ko<-as.data.frame(dcast(mr.sulf.ko, SampleID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###just dcast, nothing is being added here!
# rownames(mr.cov.sum.sulf.ko)<-mr.cov.sum.sulf.ko$SampleID
# mr.cov.sum.sulf.ko[1:4,]
#
# # sanity check
# mr.cov.sum.sulf.ko$`cysJ; sulfite reductase (NADPH) flavoprotein alpha-component [EC:1.8.1.2]`[1:4]
# head(mr.sulf.ko)
#
# #### Sulfur Heat Maps ####
# # see max & mean of summed
# max(mr.cov.sum.sulf.ko[,-1])
# min(mr.cov.sum.sulf.ko[,-1])
# max(mr.cov.sum.sulf.ko[,-1])/2
#
# # first heat map of sulfur KOs
# #heatmap(as.matrix(mr.cov.sum.sulf.ko[,-1]), scale = "none")
#
# colSums(mr.cov.sum.sulf.ko[,-1])
# #mr.cov.sum.sulf.ko2 <- mr.cov.sum.sulf.ko[,which(colSums(mr.cov.sum.sulf.ko[,-1])>10)]
#
# #heatmap(as.matrix(mr.cov.sum.sulf.ko[,-1]), scale = "none")
#
# # prep for ggplot2 heatmap -- using MR data that includes NAs so they are blocked out on heatmaps
# mr.sulf.ko[1:4,]
# mr.sulf.all<-merge(mr.sulf.ko,meta.all.scaled,by="SampleID")
# head(mr.sulf.all)
# mr.sulf.all$PlotID = factor(mr.sulf.all$SampleID, levels=unique(mr.sulf.all$SampleID[order(mr.sulf.all$Site,mr.sulf.all$SampDate)]), ordered=TRUE)
# # mr.sulf.all$SampDate<-gsub("\\."," ",mr.sulf.all$SampDate)
#
# mr.sulf.all$PathShort<-mr.sulf.all$Pathway
# mr.sulf.all$PathShort[(mr.sulf.all$PathShort) == "Dissimilatory Sulfate Redox"] <- "D.SO4 RedOx"
# mr.sulf.all$PathShort[(mr.sulf.all$PathShort) == "Assimilatory Sulfate Reduction"] <- "A.SO4 Red"
# mr.sulf.all$PathShort[(mr.sulf.all$PathShort) == "Multiple Pathways"] <- "Multi Paths"
# mr.sulf.all$PathShort[(mr.sulf.all$PathShort) == "S Disproportionation"] <- "S Disprop."
#
# mr.sulf.all$Pathway<-factor(mr.sulf.all$Pathway,levels=c("Assimilatory Sulfate Reduction","Dissimilatory Sulfate Redox","Multiple Pathways","SOX","S Disproportionation"))
# mr.sulf.all$PathShort<-factor(mr.sulf.all$PathShort,levels=c("A.SO4 Red","D.SO4 RedOx","Multi Paths","SOX","S Disprop."))
#
# mr.sulf.all$PathSpecShort<-mr.sulf.all$PathwaySpecific
# mr.sulf.all$PathSpecShort[(mr.sulf.all$PathSpecShort) == "Dissimilatory Sulfate Redox"] <- "D.SO4 RedOx"
# mr.sulf.all$PathSpecShort[(mr.sulf.all$PathSpecShort) == "Assimilatory Sulfate Reduction"] <- "A.SO4 Red"
# mr.sulf.all$PathSpecShort[(mr.sulf.all$PathSpecShort) == "Multiple Pathways"] <- "Multi"
# mr.sulf.all$PathSpecShort[(mr.sulf.all$PathSpecShort) == "Sulfur Disproportionation"] <- "Dispro"
# mr.sulf.all$PathSpecShort[(mr.sulf.all$PathSpecShort) == "Sulfide Oxidation"] <- "H2S Ox"
# mr.sulf.all$PathSpecShort[(mr.sulf.all$PathSpecShort) == "Sulfite Oxidation"] <- "SO3 Ox"
# mr.sulf.all$PathSpecShort[(mr.sulf.all$PathSpecShort) == "Thiosulfate Oxidation"] <- "S2O3 Ox"
#
# mr.sulf.all$PathwaySpecific<-factor(mr.sulf.all$PathwaySpecific,levels=c("Assimilatory Sulfate Reduction","Dissimilatory Sulfate Redox","Multiple Pathways","SOX","Sulfur Disproportionation","Sulfide Oxidation","Sulfite Oxidation","Thiosulfate Oxidation"))
# mr.sulf.all$PathSpecShort<-factor(mr.sulf.all$PathSpecShort,levels=c("A.SO4 Red","D.SO4 RedOx","Multi","Dispro","H2S Ox","SO3 Ox","S2O3 Ox"))
#
# mr.sulf.all$KO_Function.KEGG = factor(mr.sulf.all$KO_Function.KEGG, levels=unique(mr.sulf.all$KO_Function.KEGG[order(mr.sulf.all$Pathway)]), ordered=TRUE)
#
# head(mr.sulf.all)
#
# # convert all 0s to NAs so they appear gray in heat map
# mr.sulf.all$MR_SumCovPerKO_NA<-ifelse(mr.sulf.all$MR_SumCovPerKO==0,NA,mr.sulf.all$MR_SumCovPerKO)
# head(mr.sulf.all)
#
# # For heatmap color gradient
# max(mr.sulf.all$MR_SumCovPerKO, na.rm=TRUE)
# max(mr.sulf.all$MR_SumCovPerKO, na.rm=TRUE)/2
# min(mr.sulf.all$MR_SumCovPerKO, na.rm=TRUE)
#
# # Figures below
# # by SampleID
#
# sulf.hm1a<-ggplot(mr.sulf.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sulfur Metabolism in Salton Sea Dust Metagenomes",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))
#
# ggsave(sulf.hm1a,filename = "figures/MGM_Figs/Contigs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_SampID_by_Function_heatmap.png", width=18, height=27, dpi=600,create.dir = TRUE)
#
# sulf.hm1a2<-ggplot(mr.sulf.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.15) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sulfur Metabolism in Salton Sea Dust Metagenomes",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Pathway~.,scales="free_y", space = "free")
#
# ggsave(sulf.hm1a2,filename = "figures/MGM_Figs/Contigs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_SampID_by_Function_Pathway_heatmap.png", width=17, height=15, dpi=600,create.dir = TRUE)
#
# sulf.hm1a3<-ggplot(mr.sulf.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sulfur Metabolism in Salton Sea Dust Metagenomes",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~Site,scales="free_x", space = "free")
#
# ggsave(sulf.hm1a3,filename = "figures/MGM_Figs/Contigs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_SampleID_by_Function_Site_best_heatmap.png", width=25, height=27, dpi=600,create.dir = TRUE)
#
# sulf.hm1a4<-ggplot(mr.sulf.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.15) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sulfur Metabolism in Salton Sea Dust Metagenomes",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")
#
# ggsave(sulf.hm1a4,filename = "figures/MGM_Figs/Contigs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_SampID_by_Function_Pathway_heatmap2.png", width=17, height=15, dpi=600,create.dir = TRUE)
#
# sulf.hm1a5<-ggplot(mr.sulf.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sulfur Metabolism in Salton Sea Dust Metagenomes",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Pathway~Site, scales="free", space = "free")
#
# ggsave(sulf.hm1a5,filename = "figures/MGM_Figs/Contigs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_SampleID_by_Function_SampDate_Pathway_heatmap.png", width=25, height=15, dpi=600,create.dir = TRUE)
#
# sulf.hm1a6<-ggplot(mr.sulf.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sulfur Metabolism in Salton Sea Dust Metagenomes",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~Site, scales="free", space = "free")
#
# ggsave(sulf.hm1a6,filename = "figures/MGM_Figs/Contigs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_SampleID_by_Function_Site_Pathway_best_heatmap2.png", width=25, height=15, dpi=600,create.dir = TRUE)
#
# sulf.hm1a7<-ggplot(mr.sulf.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sulfur Metabolism in Salton Sea Dust Metagenomes",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathwaySpecific~Site, scales="free", space = "free")
#
# ggsave(sulf.hm1a7,filename = "figures/MGM_Figs/Contigs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_SampleID_by_Function_SampDate_PathwaySpecific_best_heatmap.png", width=25, height=20, dpi=600,create.dir = TRUE)
#
# sulf.hm1a8<-ggplot(mr.sulf.all, aes(PlotID, KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sulfur Metabolism in Salton Sea Dust Metagenomes",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathSpecShort~Site, scales="free", space = "free")
#
# ggsave(sulf.hm1a8,filename = "figures/MGM_Figs/Contigs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_SampleID_by_Function_SampDate_PathwaySpecific_best_heatmap2.png", width=25, height=20, dpi=600,create.dir = TRUE)
#
#
# # sulf.hm1c<-ggplot(mr.sulf.all, aes(interaction(SampDate,Depth_m), KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
# #   geom_tile(colour="white",size=0.25) +
# #   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sulfur Metabolism in Salton Sea Dust Metagenomes",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
# #   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
# #         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
# #         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
# #   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Pathway~.,scales="free_y", space = "free")
# #
# # ggsave(sulf.hm1c,filename = "figures/MGM_Figs/Contigs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_SampDate_Depth_by_Function_Pathway_heatmap.png", width=15, height=18, dpi=600,create.dir = TRUE)
# #
# # sulf.hm1c2<-ggplot(mr.sulf.all, aes(interaction(SampDate,Depth_m), KO_Function.KEGG, fill=MR_SumCovPerKO_NA)) +
# #   geom_tile(colour="white",size=0.25) +
# #   scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Sulfur Metabolism in Salton Sea Dust Metagenomes",subtitle="Using MR-Transformed, Gene Coverage Summed by KO",fill="mr Coverage Per KO") +
# #   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
# #         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
# #         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
# #   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")
# #
# # ggsave(sulf.hm1c2,filename = "figures/MGM_Figs/Contigs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_SampDate_Depth_by_Function_Pathway_heatmap2.png", width=15, height=18, dpi=600,create.dir = TRUE)
#### TempShock Fxns PCoA ####
## PCOA with MR normalized data first
# calculate our Euclidean distance matrix using MR data

# pull out temp shock response functions from MR normalized, summed coverages (summed coverage per KO)
tempshock.ko<-mgm.mr[,which(colnames(mgm.mr) %in% tempshock.fxns$KO_ID)] # merge MRN data w/ S fxns found in contigs from KOFamScan
tempshock.ko$SampleID<-rownames(tempshock.ko)
tempshock.ko.melt<-melt(tempshock.ko, by="SampleID")
colnames(tempshock.ko.melt)[which(names(tempshock.ko.melt) == "variable")] <- "KO_ID"
colnames(tempshock.ko.melt)[which(names(tempshock.ko.melt) == "value")] <- "MR_SumCovPerKO"
head(tempshock.ko.melt) #sanity check

mr.tempshock.ko<-merge(tempshock.ko.melt,tempshock.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(mr.tempshock.ko)
colnames(mr.tempshock.ko)[which(names(mr.tempshock.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
mr.cov.sum.tempshock.ko<-as.data.frame(dcast(mr.tempshock.ko, SampleID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###
rownames(mr.cov.sum.tempshock.ko)<-mr.cov.sum.tempshock.ko$SampleID
mr.cov.sum.tempshock.ko[1:4,1:4]

tempshock.euc.mr_dist <- dist(mr.cov.sum.tempshock.ko[,-1], method = "euclidean")

# creating our hierarcical clustering dendrogram
tempshock.euc.mr_clust <- hclust(tempshock.euc.mr_dist, method="ward.D2")

# let's make it a little nicer...
tempshock.euc.mr_dend <- as.dendrogram(tempshock.euc.mr_clust, hang=0.2)
tempshock.dend_cols <- as.character(meta.all.scaled$SampDate_Color[order.dendrogram(tempshock.euc.mr_dend)])
labels_colors(tempshock.euc.mr_dend) <- tempshock.dend_cols

plot(tempshock.euc.mr_dend, ylab="mr Euclidean Distance",cex = 0.5) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 3)
# legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c("#ef781c","#03045e","#059c3f"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# let's use our Euclidean distance matrix from before
tempshock.pcoa.mr <- pcoa(tempshock.euc.mr_dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
##save.image("data/Metagenomesssd_mr.euc.dist1_3.7.23.Rdata")

# The proportion of variances explained is in its element values$Relative_eig
tempshock.pcoa.mr$values

# extract principal coordinates
tempshock.pcoa.mr.vectors<-data.frame(tempshock.pcoa.mr$vectors)
tempshock.pcoa.mr.vectors$SampleID<-rownames(tempshock.pcoa.mr$vectors)

# merge pcoa coordinates w/ metadata
tempshock.pcoa.mr.meta<-merge(tempshock.pcoa.mr.vectors, meta.all.scaled, by.x="SampleID", by.y="SampleID")

head(tempshock.pcoa.mr.meta)

tempshock.pcoa.mr$values # pull out Relative (Relative_eig) variation % to add to axes labels

# create PCoA ggplot fig
ggplot(tempshock.pcoa.mr.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(SampDate),shape=Site), size=4)+theme_bw()+
  labs(title="PCoA: Temperature Shock Response Functions in Salton Sea Dust",subtitle="Using MR normalized, Summed Gene Coverage per KO Function",xlab="PC1 [50.41%]", ylab="PC2 [54.27%]",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=27),axis.title.y = element_text(size=27),legend.title.align=0.5, legend.title = element_text(size=27),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+ scale_shape_manual(name="Site",values = c(0,1,16,15))+
  scale_color_manual(name ="Sample Date",values=unique(tempshock.pcoa.mr.meta$SampDate_Color[order(tempshock.pcoa.mr.meta$SampDate)]))



#### UV-Damaged DNA Repair Fxns PCoA ####
## PCOA with MR normalized data first
# calculate our Euclidean distance matrix using MR data

# pull out UV-damanged DNA repair functions from MR normalized, summed coverages (summed coverage per KO)
UVrep.ko<-mgm.mr[,which(colnames(mgm.mr) %in% UVrep.fxns$KO_ID)] # merge MRN data w/ S fxns found in contigs from KOFamScan
UVrep.ko$SampleID<-rownames(UVrep.ko)
UVrep.ko.melt<-melt(UVrep.ko, by="SampleID")
colnames(UVrep.ko.melt)[which(names(UVrep.ko.melt) == "variable")] <- "KO_ID"
colnames(UVrep.ko.melt)[which(names(UVrep.ko.melt) == "value")] <- "MR_SumCovPerKO"
head(UVrep.ko.melt) #sanity check

mr.UVrep.ko<-merge(UVrep.ko.melt,UV.DNA.rep.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(mr.UVrep.ko)
colnames(mr.UVrep.ko)[which(names(mr.UVrep.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
mr.cov.sum.UVrep.ko<-as.data.frame(dcast(mr.UVrep.ko, SampleID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###
rownames(mr.cov.sum.UVrep.ko)<-mr.cov.sum.UVrep.ko$SampleID
mr.cov.sum.UVrep.ko[1:4,1:4]

UVrep.euc.mr_dist <- dist(mr.cov.sum.UVrep.ko[,-1], method = "euclidean")

# creating our hierarcical clustering dendrogram
UVrep.euc.mr_clust <- hclust(UVrep.euc.mr_dist, method="ward.D2")

# let's make it a little nicer...
UVrep.euc.mr_dend <- as.dendrogram(UVrep.euc.mr_clust, hang=0.2)
UVrep.dend_cols <- as.character(meta.all.scaled$SampDate_Color[order.dendrogram(UVrep.euc.mr_dend)])
labels_colors(UVrep.euc.mr_dend) <- UVrep.dend_cols

plot(UVrep.euc.mr_dend, ylab="mr Euclidean Distance",cex = 0.5) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 3)
# legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c("#ef781c","#03045e","#059c3f"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# let's use our Euclidean distance matrix from before
UVrep.pcoa.mr <- pcoa(UVrep.euc.mr_dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
##save.image("data/Metagenomesssd_mr.euc.dist1_3.7.23.Rdata")

# The proportion of variances explained is in its element values$Relative_eig
UVrep.pcoa.mr$values

# extract principal coordinates
UVrep.pcoa.mr.vectors<-data.frame(UVrep.pcoa.mr$vectors)
UVrep.pcoa.mr.vectors$SampleID<-rownames(UVrep.pcoa.mr$vectors)

# merge pcoa coordinates w/ metadata
UVrep.pcoa.mr.meta<-merge(UVrep.pcoa.mr.vectors, meta.all.scaled, by.x="SampleID", by.y="SampleID")

head(UVrep.pcoa.mr.meta)

UVrep.pcoa.mr$values # pull out Relative (Relative_eig) variation % to add to axes labels

# create PCoA ggplot fig
ggplot(UVrep.pcoa.mr.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(SampDate),shape=Site), size=4)+theme_bw()+
  labs(title="PCoA: UV-Damaged DNA Repair Genes in Salton Sea Dust",subtitle="Using MR normalized, Summed Gene Coverage per KO Function",xlab="PC1 [50.41%]", ylab="PC2 [54.27%]",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=27),axis.title.y = element_text(size=27),legend.title.align=0.5, legend.title = element_text(size=27),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+ scale_shape_manual(name="Site",values = c(0,1,16,15)) +
  scale_color_manual(name ="Sample Date",values=unique(UVrep.pcoa.mr.meta$SampDate_Color[order(UVrep.pcoa.mr.meta$SampDate)]))


#### Osmoprotectant Fxns PCoA ####
## PCOA with MR normalized data first
# calculate our Euclidean distance matrix using MR data

# pull out UV-damanged DNA repair functions from MR normalized, summed coverages (summed coverage per KO)
osmo.ko<-mgm.mr[,which(colnames(mgm.mr) %in% osmo.fxns$KO_ID)] # merge MRN data w/ S fxns found in contigs from KOFamScan
osmo.ko$SampleID<-rownames(osmo.ko)
osmo.ko.melt<-melt(osmo.ko, by="SampleID")
colnames(osmo.ko.melt)[which(names(osmo.ko.melt) == "variable")] <- "KO_ID"
colnames(osmo.ko.melt)[which(names(osmo.ko.melt) == "value")] <- "MR_SumCovPerKO"
head(osmo.ko.melt) #sanity check

mr.osmo.ko<-merge(osmo.ko.melt,osmo.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(mr.osmo.ko)
colnames(mr.osmo.ko)[which(names(mr.osmo.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
mr.cov.sum.osmo.ko<-as.data.frame(dcast(mr.osmo.ko, SampleID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###
rownames(mr.cov.sum.osmo.ko)<-mr.cov.sum.osmo.ko$SampleID
mr.cov.sum.osmo.ko[1:4,1:4]

osmo.euc.mr_dist <- dist(mr.cov.sum.osmo.ko[,-1], method = "euclidean")

# creating our hierarcical clustering dendrogram
osmo.euc.mr_clust <- hclust(osmo.euc.mr_dist, method="ward.D2")

# let's make it a little nicer...
osmo.euc.mr_dend <- as.dendrogram(osmo.euc.mr_clust, hang=0.2)
osmo.dend_cols <- as.character(meta.all.scaled$SampDate_Color[order.dendrogram(osmo.euc.mr_dend)])
labels_colors(osmo.euc.mr_dend) <- osmo.dend_cols

plot(osmo.euc.mr_dend, ylab="mr Euclidean Distance",cex = 0.5) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 3)
# legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c("#ef781c","#03045e","#059c3f"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# let's use our Euclidean distance matrix from before
osmo.pcoa.mr <- pcoa(osmo.euc.mr_dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
##save.image("data/Metagenomesssd_mr.euc.dist1_3.7.23.Rdata")

# The proportion of variances explained is in its element values$Relative_eig
osmo.pcoa.mr$values

# extract principal coordinates
osmo.pcoa.mr.vectors<-data.frame(osmo.pcoa.mr$vectors)
osmo.pcoa.mr.vectors$SampleID<-rownames(osmo.pcoa.mr$vectors)

# merge pcoa coordinates w/ metadata
osmo.pcoa.mr.meta<-merge(osmo.pcoa.mr.vectors, meta.all.scaled, by.x="SampleID", by.y="SampleID")

head(osmo.pcoa.mr.meta)

osmo.pcoa.mr$values # pull out Relative (Relative_eig) variation % to add to axes labels

# create PCoA ggplot fig
ggplot(osmo.pcoa.mr.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(SampDate),shape=Site), size=4)+theme_bw()+
  labs(title="PCoA: Osmoprotectant Transport Genes in Salton Sea Dust",subtitle="Using MR normalized, Summed Gene Coverage per KO Function",xlab="PC1 [50.41%]", ylab="PC2 [54.27%]",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=27),axis.title.y = element_text(size=27),legend.title.align=0.5, legend.title = element_text(size=27),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+ scale_shape_manual(name="Site",values = c(0,1,16,15))+
  scale_color_manual(name ="Sample Date",values=unique(osmo.pcoa.mr.meta$SampDate_Color[order(osmo.pcoa.mr.meta$SampDate)]))

#### ARGs PCoA ####
## PCOA with MR normalized data first
# calculate our Euclidean distance matrix using MR data

# pull out ARG functions from MR normalized, summed coverages (summed coverage per KO)
# ARG.ko<-mgm.mr[,which(colnames(mgm.mr) %in% ARG.fxns$KO_ID)] # merge MRN data w/ S fxns found in contigs from KOFamScan
# ARG.ko$SampleID<-rownames(ARG.ko)
ARG.ko.melt<-melt(ARG.ko, by="SampleID")
colnames(ARG.ko.melt)[which(names(ARG.ko.melt) == "variable")] <- "KO_ID"
colnames(ARG.ko.melt)[which(names(ARG.ko.melt) == "value")] <- "MR_SumCovPerKO"
head(ARG.ko.melt) #sanity check

mr.ARG.ko<-merge(ARG.ko.melt,ARG.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(mr.ARG.ko)
colnames(mr.ARG.ko)[which(names(mr.ARG.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
mr.cov.sum.ARG.ko<-as.data.frame(dcast(mr.ARG.ko, SampleID~KO_Function.KEGG, value.var="MR_SumCovPerKO", fun.aggregate=sum)) ###
rownames(mr.cov.sum.ARG.ko)<-mr.cov.sum.ARG.ko$SampleID
mr.cov.sum.ARG.ko[1:4,1:4]

ARG.euc.mr_dist <- dist(mr.cov.sum.ARG.ko[,-1], method = "euclidean")

# creating our hierarcical clustering dendrogram
ARG.euc.mr_clust <- hclust(ARG.euc.mr_dist, method="ward.D2")

# let's make it a little nicer...
ARG.euc.mr_dend <- as.dendrogram(ARG.euc.mr_clust, hang=0.2)
ARG.dend_cols <- as.character(meta.all.scaled$SampDate_Color[order.dendrogram(ARG.euc.mr_dend)])
labels_colors(ARG.euc.mr_dend) <- ARG.dend_cols

par(mar=c(1,1,1,1))
plot(ARG.euc.mr_dend, ylab="mr Euclidean Distance",cex = 0.5) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 3)
# legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c("#ef781c","#03045e","#059c3f"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# let's use our Euclidean distance matrix from before
ARG.pcoa.mr <- pcoa(ARG.euc.mr_dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
##save.image("data/Metagenomesssd_mr.euc.dist1_3.7.23.Rdata")

# The proportion of variances explained is in its element values$Relative_eig
ARG.pcoa.mr$values

# extract principal coordinates
ARG.pcoa.mr.vectors<-data.frame(ARG.pcoa.mr$vectors)
ARG.pcoa.mr.vectors$SampleID<-rownames(ARG.pcoa.mr$vectors)

# merge pcoa coordinates w/ metadata
ARG.pcoa.mr.meta<-merge(ARG.pcoa.mr.vectors, meta.all.scaled, by.x="SampleID", by.y="SampleID")

head(ARG.pcoa.mr.meta)

ARG.pcoa.mr$values # pull out Relative (Relative_eig) variation % to add to axes labels

# create PCoA ggplot fig
ggplot(ARG.pcoa.mr.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(SampDate),shape=Site), size=4)+theme_bw()+
  labs(title="PCoA: ARGs in Salton Sea Dust",subtitle="Using MR normalized, Summed Gene Coverage per KO Function",xlab="PC1", ylab="PC2",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=27),axis.title.y = element_text(size=27),legend.title.align=0.5, legend.title = element_text(size=27),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+ scale_shape_manual(name="Site",values = c(0,1,16,15))+
  scale_color_manual(name ="Sample Date",values=unique(ARG.pcoa.mr.meta$SampDate_Color[order(ARG.pcoa.mr.meta$SampDate)]))


#### Homogeneity of Variance - TempShock Fxns by Groups ####
## betadisper to look at homogeneity of group dispersions (aka variance) when considering multiple variables
# multivariate analogue to Levene's test of homogeneity of variances
# program finds spatial median or centroid of the group, & compare distances of group to centroid/spatial median via ANOVA

#While PERMANOVA tests differences in group means (analogous to MANOVA),
## a related test called PERMDISP can be used to evaluate homogeneity of group dispersion
#(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
## * need a distance matrix!

# calculate our Euclidean distance matrix using MR data
tempshock.euc.mr_dist

## betadisper to look at within Site variance
# first by compare dispersions by site
TS.disper<-betadisper(tempshock.euc.mr_dist, meta.all.scaled$Site)
TS.disper

permutest(TS.disper, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)

anova(TS.disper) # p = 0.7604 --> accept the Null H, spatial medians ARE NOT significantly difference across sample dates

TukeyHSD(TS.disper) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

# First make sure your data frames you're comparing are in the same exact order!!
rownames(mr.cov.sum.tempshock.ko[,-1]) %in% rownames(meta.all.scaled)

TS.pnova1<-adonis2(mr.cov.sum.tempshock.ko[,-1] ~ Site,data=meta.all.scaled,method = "euclidean",by="terms",permutations= 10000)
TS.pnova1 # p-value = 0.8397

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

TS.dist = (vegdist(mr.cov.sum.tempshock.ko[,-1], "euclidean", na.rm = TRUE)) #distance matrix using Bray's dissimilarity index for trait distribution (traits of interest only)
pair.mod1<-pairwise.adonis(TS.dist,meta.all.scaled$Site, p.adjust.m='bonferroni',perm=9999) # shows us variation for each sample to see which ones are different
pair.mod1

# Visualize dispersions
plot(TS.disper,main = "Centroids and Dispersion based on Aitchison Distance (mr Data)")
dev.off()

boxplot(TS.disper,xlab="Sample Collection Date", main = "Distance to Centroid by Category (mr Data)", sub="Based on Aitchison Distance")
dev.off()

## betadisper to look at within SampDate variance
# first by compare dispersions by site
TS.disper1<-betadisper(tempshock.euc.mr_dist, meta.all.scaled$SampDate)
TS.disper1

permutest(TS.disper1, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)

anova(TS.disper1) # p = 0.7068 --> accept the Null H, spatial medians ARE NOT significantly difference across sample dates

TukeyHSD(TS.disper1) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

# First make sure your data frames you're comparing are in the same exact order!!
rownames(mr.cov.sum.tempshock.ko[,-1]) %in% rownames(meta.all.scaled)

TS.pnova2<-adonis2(mr.cov.sum.tempshock.ko[,-1] ~ SampDate,data=meta.all.scaled,method = "euclidean",by="terms",permutations= 10000)
TS.pnova2 # p-value = 0.6371

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

TS.dist = (vegdist(mr.cov.sum.tempshock.ko[,-1], "euclidean", na.rm = TRUE)) #distance matrix using Bray's dissimilarity index for trait distribution (traits of interest only)
pair.mod2<-pairwise.adonis(TS.dist,meta.all.scaled$SampDate, p.adjust.m='bonferroni',perm=9999) # shows us variation for each sample to see which ones are different
pair.mod2

# Visualize dispersions
plot(TS.disper1,main = "Centroids and Dispersion based on Aitchison Distance (mr Data)")

boxplot(TS.disper1,xlab="Sample Collection Date", main = "Distance to Centroid by Category (mr Data)", sub="Based on Aitchison Distance")


#### Homogeneity of Variance - UV Damage DNA Repair by Groups ####
## betadisper to look at homogeneity of group dispersions (aka variance) when considering multiple variables
# multivariate analogue to Levene's test of homogeneity of variances
# program finds spatial median or centroid of the group, & compare distances of group to centroid/spatial median via ANOVA

#While PERMANOVA tests differences in group means (analogous to MANOVA),
## a related test called PERMDISP can be used to evaluate homogeneity of group dispersion
#(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
## * need a distance matrix!

# calculate our Euclidean distance matrix using MR data
UVrep.euc.mr_dist

## betadisper to look at within Site variance
# first by compare dispersions by site
UV.disper<-betadisper(UVrep.euc.mr_dist, meta.all.scaled$Site)
UV.disper

permutest(UV.disper, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)

anova(UV.disper) # p = 0.8435 --> accept the Null H, spatial medians ARE NOT significantly difference across sample dates

TukeyHSD(UV.disper) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

# First make sure your data frames you're comparing are in the same exact order!!
rownames(mr.cov.sum.UVrep.ko[,-1]) %in% rownames(meta.all.scaled)

UV.pnova1<-adonis2(mr.cov.sum.UVrep.ko[,-1] ~ Site,data=meta.all.scaled,method = "euclidean",by="terms",permutations= 10000)
UV.pnova1 # p-value = 0.6462

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

UV.dist = (vegdist(mr.cov.sum.UVrep.ko[,-1], "euclidean", na.rm = TRUE)) #distance matrix using Bray's dissimilarity index for trait distribution (traits of interest only)
pair.mod1<-pairwise.adonis(UV.dist,meta.all.scaled$Site, p.adjust.m='bonferroni',perm=9999) # shows us variation for each sample to see which ones are different
pair.mod1

# Visualize dispersions
plot(UV.disper,main = "Centroids and Dispersion based on Aitchison Distance (mr Data)")
dev.off()

boxplot(UV.disper,xlab="Sample Collection Date", main = "Distance to Centroid by Category (mr Data)", sub="Based on Aitchison Distance")
dev.off()

## betadisper to look at within SampDate variance
# first by compare dispersions by site
UV.disper1<-betadisper(UVrep.euc.mr_dist, meta.all.scaled$SampDate)
UV.disper1

permutest(UV.disper1, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)

anova(UV.disper1) # p = 0.6685 --> accept the Null H, spatial medians ARE NOT significantly difference across sample dates

TukeyHSD(UV.disper1) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

# First make sure your data frames you're comparing are in the same exact order!!
rownames(mr.cov.sum.UVrep.ko[,-1]) %in% rownames(meta.all.scaled)

UV.pnova2<-adonis2(mr.cov.sum.UVrep.ko[,-1] ~ SampDate,data=meta.all.scaled,method = "euclidean",by="terms",permutations= 10000)
UV.pnova2 # p-value = 0.3497

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

UV.dist = (vegdist(mr.cov.sum.UVrep.ko[,-1], "euclidean", na.rm = TRUE)) #distance matrix using Bray's dissimilarity index for trait distribution (traits of interest only)
pair.mod2<-pairwise.adonis(UV.dist,meta.all.scaled$SampDate, p.adjust.m='bonferroni',perm=9999) # shows us variation for each sample to see which ones are different
pair.mod2

# Visualize dispersions
plot(UV.disper1,main = "Centroids and Dispersion based on Aitchison Distance (mr Data)")

boxplot(UV.disper1,xlab="Sample Collection Date", main = "Distance to Centroid by Category (mr Data)", sub="Based on Aitchison Distance")

#### Homogeneity of Variance - Osmoprotectant Transport by Groups ####
## betadisper to look at homogeneity of group dispersions (aka variance) when considering multiple variables
# multivariate analogue to Levene's test of homogeneity of variances
# program finds spatial median or centroid of the group, & compare distances of group to centroid/spatial median via ANOVA

#While PERMANOVA tests differences in group means (analogous to MANOVA),
## a related test called PERMDISP can be used to evaluate homogeneity of group dispersion
#(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
## * need a distance matrix!

# calculate our Euclidean distance matrix using MR data
osmo.euc.mr_dist

## betadisper to look at within Site variance
# first by compare dispersions by site
osmo.disper<-betadisper(osmo.euc.mr_dist, meta.all.scaled$Site)
osmo.disper

permutest(osmo.disper, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)

anova(osmo.disper) # p = 0.7754 --> accept the Null H, spatial medians ARE NOT significantly difference across sample dates

TukeyHSD(osmo.disper) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

# First make sure your data frames you're comparing are in the same exact order!!
rownames(mr.cov.sum.osmo.ko[,-1]) %in% rownames(meta.all.scaled)

osmo.pnova1<-adonis2(mr.cov.sum.osmo.ko[,-1] ~ Site,data=meta.all.scaled,method = "euclidean",by="terms",permutations= 10000)
osmo.pnova1 # p-value = 0.8313

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

osmo.dist = (vegdist(mr.cov.sum.osmo.ko[,-1], "euclidean", na.rm = TRUE)) #distance matrix using Bray's dissimilarity index for trait distribution (traits of interest only)
pair.mod1<-pairwise.adonis(osmo.dist,meta.all.scaled$Site, p.adjust.m='bonferroni',perm=9999) # shows us variation for each sample to see which ones are different
pair.mod1

# Visualize dispersions
plot(osmo.disper,main = "Centroids and Dispersion based on Aitchison Distance (mr Data)")
dev.off()

boxplot(osmo.disper,xlab="Sample Collection Date", main = "Distance to Centroid by Category (mr Data)", sub="Based on Aitchison Distance")
dev.off()

## betadisper to look at within SampDate variance
# first by compare dispersions by site
osmo.disper1<-betadisper(osmo.euc.mr_dist, meta.all.scaled$SampDate)
osmo.disper1

permutest(osmo.disper1, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)

anova(osmo.disper1) # p = 0.296 --> accept the Null H, spatial medians ARE NOT significantly difference across sample dates

TukeyHSD(osmo.disper1) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

# First make sure your data frames you're comparing are in the same exact order!!
rownames(mr.cov.sum.osmo.ko[,-1]) %in% rownames(meta.all.scaled)

osmo.pnova2<-adonis2(mr.cov.sum.osmo.ko[,-1] ~ SampDate,data=meta.all.scaled,method = "euclidean",by="terms",permutations= 10000)
osmo.pnova2 # p-value = 0.296

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

osmo.dist = (vegdist(mr.cov.sum.osmo.ko[,-1], "euclidean", na.rm = TRUE)) #distance matrix using Bray's dissimilarity index for trait distribution (traits of interest only)
pair.mod2<-pairwise.adonis(osmo.dist,meta.all.scaled$SampDate, p.adjust.m='bonferroni',perm=9999) # shows us variation for each sample to see which ones are different
pair.mod2

# Visualize dispersions
plot(osmo.disper1,main = "Centroids and Dispersion based on Aitchison Distance (mr Data)")

boxplot(osmo.disper1,xlab="Sample Collection Date", main = "Distance to Centroid by Category (mr Data)", sub="Based on Aitchison Distance")

#### Homogeneity of Variance - ARGs by Groups ####
## betadisper to look at homogeneity of group dispersions (aka variance) when considering multiple variables
# multivariate analogue to Levene's test of homogeneity of variances
# program finds spatial median or centroid of the group, & compare distances of group to centroid/spatial median via ANOVA

#While PERMANOVA tests differences in group means (analogous to MANOVA),
## a related test called PERMDISP can be used to evaluate homogeneity of group dispersion
#(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
## * need a distance matrix!

# calculate our Euclidean distance matrix using MR data
ARG.euc.mr_dist

## betadisper to look at within Site variance
# first by compare dispersions by site
ARG.disper<-betadisper(ARG.euc.mr_dist, meta.all.scaled$Site)
ARG.disper

permutest(ARG.disper, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)

anova(ARG.disper) # p = 0.6225 --> accept the Null H, spatial medians ARE NOT significantly difference across sample dates

TukeyHSD(ARG.disper) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

# First make sure your data frames you're comparing are in the same exact order!!
rownames(mr.cov.sum.ARG.ko[,-1]) %in% rownames(meta.all.scaled)

ARG.pnova1<-adonis2(mr.cov.sum.ARG.ko[,-1] ~ Site,data=meta.all.scaled,method = "euclidean",by="terms",permutations= 10000)
ARG.pnova1 # p-value = 0.6522

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

ARG.dist = (vegdist(mr.cov.sum.ARG.ko[,-1], "euclidean", na.rm = TRUE)) #distance matrix using Bray's dissimilarity index for trait distribution (traits of interest only)
pair.mod1<-pairwise.adonis(ARG.dist,meta.all.scaled$Site, p.adjust.m='bonferroni',perm=9999) # shows us variation for each sample to see which ones are different
pair.mod1

# Visualize dispersions
plot(ARG.disper,main = "Centroids and Dispersion based on Aitchison Distance (mr Data)")
dev.off()

par(mar=c(1,1,1,1))
boxplot(ARG.disper,xlab="Sample Collection Date", main = "Distance to Centroid by Category (mr Data)", sub="Based on Aitchison Distance")
dev.off()

## betadisper to look at within SampDate variance
# first by compare dispersions by site
ARG.disper1<-betadisper(ARG.euc.mr_dist, meta.all.scaled$SampDate)
ARG.disper1

permutest(ARG.disper1, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)

anova(ARG.disper1) # p = 0.07151 --> accept the Null H, spatial medians ARE NOT significantly difference across sample dates

TukeyHSD(ARG.disper1) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

# First make sure your data frames you're comparing are in the same exact order!!
rownames(mr.cov.sum.ARG.ko[,-1]) %in% rownames(meta.all.scaled)

ARG.pnova2<-adonis2(mr.cov.sum.ARG.ko[,-1] ~ SampDate,data=meta.all.scaled,method = "euclidean",by="terms",permutations= 10000)
ARG.pnova2 # p-value = 0.1791

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

ARG.dist = (vegdist(mr.cov.sum.ARG.ko[,-1], "euclidean", na.rm = TRUE)) #distance matrix using Bray's dissimilarity index for trait distribution (traits of interest only)
pair.mod2<-pairwise.adonis(ARG.dist,meta.all.scaled$SampDate, p.adjust.m='bonferroni',perm=9999) # shows us variation for each sample to see which ones are different
pair.mod2

# Visualize dispersions
par(mar=c(1,1,1,1))
plot(ARG.disper1,main = "Centroids and Dispersion based on Aitchison Distance (mr Data)")

par(mar=c(1,1,1,1))
boxplot(ARG.disper1,xlab="Sample Collection Date", main = "Distance to Centroid by Category (mr Data)", sub="Based on Aitchison Distance")


#### Combine Fxn Summed Coverage Means to Save ####

# saving the means of the summed, scaled relative coverages for genes of interest (from across mgms)
mgm.fxn.means<-rbind(LPS.mgm.means,Osmo.mgm.means,QS.mgm.means,
                     Spor.mgm.means,TempShock.mgm.means,UVrep.mgm.means,HK.mgm.means)
mgm.fxn.means<-mgm.fxn.means[order(-mgm.fxn.means$MeanScaledCov),]

write.table(mgm.fxn.means, "data/MetagenomesKOsInContigs_MeanSummedScaledRelativeCoverage.tsv",
            sep="\t", quote=F, col.names=T,row.names=F)

#### Export Global Env for Other Scripts ####
save.image("data/Metagenomes/SSD_MGM_Fxn_BetaDiv.Rdata")
# ^ includes all data combined in object bac.dat.all, ASV table (samples are rows, ASVs are columns), meta.all.scaled, and an ASV count table (where ASVs are rows, not columns)
# Version Information
sessionInfo()