#### Set WD & Load Libraries ####
getwd() # use setwd("path/to/files") if you are not in the right directory
setwd("/Volumes/HLF_SSD/Aronson_Lab_Data/Salton_Sea/SaltonSeaWater")
suppressPackageStartupMessages({ # load packages quietly
  library(phyloseq)
  library(ggplot2)
  library(vegan)
  library(ggpubr)
  #library(scales)
  library(grid)
  library(data.table)
  library(ape)
  library(apeglm)
  library(plyr)
  library(dplyr)
  library(viridis)
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
  library(fitdistrplus)
  library(logspline)
  library(shades)
  library(ALDEx2)
  library(rstatix)
  library(devtools)
  library(decontam)
})

#### Load Data ####
load("data/Metagenomes/Analysis/mgm_analysis.Rdata") # load Rdata to global env

head(meta_scaled)
arsenic.fxns[1:4,1:4]
mgm_fxn.cov_table[1:4,1:4]

#### Relative Abundance of MGM Bins Across Samples, Depths, Seasons ####

head(mapped_meta)
# pull out total mapped reads per sample, results from BWA-mem
total_mapped_samp<-unique(data.frame(SampleID=mapped_meta$SampleID,Total_Mapped_Reads=mapped_meta$Total_Mapped_Reads))

# calculate total mapped reads per bin, and recalculate them to find total mapped reads per taxa
## this is because some bins were mapped to the same genus & species, so I am curious about these specific taxa's relative abundance in the mapped reads
gen.spec_mapped<-dcast(mapped_meta, SampleID~Genus+Species, fun.aggregate = sum, value.var="Bin_Mapped_Reads")
rownames(gen.spec_mapped)<-gen.spec_mapped$SampleID
head(gen.spec_mapped)
gen.spec_mapped2<-merge(gen.spec_mapped,total_mapped_samp, by="SampleID")
head(gen.spec_mapped2)
rownames(gen.spec_mapped2)<-gen.spec_mapped2$SampleID

# divide total mapped reads per genus&species by total mapped reads per sample to get relative abundance of reads for specific taxa across mgms
gen.sp_mapped_RA<-gen.spec_mapped2[,-c(1,ncol(gen.spec_mapped2))]/gen.spec_mapped2[,ncol(gen.spec_mapped2)]
head(gen.sp_mapped_RA)
png('figures/MGM_Figs/SSW_MGM_Genus_RelAb_Reads_4.8.23.png',width = 2200, height = 2200, res=100)
heatmap(as.matrix(t(gen.sp_mapped_RA)), scale = "none")
dev.off()

# melt down relativized ASV counts to merge with metadata
gen.sp_mapped_RA$SampleID<-rownames(gen.sp_mapped_RA)
gs_map_RA.m<-melt(gen.sp_mapped_RA)

head(gs_map_RA.m)
colnames(gs_map_RA.m)[which(names(gs_map_RA.m) == "variable")] <- "Genus_species"
colnames(gs_map_RA.m)[which(names(gs_map_RA.m) == "value")] <- "RelAb"
head(gs_map_RA.m) ## relative abundance based on sum of counts by class!

gs_map_RA.meta<-merge(gs_map_RA.m,meta_scaled, by="SampleID")
head(gs_map_RA.meta) ## relative abundance based on sum of counts by class!

# find the midpoint of RelAb
max(gs_map_RA.meta$RelAb)
min(gs_map_RA.meta$RelAb)

gs_map_RA.meta<-gs_map_RA.meta[with(gs_map_RA.meta, order(gs_map_RA.meta$SampDate, gs_map_RA.meta$Depth_m)),]
gs_map_RA.meta$SampleID <- factor(gs_map_RA.meta$SampleID, levels=unique(gs_map_RA.meta$SampleID[order(gs_map_RA.meta$SampDate, gs_map_RA.meta$Depth_m)]))

bin.hm1<-ggplot(gs_map_RA.meta, aes(gs_map_RA.meta$SampleID[order(gs_map_RA.meta$SampDate, gs_map_RA.meta$Depth_m)], Genus_species, fill= RelAb)) +geom_tile()+scale_fill_gradient2(low="skyblue",mid="white",high="orange",midpoint=0.035)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample ID", y="MAG Bin Assignments", title="Relative Abundance of MAG Bins by Sample",fill="Read Relative Abundance")+scale_x_discrete(expand = c(0,0))

ggsave(bin.hm1,filename = "figures/MGM_Figs/Heatmap_MGM_Bins_RelAb_bySample_4.10.23.png", width=12, height=10, dpi=600)


tsum1<-ggplot(mapped_meta, aes(Genus, RelAb_Map)) +
  geom_jitter(aes(color=as.numeric(as.character(Depth_m))), size=2, width=0.15, height=0) +
  scale_colour_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Genera of MAG Bins", y="Read Relative Abundance", title="MAG Bins (Genera) by Depth",color="Depth (m)")

ggsave(tsum1,filename = "figures/MGM_Figs/TaxaSummary_MGM_Bins_RelAb_Depth_4.10.23.png", width=12, height=10, dpi=600)

tsum2<-ggplot(mapped_meta, aes(Genus, RelAb_Map)) +
  geom_jitter(aes(color=as.numeric(as.character(Depth_m)),shape=SampDate), size=2, width=0.15, height=0) +
  scale_colour_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Genera of MAG Bins", y="Read Relative Abundance", title="MAG Bins (Genera) by Depth & Sample Date",color="Depth (m)", shape="Sampling Date")

ggsave(tsum2,filename = "figures/MGM_Figs/TaxaSummary_MGM_Bins_RelAb_Depth_Date_4.10.23.png", width=12, height=10, dpi=600)

#### Taxonomic Beta Diversity - Read Relative Abundance data ####
mgm.clr[1:4,1:4] # sample IDs are rows, genes are columns
mgm_fxn.counts_table[1:4,1:4] # sanity check

# check rownames of CLR & VST transformed feature count data & metadata
rownames(mgm.clr) %in% rownames(meta_scaled)

## PCOA with CLR transformed data first
# calculate our Euclidean distance matrix using CLR data
mgm.euc_dist.clr <- dist(mgm.clr, method = "euclidean")

# creating our hierarcical clustering dendrogram
mgm.euc_clust <- hclust(mgm.euc_dist.clr, method="ward.D2")

# let's make it a little nicer...
mgm.euc_dend <- as.dendrogram(mgm.euc_clust, hang=0.2)
mgm.dend_cols <- as.character(mgm_meta$SampDate_Color[order.dendrogram(mgm.euc_dend)])
labels_colors(mgm.euc_dend) <- mgm.dend_cols

plot(mgm.euc_dend, ylab="CLR Euclidean Distance",cex = 0.5) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 3)
legend("topright",legend = c("June 2021","August 2021","December 2021","April 2022"),cex=.8,col = c( "#26547c","#36ab57","#32cbff","#ff6f00"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# let's use our Euclidean distance matrix from before
mgm.pcoa.clr <- pcoa(mgm.euc_dist.clr) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
##save.image("data/ssw_clr.euc.dist1_3.7.23.Rdata")

# The proportion of variances explained is in its element values$Relative_eig
mgm.pcoa.clr$values

# extract principal coordinates
mgm.pcoa.clr.vectors<-data.frame(mgm.pcoa.clr$vectors)
mgm.pcoa.clr.vectors$SampleID<-rownames(mgm.pcoa.clr$vectors)

# merge pcoa coordinates w/ metadata
mgm.pcoa.clr.meta<-merge(mgm.pcoa.clr.vectors, mgm_meta, by.x="SampleID", by.y="SampleID")
mgm.pcoa.clr.meta$SampleMonth
mgm.pcoa.clr.meta$SampDate

head(mgm.pcoa.clr.meta)

head(mgm.pcoa.clr.meta$values) # pull out Relative (Relative_eig) variation % to add to axes labels

# create PCoA ggplot fig
pcoa1<-ggplot(mgm.pcoa.clr.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(SampDate)), size=4)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Seawater",subtitle="Using Centered-Log Ratio Feature Data",xlab="PC1 [41.14%]", ylab="PC2 [9.04%]",color="Sample Type")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Type",values=unique(mgm.pcoa.clr.meta$SampDate_Color[order(mgm.pcoa.clr.meta$SampDate)]),labels=c("June.2021"="June 2021","August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("PC1 [23.56%]") + ylab("PC2 [20.45%]")

ggsave(pcoa1,filename = "figures/MGM_Figs/SSW_MGM_pcoa_CLR_sampdate.png", width=12, height=10, dpi=600)

# sample month shape, depth color
pcoa2<-ggplot(mgm.pcoa.clr.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCoA: Metagenome Functions in Salton Seawater",subtitle="Using Centered-Log Ratio Feature Data",xlab="PC1", ylab="PC2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("June 2021","August 2021","December 2021","April 2022"),name="Sample Date") +
  xlab("PC1 [23.56%]") + ylab("PC2 [20.45%]")

ggsave(pcoa2,filename = "figures/MGM_Figs/SSW_MGM_pcoa_CLR.traits_depth.png", width=12, height=10, dpi=600)

## betadisper to look at within group variance

# first by sampling date
mgm.disper1<-betadisper(mgm.euc_dist.clr, mgm_meta$SampDate)
mgm.disper1

permutest(mgm.disper1, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#               June.2021 August.2021 December.2021 April.2022
#June.2021                 0.13000000    0.03700000      0.010
#August.2021   0.12533884                0.50300000      0.189
#December.2021 0.00835232  0.56548064                    0.038
#April.2022    0.00043402  0.20006152    0.00656846

anova(mgm.disper1) # p = 0.02656 --> reject the Null H, spatial medians ARE significantly difference across sample dates

TukeyHSD(mgm.disper1) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other

#                             diff        lwr       upr     p adj
#August.2021-June.2021      3.611143 -0.4208285 7.6431151 0.0784173 .
#December.2021-June.2021    4.511716  0.4797445 8.5436881 0.0302728 *
#April.2022-June.2021       1.576242 -2.4557301 5.6082134 0.5941200
#December.2021-August.2021  0.900573 -2.7057322 4.5068782 0.8404583
#April.2022-August.2021    -2.034902 -5.6412069 1.5714035 0.3208184
#April.2022-December.2021  -2.935475 -6.5417798 0.6708306 0.1118431

# Visualize dispersions
png('figures/MGM_Figs/SSW_MGM_pcoa_clr_betadispersion_sampledate.png',width = 700, height = 600, res=100)
plot(mgm.disper1,main = "Centroids and Dispersion based on Aitchison Distance (CLR Data)", col=colorset1$SampDate_Color)
dev.off()

png('figures/MGM_Figs/SSW_MGM_boxplot_clr_centroid_distance_sampledate.png',width = 700, height = 600, res=100)
boxplot(mgm.disper1,xlab="Sample Collection Date", main = "Distance to Centroid by Category (CLR Data)", sub="Based on Aitchison Distance", col=colorset1$SampDate_Color)
dev.off()

# What about between sampling depths?
mgm.disper2<-betadisper(mgm.euc_dist.clr, mgm_meta$Depth_m)
mgm.disper2

permutest(mgm.disper2, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons

anova(mgm.disper2) # p = 0.6044 --> accept the Null H, spatial medians are NOT significantly difference across sample dates

TukeyHSD(mgm.disper2) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other
#         diff       lwr      upr     p adj
#5-0   2.3886761 -4.507223 9.284575 0.6031569
#10-0  1.9710200 -4.924879 8.866919 0.7037811
#10-5 -0.4176561 -6.802018 5.966706 0.9809665

colfunc <- colorRampPalette(c("red", "blue"))
colfunc(3)

# Visualize dispersions
png('figures/MGM_Figs/ssw_mgm_pcoa_clr_betadispersion_depth.png',width = 700, height = 600, res=100)
plot(mgm.disper2,main = "Centroids and Dispersion based on Aitchison Distance (CLR Data)", col=colfunc(3))
dev.off()

png('figures/MGM_Figs/ssw_mgm_boxplot_clr_centroid_distance_depth.png',width = 700, height = 600, res=100)
boxplot(mgm.disper2,xlab="Sample Collection Depth", main = "Distance to Centroid by Category (CLR Data)", sub="Based on Aitchison Distance", col=colfunc(3))
dev.off()


### Export Global Env for Other Scripts ####
#save.image("data/Metagenomes/Analysis/mgm_analysis.Rdata")
# ^ includes all data combined in object bac.dat.all, ASV table (samples are rows, ASVs are columns), mgm_meta, and an ASV count table (where ASVs are rows, not columns)
# Version Information
sessionInfo()