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
  library(MoMAColors)
  library(microshades)
  library(lmtest)
})

#### Load Global Env to Import Count/ASV Tables ####
load("data/Amplicon/SSDust_16S.V3V4_W23_Data_Ready.Rdata") # save global env to Rdata file
#load("data/Amplicon/SSD_16S_CLR_EucDist_Ready.Rdata")

head(b.dust.all)
bac.ASV_table[1:4,1:4]
bac.ASV_table[(nrow(bac.ASV_table)-4):(nrow(bac.ASV_table)),(ncol(bac.ASV_table)-4):(ncol(bac.ASV_table))] # last 4 rows & cols
head(dust_meta)

#### Create Column for Gram +/- ID ####

## make column for Phyla by Gram +/- Abundances
head(b.dust.all)
b.dust.all$GramID<-ifelse((b.dust.all$Phylum == "Actinobacteria" | b.dust.all$Phylum == "Firmicutes"),"GramPos","GramNeg")
head(b.dust.all)

#### Kingdom Relative Abundance ####

# use dcast to count up ASVs within each Kingdom across all of the samples
b.kingdom_counts <- as.data.frame(dcast(b.dust.all, SampleID~Kingdom, value.var="Count", fun.aggregate=sum)) ###
head(b.kingdom_counts) # counts by kingdom per sample
dim(b.kingdom_counts)

rownames(b.kingdom_counts)<-b.kingdom_counts$SampleID
dim(b.kingdom_counts)
#b.kingdom_counts<-b.kingdom_counts[,colSums(b.kingdom_counts[,-1])>0] # drop kingdom that are not represented
dim(b.kingdom_counts) # sanity check that we dropped taxa with no hits

b.kingdom_RelAb<-data.frame(decostand(b.kingdom_counts[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.kingdom_RelAb) # sanity check to make sure the transformation worked!

b.kingdom_RelAb$SampleID<-rownames(b.kingdom_RelAb)
head(b.kingdom_RelAb)
#write.csv(b.kingdom_RelAb,"16S_Kingdom_Relative_Abundance.csv",row.names=TRUE) # good to save just in case

# melt down relativized data to merge with dust_meta
b.kingdom_m<-melt(b.kingdom_RelAb)

head(b.kingdom_m)
colnames(b.kingdom_m)[which(names(b.kingdom_m) == "variable")] <- "Kingdom"
colnames(b.kingdom_m)[which(names(b.kingdom_m) == "value")] <- "Count"
head(b.kingdom_m) ## relative abundance based on sum of counts by kingdom!

b.kingdom_RA_meta<-merge(b.kingdom_m,dust_meta, by="SampleID")
head(b.kingdom_RA_meta) ## relative abundance based on sum of counts by kingdom!
b.kingdom_RA_meta$SampleID = factor(b.kingdom_RA_meta$SampleID, levels=unique(b.kingdom_RA_meta$SampleID[order(b.kingdom_RA_meta$Site,b.kingdom_RA_meta$Seas_Coll_Year)]), ordered=TRUE)

# Barplot by SampleID
k.b1<-ggplot(b.kingdom_RA_meta, aes(x=SampleID, y=Count, fill=Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Kingdom Relative Abundance", x="SampleID", y="Relative Abundance", fill="Kingdom")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+guides(fill=guide_legend(ncol=2))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

ggsave(k.b1,filename = "figures/RelativeAbundance/Kingdom/SSD_16S_Kingdom.RA_barplot.png", width=12, height=10, dpi=600,create.dir = TRUE)

k.b1a<-ggplot(b.kingdom_RA_meta[b.kingdom_RA_meta$Kingdom=="Archaea",], aes(x=SampleID, y=Count, fill=Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Archaea Relative Abundance", x="SampleID", y="Relative Abundance",fill="Kingdom")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2))+scale_y_continuous(expand = c(0,0),limits = c(0,0.10))

ggsave(k.b1a,filename = "figures/RelativeAbundance/Kingdom/SSD_16S_Kingdom.RA_Archaea_Only_barplot.png", width=12, height=10, dpi=600,create.dir = TRUE)

head(b.kingdom_RA_meta)

# Heatmap by SampleID
k.h1<-ggplot(b.kingdom_RA_meta, aes(SampleID, Kingdom, fill= Count)) +geom_tile()+scale_fill_gradient2(low="blue3",mid="white",high="red",midpoint=0.45)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample ID", y="Microbial Kingdom", title="Microbial Kingdom & Sample Type",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0)) +
  facet_grid(. ~ SampDate, scales="free")

ggsave(k.h1,filename = "figures/RelativeAbundance/Kingdom/SSD_16S_Kingdom.RA_heatmap.png", width=12, height=10, dpi=600,create.dir = TRUE)

b.dust.all[1:4,1:4]

#### Phyla Relative Abundance ####

# use dcast to count up ASVs within each Phylum across all of the samples
b.phyla_counts <- as.data.frame(dcast(b.dust.all, SampleID~Phylum, value.var="Count", fun.aggregate=sum)) ###
head(b.phyla_counts) # counts by phyla per sample
dim(b.phyla_counts)

rownames(b.phyla_counts)<-b.phyla_counts$SampleID
dim(b.phyla_counts)
b.phyla_counts<-b.phyla_counts[,colSums(b.phyla_counts[,-1])>0] # drop phyla that are not represented
dim(b.phyla_counts) # sanity check that we dropped taxa with no hits

b.phyla_RelAb<-data.frame(decostand(b.phyla_counts[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.phyla_RelAb) # sanity check to make sure the transformation worked!

b.phyla_RelAb$SampleID<-rownames(b.phyla_RelAb)
head(b.phyla_RelAb)
#write.csv(b.phyla_RelAb,"16S_Phyla_Relative_Abundance.csv",row.names=TRUE) # good to save just in case

# melt down relativized data to merge with dust_meta
b.phyla_m<-melt(b.phyla_RelAb)

head(b.phyla_m)
colnames(b.phyla_m)[which(names(b.phyla_m) == "variable")] <- "Phylum"
colnames(b.phyla_m)[which(names(b.phyla_m) == "value")] <- "Count"
head(b.phyla_m) ## relative abundance based on sum of counts by phyla!

b.phyla_RA_meta<-merge(b.phyla_m,dust_meta, by="SampleID")
head(b.phyla_RA_meta) ## relative abundance based on sum of counts by phyla!
b.phyla_RA_meta$SampleID = factor(b.phyla_RA_meta$SampleID, levels=unique(b.phyla_RA_meta$SampleID[order(b.phyla_RA_meta$Site,b.phyla_RA_meta$Seas_Coll_Year)]), ordered=TRUE)

# Barplot by SampleID
p.b1<-ggplot(b.phyla_RA_meta, aes(x=SampleID, y=Count, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Phylum Relative Abundance", x="SampleID", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2))+scale_y_continuous(expand = c(0,0),limits = c(0,1))+
  facet_wrap(vars(Site), scales = "free")

ggsave(p.b1,filename = "figures/RelativeAbundance/Phylum/SSD_16S_Phyla.RA_barplot.png", width=12, height=10, dpi=600,create.dir = TRUE)

p.b1a<-ggplot(b.phyla_RA_meta[b.phyla_RA_meta$Count>0.05,], aes(x=SampleID, y=Count, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Phylum Relative Abundance", x="SampleID", y="Relative Abundance", fill="Phylum",subtitle="Only Taxa with Relative Abundance > 5%")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+scale_y_continuous(expand = c(0,0),limits = c(0,1)) +
  facet_wrap(vars(Site), scales = "free")

ggsave(p.b1a,filename = "figures/RelativeAbundance/Phylum/SSD_16S_Phyla.RA_5perc_barplot.png", width=12, height=10, dpi=600,create.dir = TRUE)

head(b.phyla_RA_meta)

# Heatmap by SampleID
p.h1<-ggplot(b.phyla_RA_meta, aes(SampleID, Phylum, fill= Count)) +geom_tile()+scale_fill_gradient2(low="blue3",mid="white",high="red",midpoint=0.45)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample ID", y="Microbial Phyla", title="Microbial Phyla & Sample Type",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0)) +
  facet_grid(. ~ SampDate, scales = "free")

ggsave(p.h1,filename = "figures/RelativeAbundance/Phylum/SSD_16S_Phyla.RA_heatmap.png", width=16, height=10, dpi=600,create.dir = TRUE)

b.dust.all[1:4,1:4]

#### Gram +/- Relative Abundance ####

# use dcast to count up ASVs within each Phylum across all of the samples
b.gram_counts <- as.data.frame(dcast(b.dust.all, SampleID~GramID, value.var="Count", fun.aggregate=sum)) ###
head(b.gram_counts) # counts by gram per sample
dim(b.gram_counts)

rownames(b.gram_counts)<-b.gram_counts$SampleID
dim(b.gram_counts)
#b.gram_counts<-b.gram_counts[,colSums(b.gram_counts[,-1])>0] # drop gram that are not represented
#dim(b.gram_counts) # sanity check that we dropped taxa with no hits

b.gram_RelAb<-data.frame(decostand(b.gram_counts[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.gram_RelAb) # sanity check to make sure the transformation worked!

b.gram_RelAb$SampleID<-rownames(b.gram_RelAb)
head(b.gram_RelAb)
#write.csv(b.gram_RelAb,"16S_Phyla_Relative_Abundance.csv",row.names=TRUE) # good to save just in case

# melt down relativized data to merge with dust_meta
b.gram_m<-melt(b.gram_RelAb)

head(b.gram_m)
colnames(b.gram_m)[which(names(b.gram_m) == "variable")] <- "GramID"
colnames(b.gram_m)[which(names(b.gram_m) == "value")] <- "Count"
head(b.gram_m) ## relative abundance based on sum of counts by gram!

b.gram_RA_meta<-merge(b.gram_m,dust_meta, by="SampleID")
head(b.gram_RA_meta) ## relative abundance based on sum of counts by gram!
b.gram_RA_meta$SampleID = factor(b.gram_RA_meta$SampleID, levels=unique(b.gram_RA_meta$SampleID[order(b.gram_RA_meta$Site,b.gram_RA_meta$Seas_Coll_Year)]), ordered=TRUE)

# Barplot by SampleID
p.b1<-ggplot(b.gram_RA_meta, aes(x=SampleID, y=Count, fill=GramID))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Gram +/- Relative Abundance", x="SampleID", y="Relative Abundance", fill="Gram +/-")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2))+scale_fill_manual(name ="Gram +/-",values=c("pink1","purple3"),labels=c("GramNeg"="Gram -","GramPos" ="Gram +"))+
  facet_wrap(vars(Site), scales = "free")

ggsave(p.b1,filename = "figures/RelativeAbundance/GramID/SSD_16S_GramPosorNeg_RA_barplot.png", width=12, height=10, dpi=600,create.dir = TRUE)

head(b.gram_RA_meta)

# Heatmap by SampleID
p.h1<-ggplot(b.gram_RA_meta, aes(SampleID, GramID, fill= Count)) +geom_tile()+scale_fill_gradient2(low="blue3",mid="white",high="red",midpoint=0.45)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample ID", y="Gram +/-", title="Gram +/- & Sample Type",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0)) +
  facet_grid(. ~ SampDate, scales = "free")

ggsave(p.h1,filename = "figures/RelativeAbundance/GramID/SSD_16S_Phyla.RA_heatmap.png", width=16, height=10, dpi=600,create.dir = TRUE)

b.dust.all[1:4,1:4]


#### Class Relative Abundance ####

# use dcast to count up ASVs within each Class across all of the samples
b.class_counts <- as.data.frame(dcast(b.dust.all, SampleID~Class, value.var="Count", fun.aggregate=sum)) ###
head(b.class_counts) # counts by class per sample
dim(b.class_counts)
rownames(b.class_counts)<-b.class_counts$SampleID
b.class_counts<-subset(b.class_counts, select=-c(Unknown))
dim(b.class_counts)
b.class_counts<-b.class_counts[,colSums(b.class_counts[,-1])>0] # drop classes that are not represented
dim(b.class_counts) # sanity check that we dropped taxa with no hits

b.class_RelAb<-data.frame(decostand(b.class_counts[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.class_RelAb) # sanity check to make sure the transformation worked!

b.class_RelAb$SampleID<-rownames(b.class_RelAb)
head(b.class_RelAb)
#write.csv(b.class_RelAb,"16S_class_Relative_Abundance.csv",row.names=TRUE) # good to save just in case

# melt down relativized data to merge with dust_meta
b.class_m<-melt(b.class_RelAb)

head(b.class_m)
colnames(b.class_m)[which(names(b.class_m) == "variable")] <- "Class"
colnames(b.class_m)[which(names(b.class_m) == "value")] <- "Count"
head(b.class_m) ## relative abundance based on sum of counts by class!

b.class_RA_meta<-merge(b.class_m,dust_meta, by="SampleID")
head(b.class_RA_meta) ## relative abundance based on sum of counts by class!
b.class_RA_meta$SampleID = factor(b.class_RA_meta$SampleID, levels=unique(b.class_RA_meta$SampleID[order(b.class_RA_meta$Site,b.class_RA_meta$Seas_Coll_Year)]), ordered=TRUE)

# Barplot by SampleID

c.b1<-ggplot(b.class_RA_meta, aes(x=SampleID, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Microbial Classes", x="SampleID", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5,
        legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+scale_y_continuous(expand = c(0,0),limits = c(0,1)) +
  guides(fill=guide_legend(ncol=4)) +
  facet_wrap(vars(Site), scales = "free")

ggsave(c.b1,filename = "figures/RelativeAbundance/Class/SSD_16S_Class.RA_barplot.png", width=15, height=10, dpi=600,create.dir = TRUE)

c.b2<-ggplot(b.class_RA_meta[b.class_RA_meta$Count>0.01,], aes(x=SampleID, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Microbial Classes", x="SampleID", y="Relative Abundance", fill="Class",subtitle="Only Taxa with Relative Abundance > 1%")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5,
        legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+scale_y_continuous(expand = c(0,0),limits = c(0,1)) +
  guides(fill=guide_legend(ncol=1)) +
  facet_wrap(vars(Site), scales = "free")

ggsave(c.b2,filename = "figures/RelativeAbundance/Class/SSD_16S_Class.RA_1perc_barplot.png", width=15, height=10, dpi=600,create.dir = TRUE)

c.b3<-ggplot(b.class_RA_meta[b.class_RA_meta$Count>0.05,], aes(x=SampleID, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Microbial Classes", x="SampleID", y="Relative Abundance", fill="Class",subtitle="Only Taxa with Relative Abundance > 5%")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5,
        legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+scale_y_continuous(expand = c(0,0),limits = c(0,1)) +
  guides(fill=guide_legend(ncol=1)) +
  facet_wrap(vars(Site), scales = "free")

ggsave(c.b3,filename = "figures/RelativeAbundance/Class/SSD_16S_Class.RA_5perc_barplot.png", width=15, height=10, dpi=600,create.dir = TRUE)

head(b.class_RA_meta)

# Heatmap by SampleID
c.h1<-ggplot(b.class_RA_meta, aes(SampleID, Class, fill= Count)) +geom_tile()+scale_fill_gradient2(low="blue3",mid="white",high="red",midpoint=0.3)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample ID", y="Microbial Class", title="Microbial Class & Sample Type",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0)) +
  facet_grid(. ~ SampDate,scales="free")

ggsave(c.h1,filename = "figures/RelativeAbundance/Class/SSD_16S_class.RA_heatmap.png", width=16, height=15, dpi=600,create.dir = TRUE)

b.dust.all[1:4,1:4]


#### Order Relative Abundance ####

# use dcast to count up ASVs within each Order across all of the samples
b.ord_counts <- as.data.frame(dcast(b.dust.all, SampleID~Order, value.var="Count", fun.aggregate=sum)) ###
head(b.ord_counts) # counts by class per sample
dim(b.ord_counts)
rownames(b.ord_counts)<-b.ord_counts$SampleID
b.ord_counts<-subset(b.ord_counts, select=-c(Unknown))
dim(b.ord_counts)
b.ord_counts<-b.ord_counts[,colSums(b.ord_counts[,-1])>0] # drop orders that are not represented
dim(b.ord_counts) # sanity check that we dropped taxa with no hits

b.ord_RelAb<-data.frame(decostand(b.ord_counts[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.ord_RelAb) # sanity check to make sure the transformation worked!

b.ord_RelAb$SampleID<-rownames(b.ord_RelAb)
head(b.ord_RelAb)
#write.csv(b.ord_RelAb,"16S_class_Relative_Abundance.csv",row.names=TRUE) # good to save just in case

# melt down relativized data to merge with dust_meta
b.ord_m<-melt(b.ord_RelAb)

head(b.ord_m)
colnames(b.ord_m)[which(names(b.ord_m) == "variable")] <- "Order"
colnames(b.ord_m)[which(names(b.ord_m) == "value")] <- "Count"
head(b.ord_m) ## relative abundance based on sum of counts by order!

b.ord_RA_meta<-merge(b.ord_m,dust_meta, by="SampleID")
head(b.ord_RA_meta) ## relative abundance based on sum of counts by order!
b.ord_RA_meta$SampleID = factor(b.ord_RA_meta$SampleID, levels=unique(b.ord_RA_meta$SampleID[order(b.ord_RA_meta$Site,b.ord_RA_meta$Seas_Coll_Year)]), ordered=TRUE)

# Barplot by SampleID

o.b1<-ggplot(b.ord_RA_meta, aes(x=SampleID, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Microbial Orders", x="SampleID", y="Relative Abundance", fill="Order")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5,
        legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+scale_y_continuous(expand = c(0,0),limits = c(0,1)) +
  guides(fill=guide_legend(ncol=4)) +
  facet_wrap(vars(Site), scales = "free")

ggsave(o.b1,filename = "figures/RelativeAbundance/Order/SSD_16S_Order.RA_barplot.png", width=15, height=10, dpi=600,create.dir = TRUE)

o.b2<-ggplot(b.ord_RA_meta[b.ord_RA_meta$Count>0.01,], aes(x=SampleID, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Microbial Orders", x="SampleID", y="Relative Abundance", fill="Order",subtitle="Only Taxa with Relative Abundance > 1%")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5,
        legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+scale_y_continuous(expand = c(0,0),limits = c(0,1)) +
  guides(fill=guide_legend(ncol=1)) +
  facet_wrap(vars(Site), scales = "free")

ggsave(o.b2,filename = "figures/RelativeAbundance/Order/SSD_16S_Order.RA_1perc_barplot.png", width=15, height=10, dpi=600,create.dir = TRUE)

o.b3<-ggplot(b.ord_RA_meta[b.ord_RA_meta$Count>0.05,], aes(x=SampleID, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Microbial Orders", x="SampleID", y="Relative Abundance", fill="Order",subtitle="Only Taxa with Relative Abundance > 5%")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5,
        legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+scale_y_continuous(expand = c(0,0),limits = c(0,1)) +
  guides(fill=guide_legend(ncol=1)) +
  facet_wrap(vars(Site), scales = "free")

ggsave(o.b3,filename = "figures/RelativeAbundance/Order/SSD_16S_Order.RA_5perc_barplot.png", width=15, height=10, dpi=600,create.dir = TRUE)

burkholderiales.only<-ggplot(b.ord_RA_meta[b.ord_RA_meta$Order=="Burkholderiales",], aes(x=SampleID, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Burkholderiales", x="SampleID", y="Relative Abundance", fill="Order",subtitle="Only Burkholderiales Order")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5,
        legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+scale_y_continuous(expand = c(0,0),limits = c(0,1)) +
  guides(fill=guide_legend(ncol=1)) +
  facet_wrap(vars(Site), scales = "free")

ggsave(burkholderiales.only,filename = "figures/RelativeAbundance/Order/SSD_16S_Order.RA_Burkholderiales_Only_barplot.png", width=15, height=10, dpi=600,create.dir = TRUE)

head(b.ord_RA_meta)

# Heatmap by SampleID
o.h1<-ggplot(b.ord_RA_meta, aes(SampleID, Order, fill= Count)) +geom_tile()+scale_fill_gradient2(low="blue3",mid="white",high="red",midpoint=0.3)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample ID", y="Microbial Order", title="Microbial Order & Sample Type",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0)) +
  facet_grid(. ~ SampDate,scales="free")

ggsave(o.h1,filename = "figures/RelativeAbundance/Order/SSD_16S_class.RA_heatmap.png", width=16, height=15, dpi=600,create.dir = TRUE)

b.dust.all[1:4,1:4]


#### Family Relative Abundance ####

# use dcast to count up ASVs within each Family across all of the samples
b.fam_counts <- as.data.frame(dcast(b.dust.all, SampleID~Family, value.var="Count", fun.aggregate=sum)) ###
head(b.fam_counts) # counts by fam per sample
dim(b.fam_counts)
rownames(b.fam_counts)<-b.fam_counts$SampleID
colnames(b.fam_counts)<-gsub(" ", ".",colnames(b.fam_counts))
b.fam_counts<-subset(b.fam_counts, select=-c(Unknown, Unknown.Family))
dim(b.fam_counts)
b.fam_counts<-b.fam_counts[,colSums(b.fam_counts[,-1])>0] # drop families that are not represented
dim(b.fam_counts) # sanity check that we dropped taxa with no hits

b.fam_RelAb<-data.frame(decostand(b.fam_counts[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.fam_RelAb) # sanity check to make sure the transformation worked!

b.fam_RelAb$SampleID<-rownames(b.fam_RelAb)
head(b.fam_RelAb)
#write.csv(b.fam_RelAb,"16S_fam_Relative_Abundance.csv",row.names=TRUE) # good to save just in case

# melt down relativized data to merge with dust_meta
b.fam_m<-melt(b.fam_RelAb)

head(b.fam_m)
colnames(b.fam_m)[which(names(b.fam_m) == "variable")] <- "Family"
colnames(b.fam_m)[which(names(b.fam_m) == "value")] <- "Count"
head(b.fam_m) ## relative abundance based on sum of counts by fam!
b.fam_m$Family<-gsub("^X.","",b.fam_m$Family) # get rid of leading X. in Family names
b.fam_m$Family<-gsub("\\.\\."," ",b.fam_m$Family) # get rid of .. in family name --> . is regex
b.fam_m$Family<-gsub("\\."," ",b.fam_m$Family) # get rid of . in family name --> . is regex

b.fam_RA_meta<-merge(b.fam_m,dust_meta, by="SampleID")
head(b.fam_RA_meta) ## relative abundance based on sum of counts by fam!
b.fam_RA_meta$SampleID = factor(b.fam_RA_meta$SampleID, levels=unique(b.fam_RA_meta$SampleID[order(b.fam_RA_meta$Site,b.fam_RA_meta$Seas_Coll_Year)]), ordered=TRUE)

# Barplot by SampleID

f.b1<-ggplot(b.fam_RA_meta, aes(x=SampleID, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Microbial Families", x="SampleID", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=3))+scale_y_continuous(expand = c(0,0),limits = c(0,1)) +
  facet_wrap(vars(Site), scales = "free")

ggsave(f.b1,filename = "figures/RelativeAbundance/Family/SSD_16S_fam.RA_barplot.png", width=17, height=10, dpi=600,create.dir = TRUE)

f.b1a<-ggplot(b.fam_RA_meta[b.fam_RA_meta$Count>0.01,], aes(x=SampleID, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Microbial Families", subtitle="Only Taxa with Relative Abundance > 1%",x="SampleID", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=4))+scale_y_continuous(expand = c(0,0),limits = c(0,1)) +
  facet_wrap(vars(Site), scales = "free")

ggsave(f.b1a,filename = "figures/RelativeAbundance/Family/SSD_16S_fam.RA_1perc_barplot.png", width=17, height=10, dpi=600,create.dir = TRUE)

f.b1b<-ggplot(b.fam_RA_meta[b.fam_RA_meta$Count>0.05,], aes(x=SampleID, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Microbial Families", subtitle="Only Taxa with Relative Abundance > 5%",x="SampleID", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2))+scale_y_continuous(expand = c(0,0),limits = c(0,1)) +
  facet_wrap(vars(Site), scales = "free")

ggsave(f.b1b,filename = "figures/RelativeAbundance/Family/SSD_16S_fam.RA_5perc_barplot.png", width=17, height=10, dpi=600,create.dir = TRUE)

f.b1c<-ggplot(b.fam_RA_meta[b.fam_RA_meta$Count>0.10,], aes(x=SampleID, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Microbial Families", subtitle="Only Taxa with Relative Abundance > 10%",x="SampleID", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+scale_y_continuous(expand = c(0,0),limits = c(0,1)) +
  facet_wrap(vars(Site), scales = "free")

ggsave(f.b1c,filename = "figures/RelativeAbundance/Family/SSD_16S_fam.RA_10perc_barplot.png", width=17, height=10, dpi=600,create.dir = TRUE)

head(b.fam_RA_meta)

# Heatmap by SampleID

f.h1<-ggplot(b.fam_RA_meta, aes(SampleID, Family, fill= Count)) +geom_tile()+scale_fill_gradient2(low="blue3",mid="white",high="red",midpoint=0.3)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample ID", y="Microbial Family", title="Microbial Families & Sample Type",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0)) +
  facet_grid(. ~ SampDate, scales="free")

ggsave(f.h1,filename = "figures/RelativeAbundance/Family/SSD_16S_fam.RA_heatmap.png", width=12, height=10, dpi=600,create.dir = TRUE)

# Taxonomic summary by Sample ID + Collection Period

f.ts1<-ggplot(b.fam_RA_meta[b.fam_RA_meta$Count>0.025,], aes(Family, Count)) +
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date",values=unique(b.fam_RA_meta$SampDate_Color[order(b.fam_RA_meta$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
  labs(x="Microbial Families", y="Relative Abundance", title="Salton Sea Dust: Microbial Families & Sample Date",subtitle="Includes taxa with Relative Abundance > 2.5%")+
  scale_shape_manual(values = c(7,10, 15,16)) + coord_flip()

ggsave(f.ts1,filename = "figures/RelativeAbundance/Family/SSD_16S_Family.RA_2.5perc_taxasum.png", width=20, height=23, dpi=600,create.dir = TRUE)

f.ts2<-ggplot(b.fam_RA_meta[b.fam_RA_meta$Count>0.05,], aes(Family, Count)) +
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date",values=unique(b.fam_RA_meta$SampDate_Color[order(b.fam_RA_meta$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
  labs(x="Microbial Families", y="Relative Abundance", title="Salton Sea Dust: Microbial Families & Sample Date",subtitle="Includes taxa with Relative Abundance > 5%")+
  scale_shape_manual(values = c(7,10, 15,16)) + coord_flip()

ggsave(f.ts2,filename = "figures/RelativeAbundance/Family/SSD_16S_Family.RA_5perc_taxasum.png", width=20, height=23, dpi=600,create.dir = TRUE)

b.dust.all[1:4,1:4]


# # by Family + depth
# bac.fam.dep <- as.data.frame(dcast(b.dust.all,Depth_m~Family, value.var="Count", fun.aggregate=sum)) ###
# head(bac.fam.dep) # counts by Family + sample depe
# rownames(bac.fam.dep)<-bac.fam.dep$Depth_m
# colnames(bac.fam.dep)<-gsub(" ", ".",colnames(bac.fam.dep))
# bac.fam.dep<-subset(bac.fam.dep, select=-c(Unknown, Unknown.Family))
#
# b.RA_fam.dep<-data.frame(decostand(bac.fam.dep[,-1], method="total", MARGIN=1, na.rm=TRUE))
# # relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
# rowSums(b.RA_fam.dep) # sanity check
# b.RA_fam.dep$Depth_m<-rownames(b.RA_fam.dep) # Depth_m is now a character, not a factor!
# head(b.RA_fam.dep)
#
# #melt down relativized data to merge with dust_meta
# b.fam.dep_m<-melt(b.RA_fam.dep, by="Depth_m")
#
# head(b.fam.dep_m)
# colnames(b.fam.dep_m)[which(names(b.fam.dep_m) == "variable")] <- "Family"
# colnames(b.fam.dep_m)[which(names(b.fam.dep_m) == "value")] <- "Count"
# head(b.fam.dep_m) ## relative abundance based on sum of counts by Family!
#
# #dep_meta<-unique(data.frame("Depth_m"=dust_meta$Depth_m, "Sample_Color"=dust_meta$Sample_Color))
# #p_dep_meta<-merge(dep_meta,b.fam.dep_m, by="Depth_m")
#
# # Barplot by Depth
#
# fd1<-ggplot(b.fam.dep_m, aes(x=Depth_m, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+theme_classic()+
#   labs(title = "Relative Abundance of Microbial Classes", x="Depth (m)", y="Relative Abundance", fill="Class")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
#   guides(fill=guide_legend(ncol=3))+
#   scale_y_continuous(expand = c(0,0),limits = c(0,1))
# #+ scale_x_discrete(labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))
#
# ggsave(fd1,filename = "figures/RelativeAbundance/Family/SSD_16S_fam.RA_barplot_depth.png", width=12, height=10, dpi=600,create.dir = TRUE)
#
# fd1a<-ggplot(b.fam.dep_m[b.fam.dep_m$Count>0.05,], aes(x=Depth_m, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+theme_classic()+
#   labs(title = "Relative Abundance of Microbial Classes", x="SampleID", y="Relative Abundance", fill="Class",subtitle="Only Relative Abundance > 5%")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
#   scale_y_continuous(expand = c(0,0),limits = c(0,1))
# #+ scale_x_discrete(labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))
#
# ggsave(fd1a,filename = "figures/RelativeAbundance/Family/SSD_16S_fam.RA_barplot_depth_5percent.png", width=12, height=10, dpi=600,create.dir = TRUE)
#
# # Taxonomic Summary by Depth
#
# fd2<-ggplot(b.fam.dep_m, aes(Family, Count)) +
#   geom_jitter(aes(color=as.numeric(Depth_m)), size=2, width=0.15, height=0) +
#   scale_colour_gradient2(high="blue3",low="red",midpoint=5.5,guide = guide_colourbar(reverse = TRUE)) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Microbial Family", y="Relative Abundance", title="Microbial Families & Depth",color="Depth (m)")
#
# ggsave(fd2,filename = "figures/RelativeAbundance/Family/SSD_16S_fam.RA_depth_taxasum.png", width=15, height=10, dpi=600,create.dir = TRUE)
#
# fd2a<-ggplot(b.fam.dep_m[b.fam.dep_m$Count>0.05,], aes(Family, Count)) +
#   geom_jitter(aes(color=as.numeric(Depth_m)), size=2, width=0.15, height=0) +
#   scale_colour_gradient2(low="red",high="blue3",midpoint=5.5,guide = guide_colourbar(reverse = TRUE)) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Microbial Family", y="Relative Abundance", title="Microbial Families & Depth",color="Depth (m)")
#
# ggsave(fd2a,filename = "figures/RelativeAbundance/Family/SSD_16S_fam.RA_depth_5percent_taxasum.png", width=15, height=10, dpi=600,create.dir = TRUE)
#
# # by Family + Sampling Date
# bac.fam.date <- as.data.frame(dcast(b.dust.all,SampDate~Family, value.var="Count", fun.aggregate=sum)) ###
# head(bac.fam.date) # counts by Family + sample depe
# rownames(bac.fam.date)<-bac.fam.date$SampDate
# colnames(bac.fam.date)<-gsub(" ", ".",colnames(bac.fam.date))
# bac.fam.date<-subset(bac.fam.date, select=-c(Unknown, Unknown.Family))
#
# b.RA_fam.date<-data.frame(decostand(bac.fam.date[,-1], method="total", MARGIN=1, na.rm=TRUE))
# # relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
# rowSums(b.RA_fam.date) # sanity check
# b.RA_fam.date$SampDate<-rownames(b.RA_fam.date)
# head(b.RA_fam.date)
#
# #melt down relativized data to merge with dust_meta
# b.fam.date_m<-melt(b.RA_fam.date, by="SampDate")
#
# head(b.fam.date_m)
# colnames(b.fam.date_m)[which(names(b.fam.date_m) == "variable")] <- "Family"
# colnames(b.fam.date_m)[which(names(b.fam.date_m) == "value")] <- "Count"
# head(b.fam.date_m) ## relative abundance based on sum of counts by Family!
#
# b.fam.date_m$SampDate<-factor(b.fam.date_m$SampDate, levels=c("August.2021","December.2021","April.2022"))
#
# # Barplot by Sample Date
#
# fsd1<-ggplot(b.fam.date_m, aes(x=SampDate, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+theme_classic()+
#   labs(title = "Relative Abundance of Microbial Families", x="SampleID", y="Relative Abundance", fill="Family")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
#   guides(fill=guide_legend(ncol=3))+
#   scale_y_continuous(expand = c(0,0),limits = c(0,1))+ scale_x_discrete(labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))
#
# ggsave(fsd1,filename = "figures/RelativeAbundance/Family/SSD_16S_fam.RA_barplot_sampdate.png", width=12, height=10, dpi=600,create.dir = TRUE)
#
# fsd1a<-ggplot(b.fam.date_m[b.fam.date_m$Count>0.01,], aes(x=SampDate, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+theme_classic()+
#   labs(title = "Relative Abundance of Microbial Families", x="SampleID", y="Relative Abundance", fill="Family",subtitle="Only Relative Abundance > 1%")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
#   guides(fill=guide_legend(ncol=2))+
#   scale_y_continuous(expand = c(0,0),limits = c(0,1))+ scale_x_discrete(labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))
#
# ggsave(fsd1a,filename = "figures/RelativeAbundance/Family/SSD_16S_fam.RA_barplot_sampdate_1percent.png", width=12, height=10, dpi=600,create.dir = TRUE)
#
# # Taxonomic Summary by Sample Date
#
# #b.fam.date_m2<-merge(b.fam.date_m, dust_meta, by="SampDate")
#
# colorset1 # remember which date goes with each color
#
# fsd2<-ggplot(b.fam.date_m, aes(Family, Count)) +
#   geom_jitter(aes(color=factor(SampDate)), size=2, width=0.15, height=0) +
#   scale_color_manual(name ="Sample Date", values=c("#ef781c","#03045e","#059c3f"), labels=c("August 2021","December 2021","April 2022")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=0.5,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Microbial Family", y="Relative Abundance", title="Microbial Families & Sample Date")
#
# ggsave(fsd2,filename = "figures/RelativeAbundance/Family/SSD_16S_fam.RA_date_taxasum.png", width=15, height=10, dpi=600,create.dir = TRUE)
#
# fsd3<-ggplot(b.fam.date_m[b.fam.date_m$Count>0.1,], aes(Family, Count)) +
#   geom_jitter(aes(color=factor(SampDate)), size=2, width=0.15, height=0) +
#   scale_color_manual(name ="Sample Date", values=c("#ef781c","#03045e","#059c3f"), labels=c("August 2021","December 2021","April 2022")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=0.5,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Microbial Family", y="Relative Abundance", title="Microbial Families & Sample Date")
#
# ggsave(fsd3,filename = "figures/RelativeAbundance/Family/SSD_16S_fam.RA_date_taxasum_5percent.png", width=15, height=10, dpi=600,create.dir = TRUE)


#### Genus Relative Abundance ####
# use dcast to count up ASVs within each Genus across all of the samples
head(b.dust.all)
b.dust.all.g<-subset(b.dust.all, b.dust.all$Genus!="Unknown") # drop unknown genera so they don't skew analyses
"Unknown" %in% b.dust.all.g$Genus

b.genus_counts <- as.data.frame(dcast(b.dust.all.g, SampleID~Genus, value.var="Count", fun.aggregate=sum)) ###
head(b.genus_counts) # counts by genus per sample
dim(b.genus_counts)
rownames(b.genus_counts)<-b.genus_counts$SampleID
b.genus_counts[1:4,1:4]
b.genus_counts<-b.genus_counts[,colSums(b.genus_counts[,-1])>0] # drop classes that are not represented
dim(b.genus_counts) # sanity check that we dropped taxa with no hits

b.genus_RelAb<-data.frame(decostand(b.genus_counts[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.genus_RelAb) # sanity check to make sure the transformation worked!

b.genus_RelAb$SampleID<-rownames(b.genus_RelAb)
head(b.genus_RelAb)
#write.csv(b.genus_RelAb,"16S_Genera_Relative_Abundance.csv",row.names=TRUE) # good to save just in case

# melt down relativized data to merge with dust_meta
b.genus_m<-melt(b.genus_RelAb)

head(b.genus_m)
colnames(b.genus_m)[which(names(b.genus_m) == "variable")] <- "Genus"
colnames(b.genus_m)[which(names(b.genus_m) == "value")] <- "Count"
head(b.genus_m) ## relative abundance based on sum of counts by genus!
b.genus_m$Genus<-gsub("^X.","",b.genus_m$Genus) # get rid of leading X. in Genus names
b.genus_m$Genus<-gsub("\\.\\."," ",b.genus_m$Genus) # get rid of .. in species name --> . is regex
b.genus_m$Genus<-gsub("\\."," ",b.genus_m$Genus) # get rid of . in species name --> . is regex
b.genus_m$Genus<-gsub("_"," ",b.genus_m$Genus) #
head(b.genus_m) ## relative abundance based on sum of counts by genus!

b.genus_RA_meta<-merge(b.genus_m,dust_meta, by="SampleID")
head(b.genus_RA_meta) ## relative abundance based on sum of counts by genus!
max(b.genus_RA_meta$Count)
b.genus_RA_meta$SampleID = factor(b.genus_RA_meta$SampleID, levels=unique(b.genus_RA_meta$SampleID[order(b.genus_RA_meta$Site,b.genus_RA_meta$Seas_Coll_Year)]), ordered=TRUE)
b.genus_RA_meta$Sample_Type<-"Dust"

# separate genera RelAb data by site for downstream figs
WI.gen.RA<-subset(b.genus_RA_meta,b.genus_RA_meta$Site=="WI")
BDC.gen.RA<-subset(b.genus_RA_meta,b.genus_RA_meta$Site=="BDC")
PD.gen.RA<-subset(b.genus_RA_meta,b.genus_RA_meta$Site=="PD")
DP.gen.RA<-subset(b.genus_RA_meta,b.genus_RA_meta$Site=="DP")

saveRDS(b.genus_RA_meta, file = "data/Amplicon/SSD_GenusOnly_RelativeAbundance_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

#### Genus + Species Relative Abundance ####
# use dcast to count up ASVs within each Genus across all of the samples
head(b.dust.all)
b.dust.all.g<-subset(b.dust.all, b.dust.all$Genus!="Unknown") # drop unknown genera so they don't skew analyses
"Unknown" %in% b.dust.all.g$Genus

b.gen.spec_counts <- as.data.frame(dcast(b.dust.all.g, SampleID~Genus+Species, value.var="Count", fun.aggregate=sum)) ###
head(b.gen.spec_counts) # counts by gen.spec per sample
dim(b.gen.spec_counts)
rownames(b.gen.spec_counts)<-b.gen.spec_counts$SampleID
b.gen.spec_counts[1:4,1:4]
b.gen.spec_counts<-b.gen.spec_counts[,colSums(b.gen.spec_counts[,-1])>0] # drop classes that are not represented
dim(b.gen.spec_counts) # sanity check that we dropped taxa with no hits

b.gen.spec_RelAb<-data.frame(decostand(b.gen.spec_counts[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.gen.spec_RelAb) # sanity check to make sure the transformation worked!

b.gen.spec_RelAb$SampleID<-rownames(b.gen.spec_RelAb)
head(b.gen.spec_RelAb)
#write.csv(b.gen.spec_RelAb,"16S_Genera_Relative_Abundance.csv",row.names=TRUE) # good to save just in case

# melt down relativized data to merge with dust_meta
b.gen.spec_m<-melt(b.gen.spec_RelAb)

head(b.gen.spec_m)
colnames(b.gen.spec_m)[which(names(b.gen.spec_m) == "variable")] <- "Genus_species"
colnames(b.gen.spec_m)[which(names(b.gen.spec_m) == "value")] <- "Count"
head(b.gen.spec_m) ## relative abundance based on sum of counts by gen.spec!
b.gen.spec_m$Genus_species<-gsub("^X.","",b.gen.spec_m$Genus_species) # get rid of leading X. in Genus_species names
b.gen.spec_m$Genus_species<-gsub("\\.\\."," ",b.gen.spec_m$Genus_species) # get rid of .. in species name --> . is regex
b.gen.spec_m$Genus_species<-gsub("\\."," ",b.gen.spec_m$Genus_species) # get rid of . in species name --> . is regex
b.gen.spec_m$Genus_species<-gsub("_"," ",b.gen.spec_m$Genus_species) #
head(b.gen.spec_m) ## relative abundance based on sum of counts by gen.spec!

b.gen.spec_RA_meta<-merge(b.gen.spec_m,dust_meta, by="SampleID")
head(b.gen.spec_RA_meta) ## relative abundance based on sum of counts by gen.spec!
max(b.gen.spec_RA_meta$Count)
b.gen.spec_RA_meta$SampleID = factor(b.gen.spec_RA_meta$SampleID, levels=unique(b.gen.spec_RA_meta$SampleID[order(b.gen.spec_RA_meta$Site,b.gen.spec_RA_meta$Seas_Coll_Year)]), ordered=TRUE)
b.gen.spec_RA_meta$Sample_Type<-"Dust"

# separate genera RelAb data by site for downstream figs
WI.gen.spec.RA<-subset(b.gen.spec_RA_meta,b.gen.spec_RA_meta$Site=="WI")
BDC.gen.spec.RA<-subset(b.gen.spec_RA_meta,b.gen.spec_RA_meta$Site=="BDC")
PD.gen.spec.RA<-subset(b.gen.spec_RA_meta,b.gen.spec_RA_meta$Site=="PD")
DP.gen.spec.RA<-subset(b.gen.spec_RA_meta,b.gen.spec_RA_meta$Site=="DP")

saveRDS(b.gen.spec_RA_meta, file = "data/Amplicon/SSD_GenusSpecies_RelativeAbundance_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

# Barplot by SampleID

b.gen_RAall<-ggplot(b.gen.spec_RA_meta, aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=10)) +
  facet_wrap(vars(Site), scales = "free")

ggsave(b.gen_RAall,filename = "figures/RelativeAbundance/Genus_species/SSD_16S_Genera.Spec.RA_barplot.png", width=30, height=10, dpi=600,create.dir = TRUE)

b.gen_RA0.0<-ggplot(b.gen.spec_RA_meta[b.gen.spec_RA_meta$Count>0.005,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 0.05%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=10)) +
  facet_wrap(vars(Site), scales = "free")

ggsave(b.gen_RA0.0,filename = "figures/RelativeAbundance/Genus_species/SSD_16S_Genera.Spec.RA_barplot.png", width=30, height=10, dpi=600,create.dir = TRUE)

b.gen_RA0<-ggplot(b.gen.spec_RA_meta[b.gen.spec_RA_meta$Count>0.01,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 1%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=5))

ggsave(b.gen_RA0,filename = "figures/RelativeAbundance/Genus_species/SSD_16S_Genera.Spec.RA_1perc_barplot.png", width=30, height=10, dpi=600,create.dir = TRUE)

b.gen_RA0v2<-ggplot(b.gen.spec_RA_meta[b.gen.spec_RA_meta$Count>0.01,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 1%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2))+facet_wrap(vars(Site), scales = "free")

ggsave(b.gen_RA0v2,filename = "figures/RelativeAbundance/Genus_species/SSD_16S_Genera.Spec.RA_1perc_barplot_v2.png", width=30, height=25, dpi=600,create.dir = TRUE)

b.gen_RA01<-ggplot(b.gen.spec_RA_meta[b.gen.spec_RA_meta$Count>0.02,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 2%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=5))

ggsave(b.gen_RA01,filename = "figures/RelativeAbundance/Genus_species/SSD_16S_Genera.Spec.RA_2perc_barplot.png", width=30, height=20, dpi=600,create.dir = TRUE)

b.gen_RA01v2<-ggplot(b.gen.spec_RA_meta[b.gen.spec_RA_meta$Count>0.02,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 2%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2))+
  facet_wrap(vars(Site), scales = "free")

ggsave(b.gen_RA01v2,filename = "figures/RelativeAbundance/Genus_species/SSD_16S_Genera.Spec.RA_2perc_barplot_v2.png", width=30, height=25, dpi=600,create.dir = TRUE)

b.gen_RA1a<-ggplot(b.gen.spec_RA_meta[b.gen.spec_RA_meta$Count>0.05,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 5%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2))

ggsave(b.gen_RA1a,filename = "figures/RelativeAbundance/Genus_species/SSD_16S_Genera.Spec.RA_5perc_barplot_v1.png", width=20, height=10, dpi=600,create.dir = TRUE)

b.gen_RA1v2<-ggplot(b.gen.spec_RA_meta[b.gen.spec_RA_meta$Count>0.05,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 5%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2))+facet_wrap(vars(Site), scales = "free")

ggsave(b.gen_RA1v2,filename = "figures/RelativeAbundance/Genus_species/SSD_16S_Genera.Spec.RA_5perc_barplot_v2.png", width=20, height=25, dpi=600,create.dir = TRUE)

b.gen_RA1v3<-ggplot(b.gen.spec_RA_meta[b.gen.spec_RA_meta$Count>0.05,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 5%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2))+facet_wrap(vars(Site), scales = "free")

ggsave(b.gen_RA1v3,filename = "figures/RelativeAbundance/Genus_species/SSD_16S_Genera.Spec.RA_5perc_barplot_v3.png", width=20, height=10, dpi=600,create.dir = TRUE)

b.gen_RA2v1<-ggplot(b.gen.spec_RA_meta[b.gen.spec_RA_meta$Count>0.10,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 10%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))

ggsave(b.gen_RA2v1,filename = "figures/RelativeAbundance/Genus_species/SSD_16S_Genera.Spec.RA_10perc_barplot_v1.png", width=16, height=10, dpi=600,create.dir = TRUE)

b.gen_RA2v2<-ggplot(b.gen.spec_RA_meta[b.gen.spec_RA_meta$Count>0.10,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 10%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+
  facet_wrap(vars(Site), scales = "free")

ggsave(b.gen_RA2v2,filename = "figures/RelativeAbundance/Genus_species/SSD_16S_Genera.Spec.RA_10perc_barplot_v2.png", width=16, height=25, dpi=600,create.dir = TRUE)

b.gen_RA2v3<-ggplot(b.gen.spec_RA_meta[b.gen.spec_RA_meta$Count>0.10,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 10%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+facet_wrap(vars(Site), scales = "free")

ggsave(b.gen_RA2v3,filename = "figures/RelativeAbundance/Genus_species/SSD_16S_Genera.Spec.RA_10perc_barplot_v3.png", width=16, height=10, dpi=600,create.dir = TRUE)

b.gen_RA3a<-ggplot(b.gen.spec_RA_meta[b.gen.spec_RA_meta$Count>0.15,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 15%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))

ggsave(b.gen_RA3a,filename = "figures/RelativeAbundance/Genus_species/SSD_16S_Genera.Spec.RA_15perc_barplot_v1.png", width=16, height=10, dpi=600,create.dir = TRUE)

b.gen_RA3b<-ggplot(b.gen.spec_RA_meta[b.gen.spec_RA_meta$Count>0.15,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 15%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+facet_wrap(vars(Site), scales = "free")

ggsave(b.gen_RA3b,filename = "figures/RelativeAbundance/Genus_species/SSD_16S_Genera.Spec.RA_15perc_barplot_v2.png", width=16, height=10, dpi=600,create.dir = TRUE)

b.gen_RA4<-ggplot(b.gen.spec_RA_meta[b.gen.spec_RA_meta$Count>0.25,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 25%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))

ggsave(b.gen_RA4,filename = "figures/RelativeAbundance/Genus_species/SSD_16S_Genera.Spec.RA_25perc_barplot.png", width=15, height=10, dpi=600,create.dir = TRUE)

b.gen_RA5<-ggplot(b.gen.spec_RA_meta[b.gen.spec_RA_meta$Count>0.35,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes taxa with Relative Abundance > 35%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))

ggsave(b.gen_RA5,filename = "figures/RelativeAbundance/Genus_species/SSD_16S_Genera.Spec.RA_35perc_barplot.png", width=12, height=10, dpi=600,create.dir = TRUE)

# plot just Massilia
massilia.relab<-b.gen.spec_RA_meta[grepl("Massilia",b.gen.spec_RA_meta$Genus_species),]
massilia.only<-ggplot(massilia.relab[!grepl("unknown",massilia.relab$Genus_species),], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Massilia Species", x="SampleID", y="Relative Abundance", fill="Taxa",subtitle="Only species in Massilia - excluding unknown species")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5,
        legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  guides(fill=guide_legend(ncol=1)) +
  facet_wrap(vars(Site), scales = "free")

ggsave(massilia.only,filename = "figures/RelativeAbundance/Genus_species/SSD_16S_Genera.Spec.RA_Massilia_Only_barplot.png", width=15, height=10, dpi=600,create.dir = TRUE)

# plot just Sphingomonas
sphingomonas.relab<-b.gen.spec_RA_meta[grepl("Sphingomonas",b.gen.spec_RA_meta$Genus_species),]
sphingomonas.only<-ggplot(sphingomonas.relab[!grepl("unknown",sphingomonas.relab$Genus_species),], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Sphingomonas Species", x="SampleID", y="Relative Abundance", fill="Taxa",subtitle="Only species in Sphingomonas - excluding unknown species")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5,
        legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  guides(fill=guide_legend(ncol=2)) +
  facet_wrap(vars(Site), scales = "free")

ggsave(sphingomonas.only,filename = "figures/RelativeAbundance/Genus_species/SSD_16S_Genera.Spec.RA_Sphingomonas_Only_barplot.png", width=15, height=10, dpi=600,create.dir = TRUE)

# prep for heatmap
max(b.gen.spec_RA_meta$Count)
mean(b.gen.spec_RA_meta$Count)
max(b.gen.spec_RA_meta$Count)/2 # what is the mid point of the RA here?

# Heatmap by SampleID

# g.h1<-ggplot(b.gen.spec_RA_meta[b.gen.spec_RA_meta$Count>0.01,], aes(SampleID, Genus_species, fill= Count)) +geom_tile()+scale_fill_gradient2(low="orange",mid="white",high="purple",midpoint=0.35)+
#   theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Sample ID", y="Microbial Genera", title="Microbial Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 1%",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0))
#
# ggsave(g.h1,filename = "figures/RelativeAbundance/Genus_species/SSD_16S_Genera.RA_1perc_heatmap_A.png", width=20, height=15, dpi=600,create.dir = TRUE)
#
g.h2<-ggplot(b.gen.spec_RA_meta[b.gen.spec_RA_meta$Count>0.05,], aes(SampleID, Genus_species, fill= Count)) +geom_tile()+scale_fill_gradient2(low="blue",mid="pink",high="red",midpoint=0.3)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample ID", y="Microbial Genera", title="Microbial Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 5%",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0))

ggsave(g.h2,filename = "figures/RelativeAbundance/Genus_species/SSD_16S_Genera.RA_5perc_heatmap_B.png", width=16, height=10, dpi=600,create.dir = TRUE)

b.dust.all[1:4,1:4]

# Taxonomic Summary by Sample ID + Collection Date

ggplot(b.gen.spec_RA_meta[b.gen.spec_RA_meta$Count>0.01,], aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date",values=unique(b.gen.spec_RA_meta$SampDate_Color[order(b.gen.spec_RA_meta$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Salton Sea Dust: Microbial Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 1%")+
  scale_shape_manual(values = c(7,10, 15,16)) + coord_flip()

tg0<-ggplot(b.gen.spec_RA_meta[b.gen.spec_RA_meta$Count>0.02,], aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date",values=unique(b.gen.spec_RA_meta$SampDate_Color[order(b.gen.spec_RA_meta$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Salton Sea Dust: Microbial Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 2%")+
  scale_shape_manual(values = c(7,10, 15,16)) + coord_flip()

ggsave(tg0,filename = "figures/RelativeAbundance/Genus_species/SSD_16S_Genera.RA_2perc_taxasum.png", width=20, height=23, dpi=600,create.dir = TRUE)


tg1<-ggplot(b.gen.spec_RA_meta[b.gen.spec_RA_meta$Count>0.05,], aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date",values=unique(b.gen.spec_RA_meta$SampDate_Color[order(b.gen.spec_RA_meta$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Salton Sea Dust: Microbial Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 5%")+
  scale_shape_manual(values = c(7,10, 15,16)) + coord_flip()

ggsave(tg1,filename = "figures/RelativeAbundance/Genus_species/SSD_16S_Genera.RA_5perc_taxasum.png", width=20, height=23, dpi=600,create.dir = TRUE)

tg2<-ggplot(b.gen.spec_RA_meta[b.gen.spec_RA_meta$Count>0.10,], aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date",values=unique(b.gen.spec_RA_meta$SampDate_Color[order(b.gen.spec_RA_meta$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Salton Sea Dust: Microbial Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 10%")+
  scale_shape_manual(values = c(7,10, 15,16)) + coord_flip()

ggsave(tg2,filename = "figures/RelativeAbundance/Genus_species/SSD_16S_Genera.RA_10perc_taxasum.png", width=20, height=23, dpi=600,create.dir = TRUE)

# Taxonomic summary by Sample ID + Collection Period

# tg1<-ggplot(b.gen.spec_RA_meta[b.gen.spec_RA_meta$Count>0.01,], aes(Genus_species, Count)) +
#   geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Site), size=4, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
#         axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
#   labs(x="Microbial Genera", y="Relative Abundance", title="Salton Sea Dust: Microbial Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 1%")+coord_flip()
#
# ggsave(tg1,filename = "figures/RelativeAbundance/Genus_species/SSD_16S_Genera.RA_1perc_taxasum.png", width=20, height=23, dpi=600,create.dir = TRUE)
#
# tg1a<-ggplot(b.gen.spec_RA_meta[b.gen.spec_RA_meta$Genus_species == "Massilia unknown",], aes(Genus_species, Count)) +
#   geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Site), size=4, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),
#         axis.text.x = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="", y="Relative Abundance", title="Bacterial Genus Massilia Across Samples")
#
# ggsave(tg1a,filename = "figures/RelativeAbundance/Genus_species/SSD_16S_Massilia.RA_only_taxasum.png", width=15, height=10, dpi=600,create.dir = TRUE)
#
# tg1b<-ggplot(b.gen.spec_RA_meta[b.gen.spec_RA_meta$Count>0.02,], aes(Genus_species, Count)) +
#   geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Site), size=4, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
#         axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
#   labs(x="Microbial Genera", y="Relative Abundance", title="Salton Sea Dust: Microbial Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 2%")+coord_flip()
#
# ggsave(tg1b,filename = "figures/RelativeAbundance/Genus_species/SSD_16S_Genera.RA_2perc_taxasum.png", width=15, height=23, dpi=600,create.dir = TRUE)
#
# tg1a2<-ggplot(b.gen.spec_RA_meta[b.gen.spec_RA_meta$Count>0.05,], aes(Genus_species, Count)) +
#   geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Site), size=4, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
#         axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2.5),"cm")) +
#   labs(x="Microbial Genera", y="Relative Abundance", title="Salton Sea Dust: Microbial Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 5%")+coord_flip()
#
# ggsave(tg1a2,filename = "figures/RelativeAbundance/Genus_species/SSD_16S_Genera.RA_5perc_taxasum.png", width=25, height=15, dpi=600,create.dir = TRUE)
#
# tg1b<-ggplot(b.gen.spec_RA_meta[b.gen.spec_RA_meta$Count>0.1,], aes(Genus_species, Count)) +
#   geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Site), size=4, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 3),"cm")) +
#   labs(x="Microbial Genera", y="Relative Abundance", title="Salton Sea Dust: Microbial Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 10%")+coord_flip()
#
# ggsave(tg1b,filename = "figures/RelativeAbundance/Genus_species/SSD_16S_Genera.RA_10perc_taxasum.png", width=18, height=10, dpi=600,create.dir = TRUE)
#
# tg1c<-ggplot(b.gen.spec_RA_meta[b.gen.spec_RA_meta$Count>0.15,], aes(Genus_species, Count)) +
#   geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Site), size=5, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 3),"cm")) +
#   labs(x="Microbial Genera", y="Relative Abundance", title="Salton Sea Dust: Microbial Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 15%")+coord_flip()
#
# ggsave(tg1c,filename = "figures/RelativeAbundance/Genus_species/SSD_16S_Genera.RA_15perc_taxasum.png", width=15, height=10, dpi=600,create.dir = TRUE)

#### Genus (by Kingdom) Relative Abundance ####
# use dcast to count up ASVs within each Genus across all of the samples
head(b.dust.all)
# b.dust.all.g<-subset(b.dust.all, b.dust.all$Genus!="Unknown") # drop unknown genera so they don't skew analyses
# "Unknown" %in% b.dust.all.g$Genus

b.k.genus_counts <- as.data.frame(dcast(b.dust.all.g, SampleID~Genus+Species+Kingdom, value.var="Count", fun.aggregate=sum)) ###
head(b.k.genus_counts) # counts by genus per sample
dim(b.k.genus_counts)
rownames(b.k.genus_counts)<-b.k.genus_counts$SampleID
b.k.genus_counts[1:4,1:4]
b.k.genus_counts<-b.k.genus_counts[,colSums(b.k.genus_counts[,-1])>0] # drop classes that are not represented
dim(b.k.genus_counts) # sanity check that we dropped taxa with no hits

b.k.genus_RelAb<-data.frame(decostand(b.k.genus_counts[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.k.genus_RelAb) # sanity check to make sure the transformation worked!

b.k.genus_RelAb$SampleID<-rownames(b.k.genus_RelAb)
head(b.k.genus_RelAb)
#write.csv(b.k.genus_RelAb,"16S_Genera_Relative_Abundance.csv",row.names=TRUE) # good to save just in case

# melt down relativized data to merge with dust_meta
b.k.genus_m<-melt(b.k.genus_RelAb)

head(b.k.genus_m)
colnames(b.k.genus_m)[which(names(b.k.genus_m) == "variable")] <- "Genus_species_Kingdom"
colnames(b.k.genus_m)[which(names(b.k.genus_m) == "value")] <- "Count"
head(b.k.genus_m) ## relative abundance based on sum of counts by genus!
b.k.genus_m$Genus_species_Kingdom<-gsub("^X.","",b.k.genus_m$Genus_species_Kingdom) # get rid of leading X. in Genus_species_Kingdom names
b.k.genus_m$Genus_species_Kingdom<-gsub("\\.\\."," ",b.k.genus_m$Genus_species_Kingdom) # get rid of .. in species name --> . is regex
b.k.genus_m$Genus_species_Kingdom<-gsub("\\."," ",b.k.genus_m$Genus_species_Kingdom) # get rid of . in species name --> . is regex

#b.k.genus_m$Genus_species_Kingdom<-gsub("_"," ",b.k.genus_m$Genus_species_Kingdom) #
head(b.k.genus_m) ## relative abundance based on sum of counts by genus!

b.k.genus_RA_meta<-merge(b.k.genus_m,dust_meta, by="SampleID")
head(b.k.genus_RA_meta) ## relative abundance based on sum of counts by genus!
max(b.k.genus_RA_meta$Count)
b.k.genus_RA_meta$SampleID = factor(b.k.genus_RA_meta$SampleID, levels=unique(b.k.genus_RA_meta$SampleID[order(b.k.genus_RA_meta$Site,b.k.genus_RA_meta$Seas_Coll_Year)]), ordered=TRUE)
b.k.genus_RA_meta$Sample_Type<-"Dust"

saveRDS(b.k.genus_RA_meta, file = "data/Amplicon/SSD_GenusSpecies_RelativeAbundance_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

# separate genera RelAb data by site for downstream figs
WI.k.gen.RA<-subset(b.k.genus_RA_meta,b.k.genus_RA_meta$Site=="WI")
BDC.k.gen.RA<-subset(b.k.genus_RA_meta,b.k.genus_RA_meta$Site=="BDC")
PD.k.gen.RA<-subset(b.k.genus_RA_meta,b.k.genus_RA_meta$Site=="PD")
DP.k.gen.RA<-subset(b.k.genus_RA_meta,b.k.genus_RA_meta$Site=="DP")

# Barplot by SampleID

b.k.gen_RAall<-ggplot(b.k.genus_RA_meta, aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=10))

ggsave(b.k.gen_RAall,filename = "figures/RelativeAbundance/Genus_by_Kingdom/SSD_16S_Genera.Spec.RA_barplot.png", width=30, height=10, dpi=600,create.dir = TRUE)

b.k.gen_RA0.0<-ggplot(b.k.genus_RA_meta[b.k.genus_RA_meta$Count>0.005,], aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 0.05%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=10))

ggsave(b.k.gen_RA0.0,filename = "figures/RelativeAbundance/Genus_by_Kingdom/SSD_16S_Genera.Spec.RA_barplot.png", width=30, height=10, dpi=600,create.dir = TRUE)

b.k.gen_RA0<-ggplot(b.k.genus_RA_meta[b.k.genus_RA_meta$Count>0.01,], aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 1%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=5))

ggsave(b.k.gen_RA0,filename = "figures/RelativeAbundance/Genus_by_Kingdom/SSD_16S_Genera.Spec.RA_1perc_barplot.png", width=30, height=10, dpi=600,create.dir = TRUE)

b.k.gen_RA0v2<-ggplot(b.k.genus_RA_meta[b.k.genus_RA_meta$Count>0.01,], aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 1%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2))+facet_wrap(vars(Site), scales = "free")

ggsave(b.k.gen_RA0v2,filename = "figures/RelativeAbundance/Genus_by_Kingdom/SSD_16S_Genera.Spec.RA_1perc_barplot_v2.png", width=30, height=20, dpi=600,create.dir = TRUE)

b.k.gen_RA01<-ggplot(b.k.genus_RA_meta[b.k.genus_RA_meta$Count>0.02,], aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 2%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=5))

ggsave(b.k.gen_RA01,filename = "figures/RelativeAbundance/Genus_by_Kingdom/SSD_16S_Genera.Spec.RA_2perc_barplot.png", width=30, height=10, dpi=600,create.dir = TRUE)

b.k.gen_RA01v2<-ggplot(b.k.genus_RA_meta[b.k.genus_RA_meta$Count>0.02,], aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 2%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2))+facet_wrap(vars(Site), scales = "free")

ggsave(b.k.gen_RA01v2,filename = "figures/RelativeAbundance/Genus_by_Kingdom/SSD_16S_Genera.Spec.RA_2perc_barplot_v2.png", width=30, height=20, dpi=600,create.dir = TRUE)

b.k.gen_RA1a<-ggplot(b.k.genus_RA_meta[b.k.genus_RA_meta$Count>0.05,], aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 5%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2))

ggsave(b.k.gen_RA1a,filename = "figures/RelativeAbundance/Genus_by_Kingdom/SSD_16S_Genera.Spec.RA_5perc_barplot_v1.png", width=30, height=20, dpi=600,create.dir = TRUE)

b.k.gen_RA1v2<-ggplot(b.k.genus_RA_meta[b.k.genus_RA_meta$Count>0.05,], aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 5%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2))+facet_wrap(vars(Site), scales = "free")

ggsave(b.k.gen_RA1v2,filename = "figures/RelativeAbundance/Genus_by_Kingdom/SSD_16S_Genera.Spec.RA_5perc_barplot_v2.png", width=30, height=20, dpi=600,create.dir = TRUE)

b.k.gen_RA1v3<-ggplot(b.k.genus_RA_meta[b.k.genus_RA_meta$Count>0.05,], aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 5%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2))+facet_wrap(vars(Site), scales = "free")

ggsave(b.k.gen_RA1v3,filename = "figures/RelativeAbundance/Genus_by_Kingdom/SSD_16S_Genera.Spec.RA_5perc_barplot_v3.png", width=30, height=20, dpi=600,create.dir = TRUE)

b.k.gen_RA2v1<-ggplot(b.k.genus_RA_meta[b.k.genus_RA_meta$Count>0.10,], aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 10%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))

ggsave(b.k.gen_RA2v1,filename = "figures/RelativeAbundance/Genus_by_Kingdom/SSD_16S_Genera.Spec.RA_10perc_barplot_v1.png", width=30, height=20, dpi=600,create.dir = TRUE)

b.k.gen_RA2v2<-ggplot(b.k.genus_RA_meta[b.k.genus_RA_meta$Count>0.10,], aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 10%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+facet_wrap(vars(Site), scales = "free")

ggsave(b.k.gen_RA2v2,filename = "figures/RelativeAbundance/Genus_by_Kingdom/SSD_16S_Genera.Spec.RA_10perc_barplot_v2.png", width=30, height=20, dpi=600,create.dir = TRUE)

b.k.gen_RA2v3<-ggplot(b.k.genus_RA_meta[b.k.genus_RA_meta$Count>0.10,], aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 10%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+facet_wrap(vars(Site), scales = "free")

ggsave(b.k.gen_RA2v3,filename = "figures/RelativeAbundance/Genus_by_Kingdom/SSD_16S_Genera.Spec.RA_10perc_barplot_v3.png", width=30, height=20, dpi=600,create.dir = TRUE)

b.k.gen_RA3a<-ggplot(b.k.genus_RA_meta[b.k.genus_RA_meta$Count>0.15,], aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 15%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))

ggsave(b.k.gen_RA3a,filename = "figures/RelativeAbundance/Genus_by_Kingdom/SSD_16S_Genera.Spec.RA_15perc_barplot_v1.png", width=30, height=20, dpi=600,create.dir = TRUE)

b.k.gen_RA3b<-ggplot(b.k.genus_RA_meta[b.k.genus_RA_meta$Count>0.15,], aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 15%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+facet_wrap(vars(Site), scales = "free")

ggsave(b.k.gen_RA3b,filename = "figures/RelativeAbundance/Genus_by_Kingdom/SSD_16S_Genera.Spec.RA_15perc_barplot_v2.png", width=30, height=20, dpi=600,create.dir = TRUE)

b.k.gen_RA4<-ggplot(b.k.genus_RA_meta[b.k.genus_RA_meta$Count>0.25,], aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 25%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))

ggsave(b.k.gen_RA4,filename = "figures/RelativeAbundance/Genus_by_Kingdom/SSD_16S_Genera.Spec.RA_25perc_barplot.png", width=15, height=10, dpi=600,create.dir = TRUE)

b.k.gen_RA5<-ggplot(b.k.genus_RA_meta[b.k.genus_RA_meta$Count>0.35,], aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes taxa with Relative Abundance > 35%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))

ggsave(b.k.gen_RA5,filename = "figures/RelativeAbundance/Genus_by_Kingdom/SSD_16S_Genera.Spec.RA_35perc_barplot.png", width=12, height=10, dpi=600,create.dir = TRUE)

# prep for heatmap
max(b.k.genus_RA_meta$Count)
mean(b.k.genus_RA_meta$Count)
max(b.k.genus_RA_meta$Count)/2 # what is the mid point of the RA here?

# Heatmap by SampleID

# g.h1<-ggplot(b.k.genus_RA_meta[b.k.genus_RA_meta$Count>0.01,], aes(SampleID, Genus_species_Kingdom, fill= Count)) +geom_tile()+scale_fill_gradient2(low="orange",mid="white",high="purple",midpoint=0.35)+
#   theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Sample ID", y="Microbial Genera", title="Microbial Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 1%",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0))
#
# ggsave(g.h1,filename = "figures/RelativeAbundance/Genus_by_Kingdom/SSD_16S_Genera.RA_1perc_heatmap_A.png", width=20, height=15, dpi=600,create.dir = TRUE)
#
g.h2<-ggplot(b.k.genus_RA_meta[b.k.genus_RA_meta$Count>0.05,], aes(SampleID, Genus_species_Kingdom, fill= Count)) +geom_tile()+scale_fill_gradient2(low="blue",mid="pink",high="red",midpoint=0.3)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample ID", y="Microbial Genera", title="Microbial Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 5%",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0))

ggsave(g.h2,filename = "figures/RelativeAbundance/Genus_by_Kingdom/SSD_16S_Genera.RA_5perc_heatmap_B.png", width=30, height=20, dpi=600,create.dir = TRUE)

b.dust.all[1:4,1:4]

# Taxonomic Summary by Sample ID + Collection Date

ggplot(b.k.genus_RA_meta[b.k.genus_RA_meta$Count>0.01,], aes(Genus_species_Kingdom, Count)) +
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date",values=unique(b.k.genus_RA_meta$SampDate_Color[order(b.k.genus_RA_meta$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Salton Sea Dust: Microbial Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 1%")+
  scale_shape_manual(values = c(7,10, 15,16)) + coord_flip()

tg0<-ggplot(b.k.genus_RA_meta[b.k.genus_RA_meta$Count>0.02,], aes(Genus_species_Kingdom, Count)) +
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date",values=unique(b.k.genus_RA_meta$SampDate_Color[order(b.k.genus_RA_meta$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Salton Sea Dust: Microbial Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 2%")+
  scale_shape_manual(values = c(7,10, 15,16)) + coord_flip()

ggsave(tg0,filename = "figures/RelativeAbundance/Genus_by_Kingdom/SSD_16S_Genera.RA_2perc_taxasum.png", width=20, height=23, dpi=600,create.dir = TRUE)


tg1<-ggplot(b.k.genus_RA_meta[b.k.genus_RA_meta$Count>0.05,], aes(Genus_species_Kingdom, Count)) +
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date",values=unique(b.k.genus_RA_meta$SampDate_Color[order(b.k.genus_RA_meta$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Salton Sea Dust: Microbial Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 5%")+
  scale_shape_manual(values = c(7,10, 15,16)) + coord_flip()

ggsave(tg1,filename = "figures/RelativeAbundance/Genus_by_Kingdom/SSD_16S_Genera.RA_5perc_taxasum.png", width=20, height=23, dpi=600,create.dir = TRUE)

tg2<-ggplot(b.k.genus_RA_meta[b.k.genus_RA_meta$Count>0.10,], aes(Genus_species_Kingdom, Count)) +
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date",values=unique(b.k.genus_RA_meta$SampDate_Color[order(b.k.genus_RA_meta$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Salton Sea Dust: Microbial Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 10%")+
  scale_shape_manual(values = c(7,10, 15,16)) + coord_flip()

ggsave(tg2,filename = "figures/RelativeAbundance/Genus_by_Kingdom/SSD_16S_Genera.RA_10perc_taxasum.png", width=20, height=23, dpi=600,create.dir = TRUE)

# Taxonomic summary by Sample ID + Collection Period

# tg1<-ggplot(b.k.genus_RA_meta[b.k.genus_RA_meta$Count>0.01,], aes(Genus_species_Kingdom, Count)) +
#   geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Site), size=4, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
#         axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
#   labs(x="Microbial Genera", y="Relative Abundance", title="Salton Sea Dust: Microbial Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 1%")+coord_flip()
#
# ggsave(tg1,filename = "figures/RelativeAbundance/Genus_by_Kingdom/SSD_16S_Genera.RA_1perc_taxasum.png", width=20, height=23, dpi=600,create.dir = TRUE)
#
# tg1a<-ggplot(b.k.genus_RA_meta[b.k.genus_RA_meta$Genus_species_Kingdom == "Massilia unknown",], aes(Genus_species_Kingdom, Count)) +
#   geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Site), size=4, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),
#         axis.text.x = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="", y="Relative Abundance", title="Bacterial Genus Massilia Across Samples")
#
# ggsave(tg1a,filename = "figures/RelativeAbundance/Genus_by_Kingdom/SSD_16S_Massilia.RA_only_taxasum.png", width=15, height=10, dpi=600,create.dir = TRUE)
#
# tg1b<-ggplot(b.k.genus_RA_meta[b.k.genus_RA_meta$Count>0.02,], aes(Genus_species_Kingdom, Count)) +
#   geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Site), size=4, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
#         axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
#   labs(x="Microbial Genera", y="Relative Abundance", title="Salton Sea Dust: Microbial Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 2%")+coord_flip()
#
# ggsave(tg1b,filename = "figures/RelativeAbundance/Genus_by_Kingdom/SSD_16S_Genera.RA_2perc_taxasum.png", width=15, height=23, dpi=600,create.dir = TRUE)
#
# tg1a2<-ggplot(b.k.genus_RA_meta[b.k.genus_RA_meta$Count>0.05,], aes(Genus_species_Kingdom, Count)) +
#   geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Site), size=4, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
#         axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2.5),"cm")) +
#   labs(x="Microbial Genera", y="Relative Abundance", title="Salton Sea Dust: Microbial Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 5%")+coord_flip()
#
# ggsave(tg1a2,filename = "figures/RelativeAbundance/Genus_by_Kingdom/SSD_16S_Genera.RA_5perc_taxasum.png", width=25, height=15, dpi=600,create.dir = TRUE)
#
# tg1b<-ggplot(b.k.genus_RA_meta[b.k.genus_RA_meta$Count>0.1,], aes(Genus_species_Kingdom, Count)) +
#   geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Site), size=4, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 3),"cm")) +
#   labs(x="Microbial Genera", y="Relative Abundance", title="Salton Sea Dust: Microbial Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 10%")+coord_flip()
#
# ggsave(tg1b,filename = "figures/RelativeAbundance/Genus_by_Kingdom/SSD_16S_Genera.RA_10perc_taxasum.png", width=18, height=10, dpi=600,create.dir = TRUE)
#
# tg1c<-ggplot(b.k.genus_RA_meta[b.k.genus_RA_meta$Count>0.15,], aes(Genus_species_Kingdom, Count)) +
#   geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Site), size=5, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 3),"cm")) +
#   labs(x="Microbial Genera", y="Relative Abundance", title="Salton Sea Dust: Microbial Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 15%")+coord_flip()
#
# ggsave(tg1c,filename = "figures/RelativeAbundance/Genus_by_Kingdom/SSD_16S_Genera.RA_15perc_taxasum.png", width=15, height=10, dpi=600,create.dir = TRUE)


#### Save Just Relative Abundance Data ####
save.image("data/Amplicon/SSDust_RelAb_Scaled_by_DeploymentDays_Data.Rdata")

#### Core Microbiome Analysis ####
# code for this section is from here: https://microbiome.github.io/tutorials/Core.html
# detection = detection threshold for absence/presence (strictly greater by default)
## aka detection threshold is based on relative abundance
# prevalence: how many samples have this taxa at the detection (relative abundance) threshold or higher?
b.gs.df <- as.data.frame(dcast(b.gen.spec_m, SampleID~Genus_species, value.var="Count", fun.aggregate=sum)) ###
b.gs.df[1:4,1:10] # using relative abundance of genus_species data
rownames(b.gs.df)<-b.gs.df$SampleID

b.gs.mat<-as.matrix(t(b.gs.df[,-1]))
b.gs.mat[1:4,1:4]

# relative population frequencies aka at 1% relative abundance threshold
head(prevalence(b.gs.mat, detection = 1/100, sort = TRUE),10)
# ^ output of prevalence(): For each OTU, the fraction of samples where a given OTU is detected. The output is readily given as a percentage.

# absolute population frequencies based on sample counts
head(prevalence(b.gs.mat, detection = 1/100, sort = TRUE, count = TRUE),10)

core.taxa.names<-core_members(b.gs.mat,detection=0.005,prevalence=40/100)
core.taxa.names # taxa that appear in half the samples

# With compositional (relative) abundances
det <- c(0.1, 0.5, 1, 2, 5, 10, 20)/100
det
prevalences <- seq(.25, 1, .05)
prevalences
#ggplot(d) + geom_point(aes(x, y)) + scale_x_continuous(trans="log10", limits=c(NA,1))


plot_core(b.gs.mat,
          prevalences = prevalences,
          detections = det,
          plot.type = "lineplot") +
  xlab("Relative Abundance (%)")

## a nicer core mircrobiome heatmap...

# create sequence of prevalences for heatmap, going from 0.05 to 1 in increments of 0.05
prevalences <- seq(.25, 1, .05)
prevalences

# create detection levels
det <- c(0.1, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)/100
det
#min-prevalence gets the 100th highest prevalence
core.p <- plot_core(b.gs.mat,
               plot.type = "heatmap",
               colours = brewer.pal(5, "Reds"),
               prevalences = prevalences,
               detections = det,
               min.prevalence = prevalence(b.gs.mat, sort = TRUE)[100]) +
  labs(x = "Detection Threshold\n(Relative Abundance (%))") +
  #Adjusts axis text size and legend bar height
  theme_bw() +
  theme(axis.text.y= element_text(size=14, face="italic"),
        axis.text.x.bottom=element_text(size=14),
        axis.title = element_text(size=16),
        legend.text = element_text(size=10),
        legend.title = element_text(size=14))

print(core.p)

ggsave(core.p,filename = "figures/RelativeAbundance/CoreMicrobiome/SSD_16S_CoreMicrobiomeGenera_heatmap.png", width=20, height=15, dpi=600,create.dir = TRUE)

# add prevalence labels to heatmap
core.p2<-core.p + geom_text(aes(label = round(core.p$data$Prevalence,2)), size = 6) +
  labs(title="Core Microbiome Heatmap",subtitle="Includes Overall Prevalence of Taxa") +
  theme(plot.title = element_text(size=15),plot.subtitle = element_text(size=13))
ggsave(core.p2,filename = "figures/RelativeAbundance/CoreMicrobiome/SSD_16S_CoreMicrobiomeGenera_heatmap2.png", width=22, height=15, dpi=600,create.dir = TRUE)


# pull out core microbiome taxa names
core.taxa.names<-sapply(core.p[["data"]][["Taxa"]],levels)[,1]
b.gen.spec_RA_meta[1:4,]

core.gen.meta<-b.gen.spec_RA_meta[(b.gen.spec_RA_meta$Genus_species %in% core.taxa.names),]

core.barplot1<-ggplot(core.gen.meta, aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Core Salton Sea Dust Microbiome", x="SampleID", y="Relative Abundance", subtitle="",fill="Taxa")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+
  facet_wrap(vars(Site), scales = "free")

ggsave(core.barplot1,filename = "figures/RelativeAbundance/CoreMicrobiome/SSD_16S_CoreDustMicrobiome_bySite_barplot.png", width=20, height=15, dpi=600,create.dir = TRUE)

core.barplot2<-ggplot(core.gen.meta, aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Core Salton Sea Dust Microbiome", x="SampleID", y="Relative Abundance", subtitle="",fill="Taxa")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+
  facet_wrap(vars(SampDate), scales = "free")

ggsave(core.barplot2,filename = "figures/RelativeAbundance/CoreMicrobiome/SSD_16S_CoreDustMicrobiome_bySampDate_barplot.png", width=30, height=35, dpi=600,create.dir = TRUE)

#### Find Unique Genera Per Site ####
head(dust_meta)
head(b.genus_m)

site_list<-data.frame(Site=dust_meta$Site,SampleID=dust_meta$SampleID)

b.g.RA.site<-merge(site_list, b.genus_m,by="SampleID")
head(b.g.RA.site)

# finding shared genera...
n_occur <- data.frame(table(b.g.RA.site$Genus_species)) # find frequencies of genera to see which are shared between sample types
n_occur[n_occur$Freq > 1,] # shows us which genera have a greater frequency than 2
g_shared_site<-g_site_meta[g_site_meta$Genus %in% n_occur$Var1[n_occur$Freq > 1],]
#write.csv(g_shared_site,"16S_Genera_SampleType_Shared.csv",row.names=FALSE)

g_not.shared_site<-subset(g_site_meta, !(g_site_meta$Genus %in% g_shared_site$Genus)) # subset based off of what is NOT in one dataframe from another data frame

WI.g1<-subset(b.g.RA.site,Site=="WI")
PD.g1<-subset(b.g.RA.site,Site=="PD")
BDC.g1<-subset(b.g.RA.site,Site=="BDC")
DP.g1<-subset(b.g.RA.site,Site=="DP")

# pull out taxa only in WI that's not in PD, BDC
WI.g<-subset(WI.g1, !(which(WI.g1$Genus_species %in% PD.g1$Genus_species))) # subset based off of what is NOT in one dataframe from another data frame

## Comparing Genera in Dust vs Seawater
# X Axis Breaks and Labels
lbls = paste0(as.character(c(seq(0.05, 0, -0.01), seq(0.01, 0.05, 0.01)))) # labels
brks=seq(-0.05,0.05,0.01)
g_shared_site1<-subset(g_shared_site, Site!="Soil")
g_shared_site1$Count2 <- ifelse(g_shared_site1$Site == "Seawater", -1*g_shared_site1$Count, g_shared_site1$Count)
g_shared_site1<-g_shared_site1[order(-g_shared_site1$Count2,g_shared_site1$Genus),]
g_shared_site1$GenSamp<-interaction(g_shared_site1$Genus,g_shared_site1$Site)
g_shared_site1$GenSamp<-factor(g_shared_site1$GenSamp, levels=g_shared_site1$GenSamp)
class(g_shared_site1$GenSamp)
g_shared_site1$Genus<-factor(g_shared_site1$Genus, levels=unique(g_shared_site1$Genus[sort(g_shared_site1$GenSamp)]))

share1<-ggplot(g_shared_site1, aes(x = Genus, y = -Count2, fill = Site)) +
  geom_bar(data = subset(g_shared_site1[g_shared_site1$Count2<=-0.0005,], Site == "Seawater"), stat = "identity") +
  geom_bar(data = subset(g_shared_site1[g_shared_site1$Count2>=0.0005,], Site == "Dust"), stat = "identity") +
  coord_flip()+scale_y_continuous(labels = lbls,breaks=brks)+theme_classic()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_site1$Sample_Color[rev(order(g_shared_site1$Site))]))+ylab("Relative Abundance")+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+labs(title="Microbial Genera by Sample Type", subtitle="Includes Taxa with a Relative Abundance of at least 0.05%")

ggsave(share1,filename = "figures/RelativeAbundance/SSD_16S_shared_Genera_WI.v.PD_population.pyramid.png", width=12, height=10, dpi=600,create.dir = TRUE)

#### Shared Genus Relative Abundance ####
# first let's get RelAb of taxa by site
b.g.site <- as.data.frame(dcast(b.dust.all.g, Site+SampDate~Genus+Species, value.var="Count", fun.aggregate=sum)) ###
b.g.site[1:7,1:7] # counts by genus per sample
dim(b.g.site)
rownames(b.g.site)<-interaction(b.g.site$Site, b.g.site$SampDate,sep="-")
b.g.site[1:7,1:7]
b.g.site<-b.g.site[,colSums(b.g.site[,-c(1:2)])>0] # drop genera that are not present
dim(b.g.site) # sanity check that we dropped taxa with no hits

b.g.site_RA<-data.frame(decostand(b.g.site[,-c(1:2)], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.g.site_RA) # sanity check to make sure the transformation worked!

b.g.site_RA$Site.SampDate<-rownames(b.g.site_RA)
b.g.site_RA[1:4,1:4]
#write.csv(b.genus_RelAb,"16S_Genera_Relative_Abundance.csv",row.names=TRUE) # good to save just in case

# melt down relativized data to merge with dust_meta
b.g.site_RA.m<-melt(b.g.site_RA)

head(b.g.site_RA.m)
colnames(b.g.site_RA.m)[which(names(b.g.site_RA.m) == "variable")] <- "Genus_species"
colnames(b.g.site_RA.m)[which(names(b.g.site_RA.m) == "value")] <- "Count"
head(b.g.site_RA.m) ## relative abundance based on sum of counts by genus!
b.g.site_RA.m$Genus_species<-gsub("^X.","",b.g.site_RA.m$Genus_species) # get rid of leading X. in Genus_species names
b.g.site_RA.m$Genus_species<-gsub("\\.\\."," ",b.g.site_RA.m$Genus_species) # get rid of .. in species name --> . is regex
b.g.site_RA.m$Genus_species<-gsub("\\."," ",b.g.site_RA.m$Genus_species) # get rid of . in species name --> . is regex
b.g.site_RA.m$Genus_species<-gsub("_"," ",b.g.site_RA.m$Genus_species) #
head(b.g.site_RA.m) ## relative abundance based on sum of counts by genus!
b.g.site_RA.m <- separate(data=b.g.site_RA.m, col="Site.SampDate", sep="-", into=c("Site","Samp.Date"),remove=FALSE) # use separate to create new columns (Site, SampDate) from 1 column
b.g.site_RA.m$Samp.Date<-factor(b.g.site_RA.m$Samp.Date, levels=c("July.2020","August.2020","October.2020","November.2020",
                                                                      "July.2021","August.2021","September.2021","December.2021"))
b.g.site_RA.m$Site<-factor(b.g.site_RA.m$Site, levels=c("PD","BDC","DP","WI"))
b.g.site_RA.m$Site.SampDate<-factor(b.g.site_RA.m$Site.SampDate, levels=unique(b.g.site_RA.m$Site.SampDate[order(b.g.site_RA.m$Samp.Date,b.g.site_RA.m$Site)]))
b.g.site_RA.m<-b.g.site_RA.m[b.g.site_RA.m$Count>0,]

ggplot(b.g.site_RA.m[b.g.site_RA.m$Count>0.025,], aes(Site.SampDate, Genus_species, fill= Count)) +geom_tile()+scale_fill_gradient2(low="blue",mid="pink",high="red",midpoint=0.3)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample ID", y="Microbial Genera", title="Microbial Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 5%",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0))

ggsave(g.h2,filename = "figures/RelativeAbundance/SSD_16S_Genera.RA_5perc_heatmap_B.png", width=16, height=10, dpi=600,create.dir = TRUE)

# merge dust_meta and RA data
site_meta<-unique(data.frame("Site"=dust_meta$Site, "Site_Color"=dust_meta$Site_Color))
site_meta
g_site_meta<-merge(site_meta,b.g.site_RA.m, by="Site")
g_site_meta<-subset(g_site_meta, Genus_species!="Unknown unknown")
g_site_meta<-subset(g_site_meta, Count!=0)

# finding shared genera...

n_occur <- data.frame(table(g_site_meta$Genus_species)) # find frequencies of genera to see which are shared between site
n_occur[n_occur$Freq > 1,] # shows us which genera have a greater frequency than 2
g_shared_site<-g_site_meta[g_site_meta$Genus_species %in% n_occur$Var1[n_occur$Freq > 1],]
#write.csv(g_shared_site,"16S_Genera_SampleType_Shared.csv",row.names=FALSE)

g_not.shared_site<-subset(g_site_meta, !(g_site_meta$Genus_species %in% g_shared_site$Genus_species)) # subset based off of what is NOT in one dataframe from another data frame

# sh.t1<-ggplot(g_shared_site, aes(Genus_species, Count)) +
#   geom_jitter(aes(color=factor(Site)), size=2, width=0.15, height=0) +
#   scale_color_manual(name ="Site", values=unique(g_shared_site$Site_Color[order(g_shared_site$Site)])) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13, vjust=-1),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Microbial Species", y="Relative Abundance", title="Microbial Species Shared Across Site",subtitle="Only Includes Genera Shared Across All Sites")
#
# ggsave(sh.t1,filename = "figures/RelativeAbundance/SSD_16S_shared_Genera_site_taxasum.png", width=23, height=10, dpi=600,create.dir = TRUE)

ggplot(g_shared_site[g_shared_site$Count>=0.01,], aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(Site)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Site", values=unique(g_shared_site$Site_Color[order(g_shared_site$Site)])) + theme_classic() +
  theme(axis.title.x = element_text(size=13, vjust=-1),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Microbial Species Shared Across Sites",subtitle="Only Includes Genera Shared Across All Sites")

#ggsave(sh.t1a,filename = "figures/RelativeAbundance/SSD_16S_shared_Genera_site_1perc_taxasum.png", width=23, height=10, dpi=600,create.dir = TRUE)

ggplot(g_shared_site[g_shared_site$Count>=0.01,], aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(Site)), size=3, width=0.15, height=0) +
  scale_color_manual(name ="Site", values=unique(g_shared_site$Site_Color[order(g_shared_site$Site)])) + theme_classic() +
  geom_boxplot(fill=NA, outlier.color=NA) +
  theme(axis.title.x = element_text(size=13, vjust=-1),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Microbial Species Shared Across Sites")

ggplot(g_shared_site[g_shared_site$Count>=0.0005,], aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(Site)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Site", values=unique(g_shared_site$Site_Color[order(g_shared_site$Site)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13, vjust=-1),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Microbial Species Shared Across Sites")

sh.t2<-ggplot(g_shared_site[g_shared_site$Count>=0.005,], aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(Site)), size=3, width=0.15, height=0) +
  scale_color_manual(name ="Site", values=unique(g_shared_site$Site_Color[order(g_shared_site$Site)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13, vjust=-1),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Microbial Species Shared Across Sites",subtitle="Only Includes Taxa with Relative Abundance > 0.5%")

ggsave(sh.t2,filename = "figures/RelativeAbundance/SSD_16S_shared_Genera_site_0.5perc_taxasum.png", width=23, height=10, dpi=600,create.dir = TRUE)

sh.t3<-ggplot(g_shared_site[g_shared_site$Count>=0.01,], aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(Site)), size=3, width=0.15, height=0) +
  scale_color_manual(name ="Site", values=unique(g_shared_site$Site_Color[order(g_shared_site$Site)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13, vjust=-1),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Microbial Species Shared Across Sites",subtitle="Only Includes Taxa with Relative Abundance > 1%")

ggsave(sh.t3,filename = "figures/RelativeAbundance/SSD_16S_shared_Genera_site_0.25perc_taxasum.png", width=23, height=10, dpi=600,create.dir = TRUE)

## Comparing Shared Genera in WI vs PD
# X Axis Breaks and Labels
lbls = paste0(as.character(c(seq(0.05, 0, -0.01), seq(0.01, 0.05, 0.01)))) # labels
brks=seq(-0.05,0.05,0.01)
g_shared_site1<-subset(g_shared_site, Site=="WI" | Site=="PD")
g_shared_site1$Count2 <- ifelse(g_shared_site1$Site == "WI", -1*g_shared_site1$Count, g_shared_site1$Count)
g_shared_site1<-g_shared_site1[order(-g_shared_site1$Count2,g_shared_site1$Genus_species),]
g_shared_site1$GenSamp<-interaction(g_shared_site1$Genus_species,g_shared_site1$Site)
g_shared_site1$GenSamp<-factor(g_shared_site1$GenSamp, levels=g_shared_site1$GenSamp)
class(g_shared_site1$GenSamp)
g_shared_site1$Genus_species<-factor(g_shared_site1$Genus_species, levels=unique(g_shared_site1$Genus_species[sort(g_shared_site1$GenSamp)]))

share1<-ggplot(g_shared_site1, aes(x = Genus_species, y = -Count2, fill = Site)) +
  geom_bar(data = subset(g_shared_site1[g_shared_site1$Count2<=-0.0005,], Site == "WI"), stat = "identity") +
  geom_bar(data = subset(g_shared_site1[g_shared_site1$Count2>=0.0005,], Site == "PD"), stat = "identity") +
  coord_flip()+scale_y_continuous(labels = lbls,breaks=brks)+theme_classic()+scale_fill_manual(name ="Site", values=unique(g_shared_site1$Site_Color[rev(order(g_shared_site1$Site))]))+ylab("Relative Abundance")+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Shared Microbial Genera - Wister vs. Palm Desert", subtitle="Includes Taxa with a Relative Abundance of at least 0.05%")

ggsave(share1,filename = "figures/RelativeAbundance/SSD_16S_shared_Genera_WI.v.PD_population.pyramid.png", width=12, height=10, dpi=600,create.dir = TRUE)

pp2<-ggplot(g_shared_site1[g_shared_site1$Count>=0.0005,], aes(x = reorder(Genus_species,Count), fill = Site,y = ifelse(test = Site == "PD",yes = Count, no = -Count))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(g_shared_site1$Count) * c(-1,1)) +
  coord_flip()+scale_fill_manual(name ="Site", values=unique(g_shared_site1$Site_Color[order(g_shared_site1$Site)]))+ylab("Relative Abundance")+theme_classic()+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Shared Microbial Genera - Wister vs. Palm Desert",subtitle="Includes Shared Taxa with Relative Abundance by Site of at least 0.05%")+xlab("Genus species")

ggsave(pp2,filename = "figures/RelativeAbundance/SSD_16S_shared_Genera_WI.v.PD_population.pyramid.pretty.05percent.png", width=13, height=10, dpi=600,create.dir = TRUE)

pp3<-ggplot(g_shared_site1[g_shared_site1$Count>=0.005,], aes(x = reorder(Genus_species,Count), fill = Site,y = ifelse(test = Site == "PD",yes = Count, no = -Count))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(g_shared_site1$Count) * c(-1,1)) +
  coord_flip()+scale_fill_manual(name ="Site", values=unique(g_shared_site1$Site_Color[order(g_shared_site1$Site)]))+ylab("Relative Abundance")+theme_classic()+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Shared Microbial Genera - Wister vs. Palm Desert",subtitle="Includes Shared Taxa with Relative Abundance by Site of at least 0.5%")+xlab("Genus species")

ggsave(pp3,filename = "figures/RelativeAbundance/SSD_16S_shared_Genera_WI.v.PD_0.5perc.population.pyramid.pretty.png", width=13, height=10, dpi=600,create.dir = TRUE)

pp4<-ggplot(g_shared_site1[g_shared_site1$Count>=0.010,], aes(x = reorder(Genus_species,Count), fill = Site,y = ifelse(test = Site == "PD",yes = Count, no = -Count))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(g_shared_site1$Count) * c(-1,1)) +
  coord_flip()+scale_fill_manual(name ="Site", values=unique(g_shared_site1$Site_Color[order(g_shared_site1$Site)]))+ylab("Relative Abundance")+theme_classic()+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Shared Microbial Genera - Wister vs. Palm Desert",subtitle="Includes Shared Taxa with Relative Abundance by Site of at least 1%")+xlab("Genus species")

ggsave(pp4,filename = "figures/RelativeAbundance/SSD_16S_shared_Genera_WI.v.PD_1perc_population.pyramid.pretty.png", width=13, height=10, dpi=600,create.dir = TRUE)

## Lollipop chart

lg1<-ggplot(g_shared_site1[g_shared_site1$Count>=0.005,], aes(x = reorder(Genus_species,Count),
                                                              y = ifelse(test = Site == "PD",yes = Count, no = -Count),color=g_shared_site1$Site[g_shared_site1$Count>=0.005])) +
  geom_point(stat='identity',size=3)  +
  geom_segment(aes(y = 0,
                   x = Genus_species,
                   yend = ifelse(test = Site == "PD",yes = Count, no = -Count),
                   xend = Genus_species),color = "black") +
  coord_flip()+scale_color_manual(name ="Site", values=unique(g_shared_site1$Site_Color[order(g_shared_site1$Site)]))+ylab("Relative Abundance")+theme_classic()+theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Shared Microbial Genera - Wister vs. Palm Desert",subtitle="Includes Shared Taxa with Relative Abundance by Site of at least 0.5%")+xlab("Genus species")

ggsave(lg1,filename = "figures/RelativeAbundance/SSD_16S_shared_Genera_WI.v.PD_0.5perc_lollipop_chart.png", width=13, height=10, dpi=600,create.dir = TRUE)

lg2<-ggplot(g_shared_site1[g_shared_site1$Count>=0.01,], aes(x = reorder(Genus_species,Count),
                                                              y = ifelse(test = Site == "PD",yes = Count, no = -Count),color=g_shared_site1$Site[g_shared_site1$Count>=0.01])) +
  geom_point(stat='identity',size=3)  +
  geom_segment(aes(y = 0,
                   x = Genus_species,
                   yend = ifelse(test = Site == "PD",yes = Count, no = -Count),
                   xend = Genus_species),color = "black") +
  coord_flip()+scale_color_manual(name ="Site", values=unique(g_shared_site1$Site_Color[order(g_shared_site1$Site)]))+ylab("Relative Abundance")+theme_classic()+theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Shared Microbial Genera - Wister vs. Palm Desert",subtitle="Includes Shared Taxa with Relative Abundance by Site of at least 0.5%")+xlab("Genus species")

ggsave(lg2,filename = "figures/RelativeAbundance/SSD_16S_shared_Genera_WI.v.PD_1perc_lollipop_chart.png", width=13, height=10, dpi=600,create.dir = TRUE)

## Comparing Shared Genera in WI vs BDC
# X Axis Breaks and Labels
lbls = paste0(as.character(c(seq(0.05, 0, -0.01), seq(0.01, 0.05, 0.01)))) # labels
brks=seq(-0.05,0.05,0.01)
g_shared_site2<-subset(g_shared_site, Site=="WI" | Site=="BDC")
g_shared_site2$Count2 <- ifelse(g_shared_site2$Site == "WI", -1*g_shared_site2$Count, g_shared_site2$Count)
g_shared_site2<-g_shared_site2[order(-g_shared_site2$Count2,g_shared_site2$Genus_species),]
g_shared_site2$GenSamp<-interaction(g_shared_site2$Genus_species,g_shared_site2$Site)
g_shared_site2$GenSamp<-factor(g_shared_site2$GenSamp, levels=g_shared_site2$GenSamp)
class(g_shared_site2$GenSamp)
g_shared_site2$Genus_species<-factor(g_shared_site2$Genus_species, levels=unique(g_shared_site2$Genus_species[sort(g_shared_site2$GenSamp)]))

share1<-ggplot(g_shared_site2, aes(x = Genus_species, y = -Count2, fill = Site)) +
  geom_bar(data = subset(g_shared_site2[g_shared_site2$Count2<=-0.0005,], Site == "WI"), stat = "identity") +
  geom_bar(data = subset(g_shared_site2[g_shared_site2$Count2>=0.0005,], Site == "BDC"), stat = "identity") +
  coord_flip()+scale_y_continuous(labels = lbls,breaks=brks)+theme_classic()+scale_fill_manual(name ="Site", values=unique(g_shared_site2$Site_Color[rev(order(g_shared_site2$Site))]))+ylab("Relative Abundance")+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+labs(title="Microbial Genera by Site", subtitle="Includes Taxa with a Relative Abundance of at least 0.05%")

ggsave(share1,filename = "figures/RelativeAbundance/SSD_16S_shared_Genera_WI.v.BDC_population.pyramid.png", width=12, height=10, dpi=600,create.dir = TRUE)

# pp2<-ggplot(g_shared_site2[g_shared_site2$Count>=0.0005,], aes(x = reorder(Genus_species,Count), fill = Site,y = ifelse(test = Site == "BDC",yes = Count, no = -Count))) +
#   geom_bar(stat = "identity") +
#   scale_y_continuous(labels = abs, limits = max(g_shared_site2$Count) * c(-1,1)) +
#   coord_flip()+scale_fill_manual(name ="Site", values=unique(g_shared_site2$Site_Color[order(g_shared_site2$Site)]))+ylab("Relative Abundance")+theme_classic()+
#   theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
#   labs(title="Microbial Genera by Site",subtitle="Includes Shared Taxa with Relative Abundance by Site of at least 0.05%")+xlab("Genus species")
#
# ggsave(pp2,filename = "figures/RelativeAbundance/SSD_16S_shared_Genera_WI.v.BDC_population.pyramid.pretty.05percent.png", width=13, height=10, dpi=600,create.dir = TRUE)

pp3a<-ggplot(g_shared_site2[g_shared_site2$Count>=0.005,], aes(x = reorder(Genus_species,Count), fill = Site,y = ifelse(test = Site == "BDC",yes = Count, no = -Count))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(g_shared_site2$Count) * c(-1,1)) +
  coord_flip()+scale_fill_manual(name ="Site", values=unique(g_shared_site2$Site_Color[order(g_shared_site2$Site)]))+ylab("Relative Abundance")+theme_classic()+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Shared Microbial Genera - Wister vs. Boyd Deep Canyon",subtitle="Includes Shared Taxa with Relative Abundance by Site of at least 0.5%")+xlab("Genus species")

ggsave(pp3a,filename = "figures/RelativeAbundance/SSD_16S_shared_Genera_WI.v.BDC_0.5perc.population.pyramid.pretty.png", width=13, height=10, dpi=600,create.dir = TRUE)

pp4a<-ggplot(g_shared_site2[g_shared_site2$Count>=0.010,], aes(x = reorder(Genus_species,Count), fill = Site,y = ifelse(test = Site == "BDC",yes = Count, no = -Count))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(g_shared_site2$Count) * c(-1,1)) +
  coord_flip()+scale_fill_manual(name ="Site", values=unique(g_shared_site2$Site_Color[order(g_shared_site2$Site)]))+ylab("Relative Abundance")+theme_classic()+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Shared Microbial Genera - Wister vs. Boyd Deep Canyon",subtitle="Includes Shared Taxa with Relative Abundance by Site of at least 1%")+xlab("Genus species")

ggsave(pp4a,filename = "figures/RelativeAbundance/SSD_16S_shared_Genera_WI.v.BDC_1perc_population.pyramid.pretty.png", width=13, height=10, dpi=600,create.dir = TRUE)

## Lollipop chart

lg1a<-ggplot(g_shared_site2[g_shared_site2$Count>=0.005,], aes(x = reorder(Genus_species,Count),
                                                              y = ifelse(test = Site == "BDC",yes = Count, no = -Count),color=g_shared_site2$Site[g_shared_site2$Count>=0.005])) +
  geom_point(stat='identity',size=3)  +
  geom_segment(aes(y = 0,
                   x = Genus_species,
                   yend = ifelse(test = Site == "BDC",yes = Count, no = -Count),
                   xend = Genus_species),color = "black") +
  coord_flip()+scale_color_manual(name ="Site", values=unique(g_shared_site2$Site_Color[order(g_shared_site2$Site)]))+ylab("Relative Abundance")+theme_classic()+theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Shared Microbial Genera - Wister vs. Boyd Deep Canyon",subtitle="Includes Shared Taxa with Relative Abundance by Site of at least 0.5%")+xlab("Genus species")

ggsave(lg1a,filename = "figures/RelativeAbundance/SSD_16S_shared_Genera_WI.v.BDC_0.5perc_lollipop_chart.png", width=13, height=10, dpi=600,create.dir = TRUE)

lg2a<-ggplot(g_shared_site2[g_shared_site2$Count>=0.01,], aes(x = reorder(Genus_species,Count),
                                                             y = ifelse(test = Site == "BDC",yes = Count, no = -Count),color=g_shared_site2$Site[g_shared_site2$Count>=0.01])) +
  geom_point(stat='identity',size=3)  +
  geom_segment(aes(y = 0,
                   x = Genus_species,
                   yend = ifelse(test = Site == "BDC",yes = Count, no = -Count),
                   xend = Genus_species),color = "black") +
  coord_flip()+scale_color_manual(name ="Site", values=unique(g_shared_site2$Site_Color[order(g_shared_site2$Site)]))+ylab("Relative Abundance")+theme_classic()+theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Shared Microbial Genera - Wister vs. Boyd Deep Canyon",subtitle="Includes Shared Taxa with Relative Abundance by Site of at least 0.5%")+xlab("Genus species")

ggsave(lg2a,filename = "figures/RelativeAbundance/SSD_16S_shared_Genera_WI.v.BDC_1perc_lollipop_chart.png", width=13, height=10, dpi=600,create.dir = TRUE)

g.typ.05p<-na.omit(subset(g_site_meta, Count>=0.005))

g.t.05<-ggplot(g.typ.05p, aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(Site)), size=2, width=0.15, height=0) + geom_boxplot(fill=NA, outlier.color=NA) +scale_color_manual(name ="Site", values=c(unique(g.typ.05p$Site_Color[order(g.typ.05p$Site)])),c("Seawater","Soil", "Dust")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Bacteria/Archaea & Site")

ggsave(g.t.05,filename = "figures/RelativeAbundance/SSD_16S_Gen.0.5percRA.png", width=15, height=10, dpi=600,create.dir = TRUE)

g.t.05a<-ggplot(g.typ.05p, aes(Genus_species, Count)) +
  geom_jitter(aes(color=ifelse(Count>0.01,factor(Site),"grey")), size=2, width=0.15, height=0) + geom_boxplot(fill=NA, outlier.color=NA) +scale_color_manual(name ="Site", values=c(unique(g.typ.05p$Site_Color[order(g.typ.05p$Site)]),"grey"),c("Seawater","Soil", "Dust","<1% RA")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Bacteria/Archaea & Site")

ggsave(g.t.05a,filename = "figures/RelativeAbundance/SSD_16S_Gen.0.5percRA.v2.png", width=15, height=10, dpi=600,create.dir = TRUE)

g.typ.1p<-subset(g_site_meta, Count>=0.01)
g.typ.1p<-subset(g.typ.1p, Genus_species!="Unknown")

g2<-ggplot(g.typ.1p, aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(Site)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Site", values=unique(g.typ.1p$Site_Color[order(g.typ.1p$Site)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=0.5),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera by Site", subtitle="Includes Taxa with a Relative Abundance of at least 1%")+coord_flip()

ggsave(g2,filename = "figures/RelativeAbundance/SSD_16S_Genera.RA_1percent_v1.png", width=12, height=10, dpi=600,create.dir = TRUE)

b.gen_RA.st<-ggplot(g.typ.1p, aes(x=Site, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genera by Site", x="Site", y="Relative Abundance", fill="Genus_species", subtitle="Includes Taxa with a Relative Abundance of at least 1%")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=0.5),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+guides(fill=guide_legend(ncol=1))

ggsave(b.gen_RA.st,filename = "figures/RelativeAbundance/bacterial_genera_1percent_RA_by_SampleType.png", width=12, height=10, dpi=600,create.dir = TRUE)


head(g_site_meta)
g_site_meta.no0<-subset(g_site_meta, Count!=0)

tg.h1<-ggplot(g_site_meta.no0, aes(Site, Genus, fill= Count)) +geom_tile()+scale_fill_gradient2(low="blue3", mid="white",high="red",midpoint=.025)+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.y=element_text(margin = margin(0,0)),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample Type", y="Microbial Genera", title="Microbial Genera & Sample Type",fill="Relative Abundance")+theme_classic()+scale_x_discrete(expand = c(0,0))

ggsave(tp.h1,filename = "figures/RelativeAbundance/SSD_16S_Phyla.RA_heatmap.png", width=12, height=10, dpi=600,create.dir = TRUE)

#### Look at Specific Shared Taxa Across Sites ####

head(g_site_meta)
length(g_site_meta$Site)

d.gen<-subset(g_site_meta, Site=="Dust")
sw.gen<-subset(g_site_meta, Site=="Seawater")

colnames(d.gen)[which(names(d.gen) == "Count")] <- "D.RA"
d.gen$D.RA<-as.numeric(d.gen$D.RA)
#d.gen<-subset(d.gen, D.RA>0.0001)

colnames(sw.gen)[which(names(sw.gen) == "Count")] <- "SW.RA"
sw.gen$SW.RA<-as.numeric(sw.gen$SW.RA)
#sw.gen<-subset(sw.gen, SW.RA>0.0001)

head(d.gen)
head(sw.gen)

tgen.comp<-merge(d.gen, sw.gen, by="Genus")
colnames(tgen.comp)[which(names(tgen.comp) == "Site.x")] <- "D_Samp"
colnames(tgen.comp)[which(names(tgen.comp) == "Site.y")] <- "SW_Samp"
colnames(tgen.comp)[which(names(tgen.comp) == "Sample_Color.x")] <- "D_color"
colnames(tgen.comp)[which(names(tgen.comp) == "Sample_Color.y")] <- "SW_color"
write.csv(tgen.comp,"16S_Genera_SampleType_Shared_RA.Separated.csv",row.names=FALSE)

tgen.comp$Genus[tgen.comp$D.RA==max(tgen.comp$D.RA)] # most relatively abundant genus in dust
tgen.comp$Genus[tgen.comp$SW.RA==max(tgen.comp$SW.RA)] # most relatively abundant genus in seawater
tgen.comp$Gen2 <- ifelse(tgen.comp$SW.RA>0.01 | tgen.comp$D.RA>0.01, tgen.comp$Genus, "Other")
unique(tgen.comp$Gen2) # all genera names
length(unique(tgen.comp$Gen2)) # how many colors do we need

gencol = melt(c(Bacillus="darkgreen",Blastococcus="orange",Geodermatophilus="green3", Halomonas="darkslategray4", Marinobacter="blue3", Other="black", Paracoccus="firebrick1", Roseovarius="deeppink3", Salinicoccus="cornflowerblue", Sphingomonas="purple", Spiroplasma="deepskyblue1", Truepera="darkgoldenrod2"))
head(gencol)
gencol$Gen2<-rownames(gencol)
gencol

tgen.comp<-merge(tgen.comp, gencol, by="Gen2")
head(tgen.comp)
tgen.comp$Gen_Col<-as.character(tgen.comp$Gen_Col)

tgen.comp$SumRA<-tgen.comp$SW.RA+tgen.comp$D.RA
tgen.comp<-tgen.comp[order(-tgen.comp$SumRA),]
tgen.comp$SW.RA <- -1*tgen.comp$SW.RA
tgen.comp$Genus<-factor(tgen.comp$Genus, levels=tgen.comp$Genus)
class(tgen.comp$Genus)

test1<-ggplot(tgen.comp, aes(SW.RA, D.RA, color=Gen2))+geom_point(size=2.5) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+scale_colour_manual(values=unique(tgen.comp$Gen_Col[order(tgen.comp$Gen2)]))+
  labs(x="Seawater Relative Abundance", y="Dust Relative Abundance", title="Microbial Genera & Sample Type",subtitle="Points labeled 'Other' include Genera with a Relative Abundance of < 1%")+guides(shape = guide_legend(override.aes = list(size = 3)),col=guide_legend(title="Genus"))

ggsave(test1,filename = "figures/RelativeAbundance/SSD_16S_Gen.RA_sample.type_scatterplot.png", width=12, height=10, dpi=600,create.dir = TRUE)

test2<-ggplot(tgen.comp, aes(SW.RA, D.RA, color=Genus))  +geom_point(aes(color=ifelse(SW.RA>0.01 | D.RA>0.01,Genus,"black")),size=2) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Seawater Relative Abundance", y="Dust Relative Abundance", title="Microbial Genera & Sample Type")

ggsave(test2,filename = "figures/RelativeAbundance/SSD_16S_Gen.RA_site_linear1_colortest.png", width=12, height=10, dpi=600,create.dir = TRUE)

ggplot(tgen.comp, aes(SW.RA, D.RA, color=ifelse(SW.RA>0.01 | D.RA>0.01,Genus,"black")))  +geom_point() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Seawater Relative Abundance", y="Dust Relative Abundance", title="Microbial Genera & Sample Type")

tgc1<-ggplot(tgen.comp, aes(SW.RA, D.RA, color=Genus)) +geom_point(aes(color=D.RA>0.01 | SW.RA>0.01))+ theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Seawater Relative Abundance", y="Dust Relative Abundance", title="Microbial Genera & Sample Type", subtitle="Only shows shared taxa w/ RA > 0.0001")

ggsave(tgc1,filename = "figures/RelativeAbundance/SSD_16S_Gen.RA_site_linear1.png", width=12, height=10, dpi=600,create.dir = TRUE)
## Date of note: 11/1/21 vvv
# 16S_Gen.RA_site_linear1 - cutoff is 0.0001
# 16S_Gen.RA_site_linear2 - cutoff is 0.00001
# 16S_Gen.RA_site_linear 3 - cutoff is 0.000001
# ** only 2 genera shared at 0.001 cutoff: Halomonas, Truepera

ggplot(g_shared_site, aes(x = Genus, y = Count2, fill = Site)) +
  geom_bar(data = subset(g_shared_site, Site == "Dust"), stat = "identity") +
  geom_bar(data = subset(g_shared_site, Site == "Seawater"), stat = "identity") +
  coord_flip()+scale_y_continuous(labels = lbls,breaks=brks)+theme_classic()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_site$Sample_Color[order(g_shared_site$Site)]))+ylab("Relative Abundance")+theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+labs(title="Microbial Genera by Sample Type", subtitle="Includes Taxa with a Relative Abundance of at least 0.05%")




#### Look at Unique Taxa In Sites ####
head(g_not.shared_site) # see section at 1271 to see where this df came from & how it was calculated
# ^ contains unique taxa by site

WI.uniq.b1<-subset(g_not.shared_site,Site=="WI")
WI.g.RA.meta<-subset(b.genus_RA_meta, Site=="WI")
WI.uniq.b.meta<-WI.g.RA.meta[WI.g.RA.meta$Genus_species %in% WI.uniq.b1$Genus_species,]

# Barplot by SampleID

WI.uniq.gen1<-ggplot(WI.uniq.b.meta[WI.uniq.b.meta$Count>0.0025,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genera Unique to Wister", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 0.25%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))

ggsave(WI.uniq.gen1,filename = "figures/RelativeAbundance/Wister/SSD_16S_WI_Genera.Spec.RA_0.25perc_barplot.png", width=15, height=10, dpi=600,create.dir = TRUE)

WI.uniq.gen2<-ggplot(WI.uniq.b.meta[WI.uniq.b.meta$Count>0.005,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genera Unique to Wister", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 0.5%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))

ggsave(WI.uniq.gen2,filename = "figures/RelativeAbundance/Wister/SSD_16S_WI_Genera.Spec.RA_0.5perc_barplot.png", width=15, height=10, dpi=600,create.dir = TRUE)

WI.uniq.gen3<-ggplot(WI.uniq.b.meta[WI.uniq.b.meta$Count>0.01,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genera Unique to Wister", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 1%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))

ggsave(WI.uniq.gen3,filename = "figures/RelativeAbundance/Wister/SSD_16S_WI_Genera.Spec.RA_1perc_barplot.png", width=15, height=10, dpi=600,create.dir = TRUE)

WI.uniq.gen4<-ggplot(WI.uniq.b.meta[WI.uniq.b.meta$Count>0.05,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genera Unique to Wister", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 5%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))

ggsave(WI.uniq.gen4,filename = "figures/RelativeAbundance/Wister/SSD_16S_WI_Genera.Spec.RA_5perc_barplot.png", width=16, height=10, dpi=600,create.dir = TRUE)

# Taxonomic summary by Sample ID + Collection Period

wi.uniq.tg1<-ggplot(WI.uniq.b.meta[WI.uniq.b.meta$Count>0,], aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Site), size=3, width=0.15, height=0) +
  scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera Unique to Wister by Collection",subtitle="Includes taxa with Relative Abundance > 1%")

ggsave(wi.uniq.tg1,filename = "figures/RelativeAbundance/Wister/SSD_16S_WI_Genera.RA_taxasum_1perc.png", width=23, height=10, dpi=600,create.dir = TRUE)

wi.uniq.tg1a<-ggplot(WI.uniq.b.meta[WI.uniq.b.meta$Genus_species == "Massilia unknown",], aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Site), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),
        axis.text.x = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="", y="Relative Abundance", title="Bacterial Genus Massilia Across Samples")

ggsave(wi.uniq.tg1a,filename = "figures/RelativeAbundance/Wister/SSD_16S_WI_Massilia.RA_only_taxasum.png", width=15, height=10, dpi=600,create.dir = TRUE)

wi.uniq.tg1a2<-ggplot(WI.uniq.b.meta[WI.uniq.b.meta$Count>0.05,], aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Site), size=3, width=0.15, height=0) +
  scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2.5),"cm")) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera Unique to Wister by Collection",subtitle="Includes taxa with Relative Abundance > 5%")

ggsave(wi.uniq.tg1a2,filename = "figures/RelativeAbundance/Wister/SSD_16S_WI_Genera.RA_taxasum_5perc.png", width=25, height=10, dpi=600,create.dir = TRUE)

wi.uniq.tg1b<-ggplot(WI.uniq.b.meta[WI.uniq.b.meta$Count>0.1,], aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Site), size=3, width=0.15, height=0) +
  scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 3),"cm")) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera Unique to Wister by Collection",subtitle="Includes taxa with Relative Abundance > 10%")

ggsave(wi.uniq.tg1b,filename = "figures/RelativeAbundance/Wister/SSD_16S_WI_Genera.RA_taxasum_10perc.png", width=18, height=10, dpi=600,create.dir = TRUE)

wi.uniq.tg1c<-ggplot(WI.uniq.b.meta[WI.uniq.b.meta$Count>0.15,], aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Site), size=3, width=0.15, height=0) +
  scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 3),"cm")) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera Unique to Wister by Collection",subtitle="Includes taxa with Relative Abundance > 15%")

ggsave(wi.uniq.tg1c,filename = "figures/RelativeAbundance/Wister/SSD_16S_WI_Genera.RA_taxasum_15perc.png", width=15, height=10, dpi=600,create.dir = TRUE)

#### Save Everything ####
save.image("data/Amplicon/SSDust_RelAb_DataAll.Rdata")

# Version Information
sessionInfo()