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
load("data/Amplicon/SSDust_RelAb_Scaled_by_DeploymentDays_Data.Rdata") # save global env to Rdata file

head(b.dust.all.g)
b.genus_RelAb[1:4,1:4]
head(meta.all.scaled)

rownames(meta.all.scaled)<-meta.all.scaled$SampleID

#### Subset Specific Taxa's Relative Abundance ####

# first by phylum of interest
phyOI.RelAb<-subset(b.phyla_RelAb,select=c("Actinobacteriota","Bacteroidota","Firmicutes","Proteobacteria",
                                           "SampleID"))
rownames(phyOI.RelAb)<-phyOI.RelAb$SampleID

# then by class of interest
clsOI.RelAb<-subset(b.class_RelAb,select=c("Actinobacteria","Alphaproteobacteria","Bacilli","Gammaproteobacteria",
                                           "SampleID"))
rownames(clsOI.RelAb)<-clsOI.RelAb$SampleID

# then by family of interest
famOI.RelAb<-subset(b.fam_RelAb,select=c("Oxalobacteraceae","Moraxellaceae","Paenibacillaceae","Pseudomonadaceae",
                                         "Bacillaceae","Alicyclobacillaceae","Spirosomaceae","Corynebacteriaceae",
                                         "Staphylococcaceae","SampleID"))
rownames(famOI.RelAb)<-famOI.RelAb$SampleID

# then by genus of interest
genOI.RelAb<-subset(b.genus_RelAb,select=c("Massilia","Paenibacillus","Acinetobacter","Bacillus",
                                           "Pseudomonas","Tumebacillus","Dyadobacter","Desulfovibrio","Spirosoma",
                                           "Corynebacterium","Salinicoccus","SampleID"))
rownames(genOI.RelAb)<-genOI.RelAb$SampleID

# then by CORE genera
Core.genOI.RelAb<-subset(b.genus_RelAb,select=c("Massilia","Sphingomonas","Planomicrobium","Planococcus",
                                           "Hymenobacter","Devosia","Allorhizobium.Neorhizobium.Pararhizobium.Rhizobium",
                                           "Nibribacter","Roseomonas",
                                           "Pseudomonas","Novosphingobium","Kocuria","Paracoccus","SampleID"))
rownames(Core.genOI.RelAb)<-Core.genOI.RelAb$SampleID

#### Prep Env Data for Downstream Analyses ####

# prep env data for downstream visualizations/analyses
env.vars<-meta.all.scaled[,c(4,7:8,11:12,37:46)] # subset env var data
env.vars.names<-names(meta.all.scaled[,c(4,7:8,11:12,37:46)]) # create vector list of variable names

#### Functions to Generate Scatter Plots ####

## Loop to Generate Scatter Plots
### comparing y ~ x where y = fam relative abundance and x = whatever environmental variable

# first make the scatter plot function
# scatter plot function with preset df, y var, color var, shape var
scat.plot.fxn<-function(df, x_var, y_var){
  scat.plt<-ggplot(df,aes(x=x_var, y=y_var, col=SampDate,shape=Site))+
    geom_jitter(aes(color=factor(SampDate),shape=Site), size=5, width=0.15, height=0) +
    scale_color_manual(name ="Collection Date", values=c("green2","orange","red","purple","darkgreen","darkorange3","red4","purple4"),
                       labels=c("July 2020","August 2020","October 2020","November 2020",
                                "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
    scale_shape_manual(name="Site",values = c(7,10, 15,16)) +
    theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
          axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
    labs(x = "",
         y = "",
         shape="",
         color = "R")

  return(scat.plt)
}

# Similar plot to above but x=SampDate for every plot generated:
# first here is a plot function to use that already has the parameters we want for our figure (and the SampDate var):
date.scat.plot.fxn<-function(df, y_var){
  date.scat.plt<-ggplot(df,aes(x=SampDate, y=y_var, col=SampDate,shape=Site))+
    geom_jitter(aes(color=factor(SampDate),shape=Site), size=5, width=0.15, height=0) +
    scale_color_manual(name ="Collection Date", values=c("green2","orange","red","purple","darkgreen","darkorange3","red4","purple4"),
                       labels=c("July 2020","August 2020","October 2020","November 2020",
                                "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
    scale_shape_manual(name="Site",values = c(7,10, 15,16)) +
    theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
          axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
    labs(x = "",
         y = "",
         shape="",
         color = "R")

  return(date.scat.plt)
}

#### Visualize Phyla Of Interest Rel Ab vs Env Vars ####

# merge families & scaled metadata together
phyOI.meta.all<-merge(meta.all.scaled,phyOI.RelAb,by="SampleID")
head(phyOI.meta.all)

# first let's do some basic plotting and see what looks important...
# reminder:

ggplot(data=phyOI.meta.all,aes(x=ave.relative_humidity, y=Proteobacteria, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=5,width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("green2","orange","red","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","October 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(7,10, 15,16)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Ave Relative Humidity", y="Relative Abundance", title="Proteobacteria x Average Relative Humidity",subtitle="Using Scaled Climate Data")


# then make empty list, object with famera and variable names, and plotnum variable
## plotnum will count the # of plots
phylum.scatter.plot.list<-list() # create empty list for each plot to be stored in
phylum.names<-names(phyOI.RelAb)[!names(phyOI.RelAb) %in% c("SampleID")] # pull out names of columns in df that contain "Axis" in name
env.vars.names # sanity check
plotnum <- 1 # counting our model numbers for indexes purposes in the loop

for (i in phylum.names) {
  for (j in env.vars.names){
    phyOI.scat.plot=scat.plot.fxn(phyOI.meta.all,phyOI.meta.all[,j],phyOI.meta.all[,i])
    plot.titled = phyOI.scat.plot + ggtitle(paste(as.character(i),"Relative Abundance ~",as.character(j),sep=" "),subtitle="Using Scaled Climate Data")
    phylum.scatter.plot.list[[plotnum]]=plot(plot.titled)
    ggsave(plot.titled,filename = paste("figures/RelativeAbundance/EnvVarCorrs/Phylum/SSD_",as.character(i),"~",as.character(j),"_scatterplot.png",sep=""), width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).
    #print(names(phyOI.meta.all[i]))
    #print(phyOI.meta.all[j])
    plotnum<-1+plotnum
  }
}
# check if our loop + functions worked! several examples below!
phylum.scatter.plot.list[[2]]

ggplot(data=phyOI.meta.all,aes(x=ave.air_temp, y=Actinobacteriota, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("green2","orange","red","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","October 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(7,10, 15,16)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Ave Wind Speed", y="Relative Abundance", title=paste("Actinobacteria ~ Ave Wind Speed",sep=" "),subtitle="Using Scaled Climate Data")

phylum.scatter.plot.list[[53]]

ggplot(data=phyOI.meta.all,aes(x=Developed, y=Proteobacteria, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("green2","orange","red","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","October 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(7,10, 15,16)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Developed", y="Relative Abundance", title=paste("Firmicutes ~ Developed STF",sep=" "),subtitle="Using Scaled Climate Data")

phylum.scatter.plot.list[[25]]

ggplot(data=phyOI.meta.all,aes(x=Herbaceous, y=Bacteroidota, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("green2","orange","red","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","October 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(7,10, 15,16)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Herbaceous", y="Relative Abundance", title=paste("Bacteroidota ~ Herbaceous STF",sep=" "),subtitle="Using Scaled Climate Data")

### now visualize genus by collection date via a loop!


# make empty list, object with genera and variable names, and plotnum variable
## plotnum will count the # of plots
phylum.date.plot.list<-list() # create empty list for each plot to be stored in
phylum.names<-names(phyOI.RelAb)[!names(phyOI.RelAb) %in% c("SampleID")] # pull out names of columns in df that contain "Axis" in name
#env.vars.names<-names(env.vars)
plotnum <- 1 # counting our model numbers for indexes purposes in the loop

for (i in phylum.names) {
  phylum.date.scat.plot=date.scat.plot.fxn(phyOI.meta.all,phyOI.meta.all[,i])
  date.plot.titled = phylum.date.scat.plot + ggtitle(paste(as.character(i),"Relative Abundance by Collection Date",sep=" "))
  phylum.date.plot.list[[plotnum]]=plot(date.plot.titled)
  ggsave(date.plot.titled,filename = paste("figures/RelativeAbundance/EnvVarCorrs/Phylum/SSD_",as.character(i),"_by_SampDate_scatterplot.png",sep=""), width=18, height=13, dpi=600,create.dir=TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).
  #print(names(phyOI.meta.all[i]))
  #print(phyOI.meta.all[j])
  plotnum<-1+plotnum

}

#### Correlate Phyla of Interest with Env Vars ####
head(meta.all.scaled)
head(phyOI.RelAb)

env.vars<-meta.all.scaled[,c(4,7:8,11:12,37:46)] # subset env var data

dim(env.vars) # confirming that both data frames have the same # of rows
dim(phyOI.RelAb)

rownames(env.vars) # check rownames to see if they are in the same order in both data frames
rownames(phyOI.RelAb)

# reorder data frames so they are in the same order by row (SampleID)
env.vars=env.vars[rownames(phyOI.RelAb),] ## reorder metadata to match order of CLR data

rownames(env.vars) # check rownames to see if they are in the same order in both data frames after reordering
rownames(phyOI.RelAb)

# figure out how we access p value
cor.ex<-cor.test(phyOI.RelAb$Proteobacteria, env.vars$ave.wind_direction, method="pearson", alternative="two.sided")
cor.ex$p.value
cor.ex$estimate
p.adjust(cor.ex$p.value,method="bonferroni",n=cor.ex$parameter)

multi.univar.phylum.corr.fxn<-function(dep.var.df,indep.var.df){
  # create empty lists to store stuff & model number (cornum) to keep track of models each iteration of loop in fxn
  corr_<- vector('list', ncol(dep.var.df) * ncol(indep.var.df)) # create empty list where the corr output is stored
  results_<- vector('list', ncol(dep.var.df) * ncol(indep.var.df)) # create an empty list where the corr summaries are stored
  sig.results<-vector('list', ncol(dep.var.df) * ncol(indep.var.df))
  cornum <- 1 # counting our model numbers for indexes purposes in the loop

  # run the nested loop that famerates corrs from each data frame
  ## dep.var.df[i] is dependent variable (y), indep.var.df[j] is independent variable (x) in corr
  for (i in 1:ncol(dep.var.df)){ # for each column in dep.var.df
    for (j in 1:ncol(indep.var.df)){ # for each column in indep.var.df
      # cor.test(x, y)
      corr_[[cornum]] <- cor.test(indep.var.df[,j],dep.var.df[,i], method="pearson")
      #corr_[[cornum]]$p.value<-p.adjust(corr_[[cornum]]$p.value,method="bonferroni",n=corr_[[cornum]]$parameter)
      results_[[cornum]] <-corr_[[cornum]] # save results of corr into list called results
      names(results_)[cornum]<-paste(names(dep.var.df)[i],"~",names(indep.var.df)[j]) # rename list element to contain the name of the columns used in the model

      # save only significant corrs to another list called sig.results
      ## if p-value < 0.05, save to sig.results list
      ifelse(results_[[cornum]]$p.value < 0.05, sig.results[[cornum]]<-results_[[cornum]], "Not Sig")
      names(sig.results)[cornum]<-paste(names(dep.var.df)[i],"~",names(indep.var.df)[j])
      cornum <- cornum + 1 # add 1 to modelnumber so we keep track of # of models (for indexing purposes in list)

    }
  }

  # drop all NULL elements from sig.results list so it only includes significant corrs
  sig.results[sapply(sig.results, is.null)] <- NULL

  # assign lists to global env so they are saved there are function ends
  assign("env.phylum.results.corrs", results_,envir = .GlobalEnv)
  assign("env.phylum.sig.results.corrs", sig.results,envir = .GlobalEnv)

}

multi.univar.phylum.corr.fxn(phyOI.RelAb[,!colnames(phyOI.RelAb)=="SampleID"],env.vars) # test the function!


#### Do Phyla of Interest Significantly Vary by Group ####
# phylumes are probably not normally distributed but let's just check anyway...
shapiro.test(phyOI.RelAb$Proteobacteria) # what is the p-value?
# W = 0.88324, p-value = 0.004736
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(phyOI.RelAb$Proteobacteria, col="blue") # with outliars

# visualize Q-Q plot for alpha div
# The Q-Q plot, or quantile-quantile plot, is a graphical tool to help us assess if a set of data plausibly came from some theoretical distribution such as a normal or exponential.
# For example, if we run a statistical analysis that assumes our residuals are normally distributed, we can use a normal Q-Q plot to check that assumption
# more on Q-Q plots here: https://data.library.virginia.edu/understanding-q-q-plots/
# more here too: https://grodri.github.io/glms/notes/c2s9#:~:text=8%20The%20Q%2DQ%20Plot,versus%20quantiles%20of%20a%20distribution.
qqnorm(phyOI.RelAb$Proteobacteria, pch = 1, frame = FALSE)
qqline(phyOI.RelAb$Proteobacteria, col = "red", lwd = 2)

# create columns in phyOI.meta.all for stats tests
phyOI.meta.all$SiteColYr<-interaction(phyOI.meta.all$Site,phyOI.meta.all$CollectionYear,sep=".")
phyOI.meta.all$SiteColMonth<-interaction(phyOI.meta.all$Site,phyOI.meta.all$SampleMonth,sep=".")
phyOI.meta.all$SiteSampDate<-interaction(phyOI.meta.all$Site,phyOI.meta.all$SampDate,sep=".")

# use the following statisitcal tests for variance comparisons
## Kruskal: are variances significantly different between groups
## Dunn test: which groups' variances are significant different from one another
## Fligner test: is variance homogenous aka equal across samples?

head(phyOI.meta.all)

## Proteobacteria first!!!!
# Kruskal-Wallis test is an ANOVA for non-normal data
proteobacteria.fit1<-kruskal.test(Proteobacteria ~ Site, data=phyOI.meta.all)

proteobacteria.fit1

# Instead of using Tukey test, we can use Dunn's test to see which groups significantly vary if Kruskal-Wallis test is significant
# ANOVA + Tukey for normally distributed data, Kruskal-Wallis + Dunn's test for non-normal data
rstatix::dunn_test(phyOI.meta.all, Proteobacteria ~ Site, p.adjust.method = "bonferroni", detailed = TRUE)

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(Proteobacteria ~ Site, data = phyOI.meta.all)
# Fligner-Killeen:med chi-squared = 2.7663, df = 3, p-value = 0.4291
# Which shows that the data do NOT deviate significantly from homogeneity.

# compare means
compare_means(Proteobacteria ~ Site, data=phyOI.meta.all, method="wilcox.test",p.adjust.method = "fdr")

# compare variance
compare_means(Proteobacteria ~ Site, data=phyOI.meta.all, method="kruskal.test",p.adjust.method = "bonferroni")

# what about Proteobacteria & Site + Year
proteobacteria.fit2<-kruskal.test(Proteobacteria ~ SiteColYr, data=phyOI.meta.all)

proteobacteria.fit2

# Instead of using Tukey test, we can use Dunn's test to see which groups significantly vary if Kruskal-Wallis test is significant
# ANOVA + Tukey for normally distributed data, Kruskal-Wallis + Dunn's test for non-normal data
proteobacteria.fit2.grp<-rstatix::dunn_test(phyOI.meta.all, Proteobacteria ~ SiteColYr, p.adjust.method = "bonferroni", detailed = TRUE)

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(Proteobacteria ~ SiteColYr, data = phyOI.meta.all)
# Fligner-Killeen:med chi-squared = 4.7864, df = 7, p-value = 0.686
# Which shows that the data do NOT deviate significantly from homogeneity.

proteobacteria.wilcox<-compare_means(Proteobacteria ~ SiteColYr, data=phyOI.meta.all, method="wilcox.test",p.adjust.method = "bonferroni")

compare_means(Proteobacteria ~ SiteColYr, data=phyOI.meta.all, method="kruskal.test",p.adjust.method = "bonferroni")

## Next Firmicutes!
# Kruskal-Wallis test is an ANOVA for non-normal data
firmicutes.fit1<-kruskal.test(Firmicutes ~ Site, data=phyOI.meta.all)

firmicutes.fit1
# Kruskal-Wallis chi-squared = 9.5109, df = 3, p-value = 0.02322

# Instead of using Tukey test, we can use Dunn's test to see which groups significantly vary if Kruskal-Wallis test is significant
# ANOVA + Tukey for normally distributed data, Kruskal-Wallis + Dunn's test for non-normal data
rstatix::dunn_test(phyOI.meta.all, Firmicutes ~ Site, p.adjust.method = "bonferroni", detailed = TRUE)
# DP & WI are significantly different in Firmicutes rel ab

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(Firmicutes ~ Site, data = phyOI.meta.all)
# Fligner-Killeen:med chi-squared = 5.1329, df = 3, p-value = 0.1623
# Which shows that the data do NOT deviate significantly from homogeneity.

compare_means(Firmicutes ~ Site, data=phyOI.meta.all, method="wilcox.test",p.adjust.method = "bonferroni")

compare_means(Firmicutes ~ Site, data=phyOI.meta.all, method="kruskal.test",p.adjust.method = "bonferroni")

# what about Firmicutes & Sate + Year
firmicutes.fit2<-kruskal.test(Firmicutes ~ SiteColYr, data=phyOI.meta.all)

firmicutes.fit2
#Kruskal-Wallis chi-squared = 18.914, df = 7, p-value = 0.008462

# Instead of using Tukey test, we can use Dunn's test to see which groups significantly vary if Kruskal-Wallis test is significant
# ANOVA + Tukey for normally distributed data, Kruskal-Wallis + Dunn's test for non-normal data
firmicutes.fit2.grp<-rstatix::dunn_test(phyOI.meta.all, Firmicutes ~ SiteColYr, p.adjust.method = "bonferroni", detailed = TRUE)

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(Firmicutes ~ SiteColYr, data = phyOI.meta.all)
# Fligner-Killeen:med chi-squared = 3.9144, df = 7, p-value = 0.7896
# Which shows that the data do NOT deviate significantly from homogeneity.

compare_means(Firmicutes ~ SiteColYr, data=phyOI.meta.all, method="wilcox.test",p.adjust.method = "bonferroni")

compare_means(Firmicutes ~ SiteColYr, data=phyOI.meta.all, method="kruskal.test",p.adjust.method = "bonferroni")

#### Visualize Class Of Interest Rel Ab vs Env Vars ####

# merge families & scaled metadata together
clsOI.meta.all<-merge(meta.all.scaled,clsOI.RelAb,by="SampleID")
head(clsOI.meta.all)

# first let's do some basic plotting and see what looks important...
# reminder:

ggplot(data=clsOI.meta.all,aes(x=ave.relative_humidity, y=Bacilli, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=5,width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("green2","orange","red","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","October 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(7,10, 15,16)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Ave Relative Humidity", y="Relative Abundance", title="Bacilli x Average Relative Humidity",subtitle="Using Scaled Climate Data")


# then make empty list, object with famera and variable names, and plotnum variable
## plotnum will count the # of plots
class.scatter.plot.list<-list() # create empty list for each plot to be stored in
class.names<-names(clsOI.RelAb)[!names(clsOI.RelAb) %in% c("SampleID")] # pull out names of columns in df that contain "Axis" in name
env.vars.names # sanity check
plotnum <- 1 # counting our model numbers for indexes purposes in the loop

for (i in class.names) {
  for (j in env.vars.names){
    clsOI.scat.plot=scat.plot.fxn(clsOI.meta.all,clsOI.meta.all[,j],clsOI.meta.all[,i])
    plot.titled = clsOI.scat.plot + ggtitle(paste(as.character(i),"Relative Abundance ~",as.character(j),sep=" "),subtitle="Using Scaled Climate Data")
    class.scatter.plot.list[[plotnum]]=plot(plot.titled)
    ggsave(plot.titled,filename = paste("figures/RelativeAbundance/EnvVarCorrs/Class/SSD_",as.character(i),"~",as.character(j),"_scatterplot.png",sep=""), width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).
    #print(names(clsOI.meta.all[i]))
    #print(clsOI.meta.all[j])
    plotnum<-1+plotnum
  }
}
# check if our loop + functions worked! several examples below!
class.scatter.plot.list[[2]]

ggplot(data=clsOI.meta.all,aes(x=ave.wind_speed, y=Actinobacteria, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("green2","orange","red","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","October 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(7,10, 15,16)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Ave Wind Speed", y="Relative Abundance", title=paste("Actinobacteria ~ Ave Wind Speed",sep=" "),subtitle="Using Scaled Climate Data")

class.scatter.plot.list[[53]]

ggplot(data=clsOI.meta.all,aes(x=OpenWater, y=Gammaproteobacteria, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("green2","orange","red","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","October 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(7,10, 15,16)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Open Water", y="Relative Abundance", title=paste("Gammaproteobacteria ~ Open Water STF",sep=" "),subtitle="Using Scaled Climate Data")

class.scatter.plot.list[[25]]

ggplot(data=clsOI.meta.all,aes(x=OpenWater, y=Alphaproteobacteria, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("green2","orange","red","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","October 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(7,10, 15,16)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Open Water", y="Relative Abundance", title=paste("Alphaproteobacteria ~ Open Water STF",sep=" "),subtitle="Using Scaled Climate Data")

class.scatter.plot.list[[41]]

ggplot(data=clsOI.meta.all,aes(x=SaltonSea, y=Bacilli, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("green2","orange","red","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","October 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(7,10, 15,16)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Salton Sea", y="Relative Abundance", title=paste("Bacillli ~ Salton Sea STF",sep=" "),subtitle="Using Scaled Climate Data")

### now visualize genus by collection date via a loop!


# make empty list, object with genera and variable names, and plotnum variable
## plotnum will count the # of plots
class.date.plot.list<-list() # create empty list for each plot to be stored in
class.names<-names(clsOI.RelAb)[!names(clsOI.RelAb) %in% c("SampleID")] # pull out names of columns in df that contain "Axis" in name
#env.vars.names<-names(env.vars)
plotnum <- 1 # counting our model numbers for indexes purposes in the loop

for (i in class.names) {
  class.date.scat.plot=date.scat.plot.fxn(clsOI.meta.all,clsOI.meta.all[,i])
  date.plot.titled = class.date.scat.plot + ggtitle(paste(as.character(i),"Relative Abundance by Collection Date",sep=" "))
  class.date.plot.list[[plotnum]]=plot(date.plot.titled)
  ggsave(date.plot.titled,filename = paste("figures/RelativeAbundance/EnvVarCorrs/Class/SSD_",as.character(i),"_by_SampDate_scatterplot.png",sep=""), width=18, height=13, dpi=600,create.dir=TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).
  #print(names(clsOI.meta.all[i]))
  #print(clsOI.meta.all[j])
  plotnum<-1+plotnum

}

#### Correlate Classes of Interest with Env Vars ####
head(meta.all.scaled)
head(clsOI.RelAb)

env.vars<-meta.all.scaled[,c(4,7:8,11:12,37:46)] # subset env var data

dim(env.vars) # confirming that both data frames have the same # of rows
dim(clsOI.RelAb)

rownames(env.vars) # check rownames to see if they are in the same order in both data frames
rownames(clsOI.RelAb)

# reorder data frames so they are in the same order by row (SampleID)
env.vars=env.vars[rownames(clsOI.RelAb),] ## reorder metadata to match order of CLR data

rownames(env.vars) # check rownames to see if they are in the same order in both data frames after reordering
rownames(clsOI.RelAb)

# figure out how we access p value
cor.ex<-cor.test(clsOI.RelAb$Bacilli, env.vars$ave.wind_direction, method="pearson", alternative="two.sided")
cor.ex$p.value
cor.ex$estimate
p.adjust(cor.ex$p.value,method="bonferroni",n=cor.ex$parameter)

multi.univar.class.corr.fxn<-function(dep.var.df,indep.var.df){
  # create empty lists to store stuff & model number (cornum) to keep track of models each iteration of loop in fxn
  corr_<- vector('list', ncol(dep.var.df) * ncol(indep.var.df)) # create empty list where the corr output is stored
  results_<- vector('list', ncol(dep.var.df) * ncol(indep.var.df)) # create an empty list where the corr summaries are stored
  sig.results<-vector('list', ncol(dep.var.df) * ncol(indep.var.df))
  cornum <- 1 # counting our model numbers for indexes purposes in the loop

  # run the nested loop that famerates corrs from each data frame
  ## dep.var.df[i] is dependent variable (y), indep.var.df[j] is independent variable (x) in corr
  for (i in 1:ncol(dep.var.df)){ # for each column in dep.var.df
    for (j in 1:ncol(indep.var.df)){ # for each column in indep.var.df
      # cor.test(x, y)
      corr_[[cornum]] <- cor.test(indep.var.df[,j],dep.var.df[,i], method="pearson")
      #corr_[[cornum]]$p.value<-p.adjust(corr_[[cornum]]$p.value,method="bonferroni",n=corr_[[cornum]]$parameter)
      results_[[cornum]] <-corr_[[cornum]] # save results of corr into list called results
      names(results_)[cornum]<-paste(names(dep.var.df)[i],"~",names(indep.var.df)[j]) # rename list element to contain the name of the columns used in the model

      # save only significant corrs to another list called sig.results
      ## if p-value < 0.05, save to sig.results list
      ifelse(results_[[cornum]]$p.value < 0.05, sig.results[[cornum]]<-results_[[cornum]], "Not Sig")
      names(sig.results)[cornum]<-paste(names(dep.var.df)[i],"~",names(indep.var.df)[j])
      cornum <- cornum + 1 # add 1 to modelnumber so we keep track of # of models (for indexing purposes in list)

    }
  }

  # drop all NULL elements from sig.results list so it only includes significant corrs
  sig.results[sapply(sig.results, is.null)] <- NULL

  # assign lists to global env so they are saved there are function ends
  assign("env.class.results.corrs", results_,envir = .GlobalEnv)
  assign("env.class.sig.results.corrs", sig.results,envir = .GlobalEnv)

}

multi.univar.class.corr.fxn(clsOI.RelAb[,!colnames(clsOI.RelAb)=="SampleID"],env.vars) # test the function!


#### Do Classes of Interest Significantly Vary by Group ####
# classes are probably not normally distributed but let's just check anyway...
shapiro.test(clsOI.RelAb$Bacilli) # what is the p-value?
# W = 0.83114, p-value = 0.0003982
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(clsOI.RelAb$Bacilli, col="blue") # with outliars

# visualize Q-Q plot for alpha div
# The Q-Q plot, or quantile-quantile plot, is a graphical tool to help us assess if a set of data plausibly came from some theoretical distribution such as a normal or exponential.
# For example, if we run a statistical analysis that assumes our residuals are normally distributed, we can use a normal Q-Q plot to check that assumption
# more on Q-Q plots here: https://data.library.virginia.edu/understanding-q-q-plots/
# more here too: https://grodri.github.io/glms/notes/c2s9#:~:text=8%20The%20Q%2DQ%20Plot,versus%20quantiles%20of%20a%20distribution.
qqnorm(clsOI.RelAb$Bacilli, pch = 1, frame = FALSE)
qqline(clsOI.RelAb$Bacilli, col = "red", lwd = 2)

# create columns in clsOI.meta.all for stats tests
clsOI.meta.all$SiteColYr<-interaction(clsOI.meta.all$Site,clsOI.meta.all$CollectionYear,sep=".")
clsOI.meta.all$SiteColMonth<-interaction(clsOI.meta.all$Site,clsOI.meta.all$SampleMonth,sep=".")
clsOI.meta.all$SiteSampDate<-interaction(clsOI.meta.all$Site,clsOI.meta.all$SampDate,sep=".")

# use the following statisitcal tests for variance comparisons
## Kruskal: are variances significantly different between groups
## Dunn test: which groups' variances are significant different from one another
## Fligner test: is variance homogenous aka equal across samples?

head(clsOI.meta.all)

## Bacilli first!!!!
# Kruskal-Wallis test is an ANOVA for non-normal data
bacilli.fit1<-kruskal.test(Bacilli ~ Site, data=clsOI.meta.all)

bacilli.fit1

# Instead of using Tukey test, we can use Dunn's test to see which groups significantly vary if Kruskal-Wallis test is significant
# ANOVA + Tukey for normally distributed data, Kruskal-Wallis + Dunn's test for non-normal data
rstatix::dunn_test(clsOI.meta.all, Bacilli ~ Site, p.adjust.method = "bonferroni", detailed = TRUE)

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(Bacilli ~ Site, data = clsOI.meta.all)
# Fligner-Killeen:med chi-squared = 2.7663, df = 3, p-value = 0.4291
# Which shows that the data do NOT deviate significantly from homogeneity.

# compare means
compare_means(Bacilli ~ Site, data=clsOI.meta.all, method="wilcox.test",p.adjust.method = "fdr")

# compare variance
compare_means(Bacilli ~ Site, data=clsOI.meta.all, method="kruskal.test",p.adjust.method = "bonferroni")

# what about Bacilli & Sate + Year
bacilli.fit2<-kruskal.test(Bacilli ~ SiteColYr, data=clsOI.meta.all)

bacilli.fit2

# Instead of using Tukey test, we can use Dunn's test to see which groups significantly vary if Kruskal-Wallis test is significant
# ANOVA + Tukey for normally distributed data, Kruskal-Wallis + Dunn's test for non-normal data
bacilli.fit2.grp<-rstatix::dunn_test(clsOI.meta.all, Bacilli ~ SiteColYr, p.adjust.method = "bonferroni", detailed = TRUE)

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(Bacilli ~ SiteColYr, data = clsOI.meta.all)
# Fligner-Killeen:med chi-squared = 4.7864, df = 7, p-value = 0.686
# Which shows that the data do NOT deviate significantly from homogeneity.

bacilli.wilcox<-compare_means(Bacilli ~ SiteColYr, data=clsOI.meta.all, method="wilcox.test",p.adjust.method = "bonferroni")

compare_means(Bacilli ~ SiteColYr, data=clsOI.meta.all, method="kruskal.test",p.adjust.method = "bonferroni")

## Next Gammaproteobacteria!
# Kruskal-Wallis test is an ANOVA for non-normal data
gammaproteobacteria.fit1<-kruskal.test(Gammaproteobacteria ~ Site, data=clsOI.meta.all)

gammaproteobacteria.fit1

# Instead of using Tukey test, we can use Dunn's test to see which groups significantly vary if Kruskal-Wallis test is significant
# ANOVA + Tukey for normally distributed data, Kruskal-Wallis + Dunn's test for non-normal data
rstatix::dunn_test(clsOI.meta.all, Gammaproteobacteria ~ Site, p.adjust.method = "bonferroni", detailed = TRUE)

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(Gammaproteobacteria ~ Site, data = clsOI.meta.all)
# Fligner-Killeen:med chi-squared = 2.7663, df = 3, p-value = 0.4291
# Which shows that the data do NOT deviate significantly from homogeneity.

compare_means(Gammaproteobacteria ~ Site, data=clsOI.meta.all, method="wilcox.test",p.adjust.method = "bonferroni")

compare_means(Gammaproteobacteria ~ Site, data=clsOI.meta.all, method="kruskal.test",p.adjust.method = "bonferroni")

# what about Gammaproteobacteria & Sate + Year
gammaproteobacteria.fit2<-kruskal.test(Gammaproteobacteria ~ SiteColYr, data=clsOI.meta.all)

gammaproteobacteria.fit2

# Instead of using Tukey test, we can use Dunn's test to see which groups significantly vary if Kruskal-Wallis test is significant
# ANOVA + Tukey for normally distributed data, Kruskal-Wallis + Dunn's test for non-normal data
gammaproteobacteria.fit2.grp<-rstatix::dunn_test(clsOI.meta.all, Gammaproteobacteria ~ SiteColYr, p.adjust.method = "bonferroni", detailed = TRUE)

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(Gammaproteobacteria ~ SiteColYr, data = clsOI.meta.all)
# Fligner-Killeen:med chi-squared = 2.6291, df = 7, p-value = 0.9171
# Which shows that the data do NOT deviate significantly from homogeneity.

compare_means(Gammaproteobacteria ~ SiteColYr, data=clsOI.meta.all, method="wilcox.test",p.adjust.method = "bonferroni")

compare_means(Gammaproteobacteria ~ SiteColYr, data=clsOI.meta.all, method="kruskal.test",p.adjust.method = "bonferroni")

#### Visualize Families Of Interest Rel Ab vs Env Vars ####

# merge families & scaled metadata together
famOI.meta.all<-merge(meta.all.scaled,famOI.RelAb,by="SampleID")
head(famOI.meta.all)

# first let's do some basic plotting and see what looks important...
# reminder:

ggplot(data=famOI.meta.all,aes(x=ave.air_temp, y=Oxalobacteraceae, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=5,width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("green2","orange","red","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","October 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(7,10, 15,16)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Ave Air Temp", y="Relative Abundance", title="Massilia x Average Air Temp",subtitle="Using Scaled Climate Data")


# then make empty list, object with famera and variable names, and plotnum variable
## plotnum will count the # of plots
fam.scatter.plot.list<-list() # create empty list for each plot to be stored in
fam.names<-names(famOI.RelAb)[!names(famOI.RelAb) %in% c("SampleID")] # pull out names of columns in df that contain "Axis" in name
env.vars.names # sanity check
plotnum <- 1 # counting our model numbers for indexes purposes in the loop

for (i in fam.names) {
  for (j in env.vars.names){
    famOI.scat.plot=scat.plot.fxn(famOI.meta.all,famOI.meta.all[,j],famOI.meta.all[,i])
    plot.titled = famOI.scat.plot + ggtitle(paste(as.character(i),"Relative Abundance ~",as.character(j),sep=" "),subtitle="Using Scaled Climate Data")
    fam.scatter.plot.list[[plotnum]]=plot(plot.titled)
    ggsave(plot.titled,filename = paste("figures/RelativeAbundance/EnvVarCorrs/Family/SSD_",as.character(i),"~",as.character(j),"_scatterplot.png",sep=""), width=18, height=13, dpi=600,create.dir=TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).
    #print(names(famOI.meta.all[i]))
    #print(famOI.meta.all[j])
    plotnum<-1+plotnum
  }
}
# check if our loop + functions worked! several examples below!
fam.scatter.plot.list[[2]]

ggplot(data=famOI.meta.all,aes(x=ave.wind_speed, y=Oxalobacteraceae, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("green2","orange","red","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","October 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(7,10, 15,16)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Ave Wind Speed", y="Relative Abundance", title=paste("Oxalobacteraceae ~ Ave Wind Speed",sep=" "),subtitle="Using Scaled Climate Data")

fam.scatter.plot.list[[50]]

ggplot(data=famOI.meta.all,aes(x=Forest, y=Pseudomonadaceae, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("green2","orange","red","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","October 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(7,10, 15,16)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Forest STF", y="Relative Abundance", title=paste("Pseudomonadaceae ~ Forest STF",sep=" "),subtitle="Using Scaled Climate Data")

fam.scatter.plot.list[[70]]

ggplot(data=famOI.meta.all,aes(x=Shrub, y=Bacillaceae, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("green2","orange","red","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","October 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(7,10, 15,16)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Shrub", y="Relative Abundance", title=paste("Bacillaceae ~ Shrub STF",sep=" "),subtitle="Using Scaled Climate Data")

### now visualize family by collection date via a loop!

# make empty list, object with genera and variable names, and plotnum variable
## plotnum will count the # of plots
fam.date.plot.list<-list() # create empty list for each plot to be stored in
fam.names<-names(famOI.RelAb)[!names(famOI.RelAb) %in% c("SampleID")] # pull out names of columns in df that contain "Axis" in name
#env.vars.names<-names(env.vars)
plotnum <- 1 # counting our model numbers for indexes purposes in the loop

for (i in fam.names) {
  fam.date.scat.plot=date.scat.plot.fxn(famOI.meta.all,famOI.meta.all[,i])
  date.plot.titled = fam.date.scat.plot + ggtitle(paste(as.character(i),"Relative Abundance by Collection Date",sep=" "))
  fam.date.plot.list[[plotnum]]=plot(date.plot.titled)
  ggsave(date.plot.titled,filename = paste("figures/RelativeAbundance/EnvVarCorrs/Family/SSD_",as.character(i),"_by_SampDate_scatterplot.png",sep=""), width=18, height=13, dpi=600,create.dir=TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).
  #print(names(famOI.meta.all[i]))
  #print(famOI.meta.all[j])
  plotnum<-1+plotnum

}


#### Correlate Families of Interest with Env Vars ####
head(meta.all.scaled)
head(famOI.RelAb)

env.vars<-meta.all.scaled[,c(4,7:8,11:12,37:46)] # subset env var data

dim(env.vars) # confirming that both data frames have the same # of rows
dim(famOI.RelAb)

rownames(env.vars) # check rownames to see if they are in the same order in both data frames
rownames(famOI.RelAb)

# reorder data frames so they are in the same order by row (SampleID)
env.vars=env.vars[rownames(famOI.RelAb),] ## reorder metadata to match order of CLR data

rownames(env.vars) # check rownames to see if they are in the same order in both data frames after reordering
rownames(famOI.RelAb)

# figure out how we access p value
cor.ex<-cor.test(famOI.RelAb$Bacillaceae, env.vars$ave.wind_direction, method="pearson", alternative="two.sided")
cor.ex$p.value
cor.ex$estimate
p.adjust(cor.ex$p.value,method="bonferroni",n=cor.ex$parameter)

multi.univar.fam.corr.fxn<-function(dep.var.df,indep.var.df){
  # create empty lists to store stuff & model number (cornum) to keep track of models each iteration of loop in fxn
  corr_<- vector('list', ncol(dep.var.df) * ncol(indep.var.df)) # create empty list where the corr output is stored
  results_<- vector('list', ncol(dep.var.df) * ncol(indep.var.df)) # create an empty list where the corr summaries are stored
  sig.results<-vector('list', ncol(dep.var.df) * ncol(indep.var.df))
  cornum <- 1 # counting our model numbers for indexes purposes in the loop

  # run the nested loop that famerates corrs from each data frame
  ## dep.var.df[i] is dependent variable (y), indep.var.df[j] is independent variable (x) in corr
  for (i in 1:ncol(dep.var.df)){ # for each column in dep.var.df
    for (j in 1:ncol(indep.var.df)){ # for each column in indep.var.df
      # cor.test(x, y)
      corr_[[cornum]] <- cor.test(indep.var.df[,j],dep.var.df[,i], method="pearson")
      #corr_[[cornum]]$p.value<-p.adjust(corr_[[cornum]]$p.value,method="bonferroni",n=corr_[[cornum]]$parameter)
      results_[[cornum]] <-corr_[[cornum]] # save results of corr into list called results
      names(results_)[cornum]<-paste(names(dep.var.df)[i],"~",names(indep.var.df)[j]) # rename list element to contain the name of the columns used in the model

      # save only significant corrs to another list called sig.results
      ## if p-value < 0.05, save to sig.results list
      ifelse(results_[[cornum]]$p.value < 0.05, sig.results[[cornum]]<-results_[[cornum]], "Not Sig")
      names(sig.results)[cornum]<-paste(names(dep.var.df)[i],"~",names(indep.var.df)[j])
      cornum <- cornum + 1 # add 1 to modelnumber so we keep track of # of models (for indexing purposes in list)

    }
  }

  # drop all NULL elements from sig.results list so it only includes significant corrs
  sig.results[sapply(sig.results, is.null)] <- NULL

  # assign lists to global env so they are saved there are function ends
  assign("env.fam.results.corrs", results_,envir = .GlobalEnv)
  assign("env.fam.sig.results.corrs", sig.results,envir = .GlobalEnv)

}

multi.univar.fam.corr.fxn(famOI.RelAb[,!colnames(famOI.RelAb)=="SampleID"],env.vars) # test the function!

#### Do Families of Interest Significantly Vary by Group ####
# families are probably not normally distributed but let's just check anyway...
shapiro.test(famOI.RelAb$Bacillaceae) # what is the p-value?
# p-value = 0.000402
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(famOI.RelAb$Bacillaceae, col="blue") # with outliars

# visualize Q-Q plot for alpha div
# The Q-Q plot, or quantile-quantile plot, is a graphical tool to help us assess if a set of data plausibly came from some theoretical distribution such as a normal or exponential.
# For example, if we run a statistical analysis that assumes our residuals are normally distributed, we can use a normal Q-Q plot to check that assumption
# more on Q-Q plots here: https://data.library.virginia.edu/understanding-q-q-plots/
# more here too: https://grodri.github.io/glms/notes/c2s9#:~:text=8%20The%20Q%2DQ%20Plot,versus%20quantiles%20of%20a%20distribution.
qqnorm(famOI.RelAb$Bacillaceae, pch = 1, frame = FALSE)
qqline(famOI.RelAb$Bacillaceae, col = "red", lwd = 2)

# create columns in famOI.meta.all for stats tests
famOI.meta.all$SiteColYr<-interaction(famOI.meta.all$Site,famOI.meta.all$CollectionYear,sep=".")
famOI.meta.all$SiteColMonth<-interaction(famOI.meta.all$Site,famOI.meta.all$SampleMonth,sep=".")
famOI.meta.all$SiteSampDate<-interaction(famOI.meta.all$Site,famOI.meta.all$SampDate,sep=".")

# use the following statisitcal tests for variance comparisons
## Kruskal: are variances significantly different between groups
## Dunn test: which groups' variances are significant different from one another
## Fligner test: is variance homogenous aka equal across samples?

head(famOI.meta.all)

## Bacillaceae first!!!!
# Kruskal-Wallis test is an ANOVA for non-normal data
bacillaceae.fit1<-kruskal.test(Bacillaceae ~ Site, data=famOI.meta.all)

bacillaceae.fit1

# Instead of using Tukey test, we can use Dunn's test to see which groups significantly vary if Kruskal-Wallis test is significant
# ANOVA + Tukey for normally distributed data, Kruskal-Wallis + Dunn's test for non-normal data
rstatix::dunn_test(famOI.meta.all, Bacillaceae ~ Site, p.adjust.method = "bonferroni", detailed = TRUE)

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(Bacillaceae ~ Site, data = famOI.meta.all)
# Fligner-Killeen:med chi-squared = 1.5679, df = 3, p-value = 0.6667
# Which shows that the data do NOT deviate significantly from homogeneity.

compare_means(Bacillaceae ~ Site, data=famOI.meta.all, method="wilcox.test",p.adjust.method = "bonferroni")

compare_means(Bacillaceae ~ Site, data=famOI.meta.all, method="kruskal.test",p.adjust.method = "bonferroni")

# what about Bacillaceae & Sate + Year
bacillaceae.fit2<-kruskal.test(Bacillaceae ~ SiteColYr, data=famOI.meta.all)

bacillaceae.fit2

# Instead of using Tukey test, we can use Dunn's test to see which groups significantly vary if Kruskal-Wallis test is significant
# ANOVA + Tukey for normally distributed data, Kruskal-Wallis + Dunn's test for non-normal data
bacillaceae.fit2.grp<-rstatix::dunn_test(famOI.meta.all, Bacillaceae ~ SiteColYr, p.adjust.method = "bonferroni", detailed = TRUE)

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(Bacillaceae ~ SiteColYr, data = famOI.meta.all)
# Fligner-Killeen:med chi-squared = 4.7864, df = 7, p-value = 0.686
# Which shows that the data do NOT deviate significantly from homogeneity.

compare_means(Bacillaceae ~ SiteColYr, data=famOI.meta.all, method="wilcox.test",p.adjust.method = "bonferroni")

compare_means(Bacillaceae ~ SiteColYr, data=famOI.meta.all, method="kruskal.test",p.adjust.method = "bonferroni")

## Next Oxalobacteraceae!
# Kruskal-Wallis test is an ANOVA for non-normal data
oxalobacteraceae.fit1<-kruskal.test(Oxalobacteraceae ~ Site, data=famOI.meta.all)

oxalobacteraceae.fit1

# Instead of using Tukey test, we can use Dunn's test to see which groups significantly vary if Kruskal-Wallis test is significant
# ANOVA + Tukey for normally distributed data, Kruskal-Wallis + Dunn's test for non-normal data
rstatix::dunn_test(famOI.meta.all, Oxalobacteraceae ~ Site, p.adjust.method = "bonferroni", detailed = TRUE)

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(Oxalobacteraceae ~ Site, data = famOI.meta.all)
# Fligner-Killeen:med chi-squared = 2.7663, df = 3, p-value = 0.4291
# Which shows that the data do NOT deviate significantly from homogeneity.

compare_means(Oxalobacteraceae ~ Site, data=famOI.meta.all, method="wilcox.test",p.adjust.method = "bonferroni")

compare_means(Oxalobacteraceae ~ Site, data=famOI.meta.all, method="kruskal.test",p.adjust.method = "bonferroni")

# what about Oxalobacteraceae & Sate + Year
oxalobacteraceae.fit2<-kruskal.test(Oxalobacteraceae ~ SiteColYr, data=famOI.meta.all)

oxalobacteraceae.fit2

# Instead of using Tukey test, we can use Dunn's test to see which groups significantly vary if Kruskal-Wallis test is significant
# ANOVA + Tukey for normally distributed data, Kruskal-Wallis + Dunn's test for non-normal data
oxalobacteraceae.fit2.grp<-rstatix::dunn_test(famOI.meta.all, Oxalobacteraceae ~ SiteColYr, p.adjust.method = "bonferroni", detailed = TRUE)

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(Oxalobacteraceae ~ SiteColYr, data = famOI.meta.all)
# Fligner-Killeen:med chi-squared = 2.6291, df = 7, p-value = 0.9171
# Which shows that the data do NOT deviate significantly from homogeneity.

compare_means(Oxalobacteraceae ~ SiteColYr, data=famOI.meta.all, method="wilcox.test",p.adjust.method = "bonferroni")

compare_means(Oxalobacteraceae ~ SiteColYr, data=famOI.meta.all, method="kruskal.test",p.adjust.method = "bonferroni")


#### Visualize Genera Of Interest Rel Ab vs Env Vars ####

# merge genus & scaled metadata together
genOI.meta.all<-merge(meta.all.scaled,genOI.RelAb,by="SampleID")
head(genOI.meta.all)

# first let's do some basic plotting and see what looks important...
# reminder:

ggplot(data=genOI.meta.all,aes(x=ave.wind_direction, y=Massilia, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=5,width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("green2","orange","red","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","October 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(7,10, 15,16)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Ave. Wind Direction", y="Relative Abundance", title="Massilia ~ Ave Wind Direction",subtitle="")



# then make empty list, object with genera and variable names, and plotnum variable
## plotnum will count the # of plots
gen.scatter.plot.list<-list() # create empty list for each plot to be stored in
genus.names<-names(genOI.RelAb)[!names(genOI.RelAb) %in% c("SampleID")] # pull out names of columns in df that contain "Axis" in name
env.vars.names<-names(env.vars)
plotnum <- 1 # counting our model numbers for indexes purposes in the loop

# now we are going to loop through the genus names, then the env var names, to generate scatter plots of each
## genus rel ab will be the y axis, env var (scaled) will be the x axis
for (i in genus.names) {
  for (j in env.vars.names){
    genOI.scat.plot=scat.plot.fxn(genOI.meta.all,genOI.meta.all[,j],genOI.meta.all[,i])
    plot.titled = genOI.scat.plot + ggtitle(paste(as.character(i),"Relative Abundance ~",as.character(j),sep=" "),subtitle="Using Scaled Climate Data")
    gen.scatter.plot.list[[plotnum]]=plot(plot.titled)
    ggsave(plot.titled,filename = paste("figures/RelativeAbundance/EnvVarCorrs/Genus/SSD_",as.character(i),"~",as.character(j),"_scatterplot.png",sep=""), width=18, height=13, dpi=600,create.dir=TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).
    #print(names(genOI.meta.all[i]))
    #print(genOI.meta.all[j])
    plotnum<-1+plotnum
  }
}
# check if our loop + functions worked! several examples below!
gen.scatter.plot.list[[2]]

ggplot(data=genOI.meta.all,aes(x=ave.wind_speed, y=Massilia, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("green2","orange","red","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","October 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(7,10, 15,16)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Ave Wind Speed", y="Relative Abundance", title=paste("Massilia ~ Ave Wind Speed",sep=" "),subtitle="Using Scaled Climate Data")

gen.scatter.plot.list[[50]]

ggplot(data=genOI.meta.all,aes(x=Forest, y=Bacillus, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("green2","orange","red","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","October 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(7,10, 15,16)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Ave Wind Direction", y="Relative Abundance", title=paste("Bacillus ~ Forest STF",sep=" "),subtitle="Using Scaled Climate Data")

gen.scatter.plot.list[[25]]

ggplot(data=genOI.meta.all,aes(x=OpenWater, y=Paenibacillus, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("green2","orange","red","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","October 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(7,10, 15,16)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="OpenWater", y="Relative Abundance", title=paste("Paenibacillus ~ OpenWater STF",sep=" "),subtitle="Using Scaled Climate Data")

gen.scatter.plot.list[[70]]

ggplot(data=genOI.meta.all,aes(x=Shrub, y=Pseudomonas, col=SampDate,shape=Site))+
  geom_jitter(aes(color=factor(SampDate),shape=Site), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("green2","orange","red","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","October 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  scale_shape_manual(name="Site",values = c(7,10, 15,16)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Shrub", y="Relative Abundance", title=paste("Pseudomonas ~ Shrub STF",sep=" "),subtitle="Using Scaled Climate Data")


# make empty list, object with genera and variable names, and plotnum variable
## plotnum will count the # of plots
gen.date.plot.list<-list() # create empty list for each plot to be stored in
genus.names<-names(genOI.RelAb)[!names(genOI.RelAb) %in% c("SampleID")] # pull out names of columns in df that contain "Axis" in name
#env.vars.names<-names(env.vars)
plotnum <- 1 # counting our model numbers for indexes purposes in the loop

for (i in genus.names) {
  gen.date.scat.plot=date.scat.plot.fxn(genOI.meta.all,genOI.meta.all[,i])
  date.plot.titled = gen.date.scat.plot + ggtitle(paste(as.character(i),"Relative Abundance by Collection Date",sep=" "))
  gen.date.plot.list[[plotnum]]=plot(date.plot.titled)
  ggsave(date.plot.titled,filename = paste("figures/RelativeAbundance/EnvVarCorrs/Genus/SSD_",as.character(i),"_by_SampDate_scatterplot.png",sep=""), width=18, height=13, dpi=600,create.dir=TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).
  #print(names(genOI.meta.all[i]))
  #print(genOI.meta.all[j])
  plotnum<-1+plotnum

}

#### Correlate Genera of Interest with Env Vars ####
head(meta.all.scaled)
head(genOI.RelAb)

env.vars<-meta.all.scaled[,c(4,7:8,11:12,37:46)] # subset env var data

dim(env.vars) # confirming that both data frames have the same # of rows
dim(genOI.RelAb)

rownames(env.vars) # check rownames to see if they are in the same order in both data frames
rownames(genOI.RelAb)

# reorder data frames so they are in the same order by row (SampleID)
env.vars=env.vars[rownames(genOI.RelAb),] ## reorder metadata to match order of CLR data

rownames(env.vars) # check rownames to see if they are in the same order in both data frames after reordering
rownames(genOI.RelAb)

# figure out how we access p value
cor.ex<-cor.test(genOI.RelAb$Paenibacillus, env.vars$ave.wind_direction, method="pearson", alternative="two.sided")
cor.ex$p.value
cor.ex$estimate
p.adjust(cor.ex$p.value,method="bonferroni",n=cor.ex$parameter)

multi.univar.corr.fxn<-function(dep.var.df,indep.var.df){
  # create empty lists to store stuff & model number (cornum) to keep track of models each iteration of loop in fxn
  corr_<- vector('list', ncol(dep.var.df) * ncol(indep.var.df)) # create empty list where the corr output is stored
  results_<- vector('list', ncol(dep.var.df) * ncol(indep.var.df)) # create an empty list where the corr summaries are stored
  sig.results<-vector('list', ncol(dep.var.df) * ncol(indep.var.df))
  cornum <- 1 # counting our model numbers for indexes purposes in the loop

  # run the nested loop that generates corrs from each data frame
  ## dep.var.df[i] is dependent variable (y), indep.var.df[j] is independent variable (x) in corr
  for (i in 1:ncol(dep.var.df)){ # for each column in dep.var.df
    for (j in 1:ncol(indep.var.df)){ # for each column in indep.var.df
      # cor.test(x, y)
      corr_[[cornum]] <- cor.test(indep.var.df[,j],dep.var.df[,i], method="pearson")
      corr_[[cornum]]$p.value<-p.adjust(corr_[[cornum]]$p.value,method="bonferroni",n=corr_[[cornum]]$parameter)
      results_[[cornum]] <-corr_[[cornum]] # save results of corr into list called results
      names(results_)[cornum]<-paste(names(dep.var.df)[i],"~",names(indep.var.df)[j]) # rename list element to contain the name of the columns used in the model

      # save only significant corrs to another list called sig.results
      ## if p-value < 0.05, save to sig.results list
      ifelse(results_[[cornum]]$p.value < 0.05, sig.results[[cornum]]<-results_[[cornum]], "Not Sig")
      names(sig.results)[cornum]<-paste(names(dep.var.df)[i],"~",names(indep.var.df)[j])
      cornum <- cornum + 1 # add 1 to modelnumber so we keep track of # of models (for indexing purposes in list)

    }
  }

  # drop all NULL elements from sig.results list so it only includes significant corrs
  sig.results[sapply(sig.results, is.null)] <- NULL

  # assign lists to global env so they are saved there are function ends
  assign("env.genus.results.corrs", results_,envir = .GlobalEnv)
  assign("env.genus.sig.results.corrs", sig.results,envir = .GlobalEnv)

}

multi.univar.corr.fxn(genOI.RelAb[,!colnames(genOI.RelAb)=="SampleID"],env.vars) # test the function!


#### Do Genera of Interest Significantly Vary by Group ####
# families are probably not normally distributed but let's just check anyway...
shapiro.test(genOI.RelAb$Bacillus) # what is the p-value?
# p-value = 0.0004967
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(genOI.RelAb$Bacillus, col="blue") # with outliars

# visualize Q-Q plot for alpha div
# The Q-Q plot, or quantile-quantile plot, is a graphical tool to help us assess if a set of data plausibly came from some theoretical distribution such as a normal or exponential.
# For example, if we run a statistical analysis that assumes our residuals are normally distributed, we can use a normal Q-Q plot to check that assumption
# more on Q-Q plots here: https://data.library.virginia.edu/understanding-q-q-plots/
# more here too: https://grodri.github.io/glms/notes/c2s9#:~:text=8%20The%20Q%2DQ%20Plot,versus%20quantiles%20of%20a%20distribution.
qqnorm(genOI.RelAb$Bacillus, pch = 1, frame = FALSE)
qqline(genOI.RelAb$Bacillus, col = "red", lwd = 2)

# create columns in genOI.meta.all for stats tests
genOI.meta.all$SiteColYr<-interaction(genOI.meta.all$Site,genOI.meta.all$CollectionYear,sep=".")
genOI.meta.all$SiteColMonth<-interaction(genOI.meta.all$Site,genOI.meta.all$SampleMonth,sep=".")
genOI.meta.all$SiteSampDate<-interaction(genOI.meta.all$Site,genOI.meta.all$SampDate,sep=".")

# use the following statisitcal tests for variance comparisons
## Kruskal: are variances significantly different between groups
## Dunn test: which groups' variances are significant different from one another
## Fligner test: is variance homogenous aka equal across samples?

head(genOI.meta.all)

## Bacillus first!!!!
# Kruskal-Wallis test is an ANOVA for non-normal data
massilia.fit1<-kruskal.test(Bacillus ~ Site, data=genOI.meta.all)

massilia.fit1

# Instead of using Tukey test, we can use Dunn's test to see which groups significantly vary if Kruskal-Wallis test is significant
# ANOVA + Tukey for normally distributed data, Kruskal-Wallis + Dunn's test for non-normal data
rstatix::dunn_test(genOI.meta.all, Bacillus ~ Site, p.adjust.method = "bonferroni", detailed = TRUE)

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(Bacillus ~ Site, data = genOI.meta.all)
# Fligner-Killeen:med chi-squared = 1.5679, df = 3, p-value = 0.6667
# Which shows that the data do NOT deviate significantly from homogeneity.

compare_means(Bacillus ~ Site, data=genOI.meta.all, method="wilcox.test",p.adjust.method = "bonferroni")

compare_means(Bacillus ~ Site, data=genOI.meta.all, method="kruskal.test",p.adjust.method = "bonferroni")

# what about Bacillus & Sate + Year
massilia.fit2<-kruskal.test(Bacillus ~ SiteColYr, data=genOI.meta.all)

massilia.fit2

# Instead of using Tukey test, we can use Dunn's test to see which groups significantly vary if Kruskal-Wallis test is significant
# ANOVA + Tukey for normally distributed data, Kruskal-Wallis + Dunn's test for non-normal data
massilia.fit2.grp<-rstatix::dunn_test(genOI.meta.all, Bacillus ~ SiteColYr, p.adjust.method = "bonferroni", detailed = TRUE)

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(Bacillus ~ SiteColYr, data = genOI.meta.all)
# Fligner-Killeen:med chi-squared = 4.7864, df = 7, p-value = 0.686
# Which shows that the data do NOT deviate significantly from homogeneity.

compare_means(Bacillus ~ SiteColYr, data=genOI.meta.all, method="wilcox.test",p.adjust.method = "bonferroni")

compare_means(Bacillus ~ SiteColYr, data=genOI.meta.all, method="kruskal.test",p.adjust.method = "bonferroni")

## Next Paenibacillus!
# Kruskal-Wallis test is an ANOVA for non-normal data
sphingomonas.fit1<-kruskal.test(Paenibacillus ~ Site, data=genOI.meta.all)

sphingomonas.fit1

# Instead of using Tukey test, we can use Dunn's test to see which groups significantly vary if Kruskal-Wallis test is significant
# ANOVA + Tukey for normally distributed data, Kruskal-Wallis + Dunn's test for non-normal data
rstatix::dunn_test(genOI.meta.all, Paenibacillus ~ Site, p.adjust.method = "bonferroni", detailed = TRUE)

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(Paenibacillus ~ Site, data = genOI.meta.all)
# Fligner-Killeen:med chi-squared = 3.9001, df = 3, p-value = 0.2724
# Which shows that the data do NOT deviate significantly from homogeneity.

compare_means(Paenibacillus ~ Site, data=genOI.meta.all, method="wilcox.test",p.adjust.method = "bonferroni")

compare_means(Paenibacillus ~ Site, data=genOI.meta.all, method="kruskal.test",p.adjust.method = "bonferroni")

# what about Paenibacillus & Sate + Year
sphingomonas.fit2<-kruskal.test(Paenibacillus ~ SiteColYr, data=genOI.meta.all)

sphingomonas.fit2

# Instead of using Tukey test, we can use Dunn's test to see which groups significantly vary if Kruskal-Wallis test is significant
# ANOVA + Tukey for normally distributed data, Kruskal-Wallis + Dunn's test for non-normal data
sphingomonas.fit2.grp<-rstatix::dunn_test(genOI.meta.all, Paenibacillus ~ SiteColYr, p.adjust.method = "bonferroni", detailed = TRUE)

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(Paenibacillus ~ SiteColYr, data = genOI.meta.all)
# Fligner-Killeen:med chi-squared = 2.6291, df = 7, p-value = 0.9171
# Which shows that the data do NOT deviate significantly from homogeneity.

compare_means(Paenibacillus ~ SiteColYr, data=genOI.meta.all, method="wilcox.test",p.adjust.method = "bonferroni")

compare_means(Paenibacillus ~ SiteColYr, data=genOI.meta.all, method="kruskal.test",p.adjust.method = "bonferroni")



#### Do Core Genera of Interest Significantly Vary by Group ####
# merge genus & scaled metadata together
Core.genOI.meta.all<-merge(meta.all.scaled,Core.genOI.RelAb,by="SampleID")
head(Core.genOI.meta.all)

# families are probably not normally distributed but let's just check anyway...
shapiro.test(Core.genOI.RelAb$Massilia) # what is the p-value?
# p-value = 0.0006307
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(Core.genOI.RelAb$Massilia, col="blue") # with outliars

# visualize Q-Q plot for alpha div
# The Q-Q plot, or quantile-quantile plot, is a graphical tool to help us assess if a set of data plausibly came from some theoretical distribution such as a normal or exponential.
# For example, if we run a statistical analysis that assumes our residuals are normally distributed, we can use a normal Q-Q plot to check that assumption
# more on Q-Q plots here: https://data.library.virginia.edu/understanding-q-q-plots/
# more here too: https://grodri.github.io/glms/notes/c2s9#:~:text=8%20The%20Q%2DQ%20Plot,versus%20quantiles%20of%20a%20distribution.
qqnorm(Core.genOI.RelAb$Massilia, pch = 1, frame = FALSE)
qqline(Core.genOI.RelAb$Massilia, col = "red", lwd = 2)

# create columns in Core.genOI.meta.all for stats tests
Core.genOI.meta.all$SiteColYr<-interaction(Core.genOI.meta.all$Site,Core.genOI.meta.all$CollectionYear,sep=".")
Core.genOI.meta.all$SiteColMonth<-interaction(Core.genOI.meta.all$Site,Core.genOI.meta.all$SampleMonth,sep=".")
Core.genOI.meta.all$SiteSampDate<-interaction(Core.genOI.meta.all$Site,Core.genOI.meta.all$SampDate,sep=".")

# use the following statisitcal tests for variance comparisons
## Kruskal: are variances significantly different between groups
## Dunn test: which groups' variances are significant different from one another
## Fligner test: is variance homogenous aka equal across samples?

head(Core.genOI.meta.all)

## Massilia first!!!!
# Kruskal-Wallis test is an ANOVA for non-normal data
massilia.fit1<-kruskal.test(Massilia ~ Site, data=Core.genOI.meta.all)

massilia.fit1

# Instead of using Tukey test, we can use Dunn's test to see which groups significantly vary if Kruskal-Wallis test is significant
# ANOVA + Tukey for normally distributed data, Kruskal-Wallis + Dunn's test for non-normal data
rstatix::dunn_test(Core.genOI.meta.all, Massilia ~ Site, p.adjust.method = "bonferroni", detailed = TRUE)

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(Massilia ~ Site, data = Core.genOI.meta.all)
# Fligner-Killeen:med chi-squared =0.40846, df = 3, p-value = 0.9385
# Which shows that the data do NOT deviate significantly from homogeneity.

compare_means(Massilia ~ Site, data=Core.genOI.meta.all, method="wilcox.test",p.adjust.method = "bonferroni")

compare_means(Massilia ~ Site, data=Core.genOI.meta.all, method="kruskal.test",p.adjust.method = "bonferroni")

# what about Massilia & Sate + Year
massilia.fit2<-kruskal.test(Massilia ~ SiteColYr, data=Core.genOI.meta.all)

massilia.fit2

# Instead of using Tukey test, we can use Dunn's test to see which groups significantly vary if Kruskal-Wallis test is significant
# ANOVA + Tukey for normally distributed data, Kruskal-Wallis + Dunn's test for non-normal data
massilia.fit2.grp<-rstatix::dunn_test(Core.genOI.meta.all, Massilia ~ SiteColYr, p.adjust.method = "bonferroni", detailed = TRUE)

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(Massilia ~ SiteColYr, data = Core.genOI.meta.all)
# Fligner-Killeen:med chi-squared = 4.7864, df = 7, p-value = 0.686
# Which shows that the data do NOT deviate significantly from homogeneity.

compare_means(Massilia ~ SiteColYr, data=Core.genOI.meta.all, method="wilcox.test",p.adjust.method = "bonferroni")

compare_means(Massilia ~ SiteColYr, data=Core.genOI.meta.all, method="kruskal.test",p.adjust.method = "bonferroni")

## Next Sphingomonas!
# Kruskal-Wallis test is an ANOVA for non-normal data
sphingomonas.fit1<-kruskal.test(Sphingomonas ~ Site, data=Core.genOI.meta.all)

sphingomonas.fit1

# Instead of using Tukey test, we can use Dunn's test to see which groups significantly vary if Kruskal-Wallis test is significant
# ANOVA + Tukey for normally distributed data, Kruskal-Wallis + Dunn's test for non-normal data
rstatix::dunn_test(Core.genOI.meta.all, Sphingomonas ~ Site, p.adjust.method = "bonferroni", detailed = TRUE)

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(Sphingomonas ~ Site, data = Core.genOI.meta.all)
# Fligner-Killeen:med chi-squared = 3.9001, df = 3, p-value = 0.2724
# Which shows that the data do NOT deviate significantly from homogeneity.

compare_means(Sphingomonas ~ Site, data=Core.genOI.meta.all, method="wilcox.test",p.adjust.method = "bonferroni")

compare_means(Sphingomonas ~ Site, data=Core.genOI.meta.all, method="kruskal.test",p.adjust.method = "bonferroni")

# what about Sphingomonas & Sate + Year
sphingomonas.fit2<-kruskal.test(Sphingomonas ~ SiteColYr, data=Core.genOI.meta.all)

sphingomonas.fit2

# Instead of using Tukey test, we can use Dunn's test to see which groups significantly vary if Kruskal-Wallis test is significant
# ANOVA + Tukey for normally distributed data, Kruskal-Wallis + Dunn's test for non-normal data
sphingomonas.fit2.grp<-rstatix::dunn_test(Core.genOI.meta.all, Sphingomonas ~ SiteColYr, p.adjust.method = "bonferroni", detailed = TRUE)

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(Sphingomonas ~ SiteColYr, data = Core.genOI.meta.all)
# Fligner-Killeen:med chi-squared = 2.6291, df = 7, p-value = 0.9171
# Which shows that the data do NOT deviate significantly from homogeneity.

compare_means(Sphingomonas ~ SiteColYr, data=Core.genOI.meta.all, method="wilcox.test",p.adjust.method = "bonferroni")

compare_means(Sphingomonas ~ SiteColYr, data=Core.genOI.meta.all, method="kruskal.test",p.adjust.method = "bonferroni")




#### Separate All Data by Time points ####
# create metadata df that will contain scaled chemical data
head(meta.all.scaled)

site_list<-unique(meta.all.scaled$Site) #define an array of string values
# go through metadata & create a list of data frames
## when metadata$Variable == element in site_list (aka x in this case), subset metadata by said element into elements of a list

# here the function(x) is using site_list aka x to subset metadata, when $Variable column == site_list
# Run the function so it's stored in Global Env
site_subsets<-lapply(site_list, function(x) {subset(meta.all.scaled, Site==x)})

site_subsets # sanity check1 (should see all elements in list)
site_subsets[[1]] # sanity check2 (see 1st element in list)
#rename the list elements

# name each element in list
names(site_subsets)<-site_list # * only do this if the order of names in site_list match order of the elements in site_subsets!
site_subsets$BDC # sanity check3 - should be able to pull dataframes by names rather than index now

# example of subsetting
site_subsets[[2]][1:3]
site_subsets$BDC[1:3] # should produce same ouptut as line above

site_subsets[[2]][1:2,1:2] # another example

# ^ subsetting to [[second dataframe]], [[row #, column #]]
site_subsets[[2]][[1,2]] # [[second dataframe]], [[row 1, column 2]]

# set up the function and run this to store it in our Global environment
df_specific.subset<-function(var_vec,var_subsets){
  # var_vec = vector of variable elements from specific categorical variable;
  ## e.g. vector of names from Site categorical variable (metadata sites)
  # var_subsets = list of dataframes subsetted by column$element from original dataframe;
  ## e.g. list of dataframes (each df = element of list) subsetted from metadata using vector of metadata$Site names
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
df_specific.subset(site_list, site_subsets) # used scaled metadata quantitative values

head(WI) # sanity check
WI[1:5,] # double check that our new Variable (here Site) data frames still have scaled chemical data
rownames(WI)

# matching data with user defined function -- here is the function, must run to store function in Global env
# note: metadata must have rownames for this to work!
match_dat<-function(compdata, subset_metadata){
  ### * comp data and metadata need to have row names - rownames should be Sample IDs
  subset_comp_data = pullrow<-(is.element(row.names(compdata), row.names(subset_metadata)))
  subset_comp_data=compdata[pullrow,]
  return(subset_comp_data)
}

rownames(b.genus_RelAb)
# loop through list containing each site's metadata and use match_dat to pair with CLR data
for (i in seq_along(site_subsets)){
  print(site_subsets[[i]]) # shows what is in each element within list
  #print(names(site_subsets[i]))
  new.genus.df<-match_dat(b.genus_RelAb,site_subsets[[i]])
  assign(paste0("b.gen.RA_",names(site_subsets[i])), new.genus.df,envir = .GlobalEnv)
}
# names(site_subsets[i]) --> gives us name of each element in list

# did the function work the way we wanted it to?

b.gen.RA_WI[1:4,1:4]
rownames(WI) %in% rownames(b.gen.RA_WI) # hopefully all of the rownames match, aka will get output of TRUE


#### Group Differences in Genera Relative Abundance ####
head(genOI.meta.all)

# Kruskal-Wallis test = nonparametric one-way ANOVA
kruskal.test(Bacillus ~ Site, data = genOI.meta.all)
pairwise.wilcox.test(genOI.meta.all$Bacillus, genOI.meta.all$Site, p.adjust.method = "bonf") # returns p values

# Kruskal-Wallis test = nonparametric one-way ANOVA
kruskal.test(Massilia ~ Site, data = genOI.meta.all)
pairwise.wilcox.test(genOI.meta.all$Massilia, genOI.meta.all$Site, p.adjust.method = "bonf") # returns p values

# Kruskal-Wallis test = nonparametric one-way ANOVA
kruskal.test(Massilia ~ SampDate, data = genOI.meta.all)
pairwise.wilcox.test(genOI.meta.all$Massilia, genOI.meta.all$SampDate, p.adjust.method = "bonf") # returns p values

#### Group Differences in Family Relative Abundance ####
head(famOI.meta.all)

# Kruskal-Wallis test = nonparametric one-way ANOVA
kruskal.test(Bacillaceae ~ CollectionYear, data = famOI.meta.all)
pairwise.wilcox.test(famOI.meta.all$Bacillaceae, famOI.meta.all$CollectionYear, p.adjust.method = "bonf") # returns p values

# Kruskal-Wallis test = nonparametric one-way ANOVA
kruskal.test(Oxalobacteraceae ~ Site, data = famOI.meta.all)
pairwise.wilcox.test(famOI.meta.all$Oxalobacteraceae, famOI.meta.all$Site, p.adjust.method = "bonf") # returns p values

# Kruskal-Wallis test = nonparametric one-way ANOVA
kruskal.test(Massilia ~ SampDate, data = famOI.meta.all)
pairwise.wilcox.test(famOI.meta.all$Massilia, famOI.meta.all$SampDate, p.adjust.method = "bonf") # returns p values
