#### Set WD & Load Libraries ####
getwd() # use setwd("path/to/files") if you are not in the right directory
setwd("/Volumes/HLF_SSD/Aronson_Lab_Data/Salton_Sea/SaltonSeaDust")
suppressPackageStartupMessages({ # load packages quietly
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
  #library(nationalparkcolors)
  library(shades)
  #library(ALDEx2)
  library(rstatix)
  library(devtools)
  #library(decontam)
  library(ggvegan)
  library(microbiome)
})

# NOTE: Aitchison Distance = Euclidean distance of CLR-transformed data
## CLR = center log ratio transformation

#### Load Global Env to Import Count/ASV Tables ####
load("data/Amplicon/SSDust_16S.V3V4_W23_Data_Ready.Rdata") # save global env to Rdata file
#load("data/Amplicon/SSD_16S_CLR_EucDist_Ready.Rdata")

head(b.dust.all)
bac.ASV_table[1:4,1:4]
bac.ASV_table[(nrow(bac.ASV_table)-4):(nrow(bac.ASV_table)),(ncol(bac.ASV_table)-4):(ncol(bac.ASV_table))] # last 4 rows & cols
head(dust_meta)
head(meta.all.scaled)

#### Create Centered Log-Ratio Table from ASV table ####
bac.ASV_table[1:4,1:4]
b.clr<-decostand(bac.ASV_table[,-1],method = "clr", pseudocount = 1) #CLR transformation
b.clr[1:4,1:4]
# df must have rownames are SampleIDs, columns are ASV IDs for vegan functions below

rownames(b.clr) %in% rownames(meta.all.scaled) # check if metadata and micro comp data have matching rownames

#### Look at Surface Type Frequencies by Sample for Modeling Later ####

# melt down data so that we can make a stacked barplot
STF.melt<-melt(SurfTypFreq[,c(1,3:12,15:16)],by=c("Site","SampleID","STF_Date"))
colnames(STF.melt)[which(names(STF.melt) == "variable")] <- "SurfaceType"
colnames(STF.melt)[which(names(STF.melt) == "value")] <- "Frequency"
head(STF.melt)

# create palette
colorset9 = as.data.frame(t(data.frame("BarrenLand"="peachpuff2","CropLand"="gold1","Developed"="gray","Forest"="darkgreen","Herbaceous"="limegreen",
                   "Mexico"="red1","OpenWater"="mediumblue","Others"="black","SaltonSea"="darkturquoise","Shrub"="saddlebrown")))
colorset9$SurfaceType<-rownames(colorset9)
colnames(colorset9)[which(names(colorset9) == "V1")] <- "ST_Color"
colorset9

# merge color palette and STF data together
STF.melt<-merge(STF.melt, colorset9, by="SurfaceType")
head(STF.melt)
STF.melt$SampDate_Color <- as.character(STF.melt$ST_Color)
head(SurfTypFreq)

# create factors for organizing data in plot
STF.melt$Site<-factor(STF.melt$Site,levels=c("PD","BDC","DP","WI"))
unique(STF.melt$STF_Date)
STF.melt$SampleID = factor(STF.melt$SampleID, levels=unique(STF.melt$SampleID[order(STF.melt$Site,STF.melt$STF_Date)]), ordered=TRUE)

# # plot time
# just.stfs1<-ggplot(STF.melt, aes(x=SampleID, y=Frequency, fill=SurfaceType))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
#   labs(title = "Surface Type Frequencies in Salton Sea Dust", x="SampleID", y="Frequency", subtitle="",fill="Surface Type")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
#   guides(fill=guide_legend(ncol=1))+
#   scale_fill_manual(name ="Surface Type",values=unique(STF.melt$ST_Color[order(STF.melt$SurfaceType)]),labels=c("Barren Land","Crop Land","Developed","Forest","Herbaceous","Mexico","Open Water","Others","Salton Sea","Shrub"))
# ggsave(just.stfs1,filename = "figures/SurfaceTypeFrequencies/SSD_SurfaceTypeFrequencys_barplot.png", width=10, height=10, dpi=600,create.dir = TRUE)
#
# just.stfs2<-ggplot(STF.melt, aes(x=SampleID, y=Frequency, fill=SurfaceType))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
#   labs(title = "Surface Type Frequencies in Salton Sea Dust", x="SampleID", y="Frequency", subtitle="",fill="Surface Type")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
#   guides(fill=guide_legend(ncol=1))+
#   scale_fill_manual(name ="Surface Type",values=unique(STF.melt$ST_Color[order(STF.melt$SurfaceType)]),labels=c("Barren Land","Crop Land","Developed","Forest","Herbaceous","Mexico","Open Water","Others","Salton Sea","Shrub")) +
#   facet_wrap(vars(Site), scales = "free")
# ggsave(just.stfs2,filename = "figures/SurfaceTypeFrequencies/SSD_SurfaceTypeFrequencys_bySite_barplot.png", width=10, height=10, dpi=600,create.dir = TRUE)

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

# loop through list containing each site's metadata and use match_dat to pair with CLR data
for (i in seq_along(site_subsets)){
  print(site_subsets[[i]]) # shows what is in each element within list
  #print(names(site_subsets[i]))
  new.clr.df<-match_dat(b.clr,site_subsets[[i]])
  assign(paste0("b.clr_",names(site_subsets[i])), new.clr.df,envir = .GlobalEnv)
}
# names(site_subsets[i]) --> gives us name of each element in list

# did the function work the way we wanted it to?

b.clr_WI[1:4,1:4]
rownames(WI) %in% rownames(b.clr_WI) # hopefully all of the rownames match, aka will get output of TRUE

#### Check Count Data Relationship w/ Env Variables (w/ DCA) ####
## remember, CCA assumes that our species have a unimodal relationship with our variables.
### unimodal = one maximum, think upsidedown bellcurve or something
## RDA assumes a linear relationship
## check the assumption w/ DCA
# ^ more on DCA here: https://ordination.okstate.edu/DCA.htm

# ALL data
# add pseudocount so row sums are > 0
b.clr.pseudo<-b.clr+1
b.dca = decorana(b.clr.pseudo)

#plot(b.dca) # may take too long to load, do not run unless you have to
b.dca #DCA1 axis length = 0.1192425; use RDA
## The length of first DCA axis:
## > 4 indicates heterogeneous dataset on which unimodal methods should be used (CCA),
##  < 3 indicates homogeneous dataset for which linear methods are suitable (RDA)
## between 3 and 4 both linear and unimodal methods are OK.

# BY MONTH

b.clr_WI.pseudo<-b.clr_WI+1
b.WI.dca = decorana(b.clr_WI.pseudo)
b.WI.dca #DCA1 axis length = 0.100551; use RDA

b.clr_DP.pseudo<-b.clr_DP+1
b.DP.dca = decorana(b.clr_DP.pseudo)
b.DP.dca #DCA1 axis length = 0.0917782; use RDA

b.clr_BDC.pseudo<-b.clr_BDC+1
b.BDC.dca = decorana(b.clr_BDC.pseudo)
b.BDC.dca #DCA1 axis length = 0.108458; use RDA

b.clr_PD.pseudo<-b.clr_PD+1
b.PD.dca = decorana(b.clr_PD.pseudo)
b.PD.dca #DCA1 axis length = 0.143368; use RDA

#### RDA w/ All Data ####

rownames(meta.all.scaled) %in% rownames(b.clr) # check order of DFs
head(meta.all.scaled)

rda.all.0<-rda(b.clr ~ precip_24hr_accum+ave.air_temp+ave.wind_speed+ave.relative_humidity+ave.wind_direction+BarrenLand+CropLand+Developed+Forest+Herbaceous+Mexico+OpenWater+Others+SaltonSea+Shrub,data=meta.all.scaled)

# check summary of RDA
rda.all.0
summary(rda.all.0)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.all.0) # 0.05642568
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.all.0, permutations = how(nperm=999)) # p = 0.348

## we can also do a permutation test by RDA axis
#anova(rda.all.0, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.all.0, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
# ave. wind speed, ave wind direction are sig, precip is near sig

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.all.0)
# precip_24hr_accum          ave.air_temp        ave.wind_speed ave.relative_humidity    ave.wind_direction            BarrenLand              CropLand
# 3.866125             18.532172              5.853327              6.187007              9.166089             32.556199             48.177209
# Developed                Forest            Herbaceous                Mexico             OpenWater                Others             SaltonSea
# 22.406496             30.509102              9.318751             24.567900             11.132523             30.331429             24.287988
# Shrub
# NA

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/

head(meta.all.scaled[,c(4,7:8,11:12,37:46)])
## we can use model selection instead of picking variables we think are important (by p values)
# more info on ordistep & ordiR2step here: https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
rda.all.a = ordistep(rda(b.clr ~ 1, data = meta.all.scaled[,c(4,7:8,11:12,37:46)]),
                     scope=formula(rda.all.0),
                     direction = "forward",
                     permutations = how(nperm=999))
# b.clr ~ Developed, precip 24 hr accum are sig, ave wind speed near sig

rda.all.a$anova # see significance of individual terms in model
#                               Df    AIC      F Pr(>F)
# + Developed          1 162.22 2.8670  0.003 **
# + precip_24hr_accum  1 162.00 2.0647  0.030 *

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.all.a2 = ordiR2step(rda(b.clr ~ 1, data = meta.all.scaled[,c(4,7:8,11:12,37:46)]),
                        scope=formula(rda.all.0),
                        permutations = how(nperm=999))
rda.all.a2$anova # see significance of individual terms in model
#                               R2.adj Df    AIC      F Pr(>F)
# nothing sig

# Dropping Herbaceous, Others, Forest, ave air temp, CropLand, and OpenWater due to high p values
rda.all.1<-rda(b.clr ~ precip_24hr_accum+ave.wind_speed+ave.relative_humidity+ave.wind_direction+BarrenLand+Developed+Mexico+SaltonSea+Shrub,data=meta.all.scaled)

# check summary of RDA
rda.all.1
summary(rda.all.1)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.all.1) # 0.08925873
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.all.1, permutations = how(nperm=999)) # p = 0.036

## we can also do a permutation test by RDA axis
#anova(rda.all.1, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.all.1, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
# ave. wind speed, ave wind direction are sig, precip is  sig

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.all.1)
# precip_24hr_accum        ave.wind_speed ave.relative_humidity    ave.wind_direction            BarrenLand             Developed                Mexico
# 1.988769              2.348579              3.245290              5.401027             59.266299             52.294921             10.292222
# SaltonSea                 Shrub
# 27.119313            101.123915

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/

head(meta.all.scaled[,c(4,7:8,11:12,37:46)])
## we can use model selection instead of picking variables we think are important (by p values)
# more info on ordistep & ordiR2step here: https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
rda.all.b = ordistep(rda(b.clr ~ 1, data = meta.all.scaled[,c(4,7:8,11:12,37:46)]),
                     scope=formula(rda.all.1),
                     direction = "forward",
                     permutations = how(nperm=999))
# b.clr ~ Developed, precip 24 hr accum are sig, ave wind speed near sig

rda.all.b$anova # see significance of individual terms in model
#                               Df    AIC      F Pr(>F)
# + Developed          1 162.22 2.8670  0.003 **
# + precip_24hr_accum  1 162.00 2.0647  0.030 *

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.all.b2 = ordiR2step(rda(b.clr ~ 1, data = meta.all.scaled[,c(4,7:8,11:12,37:46)]),
                        scope=formula(rda.all.1),
                        permutations = how(nperm=999))
rda.all.b2$anova # see significance of individual terms in model
#                               R2.adj Df    AIC      F Pr(>F)
# nothing sig

# dropping BarrenLand, SaltonSea, ave wind dir, ave rel humidity, Mexico, and Shrub due to high p vals and low R^2 vals
rda.all.2<-rda(b.clr ~ precip_24hr_accum+ave.wind_speed+Developed,data=meta.all.scaled)

# check summary of RDA
rda.all.2
summary(rda.all.2)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.all.2) # 0.1210451
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.all.2, permutations = how(nperm=999)) # p = 0.001

## we can also do a permutation test by RDA axis
#anova(rda.all.2, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.all.2, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
# Dveloped, ave. wind speed, precip are sig

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.all.2)
#precip_24hr_accum    ave.wind_speed         Developed
#1.014969          1.308912          1.292974

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/

head(meta.all.scaled[,c(4,7:8,11:12,37:46)])
## we can use model selection instead of picking variables we think are important (by p values)
# more info on ordistep & ordiR2step here: https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
rda.all.c = ordistep(rda(b.clr ~ 1, data = meta.all.scaled[,c(4,7:8,11:12,37:46)]),
                     scope=formula(rda.all.2),
                     direction = "forward",
                     permutations = how(nperm=999))
# b.clr ~ Developed, precip 24 hr accum are sig, ave wind speed near sig

rda.all.c$anova # see significance of individual terms in model
#                               Df    AIC      F Pr(>F)
# + Developed          1 162.22 2.8670  0.003 **
# + precip_24hr_accum  1 162.00 2.0647  0.030 *

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.all.c2 = ordiR2step(rda(b.clr ~ 1, data = meta.all.scaled[,c(4,7:8,11:12,37:46)]),
                        scope=formula(rda.all.2),
                        permutations = how(nperm=999))
rda.all.c2$anova # see significance of individual terms in model
#                               R2.adj Df    AIC      F Pr(>F)
# + Developed         0.064677  1 162.22 2.8670  0.002 **
# + precip_24hr_accum 0.101471  1 162.00 2.0647  0.034 *
#   <All variables>     0.118262



#### RDA - WI ####

rownames(WI) %in% rownames(b.clr_WI) # check order of DFs
head(WI)

# included all Surface type frequencies and it overfit the model so pulling some out...
rda.WI.0<-rda(b.clr_WI ~ precip_24hr_accum+ave.air_temp+ave.wind_speed+ave.relative_humidity+ave.wind_direction+BarrenLand+CropLand+Developed+Forest+Herbaceous+Mexico+OpenWater+Others+SaltonSea+Shrub,data=WI)

# check summary of RDA
rda.WI.0
summary(rda.WI.0)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.WI.0) # NA
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.WI.0, permutations = how(nperm=999))

## we can also do a permutation test by RDA axis
#anova(rda.WI.0, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.WI.0, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                               Df Variance      F Pr(>F)

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.WI.0)
# precip_24hr_accum          ave.air_temp        ave.wind_speed ave.relative_humidity    ave.wind_direction            BarrenLand              CropLand
# 19.07426              18.76960              44.33004              22.20618              23.12481              10.09772                    NA
# Developed                Forest            Herbaceous                Mexico             OpenWater                Others             SaltonSea
# NA                    NA                    NA                    NA                    NA                    NA                    NA
# Shrub
# NA

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(WI[,c(4,7:8,11:12,37:46)])
## we can use model selection instead of picking variables we think are important (by p values)
rda.WI.a = ordistep(rda(b.clr_WI ~ 1, data = WI[,c(4,7:8,11:12,37:46)]),
                         scope=formula(rda.WI.0),
                         direction = "forward",
                         permutations = how(nperm=999))
rda.WI.a$anova # see significance of individual terms in model
# + ave.relative_humidity  1 39.758 2.0824  0.032 *

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.WI.a2 = ordiR2step(rda(b.clr_WI ~ 1, data = WI[,c(4,7:8,11:12,37:46)]),
                            scope=formula(rda.WI.0),
                            permutations = how(nperm=999))
# too many terms


# dropped SaltonSea, Developed, Herbaceous, Mexico, Shrub, BarrenLand, precip 24 hr accum, CropLand, Forest, Others
rda.WI.1<-rda(b.clr_WI ~ ave.air_temp+ave.wind_speed+ave.relative_humidity+ave.wind_direction+OpenWater,data=WI)

# check summary of RDA
rda.WI.1
summary(rda.WI.1)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.WI.1) # -0.2450721
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.WI.1, permutations = how(nperm=999))

## we can also do a permutation test by RDA axis
#anova(rda.WI.1, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.WI.1, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                               Df Variance      F Pr(>F)

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.WI.1)
#

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(WI[,c(4,7:8,11:12,37:46)])
## we can use model selection instead of picking variables we think are important (by p values)
rda.WI.b = ordistep(rda(b.clr_WI ~ 1, data = WI[,c(4,7:8,11:12,37:46)]),
                    scope=formula(rda.WI.1),
                    direction = "forward",
                    permutations = how(nperm=999))
rda.WI.b$anova # see significance of individual terms in model
# + ave.relative_humidity  1 39.758 2.0824  0.035 *

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.WI.b2 = ordiR2step(rda(b.clr_WI ~ 1, data = WI[,c(4,7:8,11:12,37:46)]),
                       scope=formula(rda.WI.1),
                       permutations = how(nperm=999))
# too many terms

# remove ave air temp and OpenWater and ave wind speed
rda.WI.2<-rda(b.clr_WI ~ ave.relative_humidity+ave.wind_direction,data=WI)

# check summary of RDA
rda.WI.2
summary(rda.WI.2)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.WI.2) # 0.2449176
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.WI.2, permutations = how(nperm=999))

## we can also do a permutation test by RDA axis
#anova(rda.WI.2, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.WI.2, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                               Df Variance      F Pr(>F)

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.WI.2)
#

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(WI[,c(4,7:8,11:12,37:46)])
## we can use model selection instead of picking variables we think are important (by p values)
rda.WI.c = ordistep(rda(b.clr_WI ~ 1, data = WI[,c(4,7:8,11:12,37:46)]),
                    scope=formula(rda.WI.2),
                    direction = "forward",
                    permutations = how(nperm=999))
rda.WI.c$anova # see significance of individual terms in model
# + ave.relative_humidity  1 39.758 2.0824  0.035 *

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.WI.c2 = ordiR2step(rda(b.clr_WI ~ 1, data = WI[,c(4,7:8,11:12,37:46)]),
                       scope=formula(rda.WI.2),
                       permutations = how(nperm=999))
#

rda.WI.3<-rda(b.clr_WI ~ ave.relative_humidity,data=WI)

# check summary of RDA
rda.WI.3
summary(rda.WI.3)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.WI.3) # 0.1528332
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.WI.3, permutations = how(nperm=999))

## we can also do a permutation test by RDA axis
#anova(rda.WI.3, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.WI.3, by = "terms", permutations = how(nperm=999)) ### by variables

#### RDA - DP ####

rownames(DP) %in% rownames(b.clr_DP) # check order of DFs
head(DP)

# included all Surface type frequencies and it overfit the model so pulling some out...
rda.DP.0<-rda(b.clr_DP ~ precip_24hr_accum+ave.air_temp+ave.wind_speed+ave.relative_humidity+ave.wind_direction+BarrenLand+CropLand+Developed+Forest+Herbaceous+Mexico+OpenWater+Others+SaltonSea+Shrub,data=DP)

# check summary of RDA
rda.DP.0
summary(rda.DP.0)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.DP.0) # NA
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.DP.0, permutations = how(nperm=999))

## we can also do a permutation test by RDA axis
#anova(rda.DP.0, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.DP.0, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                               Df Variance      F Pr(>F)

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.DP.0)
# ave.air_temp        ave.wind_speed ave.relative_humidity    ave.wind_direction            BarrenLand              CropLand
# 18.91092              35.48104              24.87321              31.29241              11.32342              21.88877
# Developed                Forest            Herbaceous                Mexico             OpenWater                Others
# NA                    NA                    NA                    NA                    NA                    NA
# SaltonSea                 Shrub
# NA                    NA

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(DP[,c(4,7:8,11:12,37:46)])
## we can use model selection instead of picking variables we think are important (by p values)
rda.DP.a = ordistep(rda(b.clr_DP ~ 1, data = DP[,c(4,7:8,11:12,37:46)]),
                    scope=formula(rda.DP.0),
                    direction = "forward",
                    permutations = how(nperm=999))
rda.DP.a$anova # see significance of individual terms in model
# + Mexico  1 38.832 1.7083  0.022 *

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.DP.a2 = ordiR2step(rda(b.clr_DP ~ 1, data = DP[,c(4,7:8,11:12,37:46)]),
                       scope=formula(rda.DP.0),
                       permutations = how(nperm=999))
# too many terms

# check best fit model based on above results
anova(rda.DP.a, permutations = how(nperm=999)) #not significant


# drop Others, CropLand, OpenWater, SaltonSea, Forest, Developed, ave.wind speed, Shrub, ave wind direction due to high p value
rda.DP.1<-rda(b.clr_DP ~ precip_24hr_accum+ave.air_temp+ave.relative_humidity+BarrenLand+Herbaceous+Mexico,data=DP)

# check summary of RDA
rda.DP.1
summary(rda.DP.1)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.DP.1) # NA
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.DP.1, permutations = how(nperm=999))

## we can also do a permutation test by RDA axis
#anova(rda.DP.1, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.DP.1, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                               Df Variance      F Pr(>F)

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.DP.1)
#

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(DP[,c(4,7:8,11:12,37:46)])
## we can use model selection instead of picking variables we think are important (by p values)
rda.DP.b = ordistep(rda(b.clr_DP ~ 1, data = DP[,c(4,7:8,11:12,37:46)]),
                    scope=formula(rda.DP.1),
                    direction = "forward",
                    permutations = how(nperm=999))
rda.DP.b$anova # see significance of individual terms in model
# + Mexico  1 38.832 1.7083  0.022 *

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.DP.b2 = ordiR2step(rda(b.clr_DP ~ 1, data = DP[,c(4,7:8,11:12,37:46)]),
                       scope=formula(rda.DP.1),
                       permutations = how(nperm=999))
# too many terms

# check best fit model based on above results
anova(rda.DP.b, permutations = how(nperm=999)) #not significant


#### RDA - BDC ####

rownames(BDC) %in% rownames(b.clr_BDC) # check order of DFs
head(BDC)

# included all Surface type frequencies and it overfit the model so pulling some out...
rda.BDC.0<-rda(b.clr_BDC ~ precip_24hr_accum+ave.air_temp+ave.wind_speed+ave.relative_humidity+ave.wind_direction+BarrenLand+CropLand+Developed+Forest+Herbaceous+Mexico+OpenWater+Others+SaltonSea+Shrub,data=BDC)

# check summary of RDA
rda.BDC.0
summary(rda.BDC.0)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.BDC.0) # NA
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.BDC.0, permutations = how(nperm=999))

## we can also do a permutation test by RDA axis
#anova(rda.BDC.0, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.BDC.0, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                               Df Variance      F Pr(>F)

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.BDC.0)
# ave.air_temp        ave.wind_speed ave.relative_humidity    ave.wind_direction            BarrenLand              CropLand
# 18.91092              35.48104              24.87321              31.29241              11.32342              21.88877
# Developed                Forest            Herbaceous                Mexico             OpenWater                Others
# NA                    NA                    NA                    NA                    NA                    NA
# SaltonSea                 Shrub
# NA                    NA

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(BDC[,c(4,7:8,11:12,37:46)])
## we can use model selection instead of picking variables we think are important (by p values)
rda.BDC.a = ordistep(rda(b.clr_BDC ~ 1, data = BDC[,c(4,7:8,11:12,37:46)]),
                    scope=formula(rda.BDC.0),
                    direction = "forward",
                    permutations = how(nperm=999))
rda.BDC.a$anova # see significance of individual terms in model
# + OpenWater  1 63.384 2.0387  0.048 *

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.BDC.a2 = ordiR2step(rda(b.clr_BDC ~ 1, data = BDC[,c(4,7:8,11:12,37:46)]),
                       scope=formula(rda.BDC.0),
                       permutations = how(nperm=999))
# too many terms

# check best fit model based on above results
anova(rda.BDC.a, permutations = how(nperm=999)) #not significant


#### RDA - PD ####

rownames(PD) %in% rownames(b.clr_PD) # check order of DFs
head(PD)

# included all Surface type frequencies and it overfit the model so pulling some out...
rda.PD.0<-rda(b.clr_PD ~ precip_24hr_accum+ave.air_temp+ave.wind_speed+ave.relative_humidity+ave.wind_direction+BarrenLand+CropLand+Developed+Forest+Herbaceous+Mexico+OpenWater+Others+SaltonSea+Shrub,data=PD)

# check summary of RDA
rda.PD.0
summary(rda.PD.0)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.PD.0) # NA
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.PD.0, permutations = how(nperm=999))

## we can also do a permutation test by RDA axis
#anova(rda.PD.0, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.PD.0, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                               Df Variance      F Pr(>F)

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.PD.0)
# ave.air_temp        ave.wind_speed ave.relative_humidity    ave.wind_direction            BarrenLand              CropLand
# 18.91092              35.48104              24.87321              31.29241              11.32342              21.88877
# Developed                Forest            Herbaceous                Mexico             OpenWater                Others
# NA                    NA                    NA                    NA                    NA                    NA
# SaltonSea                 Shrub
# NA                    NA

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(PD[,c(4,7:8,11:12,37:46)])
## we can use model selection instead of picking variables we think are important (by p values)
rda.PD.a = ordistep(rda(b.clr_PD ~ 1, data = PD[,c(4,7:8,11:12,37:46)]),
                    scope=formula(rda.PD.0),
                    direction = "forward",
                    permutations = how(nperm=999))
rda.PD.a$anova # see significance of individual terms in model
# + OpenWater  1 63.384 2.0387  0.048 *

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.PD.a2 = ordiR2step(rda(b.clr_PD ~ 1, data = PD[,c(4,7:8,11:12,37:46)]),
                       scope=formula(rda.PD.0),
                       permutations = how(nperm=999))
# too many terms

# check best fit model based on above results
anova(rda.PD.a, permutations = how(nperm=999)) #not significant


#### Final RDAs ####
# RDA by sampling timepoint
head(meta.all.scaled)
head(b.clr)
rownames(b.clr) %in% rownames(meta.all.scaled) # sanity check 1

# all data
#rda.all.2$call # best model for all data

rda.all<-rda(b.clr ~ precip_24hr_accum + ave.wind_speed + Developed,data=meta.all.scaled)
rda.all
summary(rda.all)
RsquareAdj(rda.all) # how much variation is explained by our model? 0.1210451
anova(rda.all, permutations = how(nperm=999)) # p-value = 0.001 **
anova(rda.all, by = "terms", permutations = how(nperm=999))
#                               Df Variance      F Pr(>F)
# precip_24hr_accum  1   22.500 2.1097  0.033 *
#   ave.wind_speed     1   20.914 1.9609  0.014 *
#   Developed          1   28.238 2.6477  0.001 ***
#   Residual          24  255.965

# p value adjusted for each term
aov.rda.all.terms<-anova(rda.all, by = "terms", permutations = how(nperm=999))
p.adjust(aov.rda.all.terms$`Pr(>F)`,method="bonferroni",n=length(aov.rda.all.terms)) # adjusted pvalues
# [1] 0.120 0.056 0.008  NA

# p value adjusted for entire model
aov.rda.all<-anova(rda.all, by = NULL, permutations = how(nperm=999))
p.adjust(aov.rda.all$`Pr(>F)`,method="bonferroni",n=length(aov.rda.all)) # adjusted pvalues
# [1] 0.004   NA

# WI
#rda.WI.2$call # best model

rda.WI<-rda(b.clr_WI ~ ave.relative_humidity,data=WI)
summary(rda.WI)
RsquareAdj(rda.WI) # how much variation is explained by our model? 0.2449176
anova(rda.WI, permutations = how(nperm=999)) # p-value = 0.036
anova(rda.WI, by = "terms", permutations = how(nperm=999))
#                           Df Variance      F Pr(>F)
#
aov.rda.WI<-anova(rda.WI, by = "terms", permutations = how(nperm=999))
p.adjust(aov.rda.WI$`Pr(>F)`,method="bonferroni",n=length(aov.rda.WI)) # adjusted pvalues
# [1] 0.224 0.944   NA

# DP
#rda.DP.9$call # best model from above

rda.DP<-rda(b.clr_DP ~ ave.wind_direction + BarrenLand,data=DP)
summary(rda.DP)
RsquareAdj(rda.DP) # how much variation is explained by our model? 0.2016672
anova(rda.DP, permutations = how(nperm=999)) # p-value = 0.039
anova(rda.DP, by = "terms", permutations = how(nperm=999))
#                 Df Variance      F Pr(>F)
# ave.wind_direction  1   3497.7 2.4809  0.013 *
#   BarrenLand          1   1458.9 1.0348  0.428
# Residual            4   5639.5

aov.rda.DP<-anova(rda.DP, by = NULL, permutations = how(nperm=999))
p.adjust(aov.rda.DP$`Pr(>F)`,method="bonferroni",n=length(aov.rda.DP)) # adjusted pvalues
# 0.2

aov.rda.DP.terms<-anova(rda.DP, by = "terms", permutations = how(nperm=999))
p.adjust(aov.rda.DP.terms$`Pr(>F)`,method="bonferroni",n=length(aov.rda.DP.terms)) # adjusted pvalues
# 0.044 1.000    NA

# BDC
#rda.BDC.9$call - best model

rda.BDC<-rda(b.clr_BDC ~ BarrenLand + Developed,data=BDC)
summary(rda.BDC)
RsquareAdj(rda.BDC) # how much variation is explained by our model? 0.1373921
anova(rda.BDC, permutations = how(nperm=999)) # p-value = 0.052
anova(rda.BDC, by = "terms", permutations = how(nperm=999))
#                           Df Variance      F Pr(>F)
# BarrenLand  1   1754.3 1.8070  0.047 *
#   Developed   1   1115.2 1.1487  0.241
# Residual    4   3883.5

aov.rda.BDC<-anova(rda.BDC, by = NULL, permutations = how(nperm=999))
p.adjust(aov.rda.BDC$`Pr(>F)`,method="bonferroni",n=length(aov.rda.BDC)) # adjusted pvalues
# [1] 0.208    NA
aov.rda.BDC.terms<-anova(rda.BDC, by = "terms", permutations = how(nperm=999))
p.adjust(aov.rda.BDC.terms$`Pr(>F)`,method="bonferroni",n=length(aov.rda.BDC.terms)) # adjusted pvalues
# 0.196 1.000    NA

# best model : rda.PD.9$call
rda.PD<-rda(b.clr_PD ~ ave.wind_speed+CropLand,data=PD)
summary(rda.PD)
RsquareAdj(rda.PD) # how much variation is explained by our model? 0.1046263
anova(rda.PD, permutations = how(nperm=999)) # p-value = 0.156
anova(rda.PD, by = "terms", permutations = how(nperm=999))
#                           Df Variance      F Pr(>F)
# ave.wind_speed  1   2405.9 1.3473  0.189
# CropLand        1   2417.5 1.3538  0.247
# Residual        4   7142.8

aov.rda.PD<-anova(rda.PD, by = NULL, permutations = how(nperm=999))
p.adjust(aov.rda.PD$`Pr(>F)`,method="bonferroni",n=length(aov.rda.PD)) # adjusted pvalues
#
aov.rda.PD.terms<-anova(rda.PD, by = "terms", permutations = how(nperm=999))
p.adjust(aov.rda.PD.terms$`Pr(>F)`,method="bonferroni",n=length(aov.rda.PD.terms)) # adjusted pvalues
#

# save RDAs as R object
save.image("data/Amplicon/SSD_Amplicon_EnvDriver_RDAsOnly.Rdata")

#### Plot RDA - ALL data ####
#plot(rda.WI) # depending on how many species you have, this step may take a while
plot(rda.all, scaling = 1)
## scaling = 1 -> emphasizes relationships among sites
plot(rda.all, scaling = 2)
## scaling = 2 -> emphasizes relationships among species

# check summary of RDA
summary(rda.all)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.all) # 0.1210451 aka 12.10%
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
##anova(rda.all, permutations = how(nperm=999)) # p = 0.001, significant

png('figures/EnvDrivers/Aitchison/SSD_AllData_autoplot_rda_example.png',width = 700, height = 600, res=100)
autoplot(rda.all, arrows = TRUE,data = rda.all ,layers=c("biplot","sites"),label = FALSE, label.size = 3, shape = FALSE, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3, scale= 0)+theme_classic()
dev.off()
## FOR AUTOPLOT -> must load packagve ggvegan first

# variance partitioning of RDA
rda.all.part<-varpart(b.clr, meta.all.scaled$precip_24hr_accum, meta.all.scaled$Developed,meta.all.scaled$ave.wind_speed)
rda.all.part$part
# plot variance partitioning results
png('figures/EnvDrivers/Aitchison/SSD_AllData_RDA_VariancePartitioning.png',width = 1500, height = 1500, res=100)
plot(rda.all.part,
     Xnames = c("Accum. Precip (24 hr)", "Developed STF","Ave Wind Speed"), # name the partitions
     bg = c("steelblue3", "gray","tomato3"), alpha = 80, # colour the circles
     digits = 3, # only show 2 digits
     cex = 1.5)
dev.off()

rda.sum.all<-summary(rda.all)
rda.sum.all$sites[,1:2]
rda.sum.all$cont #cumulative proportion of variance per axis
# RDA1 = 9.97, RDA2 = 7.82

# create data frame w/ RDA axes for sites
# first check rownames of RDA & metadata, then make df
rownames(rda.sum.all$sites) %in% rownames(meta.all.scaled)
rda.axes.all1<-data.frame(RDA1=rda.sum.all$sites[,1], RDA2=rda.sum.all$sites[,2], SampleID=rownames(rda.sum.all$sites), Site=meta.all.scaled$Site, SampDate=meta.all.scaled$SampDate)

# then merge with metadata to get all category colors!
rda.axes.all<-merge(rda.axes.all1,meta.all.scaled,by=c("SampleID","Site","SampDate"))

# create data frame w/ RDA axes for variables
arrows.all<-data.frame(RDA1=rda.sum.all$biplot[,1], RDA2=rda.sum.all$biplot[,2], Label=rownames(rda.sum.all$biplot))
arrows.all
arrows.all$Label[(arrows.all$Label) == "precip_24hr_accum"] <- "Accum. Precip (24 hr)"
arrows.all$Label[(arrows.all$Label) == "ave.wind_speed"] <- "Ave Wind Speed"
arrows.all$Label[(arrows.all$Label) == "Developed"] <- "Developed STF"

rda.sum.all$cont #cumulative proportion of variance per axis
# RDA1 = 9.97, RDA2 = 7.82

rda.plot1<-ggplot(rda.axes.all, aes(x = RDA1, y = RDA2)) + geom_point(size=2) +
  geom_segment(data = arrows.all,mapping = aes(x = 0, y = 0, xend = RDA1, yend = RDA2),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.5,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.all,aes(label = Label, x = RDA1, y = RDA2, fontface="bold"))+
  coord_fixed() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))

rda.plot2<-ggplot(rda.axes.all, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=SampDate,shape=Site),size=4) +
  geom_segment(data = arrows.all,mapping = aes(x = 0, y = 0, xend = RDA1*5.5, yend = RDA2*5.5),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.8,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.all,aes(label = Label, x = RDA1*7, y = RDA2*7, fontface="bold"), size=4)+
  coord_fixed(ratio = 1, xlim = c(-8,8), ylim = c(-8,8)) + theme_classic() +
  scale_color_manual(name ="Collection Date",values=unique(rda.axes.all$SampDate_Color[order(rda.axes.all$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  scale_shape_manual(values = c(7,10, 15,16)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data") +
  xlab("RDA1 [9.97%]") + ylab("RDA2 [7.82%]")

ggsave(rda.plot2,filename = "figures/EnvDrivers/Aitchison/SSD_16S_RDA_AllData.png", width=10, height=10, dpi=600,create.dir = TRUE)


rda.plot3<-ggplot(rda.axes.all, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=SampDate,shape=Site),size=5) +
  geom_segment(data = arrows.all,mapping = aes(x = 0, y = 0, xend = RDA1*6, yend = RDA2*6),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 1,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.all,aes(label = Label, x = RDA1*8, y = RDA2*8, fontface="bold"), size=5)+
  coord_fixed(ratio = 1, xlim = c(-10,10), ylim = c(-7,7)) + theme_classic() + scale_color_manual(name ="Collection Date",values=unique(rda.axes.all$SampDate_Color[order(rda.axes.all$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  scale_shape_manual(values = c(7,10, 15,16)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data") +
  xlab("RDA1 [9.97%]") + ylab("RDA2 [7.82%]")

ggsave(rda.plot3,filename = "figures/EnvDrivers/Aitchison/SSD_16S_RDA_AllData_bigger.png", width=15, height=15, dpi=600,create.dir = TRUE)


#### Plot RDA - WI ####
#plot(rda.WI) # depending on how many species you have, this step may take a while
plot(rda.WI, scaling = 1)
## scaling = 1 -> emphasizes relationships among sites
plot(rda.WI, scaling = 2)
## scaling = 2 -> emphasizes relationships among species

# check summary of RDA
summary(rda.WI)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.WI) # 0.2449176
## ^^ use this b/c chance correlations can inflate R^2

png('figures/EnvDrivers/Aitchison/SSD_WI_autoplot_rda_example.png',width = 700, height = 600, res=100)
autoplot(rda.WI, arrows = TRUE,data = rda.WI ,layers=c("biplot","sites"),label = FALSE, label.size = 3, shape = FALSE, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3, scale= 0)+theme_classic()
dev.off()
## FOR AUTOPLOT -> must load packagve ggvegan first

# variance partitioning of RDA
rda.WI.part<-varpart(b.clr_WI, WI$ave.relative_humidity, WI$ave.wind_direction)
rda.WI.part$part
# plot variance partitioning results
png('figures/EnvDrivers/Aitchison/SSD_WI_RDA_VariancePartitioning.png',width = 900, height = 900, res=100)
plot(rda.WI.part,
     Xnames = c("Ave Rel. Humidity", "Ave Wind Dir"), # name the partitions
     bg = c("firebrick1", "purple1"), alpha = 80, # colour the circles
     digits = 3, # only show 3 digits
     cex = 1.5)
dev.off()

rda.sum.WI<-summary(rda.WI)
rda.sum.WI$sites[,1:2]
rda.sum.WI$cont #cumulative proportion of variance per axis
# RDA1=40.89%, PC1=8.77%

# create data frame w/ RDA axes for sites
rda.axes.WI<-data.frame(RDA1=rda.sum.WI$sites[,1], PC1=rda.sum.WI$sites[,2], SampleID=rownames(rda.sum.WI$sites), Site=WI$Site)

# then merge with metadata to get all category colors!
rda.axes.WI.all<-merge(rda.axes.WI,WI,by=c("SampleID","Site"))

# create data frame w/ RDA axes for variables
arrows.WI<-data.frame(RDA1=rda.sum.WI$biplot[,1], PC1=rda.sum.WI$biplot[,2], Label=rownames(rda.sum.WI$biplot))
arrows.WI
arrows.WI$Label[(arrows.WI$Label) == "ave.relative_humidity"] <- "Ave Rel Humidity"
#arrows.WI$Label[(arrows.WI$Label) == "ave.wind_direction"] <- "Ave. Wind Dir"

rda.plot5<-ggplot(rda.axes.WI.all, aes(x = RDA1, y = PC1)) + geom_point(size=2) +
  geom_segment(data = arrows.WI,mapping = aes(x = 0, y = 0, xend = RDA1, yend = PC1),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.5,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.WI,aes(label = Label, x = RDA1, y = PC1, fontface="bold"))+
  coord_fixed() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))

rda.plot6<-ggplot(rda.axes.WI.all, aes(x = RDA1, y = PC1)) + geom_point(aes(color=SampDate),size=4) +
  geom_segment(data = arrows.WI,mapping = aes(x = 0, y = 0, xend = RDA1*8, yend = PC1*8),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.8,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.WI,aes(label = Label, x = RDA1*9.85, y = PC1*9.5, fontface="bold"), size=4)+
  coord_fixed(ratio = 1,xlim = c(-10,10), ylim = c(-8,8)) + theme_classic() + scale_color_manual(name ="Collection Date",values=unique(rda.axes.WI.all$SampDate_Color[order(rda.axes.WI.all$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "September 2021", "December 2021")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Sea Dust, WI",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [40.89%]") + ylab("PC1 [8.77%]")

ggsave(rda.plot6,filename = "figures/EnvDrivers/Aitchison/SSD_16S_RDA_WI.png", width=16, height=12, dpi=600,create.dir = TRUE)

rda.plot6b<-ggplot(rda.axes.WI.all, aes(x = RDA1, y = PC1)) + geom_point(aes(color=SampDate),size=5) +
  geom_segment(data = arrows.WI,mapping = aes(x = 0, y = 0, xend = RDA1*8, yend = PC1*8),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 1,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.WI,aes(label = Label, x = RDA1*9, y = PC1*9.5, fontface="bold"), size=5)+
  coord_fixed(ratio = 1,xlim = c(-10,10), ylim = c(-8,8)) + theme_classic() + scale_color_manual(name ="Collection Date",values=unique(rda.axes.WI.all$SampDate_Color[order(rda.axes.WI.all$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "September 2021", "December 2021")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Sea Dust, WI",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [40.89%]") + ylab("PC1 [8.77%]")

ggsave(rda.plot6b,filename = "figures/EnvDrivers/Aitchison/SSD_16S_RDA_WI_bigger.png", width=15, height=15, dpi=600,create.dir = TRUE)

#### Plot RDA - DP ####
#plot(rda.DP) # depending on how many species you have, this step may take a while
plot(rda.DP, scaling = 1)
## scaling = 1 -> emphasizes relationships among sites
plot(rda.DP, scaling = 2)
## scaling = 2 -> emphasizes relationships among species

# check summary of RDA
summary(rda.DP)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.DP) # 0.2016672
## ^^ use this b/c chance correlations can inflate R^2

png('figures/EnvDrivers/Aitchison/SSD_DP_autoplot_rda_example.png',width = 700, height = 600, res=100)
autoplot(rda.DP, arrows = TRUE,data = rda.DP ,layers=c("biplot","sites"),label = FALSE, label.size = 3, shape = FALSE, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3, scale= 0)+theme_classic()
dev.off()
## FOR AUTOPLOT -> must load packagve ggvegan first

# variance partitioning of RDA
rda.DP.part<-varpart(b.clr_DP, DP$ave.wind_direction, DP$BarrenLand)
rda.DP.part$part
# plot variance partitioning results
png('figures/EnvDrivers/Aitchison/SSD_DP_RDA_VariancePartitioning.png',width = 900, height = 900, res=100)
plot(rda.DP.part,
     Xnames = c("Ave. Wind Dir.", "Barren Land STF"), # name the partitions
     bg = c("magenta", "peachpuff2"), alpha = 80, # colour the circles
     digits = 3, # only show 3 digits
     cex = 1.5)
dev.off()

rda.sum.DP<-summary(rda.DP)
rda.sum.DP$sites[,1:2]
rda.sum.DP$cont #cumulative proportion of variance per axis
# RDA1=33.25%, RDA2=13.53%

# create data frame w/ RDA axes for sites
rda.axes.DP<-data.frame(RDA1=rda.sum.DP$sites[,1], RDA2=rda.sum.DP$sites[,2], SampleID=rownames(rda.sum.DP$sites), Site=DP$Site)

# then merge with metadata to get all category colors!
rda.axes.DP.all<-merge(rda.axes.DP,DP,by=c("SampleID","Site"))

# create data frame w/ RDA axes for variables
arrows.DP<-data.frame(RDA1=rda.sum.DP$biplot[,1], RDA2=rda.sum.DP$biplot[,2], Label=rownames(rda.sum.DP$biplot))
arrows.DP$Label[(arrows.DP$Label) == "ave.wind_direction"] <- "Ave Wind Dir"
arrows.DP$Label[(arrows.DP$Label) == "BarrenLand"] <- "Barren Land STF"

rda.plot5<-ggplot(rda.axes.DP.all, aes(x = RDA1, y = RDA2)) + geom_point(size=2) +
  geom_segment(data = arrows.DP,mapping = aes(x = 0, y = 0, xend = RDA1, yend = RDA2),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.5,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.DP,aes(label = Label, x = RDA1, y = RDA2, fontface="bold"))+
  coord_fixed() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))

rda.plot6<-ggplot(rda.axes.DP.all, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=SampDate),size=4) +
  geom_segment(data = arrows.DP,mapping = aes(x = 0, y = 0, xend = RDA1*8, yend = RDA2*8),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.8,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.DP,aes(label = Label, x = RDA1*9.85, y = RDA2*9.5, fontface="bold"), size=4)+
  coord_fixed(ratio = 1) + theme_classic() + scale_color_manual(name ="Collection Date",values=unique(rda.axes.DP.all$SampDate_Color[order(rda.axes.DP.all$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "August 2021", "September 2021", "December 2021")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Sea Dust, DP",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [33.25%]") + ylab("RDA2 [13.53%]")

ggsave(rda.plot6,filename = "figures/EnvDrivers/Aitchison/SSD_16S_RDA_DP.png", width=16, height=12, dpi=600,create.dir = TRUE)

rda.plot6b<-ggplot(rda.axes.DP.all, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=SampDate),size=5) +
  geom_segment(data = arrows.DP,mapping = aes(x = 0, y = 0, xend = RDA1*8, yend = RDA2*8),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 1,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.DP,aes(label = Label, x = RDA1*9, y = RDA2*9.5, fontface="bold"), size=5)+
  coord_fixed(ratio = 1) + theme_classic() + scale_color_manual(name ="Collection Date",values=unique(rda.axes.DP.all$SampDate_Color[order(rda.axes.DP.all$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "August 2021", "September 2021", "December 2021")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Sea Dust, DP",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [33.25%]") + ylab("RDA2 [13.53%]")

ggsave(rda.plot6b,filename = "figures/EnvDrivers/Aitchison/SSD_16S_RDA_DP_bigger.png", width=15, height=15, dpi=600,create.dir = TRUE)

#### Plot RDA - BDC ####
#plot(rda.BDC) # depending on how many species you have, this step may take a while
plot(rda.BDC, scaling = 1)
## scaling = 1 -> emphasizes relationships among sites
plot(rda.BDC, scaling = 2)
## scaling = 2 -> emphasizes relationships among species

# check summary of RDA
summary(rda.BDC)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.BDC) # 0.1373921
## ^^ use this b/c chance correlations can inflate R^2

png('figures/EnvDrivers/Aitchison/SSD_BDC_autoplot_rda_example.png',width = 700, height = 600, res=100)
autoplot(rda.BDC, arrows = TRUE,data = rda.BDC ,layers=c("biplot","sites"),label = FALSE, label.size = 3, shape = FALSE, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3, scale= 0)+theme_classic()
dev.off()
## FOR AUTOPLOT -> must load packagve ggvegan first

# variance partitioning of RDA
rda.BDC.part<-varpart(b.clr_BDC, BDC$BarrenLand, BDC$Developed)
rda.BDC.part$part
# plot variance partitioning results
png('figures/EnvDrivers/Aitchison/SSD_BDC_RDA_VariancePartitioning.png',width = 900, height = 900, res=100)
plot(rda.BDC.part,
     Xnames = c("Ave. Air Temp", "Developed STF"), # name the partitions
     bg = c("peachpuff2", "lightgray"), alpha = 80, # colour the circles
     digits = 3, # only show 3 digits
     cex = 1.5)
dev.off()

rda.sum.BDC<-summary(rda.BDC)
rda.sum.BDC$sites[,1:2]
rda.sum.BDC$cont #cumulative proportion of variance per axis
# RDA1=31.19%, RDA2=11.31%

# create data frame w/ RDA axes for sites
rda.axes.BDC<-data.frame(RDA1=rda.sum.BDC$sites[,1], RDA2=rda.sum.BDC$sites[,2], SampleID=rownames(rda.sum.BDC$sites), Site=BDC$Site)

# then merge with metadata to get all category colors!
rda.axes.BDC.all<-merge(rda.axes.BDC,BDC,by=c("SampleID","Site"))

# create data frame w/ RDA axes for variables
arrows.BDC<-data.frame(RDA1=rda.sum.BDC$biplot[,1], RDA2=rda.sum.BDC$biplot[,2], Label=rownames(rda.sum.BDC$biplot))
arrows.BDC$Label[(arrows.BDC$Label) == "BarrenLand"] <- "Barren Land STF"
arrows.BDC$Label[(arrows.BDC$Label) == "Developed"] <- "Developed STF"

rda.plot5<-ggplot(rda.axes.BDC.all, aes(x = RDA1, y = RDA2)) + geom_point(size=2) +
  geom_segment(data = arrows.BDC,mapping = aes(x = 0, y = 0, xend = RDA1, yend = RDA2),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.5,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.BDC,aes(label = Label, x = RDA1, y = RDA2, fontface="bold"))+
  coord_fixed() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))

rda.plot6<-ggplot(rda.axes.BDC.all, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=SampDate),size=4) +
  geom_segment(data = arrows.BDC,mapping = aes(x = 0, y = 0, xend = RDA1*8, yend = RDA2*8),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.8,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.BDC,aes(label = Label, x = RDA1*9.85, y = RDA2*9.5, fontface="bold"), size=4)+
  coord_fixed(ratio = 1) + theme_classic() + scale_color_manual(name ="Collection Date",values=unique(rda.axes.BDC.all$SampDate_Color[order(rda.axes.BDC.all$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "September 2021", "December 2021")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Sea Dust, BDC",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [31.19%]") + ylab("RDA2 [11.31%]")

ggsave(rda.plot6,filename = "figures/EnvDrivers/Aitchison/SSD_16S_RDA_BDC.png", width=16, height=12, dpi=600,create.dir = TRUE)

rda.plot6b<-ggplot(rda.axes.BDC.all, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=SampDate),size=5) +
  geom_segment(data = arrows.BDC,mapping = aes(x = 0, y = 0, xend = RDA1*8, yend = RDA2*8),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 1,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.BDC,aes(label = Label, x = RDA1*9, y = RDA2*9.5, fontface="bold"), size=5)+
  coord_fixed(ratio = 1) + theme_classic() + scale_color_manual(name ="Collection Date",values=unique(rda.axes.BDC.all$SampDate_Color[order(rda.axes.BDC.all$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "September 2021", "December 2021")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Sea Dust, BDC",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [31.19%]") + ylab("RDA2 [11.31%]")

ggsave(rda.plot6b,filename = "figures/EnvDrivers/Aitchison/SSD_16S_RDA_BDC_bigger.png", width=15, height=15, dpi=600,create.dir = TRUE)

#### Plot RDA - PD ####
#plot(rda.PD) # depending on how many species you have, this step may take a while
plot(rda.PD, scaling = 1)
## scaling = 1 -> emphasizes relationships among sites
plot(rda.PD, scaling = 2)
## scaling = 2 -> emphasizes relationships among species

# check summary of RDA
summary(rda.PD)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.PD) # 0.1046263
## ^^ use this b/c chance correlations can inflate R^2

png('figures/EnvDrivers/Aitchison/SSD_PD_autoplot_rda_example.png',width = 700, height = 600, res=100)
autoplot(rda.PD, arrows = TRUE,data = rda.PD ,layers=c("biplot","sites"),label = FALSE, label.size = 3, shape = FALSE, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3, scale= 0)+theme_classic()
dev.off()
## FOR AUTOPLOT -> must load packagve ggvegan first

# variance partitioning of RDA
rda.PD.part<-varpart(b.clr_PD, PD$ave.wind_speed, PD$CropLand)
rda.PD.part$part
# plot variance partitioning results
png('figures/EnvDrivers/Aitchison/SSD_PD_RDA_VariancePartitioning.png',width = 900, height = 900, res=100)
plot(rda.PD.part,
     Xnames = c("Ave. Wind Speed", "Crop Land STF"), # name the partitions
     bg = c("blue1", "gold1"), alpha = 80, # colour the circles
     digits = 3, # only show 3 digits
     cex = 1.5)
dev.off()

rda.sum.PD<-summary(rda.PD)
rda.sum.PD$sites[,1:2]
rda.sum.PD$cont #cumulative proportion of variance per axis
# RDA1=20.89%, RDA2=19.42%

# create data frame w/ RDA axes for sites
rda.axes.PD<-data.frame(RDA1=rda.sum.PD$sites[,1], RDA2=rda.sum.PD$sites[,2], SampleID=rownames(rda.sum.PD$sites), Site=PD$Site)

# then merge with metadata to get all category colors!
rda.axes.PD.all<-merge(rda.axes.PD,PD,by=c("SampleID","Site"))

# create data frame w/ RDA axes for variables
arrows.PD<-data.frame(RDA1=rda.sum.PD$biplot[,1], RDA2=rda.sum.PD$biplot[,2], Label=rownames(rda.sum.PD$biplot))
arrows.PD$Label[(arrows.PD$Label) == "ave.wind_speed"] <- "Ave Wind Speed"
arrows.PD$Label[(arrows.PD$Label) == "CropLand"] <- "Crop Land STF"

rda.plot5<-ggplot(rda.axes.PD.all, aes(x = RDA1, y = RDA2)) + geom_point(size=2) +
  geom_segment(data = arrows.PD,mapping = aes(x = 0, y = 0, xend = RDA1, yend = RDA2),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.5,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.PD,aes(label = Label, x = RDA1, y = RDA2, fontface="bold"))+
  coord_fixed() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))

rda.plot6<-ggplot(rda.axes.PD.all, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=SampDate),size=4) +
  geom_segment(data = arrows.PD,mapping = aes(x = 0, y = 0, xend = RDA1*8, yend = RDA2*8),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.8,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.PD,aes(label = Label, x = RDA1*9.85, y = RDA2*9.5, fontface="bold"), size=4)+
  coord_fixed(ratio = 1) + theme_classic() + scale_color_manual(name ="Collection Date",values=unique(rda.axes.PD.all$SampDate_Color[order(rda.axes.PD.all$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "September 2021", "December 2021")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Sea Dust, PD",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [20.89%]") + ylab("RDA2 [19.42%]")

ggsave(rda.plot6,filename = "figures/EnvDrivers/Aitchison/SSD_16S_RDA_PD.png", width=16, height=12, dpi=600,create.dir = TRUE)

rda.plot6b<-ggplot(rda.axes.PD.all, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=SampDate),size=5) +
  geom_segment(data = arrows.PD,mapping = aes(x = 0, y = 0, xend = RDA1*8, yend = RDA2*8),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 1,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.PD,aes(label = Label, x = RDA1*9, y = RDA2*9.5, fontface="bold"), size=5)+
  coord_fixed(ratio = 1) + theme_classic() + scale_color_manual(name ="Collection Date",values=unique(rda.axes.PD.all$SampDate_Color[order(rda.axes.PD.all$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "September 2021", "December 2021")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Sea Dust, PD",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [20.89%]") + ylab("RDA2 [19.42%]")

ggsave(rda.plot6b,filename = "figures/EnvDrivers/Aitchison/SSD_16S_RDA_PD_bigger.png", width=15, height=15, dpi=600,create.dir = TRUE)

#### Save Progress ####

save.image("data/SSD_Amplicon_EnvDriver.Rdata")
