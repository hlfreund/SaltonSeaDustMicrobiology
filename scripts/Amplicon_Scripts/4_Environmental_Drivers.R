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
  library(nationalparkcolors)
  library(shades)
  library(ALDEx2)
  library(rstatix)
  library(devtools)
  library(decontam)
  library(ggvegan)
  library(microbiome)
})

#### Load Global Env to Import Count/ASV Tables ####
load("data/Amplicon/SSDust_16S.V3V4_W23_Data_Ready.Rdata") # save global env to Rdata file
#load("data/Amplicon/SSD_16S_CLR_EucDist_Ready.Rdata")

head(b.dust.all)
bac.ASV_table[1:4,1:4]
bac.ASV_table[(nrow(bac.ASV_table)-4):(nrow(bac.ASV_table)),(ncol(bac.ASV_table)-4):(ncol(bac.ASV_table))] # last 4 rows & cols
head(dust_meta)
head(dust.meta.surf)

#### Create Centered Log-Ratio Table from ASV table ####
bac.ASV_table[1:4,1:4]
b.clr<-decostand(bac.ASV_table[,-1],method = "clr", pseudocount = 1) #CLR transformation
b.clr[1:4,1:4]
# df must have rownames are SampleIDs, columns are ASV IDs for vegan functions below\

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

# plot time
ggplot(STF.melt, aes(x=SampleID, y=Frequency, fill=SurfaceType))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Surface Type Frequencies in Salton Sea Dust", x="SampleID", y="Frequency", subtitle="",fill="Surface Type")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(name ="Surface Type",values=unique(STF.melt$ST_Color[order(STF.melt$SurfaceType)]),labels=c("Barren Land","Crop Land","Developed","Forest","Herbaceous","Mexico","Open Water","Others","Salton Sea","Shrub"))

ggplot(STF.melt, aes(x=SampleID, y=Frequency, fill=SurfaceType))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Surface Type Frequencies in Salton Sea Dust", x="SampleID", y="Frequency", subtitle="",fill="Surface Type")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(name ="Surface Type",values=unique(STF.melt$ST_Color[order(STF.melt$SurfaceType)]),labels=c("Barren Land","Crop Land","Developed","Forest","Herbaceous","Mexico","Open Water","Others","Salton Sea","Shrub")) +
  facet_wrap(vars(Site), scales = "free")

#### Separate All Data by Time points ####
# create metadata df that will contain scaled chemical data
head(dust.meta.surf)

site_list<-unique(dust.meta.surf$Site) #define an array of string values
# go through metadata & create a list of data frames
## when metadata$Variable == element in site_list (aka x in this case), subset metadata by said element into elements of a list

# here the function(x) is using site_list aka x to subset metadata, when $Variable column == site_list
# Run the function so it's stored in Global Env
site_subsets<-lapply(site_list, function(x) {subset(dust.meta.surf, Site==x)})

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
match_dat<-function(compdata, subset_metadata){
  subset_comp_data = pullrow<-(is.element(row.names(compdata), row.names(subset_metadata)))
  ### * comp data and metadata need to have row names - rownames should be Sample IDs
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
b.dca #DCA1 axis length = 0.47172; use RDA
## The length of first DCA axis:
## > 4 indicates heterogeneous dataset on which unimodal methods should be used (CCA),
##  < 3 indicates homogeneous dataset for which linear methods are suitable (RDA)
## between 3 and 4 both linear and unimodal methods are OK.

# BY MONTH

b.clr_WI.pseudo<-b.clr_WI+1
b.WI.dca = decorana(b.clr_WI.pseudo)
b.WI.dca #DCA1 axis length = 0.54543; use RDA

b.clr_DP.pseudo<-b.clr_DP+1
b.DP.dca = decorana(b.clr_DP.pseudo)
b.DP.dca #DCA1 axis length = 0.64485; use RDA

b.clr_BDC.pseudo<-b.clr_BDC+1
b.BDC.dca = decorana(b.clr_BDC.pseudo)
b.BDC.dca #DCA1 axis length = 0.46958; use RDA

b.clr_PD.pseudo<-b.clr_PD+1
b.PD.dca = decorana(b.clr_PD.pseudo)
b.PD.dca #DCA1 axis length = 0.52948; use RDA

#### RDA w/ All Data ####

rownames(dust.meta.surf) %in% rownames(b.clr) # check order of DFs
head(dust.meta.surf)

rda.all.0<-rda(b.clr ~ BarrenLand+CropLand+Developed+Forest+Herbaceous+Mexico+OpenWater+Others+SaltonSea+Shrub,data=dust.meta.surf)

# check summary of RDA
rda.all.0
summary(rda.all.0)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.all.0) # -0.01712694 -- bad model!
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.all.0, permutations = how(nperm=999)) # p = 0.607

## we can also do a permutation test by RDA axis
#anova(rda.all.0, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.all.0, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#        Df Variance      F Pr(>F)
#Developed   1    567.1 1.5630  0.065 .

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.all.0)
# BarrenLand   CropLand  Developed     Forest Herbaceous     Mexico  OpenWater     Others  SaltonSea
# 19.621677  22.242389   2.211668   8.933791   4.592905   9.905607   2.623602  19.000430  13.792275
# Shrub
# NA

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/

head(dust.meta.surf[,c(23:32)])
## we can use model selection instead of picking variables we think are important (by p values)
# more info on ordistep & ordiR2step here: https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
rda.all.a = ordistep(rda(b.clr ~ 1, data = dust.meta.surf[,c(23:32)]),
                     scope=formula(rda.all.0),
                     direction = "forward",
                     permutations = how(nperm=999))
# b.clr ~ Developed = best model
rda.all.a$anova # see significance of individual terms in model
#                               Df    AIC      F Pr(>F)

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.all.a2 = ordiR2step(rda(b.clr ~ 1, data = dust.meta.surf[,c(23:32)]),
                        scope=formula(rda.all.0),
                        permutations = how(nperm=999))
rda.all.a2$anova # see significance of individual terms in model
#                               R2.adj Df    AIC      F Pr(>F)

# check best fit model based on above results
anova(rda.all.a, permutations = how(nperm=999)) # p =  0.001, significant

# Let's check again removing Shrub
rda.all1<-rda(b.clr ~ BarrenLand+CropLand+Developed+Herbaceous+Others+SaltonSea+Mexico,data=dust.meta.surf)
summary(rda.all1)
RsquareAdj(rda.all1) # how much variation is explained by our model? %
anova(rda.all1, by = "terms", permutations = how(nperm=999)) ### by variables
# Developed  near sig

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.all1)
# BarrenLand   CropLand  Developed Herbaceous     Others  SaltonSea     Mexico
# 12.256833  13.818835   2.180833   3.817406   7.054276  13.679534   9.664598

head(dust.meta.surf[,c(23:32)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.all.b1 = ordistep(rda(b.clr ~ 1, data = dust.meta.surf[,c(23:31)]),
                      scope=formula(rda.all1),
                      direction = "forward",
                      permutations = how(nperm=999))
# b.clr ~  Developed = best model
rda.all.b1$anova # see significance of individual terms in model
#                               Df    AIC      F Pr(>F)
# b.clr ~ Developed

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.all.b2 = ordiR2step(rda(b.clr ~ 1, data = dust.meta.surf[,c(23:31)]),
                        scope=formula(rda.all1),
                        permutations = how(nperm=999))
# b.clr ~  = best model
rda.all.b2$anova # see significance of individual terms in model
#                               R2.adj Df    AIC      F Pr(>F)
# nothing significant

# check best fit model based on above results
anova(rda.all.b1, permutations = how(nperm=999)) # p =  0.001, significant

# compare model fits to each other
anova(rda.all.0, rda.all.b1)

# now let's remove Shrub + Open Water + Others
rda.all2<-rda(b.clr ~  BarrenLand+CropLand+Developed+Herbaceous+Forest+SaltonSea+Mexico,data=dust.meta.surf)
summary(rda.all2)
RsquareAdj(rda.all2) # how much variation is explained by our model?
anova(rda.all2, by = "terms", permutations = how(nperm=999)) ### by variables
# Developed is sig

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.all2)
#BarrenLand   CropLand  Developed Herbaceous     Forest  SaltonSea     Mexico
#16.413506   9.768607   2.172474   3.267892   3.849897  12.015852   9.617201

head(dust.meta.surf[,c(23:32)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.all.c1 = ordistep(rda(b.clr ~ 1, data = dust.meta.surf[,c(23:28,31)]),
                      scope=formula(rda.all2),
                      direction = "forward",
                      permutations = how(nperm=999))
# b.clr ~ Developed
rda.all.c1$anova # see significance of individual terms in model
#                               Df    AIC      F Pr(>F)
# Developed

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.all.c2 = ordiR2step(rda(b.clr ~ 1, data = dust.meta.surf[,c(23:28,31)]),
                        scope=formula(rda.all2),
                        permutations = how(nperm=999))
# Developed
rda.all.c2$anova # see significance of individual terms in model

# messed around with a bunch of variables and Developed is always significant, and Salton Sea doens't contribute much to the variation

rda.all3<-rda(b.clr ~  Developed+SaltonSea,data=dust.meta.surf)
summary(rda.all3)
RsquareAdj(rda.all3) # how much variation is explained by our model?
anova(rda.all3, by = "terms", permutations = how(nperm=999)) ### by variables
# Developed near sig

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.all3)
# BarrenLand   CropLand  Developed  SaltonSea      Shrub
# 43.17708   14.77166   43.34301   30.04297  122.83459

head(dust.meta.surf[,c(23:32)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.all.d1 = ordistep(rda(b.clr ~ 1, data = dust.meta.surf[,c(25,31)]),
                      scope=formula(rda.all3),
                      direction = "forward",
                      permutations = how(nperm=999))
# b.clr ~ Developed
rda.all.d1$anova # see significance of individual terms in model
#                               Df    AIC      F Pr(>F)

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.all.d2 = ordiR2step(rda(b.clr ~ 1, data = dust.meta.surf[,c(25,31)]),
                        scope=formula(rda.all3),
                        permutations = how(nperm=999))
# b.clr ~ Dissolved_OrganicMatter_RFU + Temp_DegC  + DO_%Local = best model
rda.all.d2$anova # see significance of individual terms in model

# # now let's remove Mexico + Open Water + Others + Herbaceous + Shrub
# rda.all4<-rda(b.clr ~  BarrenLand+CropLand+Developed+SaltonSea,data=dust.meta.surf)
# summary(rda.all4)
# RsquareAdj(rda.all4) # how much variation is explained by our model?
# anova(rda.all4, by = "terms", permutations = how(nperm=999)) ### by variables
# # Developed near sig
#
# ## this will help us interpret our RDA and we can see some variable are not significant
# vif.cca(rda.all4)
# # BarrenLand   CropLand  Developed  SaltonSea      Shrub
# # 43.17708   14.77166   43.34301   30.04297  122.83459
#
# head(dust.meta.surf[,c(23:32)])
# ## we can use model selection instead of picking variables we think are important -- based on p values
# rda.all.e1 = ordistep(rda(b.clr ~ 1, data = dust.meta.surf[,c(23:25,27:28,31:32)]),
#                       scope=formula(rda.all4),
#                       direction = "forward",
#                       permutations = how(nperm=999))
# # b.clr ~ Developed
# rda.all.e1$anova # see significance of individual terms in model
# #                               Df    AIC      F Pr(>F)
#
# # Can also use model selection to pick variables by which ones increase variation (R^2)
# rda.all.e2 = ordiR2step(rda(b.clr ~ 1, data = dust.meta.surf[,c(23:25,27:28,31:32)]),
#                         scope=formula(rda.all4),
#                         permutations = how(nperm=999))
# # b.clr ~ Dissolved_OrganicMatter_RFU + Temp_DegC  + DO_%Local = best model
# rda.all.e2$anova # see significance of individual terms in model

#### RDA - Wister ####

rownames(WI) %in% rownames(b.clr_WI) # check order of DFs
head(WI)

# included all Surface type frequencies and it overfit the model so pulling some out...
## dropping forest, open waters, others, and shrub first (after some experimenting)
rda.WI.0<-rda(b.clr_WI ~ BarrenLand+CropLand+Developed+Forest+Herbaceous+Mexico+OpenWater+Others+SaltonSea+Shrub,data=WI)

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
# BarrenLand   CropLand  Developed     Forest Herbaceous     Mexico  OpenWater     Others  SaltonSea      Shrub
# 46.92171   17.93014   14.81736   23.98465  127.78411   32.43178         NA         NA         NA         NA

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(WI)
## we can use model selection instead of picking variables we think are important (by p values)
rda.WI.a = ordistep(rda(b.clr_WI ~ 1, data = WI[,c(23:32)]),
                         scope=formula(rda.WI.0),
                         direction = "forward",
                         permutations = how(nperm=999))
rda.WI.a$anova # see significance of individual terms in model
# OpenWater near sig but appears in 1 sample...

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.WI.a2 = ordiR2step(rda(b.clr_WI ~ 1, data = WI[,c(23:32)]),
                            scope=formula(rda.WI.0),
                            permutations = how(nperm=999))
# too many terms

# check best fit model based on above results
anova(rda.WI.a, permutations = how(nperm=999)) #not significant
#anova(rda.WI.a2, permutations = how(nperm=999)) # p =  0.001, significant

# Dropping Open Water. because it's in 1 sample, dropping Forest because neglible values across samples, & Developed and Others because it's evenly distributed across samples
rda.WI.1<-rda(b.clr_WI ~ BarrenLand+CropLand+Herbaceous+Mexico+SaltonSea+Shrub,data=WI)
summary(rda.WI.1)
RsquareAdj(rda.WI.1) # how much variation is explained by our model?
anova(rda.WI.1, by = "terms", permutations = how(nperm=999)) ### by variables
#                             Df Variance      F Pr(>F)

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.WI.1)
#  BarrenLand   CropLand Herbaceous     Mexico  SaltonSea      Shrub
# 5016.8130   989.6582   336.7180   262.7472  2793.3171 28085.3330

head(WI)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.WI.b1 = ordistep(rda(b.clr_WI ~ 1, data = WI[,c(23:24,27:28,31:32)]),
                          scope=formula(rda.WI.1),
                          direction = "forward",
                          permutations = how(nperm=999))
rda.WI.b1$anova
#
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.WI.b2 = ordiR2step(rda(b.clr_WI ~ 1, data = WI[,c(23:24,27:28,31:32)]),
                            scope=formula(rda.WI.1),
                            permutations = how(nperm=999))
# too many variables

rda.WI.b2$anova
# check best fit model based on above results
anova(rda.WI.b1, permutations = how(nperm=999)) # p =  0.001, significant

anova(rda.WI.0, rda.WI.1) #

# dropping Shrub because much higher VIF...
rda.WI.2<-rda(b.clr_WI ~ BarrenLand+CropLand+Herbaceous+Mexico+SaltonSea,data=WI)
summary(rda.WI.2)
RsquareAdj(rda.WI.2) # how much variation is explained by our model? 0.04954029
anova(rda.WI.2, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant

anova(rda.WI.2, by=NULL,permutations = how(nperm=999)) # p =  0.472

vif.cca(rda.WI.2)
#BarrenLand   CropLand Herbaceous     Mexico  SaltonSea
#17.48460   13.77561   39.03967   30.43097   30.03854

# check if ORP & Sulfide are significantly correlated in August, which they are [strong, sig negative corr]
#cor.test(dust.meta.surf[metadata$SampDate=="WI",]$Sulfide_microM, dust.meta.surf[metadata$SampDate=="WI",]$ORP_mV, method="pearson") # ******

head(WI)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.WI.c1 = ordistep(rda(b.clr_WI ~ 1, data = WI[,c(23:24,27:28,31)]),
                          scope=formula(rda.WI.2),
                          direction = "forward",
                          permutations = how(nperm=999))
#
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.WI.c2 = ordiR2step(rda(b.clr_WI ~ 1, data = WI[,c(23:24,27:28,31)]),
                            scope=formula(rda.WI.2),
                            permutations = how(nperm=999))
# BarrenLand?/ not sig though

# check best fit model based on above results
anova(rda.WI.c1, permutations = how(nperm=999)) #


rda.WI.3<-rda(b.clr_WI ~ BarrenLand+CropLand+Mexico+SaltonSea,data=WI)
summary(rda.WI.3)
RsquareAdj(rda.WI.3) # how much variation is explained by our model? 0.07253236%
anova(rda.WI.3, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)

anova(rda.WI.3, by=NULL,permutations = how(nperm=999)) # p =  0.877

vif.cca(rda.WI.3)
#BarrenLand   CropLand     Mexico  SaltonSea
#8.786964  13.247242  30.220219  23.164521

head(WI)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.WI.d1 = ordistep(rda(b.clr_WI ~ 1, data = WI[,c(23:24,28,31)]),
                          scope=formula(rda.WI.3),
                          direction = "forward",
                          permutations = how(nperm=999))
# nothing significant
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.WI.d2 = ordiR2step(rda(b.clr_WI ~ 1, data = WI[,c(23:24,28,31)]),
                            scope=formula(rda.WI.3),
                            permutations = how(nperm=999))
# nothing sig

rda.WI.4<-rda(b.clr_WI ~ BarrenLand+CropLand+SaltonSea,data=WI)
summary(rda.WI.4)
RsquareAdj(rda.WI.4) # how much variation is explained by our model? -0.2086308
anova(rda.WI.4, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)

anova(rda.WI.4, by=NULL,permutations = how(nperm=999)) # p =  0.909

vif.cca(rda.WI.4)
#BarrenLand   CropLand  SaltonSea
#4.954979   1.769520   6.515210

head(WI)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.WI.e1 = ordistep(rda(b.clr_WI ~ 1, data = WI[,c(23:24,31)]),
                     scope=formula(rda.WI.4),
                     direction = "forward",
                     permutations = how(nperm=999))
# nothing significant
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.WI.e2 = ordiR2step(rda(b.clr_WI ~ 1, data = WI[,c(23:24,31)]),
                       scope=formula(rda.WI.4),
                       permutations = how(nperm=999))
# nothing sig


rda.WI.5<-rda(b.clr_WI ~ BarrenLand+SaltonSea,data=WI)
summary(rda.WI.5)
RsquareAdj(rda.WI.5) # how much variation is explained by our model? -0.02387236
anova(rda.WI.5, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)

anova(rda.WI.5, by=NULL,permutations = how(nperm=999)) # p =  0.909

vif.cca(rda.WI.5)
#BarrenLand  SaltonSea
#4.434495   4.434495

head(WI)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.WI.f1 = ordistep(rda(b.clr_WI ~ 1, data = WI[,c(23,31)]),
                     scope=formula(rda.WI.5),
                     direction = "forward",
                     permutations = how(nperm=999))
# nothing significant
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.WI.f2 = ordiR2step(rda(b.clr_WI ~ 1, data = WI[,c(23,31)]),
                       scope=formula(rda.WI.5),
                       permutations = how(nperm=999))
# nothing sig

anova(rda(b.clr_WI ~ BarrenLand,data=WI), by=NULL,permutations = how(nperm=999)) # p =  0.386


#### RDA - DP ####

rownames(DP) %in% rownames(b.clr_DP) # check order of DFs
head(DP)

rda.DP.0<-rda(b.clr_DP ~ BarrenLand+CropLand+Developed+Forest+Herbaceous+Mexico+OpenWater+Others+SaltonSea+Shrub,data=DP)

# check summary of RDA
rda.DP.0
summary(rda.DP.0)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.DP.0) # 1.5%
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.DP.0, permutations = how(nperm=999)) # not significant

## we can also do a permutation test by RDA axis
#anova(rda.DP.0, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.DP.0, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)
#

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.DP.0)
# BarrenLand    CropLand   Developed      Forest  Herbaceous      Mexico   OpenWater      Others   SaltonSea       Shrub
# 832.929159  180.295453    8.252567   21.344448 1532.560838   11.005239          NA          NA          NA          NA

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(DP)
## we can use model selection instead of picking variables we think are important (by p values)
rda.DP.a = ordistep(rda(b.clr_DP ~ 1, data = DP[,c(23:32)]),
                         scope=formula(rda.DP.0),
                         direction = "forward",
                         permutations = how(nperm=999))
# b.clr_DP ~ BarrenLand   - best model
rda.DP.a$anova # see significance of individual terms in model

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.DP.a2 = ordiR2step(rda(b.clr_DP ~ 1, data = DP[,c(23:32)]),
                            scope=formula(rda.DP.0),
                            permutations = how(nperm=999))
# nothing sig

# check best fit model based on above results
anova(rda.DP.a, permutations = how(nperm=999))
#anova(rda.DP.a2, permutations = how(nperm=999)) # not significant

# Let's get rid of Others and OpenWater
rda.DP.1<-rda(b.clr_DP ~ BarrenLand+CropLand+Developed+Forest+Herbaceous+Mexico+SaltonSea+Shrub,data=DP)
summary(rda.DP.1)
RsquareAdj(rda.DP.1) # how much variation is explained by our model? 4.11%
anova(rda.DP.1, by = "terms", permutations = how(nperm=999)) ### by variables
#                            Df Variance      F Pr(>F)


## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.DP.1)
#BarrenLand    CropLand   Developed      Forest  Herbaceous      Mexico   SaltonSea       Shrub
# 832.929159  180.295453    8.252567   21.344448 1532.560838   11.005239          NA          NA

head(DP)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.DP.b1 = ordistep(rda(b.clr_DP ~ 1, data = DP[,c(23:28,31:32)]),
                          scope=formula(rda.DP.1),
                          direction = "forward",
                          permutations = how(nperm=999))
#b.clr_DP ~ BarrenLand

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.DP.b2 = ordiR2step(rda(b.clr_DP ~ 1, data = DP[,c(23:28,31:32)]),
                            scope=formula(rda.DP.1),
                            permutations = how(nperm=999))
# too many vars

# check best fit model based on above results
anova(rda.DP.b1, permutations = how(nperm=999))

# dropping Herbaceous (super high VIF) and Mexico
rda.DP.2<-rda(b.clr_DP ~ BarrenLand+CropLand+Developed+Forest+SaltonSea+Shrub,data=DP)
summary(rda.DP.2)
RsquareAdj(rda.DP.2) # how much variation is explained by our model? NA
anova(rda.DP.2, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                           Df Variance      F Pr(>F)
#

vif.cca(rda.DP.2)
#BarrenLand    CropLand   Developed      Forest   SaltonSea       Shrub
#335.775690   52.437493   23.362395    3.631513  412.856500 1846.649507

head(DP)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.DP.c1 = ordistep(rda(b.clr_DP ~ 1, data = DP[,c(23:26,31:32)]),
                          scope=formula(rda.DP.2),
                          direction = "forward",
                          permutations = how(nperm=999))
# b.clr_DP ~ BarrenLand

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.DP.c2 = ordiR2step(rda(b.clr_DP ~ 1, data = DP[,c(23:26,31:32)]),
                            scope=formula(rda.DP.2),
                            permutations = how(nperm=999))
# too many terms

# check best fit model based on above result
anova(rda.DP.0, rda.DP.2) # no significant difference

# drop Forest because contribution in surface type freqs is neglible
rda.DP.3<-rda(b.clr_DP ~ BarrenLand+CropLand+Developed+SaltonSea+Shrub,data=DP)
summary(rda.DP.3)
RsquareAdj(rda.DP.3) # how much variation is explained by our model? -0.309622
anova(rda.DP.3, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                           Df Variance      F Pr(>F)

vif.cca(rda.DP.3)
#BarrenLand   CropLand  Developed  SaltonSea      Shrub
#278.52299   52.10296   12.69016  276.88558 1491.04317

head(DP)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.DP.d1 = ordistep(rda(b.clr_DP ~ 1, data = DP[,c(23:25,31:32)]),
                          scope=formula(rda.DP.3),
                          direction = "forward",
                          permutations = how(nperm=999))
# b.clr_DP ~ BarrenLand

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.DP.d2 = ordiR2step(rda(b.clr_DP ~ 1, data = DP[,c(23:25,31:32)]),
                            scope=formula(rda.DP.3),
                            permutations = how(nperm=999))
# nothing significant

# check best fit model based on above result
anova(rda.DP.0, rda.DP.3) # NaN result
anova(rda.DP.2, rda.DP.3) # NaN result

# drop Developed, lowest R^2 contribution
rda.DP.4<-rda(b.clr_DP ~ BarrenLand+CropLand+SaltonSea+Shrub,data=DP)
summary(rda.DP.4)
RsquareAdj(rda.DP.4) # how much variation is explained by our model? 0.1630871
anova(rda.DP.4, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                           Df Variance      F Pr(>F)
#BarrenLand  1   3463.9 2.3437  0.040 *

vif.cca(rda.DP.4)
#BarrenLand   CropLand  SaltonSea      Shrub
#53.46721   29.95787   39.02991  221.31277

head(DP)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.DP.e1 = ordistep(rda(b.clr_DP ~ 1, data = DP[,c(23:24,31:32)]),
                          scope=formula(rda.DP.4),
                          direction = "forward",
                          permutations = how(nperm=999))
# b.clr_DP ~ BarrenLand

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.DP.e2 = ordiR2step(rda(b.clr_DP ~ 1, data = DP[,c(23:24,31:32)]),
                            scope=formula(rda.DP.4),
                            permutations = how(nperm=999))
# b.clr_DP ~ BarrenLand

# forest contribution in DP sits across time is so small so we are dropping this variable
rda.DP.5<-rda(b.clr_DP ~ BarrenLand+SaltonSea+Shrub,data=DP)
summary(rda.DP.5)
RsquareAdj(rda.DP.5) # how much variation is explained by our model? 0.291736
anova(rda.DP.5, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                           Df Variance      F Pr(>F)
#BarrenLand  1   3463.9 2.7694  0.015 *

vif.cca(rda.DP.5)
#BarrenLand  SaltonSea      Shrub
#37.53542   38.63878  128.38471

head(DP)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.DP.f1 = ordistep(rda(b.clr_DP ~ 1, data = DP[,c(23,31:32)]),
                     scope=formula(rda.DP.5),
                     direction = "forward",
                     permutations = how(nperm=999))
# b.clr_DP ~ BarrenLand

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.DP.f2 = ordiR2step(rda(b.clr_DP ~ 1, data = DP[,c(23,31:32)]),
                       scope=formula(rda.DP.5),
                       permutations = how(nperm=999))
# b.clr_DP ~ BarrenLand

rda.DP.6<-rda(b.clr_DP ~ BarrenLand+SaltonSea,data=DP)
summary(rda.DP.6)
RsquareAdj(rda.DP.6) # how much variation is explained by our model? 0.291736
anova(rda.DP.6, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                           Df Variance      F Pr(>F)
#BarrenLand  1   3463.9 2.7694  0.015 *

vif.cca(rda.DP.6)
#BarrenLand  SaltonSea
#2.937236   2.937236

head(DP)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.DP.g1 = ordistep(rda(b.clr_DP ~ 1, data = DP[,c(23,31)]),
                     scope=formula(rda.DP.6),
                     direction = "forward",
                     permutations = how(nperm=999))
# b.clr_DP ~ BarrenLand

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.DP.g2 = ordiR2step(rda(b.clr_DP ~ 1, data = DP[,c(23,31)]),
                       scope=formula(rda.DP.6),
                       permutations = how(nperm=999))
# nothing significant

#### RDA - BDC ####

rownames(BDC) %in% rownames(b.clr_BDC) # check order of DFs
head(BDC)

rda.BDC.0<-rda(b.clr_BDC ~ BarrenLand+CropLand+Developed+Forest+Herbaceous+Mexico+OpenWater+Others+SaltonSea+Shrub,data=BDC)

# check summary of RDA
rda.BDC.0
summary(rda.BDC.0)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.BDC.0) # NA
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.BDC.0, permutations = how(nperm=999)) # not significant

## we can also do a permutation test by RDA axis
#anova(rda.BDC.0, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.BDC.0, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#   nothing sig

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.BDC.0)
#BarrenLand   CropLand  Developed     Forest Herbaceous     Mexico  OpenWater     Others  SaltonSea      Shrub
#759.56624   20.32010   27.46696   32.12073  837.88399  100.97045         NA         NA         NA         NA

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(BDC)
## we can use model selection instead of picking variables we think are important (by p values)
rda.BDC.a = ordistep(rda(b.clr_BDC ~ 1, data = BDC[,c(23:32)]),
                         scope=formula(rda.BDC.0),
                         direction = "forward",
                         permutations = how(nperm=999))
# b.clr_BDC ~ BarrenLand, Developed both near sig
rda.BDC.a$anova # see significance of individual terms in model

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.BDC.a2 = ordiR2step(rda(b.clr_BDC ~ 1, data = BDC[,c(23:32)]),
                            scope=formula(rda.BDC.0),
                            permutations = how(nperm=999))
# too many terms

# check best fit model based on above results
#anova(rda.BDC.a, permutations = how(nperm=999)) # p =  0.036, significant

# drop Mexico, drop Others, drop SaltonSea
rda.BDC.1<-rda(b.clr_BDC ~ BarrenLand+CropLand+Developed+Forest+Herbaceous+OpenWater+Shrub,data=BDC)
summary(rda.BDC.1)
RsquareAdj(rda.BDC.1) # how much variation is explained by our model? NA
anova(rda.BDC.1, by = "terms", permutations = how(nperm=999)) ### by variables
#  nothing significant

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.BDC.1)
# BarrenLand   CropLand  Developed     Forest Herbaceous  OpenWater      Shrub
# 35.714635  86.465824   3.001701 728.362534 621.868135  74.052763         NA

head(BDC)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.BDC.b1 = ordistep(rda(b.clr_BDC ~ 1, data = BDC[,c(23:27,29,32)]),
                          scope=formula(rda.BDC.1),
                          direction = "forward",
                          permutations = how(nperm=999))
# b.clr_BDC ~ BarrenLand, Developed, and Forest are near sig
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.BDC.b2 = ordiR2step(rda(b.clr_BDC ~ 1, data = BDC[,c(23:27,29,32)]),
                            scope=formula(rda.BDC.1),
                            permutations = how(nperm=999))
# too many vars

# check best fit model based on above results
# anova(rda.BDC.b1, permutations = how(nperm=999))
#
# anova(rda.BDC.0, rda.BDC.1) # no significant difference

# drop CropLand and Shrub because these vars do not really change acros samples within BDC
rda.BDC.2<-rda(b.clr_BDC ~ BarrenLand+Developed+Forest+Herbaceous+OpenWater,data=BDC)
summary(rda.BDC.2)
RsquareAdj(rda.BDC.2) # how much variation is explained by our model? 0.3075986
anova(rda.BDC.2, by = "terms", permutations = how(nperm=999)) ### by variables
# Df Variance      F Pr(>F)
# BarrenLand  1  1754.33 2.2512  0.033 *
anova(rda.BDC.2, by = NULL, permutations = how(nperm=999)) ### model not sig

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.BDC.2)
# BarrenLand  Developed     Forest Herbaceous  OpenWater
# 27.306774   2.932190  22.858615  73.237064   2.529495

head(BDC)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.BDC.c1 = ordistep(rda(b.clr_BDC ~ 1, data = BDC[,c(23,25:27,29)]),
                          scope=formula(rda.BDC.2),
                          direction = "forward",
                          permutations = how(nperm=999))
# b.clr_BDC ~ BarrenLand, Developed, Forest all near sig
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.BDC.c2 = ordiR2step(rda(b.clr_BDC ~ 1, data = BDC[,c(23,25:27,29)]),
                            scope=formula(rda.BDC.2),
                            permutations = how(nperm=999))

# BarrenLand near sig
anova(rda.BDC.c1, permutations = how(nperm=999)) # 0.04

anova(rda.BDC.0, rda.BDC.2) # no significant difference

# drop Open water
rda.BDC.3<-rda(b.clr_BDC ~ BarrenLand+Developed+Forest+Herbaceous,data=BDC)
summary(rda.BDC.3)
RsquareAdj(rda.BDC.3) # how much variation is explained by our model? 0.160449
anova(rda.BDC.3, by = "terms", permutations = how(nperm=999)) ### by variables
#               Df Variance      F Pr(>F)
# Barren Land

anova(rda.BDC.3, by = NULL, permutations = how(nperm=999)) # not sig

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.BDC.3)
# BarrenLand  Developed     Forest Herbaceous
# 26.46776    2.66972   19.28559   61.38630

head(BDC)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.BDC.e1 = ordistep(rda(b.clr_BDC ~ 1, data = BDC[,c(23,25:27)]),
                          scope=formula(rda.BDC.3),
                          direction = "forward",
                          permutations = how(nperm=999))
# Barren Land, Developed, and Forest are all near sig
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.BDC.e2 = ordiR2step(rda(b.clr_BDC ~ 1, data = BDC[,c(23,25:27)]),
                            scope=formula(rda.BDC.3),
                            permutations = how(nperm=999))
# Barren land near sig

# drop Herbaceous
rda.BDC.4<-rda(b.clr_BDC ~ BarrenLand+Developed+Forest,data=BDC)
summary(rda.BDC.4)
RsquareAdj(rda.BDC.4) # how much variation is explained by our model? 0.2037182
anova(rda.BDC.4, by = "terms", permutations = how(nperm=999)) ### by variables
#               Df Variance      F Pr(>F)
# BarrenLand  1   1754.3 1.9575  0.045 *

anova(rda.BDC.4, by = NULL, permutations = how(nperm=999)) # 0.089 near sig

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.BDC.4)
# BarrenLand  Developed     Forest
# 5.699074   2.533920   4.758348

head(BDC)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.BDC.f1 = ordistep(rda(b.clr_BDC ~ 1, data = BDC[,c(23,25:26)]),
                      scope=formula(rda.BDC.4),
                      direction = "forward",
                      permutations = how(nperm=999))
# Barren Land, Developed, and Forest are all near sig
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.BDC.f2 = ordiR2step(rda(b.clr_BDC ~ 1, data = BDC[,c(23,25:26)]),
                        scope=formula(rda.BDC.4),
                        permutations = how(nperm=999))
# Barren land near sig

# drop Herbaceous
rda.BDC.5<-rda(b.clr_BDC ~ BarrenLand+Developed,data=BDC)
summary(rda.BDC.5)
RsquareAdj(rda.BDC.5) # how much variation is explained by our model? 0.1373921
anova(rda.BDC.5, by = "terms", permutations = how(nperm=999)) ### by variables
#               Df Variance      F Pr(>F)
# BarrenLand  1   1754.3 1.8070  0.036 *

anova(rda.BDC.5, by = NULL, permutations = how(nperm=999)) # 0.066 near sig

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.BDC.5)
# BarrenLand  Developed
# 1.214723   1.214723

head(BDC)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.BDC.g1 = ordistep(rda(b.clr_BDC ~ 1, data = BDC[,c(23,25)]),
                      scope=formula(rda.BDC.5),
                      direction = "forward",
                      permutations = how(nperm=999))
# Barren Land, Developed,  are all near sig
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.BDC.g2 = ordiR2step(rda(b.clr_BDC ~ 1, data = BDC[,c(23,25)]),
                        scope=formula(rda.BDC.5),
                        permutations = how(nperm=999))
# Barren land near sig

# drop Herbaceous
rda.BDC.6<-rda(b.clr_BDC ~ BarrenLand,data=BDC)
summary(rda.BDC.6)
RsquareAdj(rda.BDC.6) # how much variation is explained by our model? 0.1117436
anova(rda.BDC.6, by = "terms", permutations = how(nperm=999)) ### by variables
#               Df Variance      F Pr(>F)
# BarrenLand  1   1754.3 1.7548  0.063 .

anova(rda.BDC.6, by = NULL, permutations = how(nperm=999)) # 0.064 near sig

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.BDC.6)
# BarrenLand
# 1

head(BDC)
# ## we can use model selection instead of picking variables we think are important -- based on p values
# rda.BDC.h1 = ordistep(rda(b.clr_BDC ~ 1, data = BDC[,c(23)]),
#                       scope=formula(rda.BDC.6),
#                       direction = "forward",
#                       permutations = how(nperm=999))
# # Barren Land, Developed,  are all near sig
# # Can also use model selection to pick variables by which ones increase variation (R^2)
# rda.BDC.h2 = ordiR2step(rda(b.clr_BDC ~ 1, data = BDC[,c(23)]),
#                         scope=formula(rda.BDC.6),
#                         permutations = how(nperm=999))
# # Barren land near sig

#### RDA - PD ####

rownames(PD) %in% rownames(b.clr_PD) # check order of DFs
head(PD)

rda.PD.0<-rda(b.clr_PD ~ BarrenLand+CropLand+Developed+Forest+Herbaceous+Mexico+OpenWater+Others+SaltonSea+Shrub,data=PD)

# check summary of RDA
rda.PD.0
summary(rda.PD.0)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.PD.0) # NA
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.PD.0, permutations = how(nperm=999)) # not significant

## we can also do a permutation test by RDA axis
#anova(rda.PD.0, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.PD.0, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#   nothing sig

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.PD.0)
# BarrenLand   CropLand  Developed     Forest Herbaceous     Mexico  OpenWater     Others  SaltonSea      Shrub
# 1748.73903   27.53067  554.79277  151.34533  242.54621   25.64990         NA         NA         NA         NA

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(PD)
## we can use model selection instead of picking variables we think are important (by p values)
rda.PD.a = ordistep(rda(b.clr_PD ~ 1, data = PD[,c(23:32)]),
                     scope=formula(rda.PD.0),
                     direction = "forward",
                     permutations = how(nperm=999))
# b.clr_PD ~
rda.PD.a$anova # see significance of individual terms in model

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.PD.a2 = ordiR2step(rda(b.clr_PD ~ 1, data = PD[,c(23:32)]),
                        scope=formula(rda.PD.0),
                        permutations = how(nperm=999))
# nothing

# check best fit model based on above results
#anova(rda.PD.a, permutations = how(nperm=999)) # p =  0.036, significant

# dropping STFs that are really neglible - Salton Sea, Mexico, Open Water, Others
rda.PD.1<-rda(b.clr_PD ~ BarrenLand+CropLand+Developed+Forest+Herbaceous+Shrub,data=PD)
summary(rda.PD.1)
RsquareAdj(rda.PD.1) # how much variation is explained by our model? NA
anova(rda.PD.1, by = "terms", permutations = how(nperm=999)) ### by variables
#  nothing significant

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.PD.1)
#   BarrenLand     CropLand    Developed       Forest   Herbaceous        Shrub
# 32885968.960     1434.228 44921912.490    50677.350  4923562.060 12025114.331

head(PD)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.PD.b1 = ordistep(rda(b.clr_PD ~ 1, data = PD[,c(23:27,32)]),
                      scope=formula(rda.PD.1),
                      direction = "forward",
                      permutations = how(nperm=999))
#
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.PD.b2 = ordiR2step(rda(b.clr_PD ~ 1, data = PD[,c(23:27,32)]),
                        scope=formula(rda.PD.1),
                        permutations = how(nperm=999))
# too many vars

# check best fit model based on above results
anova(rda.PD.b1, permutations = how(nperm=999))

anova(rda.PD.0, rda.PD.1) # no significant difference

rda.PD.2<-rda(b.clr_PD ~ DO_Percent_Local+Dissolved_OrganicMatter_RFU+Sulfate_milliM,data=PD)
summary(rda.PD.2)
RsquareAdj(rda.PD.2) # how much variation is explained by our model? 0.47%
anova(rda.PD.2, by = "terms", permutations = how(nperm=999)) ### by variables
# ^ nothing significant
anova(rda.PD.2, by = NULL, permutations = how(nperm=999)) ### model not sig

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.PD.2)
#DO_Percent_Local Dissolved_OrganicMatter_RFU              Sulfate_milliM
#12.670210                   10.431717                    1.770423

head(PD)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.PD.c1 = ordistep(rda(b.clr_PD ~ 1, data = PD[,c(8,14:15)]),
                      scope=formula(rda.PD.2),
                      direction = "forward",
                      permutations = how(nperm=999))
# b.clr_PD ~ DOM
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.PD.c2 = ordiR2step(rda(b.clr_PD ~ 1, data = PD[,c(8,14:15)]),
                        scope=formula(rda.PD.2),
                        permutations = how(nperm=999))

# nothing
anova(rda.PD.c1, permutations = how(nperm=999)) # 0.04

anova(rda.PD.0, rda.PD.2) # no significant difference

# dropping Herbaceous + Shrub because its relatively consistent across samples
rda.PD.3<-rda(b.clr_PD ~ BarrenLand+CropLand+Developed+Forest,data=PD)
summary(rda.PD.3)
RsquareAdj(rda.PD.3) # how much variation is explained by our model? 0.3671824
anova(rda.PD.3, by = "terms", permutations = how(nperm=999)) ### by variables
#               Df Variance      F Pr(>F)
# cropland and develpoed near sig

anova(rda.PD.3, by = NULL, permutations = how(nperm=999)) #

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.PD.3)
# BarrenLand   CropLand  Developed     Forest
# 39.291988   8.753293  63.736732  30.542590

head(PD)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.PD.d1 = ordistep(rda(b.clr_PD ~ 1, data = PD[,c(23:26)]),
                      scope=formula(rda.PD.3),
                      direction = "forward",
                      permutations = how(nperm=999))
#
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.PD.d2 = ordiR2step(rda(b.clr_PD ~ 1, data = PD[,c(23:26)]),
                        scope=formula(rda.PD.3),
                        permutations = how(nperm=999))
# b.clr_PD ~ cropland? not sig though

# dropping developed due to high vif
rda.PD.4<-rda(b.clr_PD ~ BarrenLand+CropLand+Developed,data=PD)
summary(rda.PD.4)
RsquareAdj(rda.PD.4) # how much variation is explained by our model? 0.1452416
anova(rda.PD.4, by = "terms", permutations = how(nperm=999)) ### by variables
#               Df Variance      F Pr(>F)

anova(rda.PD.4, by = NULL, permutations = how(nperm=999)) #

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.PD.4)
# BarrenLand   CropLand  Developed
# 33.217752   8.563484  56.207852

head(PD)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.PD.e1 = ordistep(rda(b.clr_PD ~ 1, data = PD[,c(23:25)]),
                     scope=formula(rda.PD.4),
                     direction = "forward",
                     permutations = how(nperm=999))
#
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.PD.e2 = ordiR2step(rda(b.clr_PD ~ 1, data = PD[,c(23:25)]),
                       scope=formula(rda.PD.4),
                       permutations = how(nperm=999))
# b.clr_PD ~ cropland? not sig though

# dropping developed due to high vif
rda.PD.5<-rda(b.clr_PD ~ BarrenLand+CropLand,data=PD)
summary(rda.PD.5)
RsquareAdj(rda.PD.5) # how much variation is explained by our model? 0.05761236
anova(rda.PD.5, by = "terms", permutations = how(nperm=999)) ### by variables
#               Df Variance      F Pr(>F)

anova(rda.PD.5, by = NULL, permutations = how(nperm=999)) #

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.PD.5)
# Developed  CropLand
# 2.382466  2.382466

head(PD)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.PD.f1 = ordistep(rda(b.clr_PD ~ 1, data = PD[,c(23:24)]),
                     scope=formula(rda.PD.5),
                     direction = "forward",
                     permutations = how(nperm=999))
#
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.PD.f2 = ordiR2step(rda(b.clr_PD ~ 1, data = PD[,c(23:24)]),
                       scope=formula(rda.PD.5),
                       permutations = how(nperm=999))
# b.clr_PD ~ cropland? not sig though...



#### Final RDAs ####
# RDA by sampling timepoint
head(dust.meta.surf)
head(b.clr)
rownames(b.clr) %in% rownames(dust.meta.surf) # sanity check 1

# all data
#rda.all2$call # best model for all data

rda.all<-rda(b.clr ~ Temp_DegC + Dissolved_OrganicMatter_RFU + DO_Percent_Local,data=dust.meta.surf)
rda.all
summary(rda.all)
RsquareAdj(rda.all) # how much variation is explained by our model? 49.56% variation
anova(rda.all, permutations = how(nperm=999)) # p-value = 0.001
anova(rda.all, by = "terms", permutations = how(nperm=999))
#                               Df Variance      F Pr(>F)
# Temp_DegC                    1   335.51 13.8841  0.001 ***
# Dissolved_OrganicMatter_RFU  1   169.84  7.0285  0.001 ***
# DO_Percent_Local             1   113.20  4.6845  0.001 ***
#Residual                    43   992.83
aov.rda.all<-anova(rda.all, by = "terms", permutations = how(nperm=999))
p.adjust(aov.rda.all$`Pr(>F)`,method="bonferroni",n=3) # adjusted pvalues

# WI
#rda.WI.4$call # best model

rda.WI<-rda(b.clr_WI ~ Dissolved_OrganicMatter_RFU+Sulfide_microM,data=WI)
summary(rda.WI)
RsquareAdj(rda.WI) # how much variation is explained by our model? 14.28%
anova(rda.WI, permutations = how(nperm=999)) # p-value = 0.005 **
anova(rda.WI, by = "terms", permutations = how(nperm=999))
#                           Df Variance      F Pr(>F)
#Dissolved_OrganicMatter_RFU  1   157.93 1.9302  0.001 ***
#Sulfide_microM               1   101.16 1.2363  0.148
#Residual                     5   409.10
aov.rda.aug<-anova(rda.WI, by = "terms", permutations = how(nperm=999))
p.adjust(aov.rda.aug$`Pr(>F)`,method="bonferroni",n=3) # adjusted pvalues

# DP
#rda.DP.2$call # best model from above

rda.DP<-rda(b.clr_DP ~ ORP_mV + Sulfate_milliM,data=DP)
summary(rda.DP)
RsquareAdj(rda.DP) # how much variation is explained by our model? 5.3%
anova(rda.DP, permutations = how(nperm=999)) # p-value = 0.005
anova(rda.DP, by = "terms", permutations = how(nperm=999))
#                 Df Variance      F Pr(>F)
# ORP_mV          1   77.391 1.3240  0.003 **
# Sulfate_milliM  1   62.509 1.0694  0.192
# Residual        5  292.258
aov.rda.dec<-anova(rda.DP, by = "terms", permutations = how(nperm=999))
p.adjust(aov.rda.dec$`Pr(>F)`,method="bonferroni",n=3) # adjusted pvalues

# BDC
#rda.BDC.3$call  #best mode

rda.BDC<-rda(b.clr_BDC ~ Dissolved_OrganicMatter_RFU + Sulfate_milliM,data=BDC)
summary(rda.BDC)
RsquareAdj(rda.BDC) # how much variation is explained by our model? 2.61%
anova(rda.BDC, permutations = how(nperm=999)) # p-value = 0.039
anova(rda.BDC, by = "terms", permutations = how(nperm=999))
#                           Df Variance      F Pr(>F)
# Dissolved_OrganicMatter_RFU  1   60.077 1.1195  0.011 *
# Sulfate_milliM               1   57.312 1.0680  0.092 .
# Residual                     5  268.313
aov.rda.apr<-anova(rda.BDC, by = "terms", permutations = how(nperm=999))
p.adjust(aov.rda.apr$`Pr(>F)`,method="bonferroni",n=3) # adjusted pvalues

# save RDAs as R object
save.image("data/SSW_Amplicon_EnvDriver_RDAsOnly.Rdata")

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
RsquareAdj(rda.all) # 49.56%
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
##anova(rda.all, permutations = how(nperm=999)) # p = 0.001, significant

png('figures/EnvDrivers/SSW_AllData_autoplot_rda_example.png',width = 700, height = 600, res=100)
autoplot(rda.all, arrows = TRUE,data = rda.all ,layers=c("biplot","sites"),label = FALSE, label.size = 3, shape = FALSE, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3, scale= 0)+theme_classic()
dev.off()
## FOR AUTOPLOT -> must load packagve ggvegan first


# variance partitioning of RDA
rda.all.part<-varpart(b.clr, dust.meta.surf$Temp_DegC, dust.meta.surf$Dissolved_OrganicMatter_RFU,dust.meta.surf$DO_Percent_Local)
rda.all.part$part
# plot variance partitioning results
png('figures/EnvDrivers/SSW_AllData_RDA_VariancePartitioning.png',width = 900, height = 900, res=100)
plot(rda.all.part,
     Xnames = c("Temp (C)", "DOM (RFU)","DO%"), # name the partitions
     bg = c("#ef476f", "#ffbe0b","skyblue"), alpha = 80, # colour the circles
     digits = 3, # only show 2 digits
     cex = 1.5)
dev.off()

rda.sum.all<-summary(rda.all)
rda.sum.all$sites[,1:2]
rda.sum.all$cont #cumulative proportion of variance per axis
# RDA1 = 30.8, RDA2 = 23.65

# create data frame w/ RDA axes for sites
# first check rownames of RDA & metadata, then make df
rownames(rda.sum.all$sites) %in% rownames(dust.meta.surf)
rda.axes.all<-data.frame(RDA1=rda.sum.all$sites[,1], RDA2=rda.sum.all$sites[,2], SampleID=rownames(rda.sum.all$sites), Depth_m=dust.meta.surf$Depth_m, SampDate=dust.meta.surf$SampDate)

# create data frame w/ RDA axes for variables
arrows.all<-data.frame(RDA1=rda.sum.all$biplot[,1], RDA2=rda.sum.all$biplot[,2], Label=rownames(rda.sum.all$biplot))
#arrows.all$Label[(arrows.all$Label) == "ORP_mV"] <- "ORP (mV)"
arrows.all$Label[(arrows.all$Label) == "Dissolved_OrganicMatter_RFU"] <- "DOM (RFU)"
arrows.all$Label[(arrows.all$Label) == "DO_Percent_Local"] <- "DO %"
arrows.all$Label[(arrows.all$Label) == "Temp_DegC"] <- "Temp (C)"

rda.sum.all$cont #cumulative proportion of variance per axis
# RDA1=30.8%, RDA2=23.65%

rda.plot1<-ggplot(rda.axes.all, aes(x = RDA1, y = RDA2)) + geom_point(size=2) +
  geom_segment(data = arrows.all,mapping = aes(x = 0, y = 0, xend = RDA1, yend = RDA2),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.5,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.all,aes(label = Label, x = RDA1, y = RDA2, fontface="bold"))+
  coord_fixed() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))

rda.plot2<-ggplot(rda.axes.all, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=as.numeric(as.character(Depth_m)),shape=SampDate),size=4) +
  geom_segment(data = arrows.all,mapping = aes(x = 0, y = 0, xend = RDA1*5.5, yend = RDA2*5.5),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.8,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.all,aes(label = Label, x = RDA1*7, y = RDA2*7, fontface="bold"), size=4)+
  coord_fixed(ratio = 1, xlim = c(-8,8), ylim = c(-8,8)) + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  scale_shape_discrete(labels=c("WI","DP","BDC"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Seawater",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [30.80%]") + ylab("RDA2 [23.65%]")

ggsave(rda.plot2,filename = "figures/EnvDrivers/SSW_16S_RDA_AllData.png", width=10, height=10, dpi=600)


rda.plot3<-ggplot(rda.axes.all, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=as.numeric(as.character(Depth_m)),shape=SampDate),size=5) +
  geom_segment(data = arrows.all,mapping = aes(x = 0, y = 0, xend = RDA1*6, yend = RDA2*6),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 1,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.all,aes(label = Label, x = RDA1*8, y = RDA2*8, fontface="bold"), size=5)+
  coord_fixed(ratio = 1, xlim = c(-10,10), ylim = c(-10,10)) + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  scale_shape_discrete(labels=c("WI","DP","BDC"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Seawater",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [30.80%]") + ylab("RDA2 [23.65%]")

ggsave(rda.plot3,filename = "figures/EnvDrivers/SSW_16S_RDA_AllData_bigger.png", width=15, height=15, dpi=600)

#### Plot RDA - Aug 2021 ####
#plot(rda.WI) # depending on how many species you have, this step may take a while
plot(rda.WI, scaling = 1)
## scaling = 1 -> emphasizes relationships among sites
plot(rda.WI, scaling = 2)
## scaling = 2 -> emphasizes relationships among species

# check summary of RDA
summary(rda.WI)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.WI) # 14.28%
## ^^ use this b/c chance correlations can inflate R^2

png('figures/EnvDrivers/SSW_WI_autoplot_rda_example.png',width = 700, height = 600, res=100)
autoplot(rda.WI, arrows = TRUE,data = rda.WI ,layers=c("biplot","sites"),label = FALSE, label.size = 3, shape = FALSE, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3, scale= 0)+theme_classic()
dev.off()
## FOR AUTOPLOT -> must load packagve ggvegan first

# variance partitioning of RDA
rda.WI.part<-varpart(b.clr_WI, WI$Dissolved_OrganicMatter_RFU, WI$Sulfide_microM)
rda.WI.part$part
# plot variance partitioning results
png('figures/EnvDrivers/SSW_WI_RDA_VariancePartitioning.png',width = 900, height = 900, res=100)
plot(rda.WI.part,
     Xnames = c("DOM (RFU)", "Sulfide (microM)"), # name the partitions
     bg = c("#ffbe0b", "darkgreen"), alpha = 80, # colour the circles
     digits = 3, # only show 3 digits
     cex = 1.5)
dev.off()

rda.sum.WI<-summary(rda.WI)
rda.sum.WI$sites[,1:2]
rda.sum.WI$cont #cumulative proportion of variance per axis
# RDA1=26.71%, RDA2=12.06%

# create data frame w/ RDA axes for sites
rda.axes.WI<-data.frame(RDA1=rda.sum.WI$sites[,1], RDA2=rda.sum.WI$sites[,2], SampleID=rownames(rda.sum.WI$sites), Depth_m=WI$Depth_m)

# create data frame w/ RDA axes for variables
arrows.WI<-data.frame(RDA1=rda.sum.WI$biplot[,1], RDA2=rda.sum.WI$biplot[,2], Label=rownames(rda.sum.WI$biplot))
#arrows.WI$Label[(arrows.WI$Label) == "ORP_mV"] <- "ORP (mV)"
arrows.WI$Label[(arrows.WI$Label) == "Dissolved_OrganicMatter_RFU"] <- "DOM (RFU)"
#arrows.WI$Label[(arrows.WI$Label) == "Sulfate_milliM"] <- "Sulfate (milliM)"
arrows.WI$Label[(arrows.WI$Label) == "Sulfide_microM"] <- "Sulfide (microM)"

rda.plot5<-ggplot(rda.axes.WI, aes(x = RDA1, y = RDA2)) + geom_point(size=2) +
  geom_segment(data = arrows.WI,mapping = aes(x = 0, y = 0, xend = RDA1, yend = RDA2),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.5,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.WI,aes(label = Label, x = RDA1, y = RDA2, fontface="bold"))+
  coord_fixed() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))

rda.plot6<-ggplot(rda.axes.WI, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=as.numeric(as.character(Depth_m))),size=4) +
  geom_segment(data = arrows.WI,mapping = aes(x = 0, y = 0, xend = RDA1*8, yend = RDA2*8),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.8,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.WI,aes(label = Label, x = RDA1*9.85, y = RDA2*9.5, fontface="bold"), size=4)+
  coord_fixed(ratio = 1, xlim = c(-5,15), ylim = c(-10,10)) + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Seawater, WI",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [26.71%]") + ylab("RDA2 [12.06%]")

ggsave(rda.plot6,filename = "figures/EnvDrivers/SSW_16S_RDA_WI.png", width=16, height=12, dpi=600)

rda.plot6b<-ggplot(rda.axes.WI, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=as.numeric(as.character(Depth_m))),size=5) +
  geom_segment(data = arrows.WI,mapping = aes(x = 0, y = 0, xend = RDA1*8, yend = RDA2*8),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 1,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.WI,aes(label = Label, x = RDA1*9, y = RDA2*9.5, fontface="bold"), size=5)+
  coord_fixed(ratio = 1, xlim = c(-10,10), ylim = c(-10,10)) + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Seawater, WI",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [26.71%]") + ylab("RDA2 [12.06%]")

ggsave(rda.plot6b,filename = "figures/EnvDrivers/SSW_16S_RDA_WI_bigger.png", width=15, height=15, dpi=600)

#### Plot RDA - Dec 2021 ####
#plot(rda.DP) # depending on how many species you have, this step may take a while
plot(rda.DP, scaling = 1)
## scaling = 1 -> emphasizes relationships among sites
plot(rda.DP, scaling = 2)
## scaling = 2 -> emphasizes relationships among species

# check summary of RDA
summary(rda.DP)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.DP) # 0.0532124
## ^^ use this b/c chance correlations can inflate R^2

png('figures/EnvDrivers/SSW_DP_autoplot_rda_example.png',width = 700, height = 600, res=100)
autoplot(rda.DP, arrows = TRUE,data = rda.DP ,layers=c("biplot","sites"),label = FALSE, label.size = 3, shape = FALSE, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3, scale= 0)+theme_classic()
dev.off()
## FOR AUTOPLOT -> must load packagve ggvegan first

# variance partitioning of RDA
rda.DP.part<-varpart(b.clr_DP, DP$ORP_mV, DP$Sulfate_milliM)
rda.DP.part$part
# plot variance partitioning results
png('figures/EnvDrivers/SSW_DP_RDA_VariancePartitioning.png',width = 900, height = 900, res=100)
plot(rda.DP.part,
     Xnames = c("ORP (mV)", "Sulfate (milliM)"), # name the partitions
     bg = c("#3a0ca3", "#8ac926"), alpha = 80, # colour the circles
     digits = 3, # only show 3 digits
     cex = 1.5)
dev.off()

rda.sum.DP<-summary(rda.DP)
rda.sum.DP$sites[,1:2]
rda.sum.DP$cont # cumulative proportion of variation per axis
# RDA1 = 18.19, RDA2 = 14.18

# create data frame w/ RDA axes for sites
rda.axes.DP<-data.frame(RDA1=rda.sum.DP$sites[,1], RDA2=rda.sum.DP$sites[,2], SampleID=rownames(rda.sum.DP$sites), Depth_m=DP$Depth_m)

# create data frame w/ RDA axes for variables
arrows.DP<-data.frame(RDA1=rda.sum.DP$biplot[,1], RDA2=rda.sum.DP$biplot[,2], Label=rownames(rda.sum.DP$biplot))
arrows.DP$Label[(arrows.DP$Label) == "ORP_mV"] <- "ORP (mV)"
arrows.DP$Label[(arrows.DP$Label) == "Sulfate_milliM"] <- "Sulfate (milliM)"
#arrows.DP$Label[(arrows.DP$Label) == "DO_Percent_Local"] <- "DO%"

rda.plot7<-ggplot(rda.axes.DP, aes(x = RDA1, y = RDA2)) + geom_point(size=2) +
  geom_segment(data = arrows.DP,mapping = aes(x = 0, y = 0, xend = RDA1, yend = RDA2),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.5,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.DP,aes(label = Label, x = RDA1, y = RDA2, fontface="bold"))+
  coord_fixed() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))

rda.plot8<-ggplot(rda.axes.DP, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=as.numeric(as.character(Depth_m))),size=4) +
  geom_segment(data = arrows.DP,mapping = aes(x = 0, y = 0, xend = RDA1*9, yend = RDA2*9),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.8,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.DP,aes(label = Label, x = RDA1*10.5, y = RDA2*10, fontface="bold"), size=4)+
  coord_fixed(ratio = 1, xlim = c(-10,11), ylim = c(-10,10)) + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Seawater",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [18.19%]") + ylab("RDA2 [14.18%]")

ggsave(rda.plot8,filename = "figures/EnvDrivers/SSW_16S_RDA_DP.png", width=15, height=12, dpi=600)

rda.plot8b<-ggplot(rda.axes.DP, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=as.numeric(as.character(Depth_m))),size=5) +
  geom_segment(data = arrows.DP,mapping = aes(x = 0, y = 0, xend = RDA1*9, yend = RDA2*9),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 1,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.DP,aes(label = Label, x = RDA1*10.5, y = RDA2*10.5, fontface="bold"), size=5)+
  coord_fixed(ratio = 1, xlim = c(-10,11), ylim = c(-10,10)) + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Seawater, DP",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [18.19%]") + ylab("RDA2 [14.18%]")

ggsave(rda.plot8b,filename = "figures/EnvDrivers/SSW_16S_RDA_DP_bigger.png", width=15, height=15, dpi=600)

#### Plot RDA - Apr 2022 ####
#plot(rda.DP) # depending on how many species you have, this step may take a while
plot(rda.BDC, scaling = 1)
## scaling = 1 -> emphasizes relationships among sites
plot(rda.BDC, scaling = 2)
## scaling = 2 -> emphasizes relationships among species

# check summary of RDA
summary(rda.BDC)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.BDC) # 2.61%
## ^^ use this b/c chance correlations can inflate R^2

png('figures/EnvDrivers/SSW_BDC_autoplot_rda_example.png',width = 700, height = 600, res=100)
autoplot(rda.BDC, arrows = TRUE,data = rda.BDC ,layers=c("biplot","sites"),label = FALSE, label.size = 3, shape = FALSE, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3, scale= 0)+theme_classic()
dev.off()
## FOR AUTOPLOT -> must load packagve ggvegan first

# variance partitioning of RDA
rda.BDC.part<-varpart(b.clr_BDC, BDC$Dissolved_OrganicMatter_RFU, BDC$Sulfate_milliM)
rda.BDC.part$part
# plot variance partitioning results
png('figures/EnvDrivers/SSW_BDC_RDA_VariancePartitioning.png',width = 900, height = 900, res=100)
plot(rda.BDC.part,
     Xnames = c("DOM (RFU)", "Sulfate (milliM)"), # name the partitions
     bg = c("#ffbe0b", "#8ac926"), alpha = 80, # colour the circles
     digits = 3, # only show 3 digits
     cex = 1.5)
dev.off()

rda.sum.BDC<-summary(rda.BDC)
rda.sum.BDC$sites[,1:2]
rda.sum.BDC$cont
# RDA1 = 16.04, RDA2 = 14.40

# create data frame w/ RDA axes for sites
rda.axes.BDC<-data.frame(RDA1=rda.sum.BDC$sites[,1], RDA2=rda.sum.BDC$sites[,2], SampleID=rownames(rda.sum.BDC$sites), Depth_m=BDC$Depth_m)

# create data frame w/ RDA axes for variables
arrows.BDC<-data.frame(RDA1=rda.sum.BDC$biplot[,1], RDA2=rda.sum.BDC$biplot[,2], Label=rownames(rda.sum.BDC$biplot))
arrows.BDC$Label[(arrows.BDC$Label) == "Dissolved_OrganicMatter_RFU"] <- "DOM (RFU)"
arrows.BDC$Label[(arrows.BDC$Label) == "Sulfate_milliM"] <- "Sulfate (milliM)"
#arrows.BDC$Label[(arrows.BDC$Label) == "DO_Percent_Local"] <- "DO%"
#arrows.BDC$Label[(arrows.BDC$Label) == "Temp_DegC"] <- "Temp (C)"

rda.plot9<-ggplot(rda.axes.BDC, aes(x = RDA1, y = RDA2)) + geom_point(size=2) +
  geom_segment(data = arrows.BDC,mapping = aes(x = 0, y = 0, xend = RDA1, yend = RDA2),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.5,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.BDC,aes(label = Label, x = RDA1, y = RDA2, fontface="bold"))+
  coord_fixed() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))

rda.plot10<-ggplot(rda.axes.BDC, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=as.numeric(as.character(Depth_m))),size=4) +
  geom_segment(data = arrows.BDC,mapping = aes(x = 0, y = 0, xend = RDA1*8, yend = RDA2*8),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.8,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.BDC,aes(label = Label, x = RDA1*9, y = RDA2*9, fontface="bold"), size=4)+
  coord_fixed(ratio = 1, xlim = c(-10,10), ylim = c(-10,10)) + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Seawater, BDC",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [16.04%]") + ylab("RDA2 [14.40%]")

ggsave(rda.plot10,filename = "figures/EnvDrivers/SSW_16S_RDA_April2022.png", width=15, height=12, dpi=600)

rda.plot10b<-ggplot(rda.axes.BDC, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=as.numeric(as.character(Depth_m))),size=5) +
  geom_segment(data = arrows.BDC,mapping = aes(x = 0, y = 0, xend = RDA1*9, yend = RDA2*9),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 1,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.BDC,aes(label = Label, x = RDA1*10, y = RDA2*10, fontface="bold"), size=5)+
  coord_fixed(ratio = 1, xlim = c(-10,10), ylim = c(-10,10)) + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Seawater, DP",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [16.04%]") + ylab("RDA2 [14.40]")

ggsave(rda.plot10b,filename = "figures/EnvDrivers/SSW_16S_RDA_April2022_bigger.png", width=15, height=15, dpi=600)

#### Save Progress ####

save.image("data/SSW_Amplicon_EnvDriver.Rdata")
