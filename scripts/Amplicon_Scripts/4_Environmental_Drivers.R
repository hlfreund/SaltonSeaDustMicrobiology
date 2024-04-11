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

rownames(meta.all.scaled) %in% rownames(b.clr) # check order of DFs
head(meta.all.scaled)

rda.all.0<-rda(b.clr ~ ave.air_temp+ave.wind_speed+ave.wind_gust+ave.wind_direction+BarrenLand+CropLand+Developed+Forest+Herbaceous+Mexico+OpenWater+Others+SaltonSea+Shrub,data=meta.all.scaled)

# check summary of RDA
rda.all.0
summary(rda.all.0)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.all.0) # -0.01438676 -- bad model!
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
#ave.wind_speed      1    519.7 1.4361  0.099 .
#ave.wind_gust       1    639.1 1.7661  0.039 *

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.all.0)
# ave.air_temp     ave.wind_speed      ave.wind_gust ave.wind_direction         BarrenLand           CropLand
# 9.103113         120.337131          79.913952           6.467531          53.587445          27.364864
# Developed             Forest         Herbaceous             Mexico          OpenWater             Others
# 14.033592          18.422991          20.376398          20.427801           7.477958          31.848762
# SaltonSea              Shrub
# 17.653942                 NA

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/

head(meta.all.scaled[,c(5:7,10,33:42)])
## we can use model selection instead of picking variables we think are important (by p values)
# more info on ordistep & ordiR2step here: https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
rda.all.a = ordistep(rda(b.clr ~ 1, data = meta.all.scaled[,c(5:7,10,33:42)]),
                     scope=formula(rda.all.0),
                     direction = "forward",
                     permutations = how(nperm=999))
# b.clr ~ ave.wind_gust + ave.air_temp
rda.all.a$anova # see significance of individual terms in model
#                               Df    AIC      F Pr(>F)

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.all.a2 = ordiR2step(rda(b.clr ~ 1, data = meta.all.scaled[,c(5:7,10,33:42)]),
                        scope=formula(rda.all.0),
                        permutations = how(nperm=999))
rda.all.a2$anova # see significance of individual terms in model
#                               R2.adj Df    AIC      F Pr(>F)
# nothing sig

# check best fit model based on above results
anova(rda.all.a, permutations = how(nperm=999)) # p =  0.007, significant

# Let's check again removing OpenWater & Others since they overall have little contribution
rda.all1<-rda(b.clr ~ ave.air_temp+ave.wind_speed+ave.wind_gust+ave.wind_direction+BarrenLand+CropLand+Developed+Herbaceous+SaltonSea+Mexico+Shrub+Forest,data=meta.all.scaled)
summary(rda.all1)
RsquareAdj(rda.all1) # how much variation is explained by our model? % -0.004701317
anova(rda.all1, by = "terms", permutations = how(nperm=999)) ### by variables
#                   Df Variance      F Pr(>F)
# ave.wind_speed      1    519.7 1.4559  0.081 .
# ave.wind_gust       1    639.1 1.7903  0.042 *

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.all1)
# ave.air_temp     ave.wind_speed      ave.wind_gust ave.wind_direction         BarrenLand           CropLand
# 5.077477         120.035138          79.141682           4.890668        1284.095777         108.157727
# Developed         Herbaceous          SaltonSea             Mexico              Shrub             Forest
# 780.066017          49.286151         374.765973          21.207899        1910.187696          31.580231

head(meta.all.scaled[,c(5:7,10,33:38,41:42)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.all.b1 = ordistep(rda(b.clr ~ 1, data = meta.all.scaled[,c(5:7,10,33:38,41:42)]),
                      scope=formula(rda.all1),
                      direction = "forward",
                      permutations = how(nperm=999))
# b.clr ~  Developed = best model
rda.all.b1$anova # see significance of individual terms in model
#                               Df    AIC      F Pr(>F)
# b.clr ~ ave.wind_gust + ave.air_temp

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.all.b2 = ordiR2step(rda(b.clr ~ 1, data = meta.all.scaled[,c(5:7,10,33:38,41:42)]),
                        scope=formula(rda.all1),
                        permutations = how(nperm=999))
# nothing
rda.all.b2$anova # see significance of individual terms in model
#                               R2.adj Df    AIC      F Pr(>F)
# nothing significant

# check best fit model based on above results
anova(rda.all.b1, permutations = how(nperm=999)) # p =  0.007, significant

# compare model fits to each other
anova(rda.all.0, rda.all.b1)

# now let's remove Shrub + Open Water + Others
rda.all2<-rda(b.clr ~ ave.air_temp+ave.wind_speed+ave.wind_gust+ave.wind_direction+BarrenLand+CropLand+Developed+Herbaceous+Forest+SaltonSea+Mexico,data=meta.all.scaled)
summary(rda.all2)
RsquareAdj(rda.all2) # how much variation is explained by our model? # 0.005726867
anova(rda.all2, by = "terms", permutations = how(nperm=999)) ### by variables
# ave.wind_speed      1    519.7 1.4652  0.090 .
# ave.wind_gust       1    639.1 1.8018  0.038 *
# Developed           1    491.8 1.3865  0.096 .

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.all2)
#ave.air_temp     ave.wind_speed      ave.wind_gust ave.wind_direction         BarrenLand           CropLand          Developed
# 3.692272         119.215337          78.545885           4.865429          49.363743          12.921990          12.579899
# Herbaceous             Forest          SaltonSea             Mexico
# 17.308532          11.547484          13.962061          15.298904

head(meta.all.scaled[,c(5:7,10,33:38,41)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.all.c1 = ordistep(rda(b.clr ~ 1, data = meta.all.scaled[,c(5:7,10,33:38,41)]),
                      scope=formula(rda.all2),
                      direction = "forward",
                      permutations = how(nperm=999))
# b.clr ~ ave.wind_gust + ave.air_temp
rda.all.c1$anova # see significance of individual terms in model
#                               Df    AIC      F Pr(>F)

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.all.c2 = ordiR2step(rda(b.clr ~ 1, data = meta.all.scaled[,c(5:7,10,33:38,41)]),
                        scope=formula(rda.all2),
                        permutations = how(nperm=999))
#
rda.all.c2$anova # see significance of individual terms in model

# messed around with a bunch of variables and Developed is always significant, and Salton Sea doens't contribute much to the variation

# Remove Forest & Mexico from model
rda.all3<-rda(b.clr ~  ave.air_temp+ave.wind_speed+ave.wind_gust+ave.wind_direction+BarrenLand+CropLand+Developed+Herbaceous+SaltonSea,data=meta.all.scaled)
summary(rda.all3)
RsquareAdj(rda.all3) # how much variation is explained by our model? 0.01872894
anova(rda.all3, by = "terms", permutations = how(nperm=999)) ### by variables
# ave.wind_speed      1    519.7 1.4846  0.064 .
# ave.wind_gust       1    639.1 1.8257  0.032 *

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.all3)
#ave.air_temp     ave.wind_speed      ave.wind_gust ave.wind_direction         BarrenLand           CropLand          Developed
# 2.905281          92.917279          63.741495           3.483810          29.758033           5.158057           9.875250
# Herbaceous          SaltonSea
# 12.380115          11.130141

head(meta.all.scaled[,c(5:7,10,33:35,37,41)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.all.d1 = ordistep(rda(b.clr ~ 1, data = meta.all.scaled[,c(5:7,10,33:35,37,41)]),
                      scope=formula(rda.all3),
                      direction = "forward",
                      permutations = how(nperm=999))
# b.clr ~ ave.wind_gust + ave.air_temp
rda.all.d1$anova # see significance of individual terms in model
#                               Df    AIC      F Pr(>F)

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.all.d2 = ordiR2step(rda(b.clr ~ 1, data = meta.all.scaled[,c(5:7,10,33:35,37,41)]),
                        scope=formula(rda.all3),
                        permutations = how(nperm=999))
# b.clr ~ Dissolved_OrganicMatter_RFU + Temp_DegC  + DO_%Local = best model
rda.all.d2$anova # see significance of individual terms in model

# Remove Herbaceous & BarrenLand
rda.all4<-rda(b.clr ~  ave.air_temp+ave.wind_speed+ave.wind_gust+ave.wind_direction+CropLand+Developed+SaltonSea,data=meta.all.scaled)
summary(rda.all4)
RsquareAdj(rda.all4) # how much variation is explained by our model? 0.02960742
anova(rda.all4, by = "terms", permutations = how(nperm=999)) ### by variables
#ave.air_temp        1    471.3 1.3613  0.095 .
# ave.wind_speed      1    519.7 1.5013  0.074 .
#ave.wind_gust       1    639.1 1.8461  0.025 *

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.all4)
#ave.air_temp     ave.wind_speed      ave.wind_gust ave.wind_direction           CropLand          Developed          SaltonSea
# 2.475890          25.171008          24.561732           2.676141           4.904833           3.747255           5.554783

head(meta.all.scaled[,c(5:7,10,34:35,37,41)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.all.d1 = ordistep(rda(b.clr ~ 1, data = meta.all.scaled[,c(5:7,10,34:35,37,41)]),
                      scope=formula(rda.all4),
                      direction = "forward",
                      permutations = how(nperm=999))
# b.clr ~ ave.wind_gust + ave.air_temp + Developed slightly
rda.all.d1$anova # see significance of individual terms in model
#                               Df    AIC      F Pr(>F)

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.all.d2 = ordiR2step(rda(b.clr ~ 1, data = meta.all.scaled[,c(5:7,10,34:35,37,41)]),
                        scope=formula(rda.all4),
                        permutations = how(nperm=999))
#nothing
rda.all.d2$anova # see significance of individual terms in model

# Remove CropLand
rda.all5<-rda(b.clr ~  ave.air_temp+ave.wind_speed+ave.wind_gust+ave.wind_direction+Developed+SaltonSea,data=meta.all.scaled)
summary(rda.all5)
RsquareAdj(rda.all5) # how much variation is explained by our model? 0.05818993
anova(rda.all5, by = "terms", permutations = how(nperm=999)) ### by variables
#ave.air_temp        1    471.3 1.3613  0.095 .
# ave.wind_speed      1    519.7 1.5013  0.074 .
#ave.wind_gust       1    639.1 1.8461  0.025 *

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.all5)
#ave.air_temp     ave.wind_speed      ave.wind_gust ave.wind_direction          Developed          SaltonSea
#2.393376          18.716710          19.579387           2.558826           3.743855           4.396717

head(meta.all.scaled[,c(5:7,10,35,41)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.all.e1 = ordistep(rda(b.clr ~ 1, data = meta.all.scaled[,c(5:7,10,35,41)]),
                      scope=formula(rda.all5),
                      direction = "forward",
                      permutations = how(nperm=999))
# b.clr ~ ave.wind_gust + ave.air_temp
rda.all.e1$anova # see significance of individual terms in model
#                               Df    AIC      F Pr(>F)

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.all.e2 = ordiR2step(rda(b.clr ~ 1, data = meta.all.scaled[,c(5:7,10,35,41)]),
                        scope=formula(rda.all5),
                        permutations = how(nperm=999))
#nothing
rda.all.e2$anova # see significance of individual terms in model

# drop SaltonSea + Wind Direction
rda.all6<-rda(b.clr ~  ave.air_temp+ave.wind_speed+ave.wind_gust+Developed,data=meta.all.scaled)
summary(rda.all6)
RsquareAdj(rda.all6) # how much variation is explained by our model? 0.07339442
anova(rda.all6, by = "terms", permutations = how(nperm=999)) ### by variables
#ave.air_temp    1    471.3 1.4257  0.097 .
# ave.wind_speed  1    519.7 1.5722  0.050 *
#   ave.wind_gust   1    639.1 1.9334  0.016 *
#   Developed       1    399.1 1.2074  0.194

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.all6)
#ave.air_temp ave.wind_speed  ave.wind_gust      Developed
#1.126212       9.471902      14.179739       3.148465

head(meta.all.scaled[,c(5:7,35)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.all.f1 = ordistep(rda(b.clr ~ 1, data = meta.all.scaled[,c(5:7,35)]),
                      scope=formula(rda.all6),
                      direction = "forward",
                      permutations = how(nperm=999))
# b.clr ~ ave.wind_gust + ave.air_temp
rda.all.f1$anova # see significance of individual terms in model
#                               Df    AIC      F Pr(>F)

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.all.f2 = ordiR2step(rda(b.clr ~ 1, data = meta.all.scaled[,c(5:7,35)]),
                        scope=formula(rda.all6),
                        permutations = how(nperm=999))
#nothing
rda.all.f2$anova # see significance of individual terms in model

# drop wind speed
rda.all7<-rda(b.clr ~  ave.air_temp+ave.wind_gust+Developed,data=meta.all.scaled)
summary(rda.all7)
RsquareAdj(rda.all7) # how much variation is explained by our model? 0.07076976
anova(rda.all7, by = "terms", permutations = how(nperm=999)) ### by variables
#ave.air_temp   1    471.3 1.4216  0.098 .
# ave.wind_gust  1    774.9 2.3375  0.009 **
#   Developed      1    430.0 1.2972  0.133
anova(rda.all7, by = NULL, permutations = how(nperm=999)) ### by model --> p = 0.005 **

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.all7)
#ave.air_temp ave.wind_gust     Developed
#1.115160      1.681600      1.704302

head(meta.all.scaled[,c(5,7,35)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.all.g1 = ordistep(rda(b.clr ~ 1, data = meta.all.scaled[,c(5,7,35)]),
                      scope=formula(rda.all7),
                      direction = "forward",
                      permutations = how(nperm=999))
# b.clr ~ ave.wind_gust + ave.air_temp
rda.all.g1$anova # see significance of individual terms in model
#                               Df    AIC      F Pr(>F)

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.all.g2 = ordiR2step(rda(b.clr ~ 1, data = meta.all.scaled[,c(5,7,35)]),
                        scope=formula(rda.all7),
                        permutations = how(nperm=999))
#nothing
rda.all.g2$anova # see significance of individual terms in model

#### RDA - WI ####

rownames(WI) %in% rownames(b.clr_WI) # check order of DFs
head(WI)

# included all Surface type frequencies and it overfit the model so pulling some out...
## dropping forest, open waters, others, and shrub first (after some experimenting)
rda.WI.0<-rda(b.clr_WI ~ ave.air_temp+ave.wind_speed+ave.wind_gust+ave.wind_direction+BarrenLand+CropLand+Developed+Forest+Herbaceous+Mexico+OpenWater+Others+SaltonSea+Shrub,data=WI)

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
# ave.air_temp     ave.wind_speed      ave.wind_gust ave.wind_direction         BarrenLand           CropLand          Developed
# 172.99090          837.17854          669.64301          254.24173           17.10031           32.75442                 NA
# Forest         Herbaceous             Mexico          OpenWater             Others          SaltonSea              Shrub
# NA                 NA                 NA                 NA                 NA                 NA                 NA

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(WI[,c(5:7,10,33:42)])
## we can use model selection instead of picking variables we think are important (by p values)
rda.WI.a = ordistep(rda(b.clr_WI ~ 1, data = WI[,c(5:7,10,33:42)]),
                         scope=formula(rda.WI.0),
                         direction = "forward",
                         permutations = how(nperm=999))
rda.WI.a$anova # see significance of individual terms in model
# + OpenWater           1 63.384 2.0387  0.059 .
# + ave.air_temp        1 63.882 1.5552  0.090 .

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.WI.a2 = ordiR2step(rda(b.clr_WI ~ 1, data = WI[,c(5:7,10,33:42)]),
                            scope=formula(rda.WI.0),
                            permutations = how(nperm=999))
# too many terms

# check best fit model based on above results
anova(rda.WI.a, permutations = how(nperm=999)) #not significant
#anova(rda.WI.a2, permutations = how(nperm=999)) # p =  0.001, significant

# Dropping Open Water. because it's in 1 sample, dropping Forest because neglible values across samples
rda.WI.1<-rda(b.clr_WI ~ ave.air_temp+ave.wind_speed+ave.wind_gust+ave.wind_direction+BarrenLand+CropLand+Developed+Herbaceous+Mexico+Others+SaltonSea+Shrub,data=WI)
summary(rda.WI.1)
RsquareAdj(rda.WI.1) # how much variation is explained by our model?
anova(rda.WI.1, by = "terms", permutations = how(nperm=999)) ### by variables
#                             Df Variance      F Pr(>F)

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.WI.1)
# ave.air_temp     ave.wind_speed      ave.wind_gust ave.wind_direction         BarrenLand           CropLand          Developed
# 172.99090          837.17854          669.64301          254.24173           17.10031           32.75442                 NA
# Herbaceous             Mexico             Others          SaltonSea              Shrub
# NA                 NA                 NA                 NA                 NA

head(WI[,c(5:7,10,33:42)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.WI.b1 = ordistep(rda(b.clr_WI ~ 1, data = WI[,c(5:7,10,33:35,37:38,40:42)]),
                          scope=formula(rda.WI.1),
                          direction = "forward",
                          permutations = how(nperm=999))
#                      Df    AIC      F Pr(>F)
#+ ave.air_temp        1 63.882 1.5552  0.083 .
rda.WI.b1$anova
#
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.WI.b2 = ordiR2step(rda(b.clr_WI ~ 1, data = WI[,c(5:7,10,33:35,37:38,40:42)]),
                            scope=formula(rda.WI.1),
                            permutations = how(nperm=999))
# too many variables

rda.WI.b2$anova
# check best fit model based on above results
anova(rda.WI.b1, permutations = how(nperm=999)) #

anova(rda.WI.0, rda.WI.1) #

# dropping Others
rda.WI.2<-rda(b.clr_WI ~ ave.air_temp+ave.wind_speed+ave.wind_gust+ave.wind_direction+BarrenLand+CropLand+Developed+Herbaceous+Mexico+SaltonSea+Shrub,data=WI)
summary(rda.WI.2)
RsquareAdj(rda.WI.2) # how much variation is explained by our model? NA
anova(rda.WI.2, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant

anova(rda.WI.2, by=NULL,permutations = how(nperm=999)) # p =  0.472

vif.cca(rda.WI.2)
#ave.air_temp     ave.wind_speed      ave.wind_gust ave.wind_direction         BarrenLand           CropLand          Developed
# 172.99090          837.17854          669.64301          254.24173           17.10031           32.75442                 NA
# Herbaceous             Mexico          SaltonSea              Shrub
# NA                 NA                 NA                 NA

# check if ORP & Sulfide are significantly correlated in August, which they are [strong, sig negative corr]
#cor.test(meta.all.scaled[metadata$SampDate=="WI",]$Sulfide_microM, meta.all.scaled[metadata$SampDate=="WI",]$ORP_mV, method="pearson") # ******

head(WI[,c(5:7,10,33:42)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.WI.c1 = ordistep(rda(b.clr_WI ~ 1, data = WI[,c(5:7,10,33:35,37:38,41:42)]),
                          scope=formula(rda.WI.2),
                          direction = "forward",
                          permutations = how(nperm=999))
# air temp near sig
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.WI.c2 = ordiR2step(rda(b.clr_WI ~ 1, data = WI[,c(5:7,10,33:35,37:38,41:42)]),
                            scope=formula(rda.WI.2),
                            permutations = how(nperm=999))
# too many terms

# check best fit model based on above results
anova(rda.WI.c1, permutations = how(nperm=999)) #

# dropped Developed because it's consistent across samples
rda.WI.3<-rda(b.clr_WI ~ ave.air_temp+ave.wind_speed+ave.wind_gust+ave.wind_direction+BarrenLand+CropLand+Herbaceous+Mexico+SaltonSea+Shrub,data=WI)
summary(rda.WI.3)
RsquareAdj(rda.WI.3) # how much variation is explained by our model? NA%
anova(rda.WI.3, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)

anova(rda.WI.3, by=NULL,permutations = how(nperm=999)) # p =  0.877

vif.cca(rda.WI.3)
# ave.air_temp     ave.wind_speed      ave.wind_gust ave.wind_direction         BarrenLand           CropLand         Herbaceous
# 172.99090          837.17854          669.64301          254.24173           17.10031           32.75442                 NA
# Mexico          SaltonSea              Shrub
# NA                 NA                 NA

head(WI[,c(5:7,10,33:42)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.WI.d1 = ordistep(rda(b.clr_WI ~ 1, data = WI[,c(5:7,10,33:34,37:38,41:42)]),
                          scope=formula(rda.WI.3),
                          direction = "forward",
                          permutations = how(nperm=999))
# ave air temp
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.WI.d2 = ordiR2step(rda(b.clr_WI ~ 1, data = WI[,c(5:7,10,33:34,37:38,41:42)]),
                            scope=formula(rda.WI.3),
                            permutations = how(nperm=999))
# too many vars

# drop CropLand because of last ordistep results (had highest p value)
rda.WI.4<-rda(b.clr_WI ~ ave.air_temp+ave.wind_speed+ave.wind_gust+ave.wind_direction+BarrenLand+Herbaceous+Mexico+SaltonSea+Shrub,data=WI)
summary(rda.WI.4)
RsquareAdj(rda.WI.4) # how much variation is explained by our model? NA
anova(rda.WI.4, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)
# ave air
anova(rda.WI.4, by=NULL,permutations = how(nperm=999)) # p =  0.909

vif.cca(rda.WI.4)
# ave.air_temp     ave.wind_speed      ave.wind_gust ave.wind_direction         BarrenLand         Herbaceous
# 13.28082          729.59764          580.61867           67.61179           59.38119           51.41826
# Mexico          SaltonSea              Shrub
# NA                 NA                 NA

head(WI[,c(5:7,10,33:42)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.WI.e1 = ordistep(rda(b.clr_WI ~ 1, data = WI[,c(5:7,10,33,37:38,41:42)]),
                     scope=formula(rda.WI.4),
                     direction = "forward",
                     permutations = how(nperm=999))
# ave air temp is near sig
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.WI.e2 = ordiR2step(rda(b.clr_WI ~ 1, data = WI[,c(5:7,10,33,37:38,41:42)]),
                       scope=formula(rda.WI.4),
                       permutations = how(nperm=999))
# too many vars

# dropped Mexico
rda.WI.5<-rda(b.clr_WI ~ ave.air_temp+ave.wind_speed+ave.wind_gust+ave.wind_direction+BarrenLand+Herbaceous+SaltonSea+Shrub,data=WI)
summary(rda.WI.5)
RsquareAdj(rda.WI.5) # how much variation is explained by our model? -0.5715959
anova(rda.WI.5, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)

anova(rda.WI.5, by=NULL,permutations = how(nperm=999)) # p

vif.cca(rda.WI.5)
# ave.air_temp     ave.wind_speed      ave.wind_gust ave.wind_direction         BarrenLand         Herbaceous
# 13.28082          729.59764          580.61867           67.61179           59.38119           51.41826
# SaltonSea              Shrub
# NA                 NA

head(WI[,c(5:7,10,33:42)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.WI.f1 = ordistep(rda(b.clr_WI ~ 1, data = WI[,c(5:7,10,33,37,41:42)]),
                     scope=formula(rda.WI.5),
                     direction = "forward",
                     permutations = how(nperm=999))
# ave air temp near sig
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.WI.f2 = ordiR2step(rda(b.clr_WI ~ 1, data = WI[,c(5:7,10,33,37,41:42)]),
                       scope=formula(rda.WI.5),
                       permutations = how(nperm=999))
# nothing sig

# dropped salton sea & wind direction
rda.WI.6<-rda(b.clr_WI ~ ave.air_temp+ave.wind_gust+Shrub,data=WI)
summary(rda.WI.6)
RsquareAdj(rda.WI.6) # how much variation is explained by our model? 0.2284503
anova(rda.WI.6, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)
# ave.air_temp   1   1883.6 1.8450  0.082 .

anova(rda.WI.6, by=NULL,permutations = how(nperm=999)) # p = 0.083 .

vif.cca(rda.WI.6)
#ave.air_temp ave.wind_gust         Shrub
#7.000397      5.243932     15.156672

head(WI[,c(5:7,10,33:42)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.WI.e1 = ordistep(rda(b.clr_WI ~ 1, data = WI[,c(5,7,42)]),
                     scope=formula(rda.WI.6),
                     direction = "forward",
                     permutations = how(nperm=999))
# ave air temp near sig
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.WI.e2 = ordiR2step(rda(b.clr_WI ~ 1, data = WI[,c(5,7,42)]),
                       scope=formula(rda.WI.6),
                       permutations = how(nperm=999))
# ave air temp near sig


rda.WI.7<-rda(b.clr_WI ~ ave.air_temp,data=WI)
summary(rda.WI.7)
RsquareAdj(rda.WI.7) # how much variation is explained by our model? 0.08469961
anova(rda.WI.7, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)
# ave.air_temp  1   1883.6 1.5552  0.088 .

anova(rda.WI.7, by=NULL,permutations = how(nperm=999)) # p = 0.089 .

head(WI[,c(5:7,10,33:42)])
## we can use model selection instead of picking variables we think are important -- based on p values
# rda.WI.f1 = ordistep(rda(b.clr_WI ~ 1, data = WI[,c(5)]),
#                      scope=formula(rda.WI.7),
#                      direction = "forward",
#                      permutations = how(nperm=999))
# # ave air temp near sig
# # Can also use model selection to pick variables by which ones increase variation (R^2)
# rda.WI.f2 = ordiR2step(rda(b.clr_WI ~ 1, data = WI[,c(5,7,42)]),
#                        scope=formula(rda.WI.7),
#                        permutations = how(nperm=999))
# # ave air temp near sig

#### RDA - DP ####

rownames(DP) %in% rownames(b.clr_DP) # check order of DFs
head(DP)

# included all Surface type frequencies and it overfit the model so pulling some out...
## dropping forest, open waters, others, and shrub first (after some experimenting)
rda.DP.0<-rda(b.clr_DP ~ ave.air_temp+ave.wind_speed+ave.wind_gust+ave.wind_direction+BarrenLand+CropLand+Developed+Forest+Herbaceous+Mexico+OpenWater+Others+SaltonSea+Shrub,data=DP)

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
# ave.air_temp     ave.wind_speed      ave.wind_gust ave.wind_direction         BarrenLand           CropLand          Developed
# 213.6287         10192.5296         17757.5822          9986.9425         20655.3788           462.7482                 NA
# Forest         Herbaceous             Mexico          OpenWater             Others          SaltonSea              Shrub
# NA                 NA                 NA                 NA                 NA                 NA                 NA

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(DP[,c(5:7,10,33:42)])
## we can use model selection instead of picking variables we think are important (by p values)
rda.DP.a = ordistep(rda(b.clr_DP ~ 1, data = DP[,c(5:7,10,33:42)]),
                    scope=formula(rda.DP.0),
                    direction = "forward",
                    permutations = how(nperm=999))
rda.DP.a$anova # see significance of individual terms in model
# Df    AIC      F Pr(>F)
# + BarrenLand  1 65.028 2.4284  0.002 **

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.DP.a2 = ordiR2step(rda(b.clr_DP ~ 1, data = DP[,c(5:7,10,33:42)]),
                       scope=formula(rda.DP.0),
                       permutations = how(nperm=999))
# too many terms

# check best fit model based on above results
anova(rda.DP.a, permutations = how(nperm=999)) #not significant
#anova(rda.DP.a2, permutations = how(nperm=999)) # p =  0.001, significant

# Dropping Open Water & Forest & Others because neglible values across samples
rda.DP.1<-rda(b.clr_DP ~ ave.air_temp+ave.wind_speed+ave.wind_gust+ave.wind_direction+BarrenLand+CropLand+Developed+Herbaceous+Mexico+SaltonSea+Shrub,data=DP)
summary(rda.DP.1)
RsquareAdj(rda.DP.1) # how much variation is explained by our model?
anova(rda.DP.1, by = "terms", permutations = how(nperm=999)) ### by variables
#                             Df Variance      F Pr(>F)

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.DP.1)
# ave.air_temp     ave.wind_speed      ave.wind_gust ave.wind_direction         BarrenLand           CropLand          Developed
# 213.6287         10192.5296         17757.5822          9986.9425         20655.3788           462.7482                 NA
# Herbaceous             Mexico             Others          SaltonSea              Shrub
# NA                 NA                 NA                 NA                 NA

head(DP[,c(5:7,10,33:42)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.DP.b1 = ordistep(rda(b.clr_DP ~ 1, data = DP[,c(5:7,10,33:35,37:38,41:42)]),
                     scope=formula(rda.DP.1),
                     direction = "forward",
                     permutations = how(nperm=999))
#                      Df    AIC      F Pr(>F)
#+ ave.air_temp        1 63.882 1.5552  0.083 .
rda.DP.b1$anova
#             Df    AIC      F Pr(>F)
# + BarrenLand  1 65.028 2.4284  0.001 ***

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.DP.b2 = ordiR2step(rda(b.clr_DP ~ 1, data = DP[,c(5:7,10,33:35,37:38,41:42)]),
                       scope=formula(rda.DP.1),
                       permutations = how(nperm=999))
# too many variables

rda.DP.b2$anova
# check best fit model based on above results
anova(rda.DP.b1, permutations = how(nperm=999)) # p =  0.001, significant

anova(rda.DP.0, rda.DP.1) #

# dropping Mexico
rda.DP.2<-rda(b.clr_DP ~ ave.air_temp+ave.wind_speed+ave.wind_gust+ave.wind_direction+BarrenLand+CropLand+Developed+Herbaceous+SaltonSea+Shrub,data=DP)
summary(rda.DP.2)
RsquareAdj(rda.DP.2) # how much variation is explained by our model? NA
anova(rda.DP.2, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant

anova(rda.DP.2, by=NULL,permutations = how(nperm=999)) # p =  0.472

vif.cca(rda.DP.2)
# ave.air_temp     ave.wind_speed      ave.wind_gust ave.wind_direction         BarrenLand           CropLand          Developed
# 213.6287         10192.5296         17757.5822          9986.9425         20655.3788           462.7482                 NA
# Herbaceous          SaltonSea              Shrub
# NA                 NA                 NA

# check if ORP & Sulfide are significantly correlated in August, which they are [strong, sig negative corr]
#cor.test(meta.all.scaled[metadata$SampDate=="DP",]$Sulfide_microM, meta.all.scaled[metadata$SampDate=="DP",]$ORP_mV, method="pearson") # ******

head(DP[,c(5:7,10,33:42)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.DP.c1 = ordistep(rda(b.clr_DP ~ 1, data = DP[,c(5:7,10,33:35,37,41:42)]),
                     scope=formula(rda.DP.2),
                     direction = "forward",
                     permutations = how(nperm=999))
# BarrenLand near sig
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.DP.c2 = ordiR2step(rda(b.clr_DP ~ 1, data = DP[,c(5:7,10,33:35,37,41:42)]),
                       scope=formula(rda.DP.2),
                       permutations = how(nperm=999))
# too many terms

# check best fit model based on above results
anova(rda.DP.c1, permutations = how(nperm=999)) #

# dropped Developed b/c looks quite consistent across DP samples, dropping Herbaceous too...
rda.DP.3<-rda(b.clr_DP ~ ave.air_temp+ave.wind_speed+ave.wind_gust+ave.wind_direction+BarrenLand+CropLand+SaltonSea+Shrub,data=DP)
summary(rda.DP.3)
RsquareAdj(rda.DP.3) # how much variation is explained by our model? NA%
anova(rda.DP.3, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)

anova(rda.DP.3, by=NULL,permutations = how(nperm=999)) # p =  0.877

vif.cca(rda.DP.3)
# ave.air_temp     ave.wind_speed      ave.wind_gust ave.wind_direction         BarrenLand           CropLand         Herbaceous
# 213.6287         10192.5296         17757.5822          9986.9425         20655.3788           462.7482                 NA
# SaltonSea              Shrub
# NA                 NA

head(DP[,c(5:7,10,33:42)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.DP.d1 = ordistep(rda(b.clr_DP ~ 1, data = DP[,c(5:7,10,33:34,41:42)]),
                     scope=formula(rda.DP.3),
                     direction = "forward",
                     permutations = how(nperm=999))
# BarrenLand
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.DP.d2 = ordiR2step(rda(b.clr_DP ~ 1, data = DP[,c(5:7,10,33:34,41:42)]),
                       scope=formula(rda.DP.3),
                       permutations = how(nperm=999))
# nothing sig

# drop CropLand
rda.DP.4<-rda(b.clr_DP ~ ave.air_temp+ave.wind_speed+ave.wind_gust+ave.wind_direction+BarrenLand+SaltonSea+Shrub,data=DP)
summary(rda.DP.4)
RsquareAdj(rda.DP.4) # how much variation is explained by our model? NA
anova(rda.DP.4, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)
# ave air
anova(rda.DP.4, by=NULL,permutations = how(nperm=999)) # p =  0.909

vif.cca(rda.DP.4)
# ave.air_temp     ave.wind_speed      ave.wind_gust ave.wind_direction         BarrenLand          SaltonSea              Shrub
# 1211.4525         29657.8747         38664.0469          8130.9335          7327.3150           169.8904                 NA

head(DP[,c(5:7,10,33:42)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.DP.e1 = ordistep(rda(b.clr_DP ~ 1, data = DP[,c(5:7,10,33,41:42)]),
                     scope=formula(rda.DP.4),
                     direction = "forward",
                     permutations = how(nperm=999))
# BarrenLand
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.DP.e2 = ordiR2step(rda(b.clr_DP ~ 1, data = DP[,c(5:7,10,33,41:42)]),
                       scope=formula(rda.DP.4),
                       permutations = how(nperm=999))
# nothing sig

# dropped ave.wind_gust
rda.DP.5<-rda(b.clr_DP ~ ave.air_temp+ave.wind_speed+ave.wind_direction+BarrenLand+Shrub+SaltonSea,data=DP)
summary(rda.DP.5)
RsquareAdj(rda.DP.5) # how much variation is explained by our model? -0.5715959
anova(rda.DP.5, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)

anova(rda.DP.5, by=NULL,permutations = how(nperm=999)) # p

vif.cca(rda.DP.5)
#ave.air_temp     ave.wind_speed ave.wind_direction         BarrenLand              Shrub          SaltonSea
#535.9567           101.2455           843.4096          6776.6750          1715.7142           249.6605

head(DP[,c(5:7,10,33:42)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.DP.f1 = ordistep(rda(b.clr_DP ~ 1, data = DP[,c(5:6,10,33,41:42)]),
                     scope=formula(rda.DP.5),
                     direction = "forward",
                     permutations = how(nperm=999))
# BarrenLand
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.DP.f2 = ordiR2step(rda(b.clr_DP ~ 1, data = DP[,c(5:6,10,33,41:42)]),
                       scope=formula(rda.DP.5),
                       permutations = how(nperm=999))
# nothing sig

# dropped ave.wind speed because it seems to be the least important in the model assembly with ordistep
rda.DP.6<-rda(b.clr_DP ~ ave.air_temp+ave.wind_direction+BarrenLand+Shrub+SaltonSea,data=DP)
summary(rda.DP.6)
RsquareAdj(rda.DP.6) # how much variation is explained by our model? 0.2284503
anova(rda.DP.6, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)
# ave.air_temp   1   1883.6 1.8450  0.082 .

anova(rda.DP.6, by=NULL,permutations = how(nperm=999))

vif.cca(rda.DP.6)
#ave.air_temp ave.wind_direction         BarrenLand              Shrub          SaltonSea
# 47.93774           22.45774          272.13958          315.93436           86.55458

head(DP[,c(5:7,10,33:42)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.DP.e1 = ordistep(rda(b.clr_DP ~ 1, data = DP[,c(5,10,33,41:42)]),
                     scope=formula(rda.DP.6),
                     direction = "forward",
                     permutations = how(nperm=999))
# BarrenLand
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.DP.e2 = ordiR2step(rda(b.clr_DP ~ 1, data = DP[,c(5,10,33,41:42)]),
                       scope=formula(rda.DP.6),
                       permutations = how(nperm=999))

# dropping Salton Sea + SHrub
rda.DP.7<-rda(b.clr_DP ~ ave.air_temp+ave.wind_direction+BarrenLand,data=DP)
summary(rda.DP.7)
RsquareAdj(rda.DP.7) # how much variation is explained by our model? 0.1565372
anova(rda.DP.7, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)
# ave.air_temp   1   3209.4 2.1546  0.045 *

anova(rda.DP.7, by=NULL,permutations = how(nperm=999)) # p = 0.089 .

vif.cca(rda.DP.7)
# ave.air_temp ave.wind_direction         BarrenLand
#21.71411           14.11877           25.57201

head(DP[,c(5:7,10,33:42)])
# we can use model selection instead of picking variables we think are important -- based on p values
rda.DP.f1 = ordistep(rda(b.clr_DP ~ 1, data = DP[,c(5,10,33)]),
                     scope=formula(rda.DP.7),
                     direction = "forward",
                     permutations = how(nperm=999))
# BarrenLand
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.DP.f2 = ordiR2step(rda(b.clr_DP ~ 1, data = DP[,c(5,10,33)]),
                       scope=formula(rda.DP.7),
                       permutations = how(nperm=999))

#dropping ave. air temp
rda.DP.8<-rda(b.clr_DP ~ ave.wind_direction+BarrenLand,data=DP)
summary(rda.DP.8)
RsquareAdj(rda.DP.8) # how much variation is explained by our model? 0.2016672
anova(rda.DP.8, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)
# ave.air_temp   1   3497.7 2.4809  0.015 *

anova(rda.DP.8, by=NULL,permutations = how(nperm=999)) # p = 0.045 *

vif.cca(rda.DP.8)
# ave.wind_direction         BarrenLand
# 13.03326           13.03326

head(DP[,c(5:7,10,33:42)])
# we can use model selection instead of picking variables we think are important -- based on p values
rda.DP.g1 = ordistep(rda(b.clr_DP ~ 1, data = DP[,c(10,33)]),
                     scope=formula(rda.DP.8),
                     direction = "forward",
                     permutations = how(nperm=999))
# BarrenLand
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.DP.g2 = ordiR2step(rda(b.clr_DP ~ 1, data = DP[,c(10,33)]),
                       scope=formula(rda.DP.8),
                       permutations = how(nperm=999))
# ave. wind_direction

# Final DP RDA: rda.DP.8

#### RDA - BDC ####

rownames(BDC) %in% rownames(b.clr_BDC) # check order of DFs
head(BDC)

# included all Surface type frequencies and it overfit the model so pulling some out...
## dropping forest, open waters, others, and shrub first (after some experimenting)
rda.BDC.0<-rda(b.clr_BDC ~ ave.air_temp+ave.wind_speed+ave.wind_gust+ave.wind_direction+BarrenLand+CropLand+Developed+Forest+Herbaceous+Mexico+OpenWater+Others+SaltonSea+Shrub,data=BDC)

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
# ave.air_temp     ave.wind_speed      ave.wind_gust ave.wind_direction         BarrenLand           CropLand
# 3.217787          11.231785          13.732170          11.154667           4.537670           8.753092
# Developed             Forest         Herbaceous             Mexico          OpenWater             Others
# NA                 NA                 NA                 NA                 NA                 NA
# SaltonSea              Shrub
# NA                 NA

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(BDC[,c(5:7,10,33:42)])
## we can use model selection instead of picking variables we think are important (by p values)
rda.BDC.a = ordistep(rda(b.clr_BDC ~ 1, data = BDC[,c(5:7,10,33:42)]),
                    scope=formula(rda.BDC.0),
                    direction = "forward",
                    permutations = how(nperm=999))
rda.BDC.a$anova # see significance of individual terms in model
# Df    AIC      F Pr(>F)
# + BarrenLand  1 62.539 1.7548  0.049 *

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.BDC.a2 = ordiR2step(rda(b.clr_BDC ~ 1, data = BDC[,c(5:7,10,33:42)]),
                       scope=formula(rda.BDC.0),
                       permutations = how(nperm=999))
# too many terms

# check best fit model based on above results
anova(rda.BDC.a, permutations = how(nperm=999)) #not significant
#anova(rda.BDC.a2, permutations = how(nperm=999)) # p =  0.001, significant

# Dropping SaltonSea and Mexico and Others
rda.BDC.1<-rda(b.clr_BDC ~ ave.air_temp+ave.wind_speed+ave.wind_gust+ave.wind_direction+BarrenLand+CropLand+Developed+Herbaceous+Forest+OpenWater+Shrub,data=BDC)
summary(rda.BDC.1)
RsquareAdj(rda.BDC.1) # how much variation is explained by our model?
anova(rda.BDC.1, by = "terms", permutations = how(nperm=999)) ### by variables
#                             Df Variance      F Pr(>F)

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.BDC.1)
#  ave.air_temp     ave.wind_speed      ave.wind_gust ave.wind_direction         BarrenLand           CropLand
# 3.217787          11.231785          13.732170          11.154667           4.537670           8.753092
# Developed         Herbaceous             Forest            OpenWater              Shrub
# NA                 NA                 NA                      NA                 NA

head(BDC[,c(5:7,10,33:42)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.BDC.b1 = ordistep(rda(b.clr_BDC ~ 1, data = BDC[,c(5:7,10,33:37,39,42)]),
                     scope=formula(rda.BDC.1),
                     direction = "forward",
                     permutations = how(nperm=999))
#                      Df    AIC      F Pr(>F)
#+ BarrenLand          1 62.539 1.7548  0.056 .
# + Developed           1 62.645 1.6538  0.082 .
# + ave.air_temp        1 62.680 1.6207  0.088 .
# + Forest              1 62.785 1.5223  0.090 .

rda.BDC.b1$anova


# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.BDC.b2 = ordiR2step(rda(b.clr_BDC ~ 1, data = BDC[,c(5:7,10,33:37,39,42)]),
                       scope=formula(rda.BDC.1),
                       permutations = how(nperm=999))
# too many variables

rda.BDC.b2$anova
# check best fit model based on above results
anova(rda.BDC.b1, permutations = how(nperm=999)) # p =  0.001, significant

anova(rda.BDC.0, rda.BDC.1) #

# dropping OpenWater
rda.BDC.2<-rda(b.clr_BDC ~ ave.air_temp+ave.wind_speed+ave.wind_gust+ave.wind_direction+BarrenLand+CropLand+Developed+Herbaceous+Forest+Shrub,data=BDC)
summary(rda.BDC.2)
RsquareAdj(rda.BDC.2) # how much variation is explained by our model? NA
anova(rda.BDC.2, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant

anova(rda.BDC.2, by=NULL,permutations = how(nperm=999)) # p =  0.472

vif.cca(rda.BDC.2)
# ave.air_temp     ave.wind_speed      ave.wind_gust ave.wind_direction         BarrenLand           CropLand
# 3.217787          11.231785          13.732170          11.154667           4.537670           8.753092
# Developed         Herbaceous             Forest              Shrub
# NA                 NA                 NA                 NA

# check if ORP & Sulfide are significantly correlated in August, which they are [strong, sig negative corr]
#cor.test(meta.all.scaled[metadata$SampDate=="BDC",]$Sulfide_microM, meta.all.scaled[metadata$SampDate=="BDC",]$ORP_mV, method="pearson") # ******

head(BDC[,c(5:7,10,33:42)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.BDC.c1 = ordistep(rda(b.clr_BDC ~ 1, data = BDC[,c(5:7,10,33:37,42)]),
                     scope=formula(rda.BDC.2),
                     direction = "forward",
                     permutations = how(nperm=999))
# + BarrenLand          1 62.539 1.7548  0.053 .
# + ave.air_temp        1 62.680 1.6207  0.073 .
# + Developed           1 62.645 1.6538  0.085 .
# + Forest              1 62.785 1.5223  0.086 .
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.BDC.c2 = ordiR2step(rda(b.clr_BDC ~ 1, data = BDC[,c(5:7,10,33:37,42)]),
                       scope=formula(rda.BDC.2),
                       permutations = how(nperm=999))
# too many terms

# check best fit model based on above results
anova(rda.BDC.c1, permutations = how(nperm=999)) #

# dropped Shrub b/c looks quite consistent across BDC samples
rda.BDC.3<-rda(b.clr_BDC ~ ave.air_temp+ave.wind_speed+ave.wind_gust+ave.wind_direction+BarrenLand+CropLand+Developed+Herbaceous+Forest,data=BDC)
summary(rda.BDC.3)
RsquareAdj(rda.BDC.3) # how much variation is explained by our model? NA%
anova(rda.BDC.3, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)

anova(rda.BDC.3, by=NULL,permutations = how(nperm=999))

vif.cca(rda.BDC.3)
#ave.air_temp     ave.wind_speed      ave.wind_gust ave.wind_direction         BarrenLand           CropLand
# 3.217787          11.231785          13.732170          11.154667           4.537670           8.753092
# Developed         Herbaceous             Forest
# NA                 NA                 NA

head(BDC[,c(5:7,10,33:42)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.BDC.d1 = ordistep(rda(b.clr_BDC ~ 1, data = BDC[,c(5:7,10,33:37)]),
                     scope=formula(rda.BDC.3),
                     direction = "forward",
                     permutations = how(nperm=999))
#                     Df    AIC      F Pr(>F)
# + BarrenLand          1 62.539 1.7548  0.058 .
# + Developed           1 62.645 1.6538  0.079 .
# + ave.air_temp        1 62.680 1.6207  0.082 .

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.BDC.d2 = ordiR2step(rda(b.clr_BDC ~ 1, data = BDC[,c(5:7,10,33:37)]),
                       scope=formula(rda.BDC.3),
                       permutations = how(nperm=999))
# nothing sig

# drop CropLand & Herbaceous - also quite consistent across samples
rda.BDC.4<-rda(b.clr_BDC ~ ave.air_temp+ave.wind_speed+ave.wind_gust+ave.wind_direction+BarrenLand+Developed+Forest,data=BDC)
summary(rda.BDC.4)
RsquareAdj(rda.BDC.4) # how much variation is explained by our model? NA
anova(rda.BDC.4, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)
# ave air
anova(rda.BDC.4, by=NULL,permutations = how(nperm=999)) # p =  0.909

vif.cca(rda.BDC.4)
# ave.air_temp     ave.wind_speed      ave.wind_gust ave.wind_direction         BarrenLand          Developed
# 4.685982          13.989703          16.439060           7.179747           4.810692           5.232373
# Forest
# NA

head(BDC[,c(5:7,10,33:42)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.BDC.e1 = ordistep(rda(b.clr_BDC ~ 1, data = BDC[,c(5:7,10,33,35:36)]),
                     scope=formula(rda.BDC.4),
                     direction = "forward",
                     permutations = how(nperm=999))
#                      Df    AIC      F Pr(>F)
# + BarrenLand          1 62.539 1.7548  0.059 .
# + Developed           1 62.645 1.6538  0.067 .
# + ave.air_temp        1 62.680 1.6207  0.078 .
# + Forest              1 62.785 1.5223  0.100 .

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.BDC.e2 = ordiR2step(rda(b.clr_BDC ~ 1, data = BDC[,c(5:7,10,33,35:36)]),
                       scope=formula(rda.BDC.4),
                       permutations = how(nperm=999))
# nothing sig, too many vars

# dropped ave.wind gust, has highest VIF
rda.BDC.5<-rda(b.clr_BDC ~ ave.air_temp+ave.wind_speed+ave.wind_direction+BarrenLand+Developed+Forest,data=BDC)
summary(rda.BDC.5)
RsquareAdj(rda.BDC.5) # how much variation is explained by our model?
anova(rda.BDC.5, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)

anova(rda.BDC.5, by=NULL,permutations = how(nperm=999)) # p = 0.38

vif.cca(rda.BDC.5)
# ave.air_temp     ave.wind_speed ave.wind_direction         BarrenLand          Developed             Forest
# 227.374922           7.814828         783.093416         765.128298         347.296027        3327.342563

head(BDC[,c(5:7,10,33:42)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.BDC.f1 = ordistep(rda(b.clr_BDC ~ 1, data = BDC[,c(5:6,10,33,35:36)]),
                     scope=formula(rda.BDC.5),
                     direction = "forward",
                     permutations = how(nperm=999))
#                      Df    AIC      F Pr(>F)
# + BarrenLand          1 62.539 1.7548  0.058 .
# + Developed           1 62.645 1.6538  0.066 .
# + ave.air_temp        1 62.680 1.6207  0.085 .
# + Forest              1 62.785 1.5223  0.099 .

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.BDC.f2 = ordiR2step(rda(b.clr_BDC ~ 1, data = BDC[,c(5:6,10,33,35:36)]),
                       scope=formula(rda.BDC.5),
                       permutations = how(nperm=999))
# too many vars

# dropped ave.wind direction
rda.BDC.6<-rda(b.clr_BDC ~ ave.air_temp+ave.wind_speed+BarrenLand+Developed+Forest,data=BDC)
summary(rda.BDC.6)
RsquareAdj(rda.BDC.6) # how much variation is explained by our model? 0.06367732
anova(rda.BDC.6, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)
# nothing

anova(rda.BDC.6, by=NULL,permutations = how(nperm=999))

vif.cca(rda.BDC.6)
#ave.air_temp ave.wind_speed     BarrenLand      Developed         Forest
#6.384892       4.985955      16.331225       2.800031      28.808080

head(BDC[,c(5:7,10,33:42)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.BDC.e1 = ordistep(rda(b.clr_BDC ~ 1, data = BDC[,c(5:6,33,35:36)]),
                     scope=formula(rda.BDC.6),
                     direction = "forward",
                     permutations = how(nperm=999))
#                  Df    AIC      F Pr(>F)
# + BarrenLand      1 62.539 1.7548  0.054 .
# + Developed       1 62.645 1.6538  0.061 .
# + ave.air_temp    1 62.680 1.6207  0.087 .

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.BDC.e2 = ordiR2step(rda(b.clr_BDC ~ 1, data = BDC[,c(5:6,33,35:36)]),
                       scope=formula(rda.BDC.6),
                       permutations = how(nperm=999))

# dropping ave.wind_speed
rda.BDC.7<-rda(b.clr_BDC ~ ave.air_temp+BarrenLand+Developed+Forest,data=BDC)
summary(rda.BDC.7)
RsquareAdj(rda.BDC.7) # how much variation is explained by our model? 0.1694793
anova(rda.BDC.7, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)
#ave.air_temp  1   1653.1 1.7685  0.080 .

anova(rda.BDC.7, by=NULL,permutations = how(nperm=999)) # p = 0.196

vif.cca(rda.BDC.7)
# ave.air_temp   BarrenLand    Developed       Forest
# 2.844049     5.871025     2.539155     7.509786

head(BDC[,c(5:7,10,33:42)])
# we can use model selection instead of picking variables we think are important -- based on p values
rda.BDC.f1 = ordistep(rda(b.clr_BDC ~ 1, data = BDC[,c(5,33,35:36)]),
                     scope=formula(rda.BDC.7),
                     direction = "forward",
                     permutations = how(nperm=999))
# + BarrenLand    1 62.539 1.7548  0.066 .
# + Developed     1 62.645 1.6538  0.071 .
# + ave.air_temp  1 62.680 1.6207  0.097 .

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.BDC.f2 = ordiR2step(rda(b.clr_BDC ~ 1, data = BDC[,c(5,33,35:36)]),
                       scope=formula(rda.BDC.7),
                       permutations = how(nperm=999))
#+ BarrenLand  1 62.539 1.7548  0.062 .

#dropping Forest
rda.BDC.8<-rda(b.clr_BDC ~ ave.air_temp+BarrenLand+Developed,data=BDC)
summary(rda.BDC.8)
RsquareAdj(rda.BDC.8) # how much variation is explained by our model? 0.1956148
anova(rda.BDC.8, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)
# ave.air_temp   1   1653.1 1.8259  0.069 .

anova(rda.BDC.8, by=NULL,permutations = how(nperm=999)) # p = 0.056 .

vif.cca(rda.BDC.8)
# ave.air_temp   BarrenLand    Developed
# 1.802045     2.120474     1.621290

head(BDC[,c(5:7,10,33:42)])
# we can use model selection instead of picking variables we think are important -- based on p values
rda.BDC.g1 = ordistep(rda(b.clr_BDC ~ 1, data = BDC[,c(5,33,35)]),
                     scope=formula(rda.BDC.8),
                     direction = "forward",
                     permutations = how(nperm=999))
# + BarrenLand    1 62.539 1.7548  0.054 .
# + Developed     1 62.645 1.6538  0.063 .
# + ave.air_temp  1 62.680 1.6207  0.083 .

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.BDC.g2 = ordiR2step(rda(b.clr_BDC ~ 1, data = BDC[,c(5,33,35)]),
                       scope=formula(rda.BDC.8),
                       permutations = how(nperm=999))
# + BarrenLand  1 62.539 1.7548  0.051 .

# Dropped Developed
rda.BDC.9<-rda(b.clr_BDC ~ ave.air_temp+BarrenLand,data=BDC)
summary(rda.BDC.9)
RsquareAdj(rda.BDC.9) # how much variation is explained by our model? 0.1591731
anova(rda.BDC.9, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)
# ave.air_temp   1 1653.1 1.7468  0.067 .

anova(rda.BDC.9, by=NULL,permutations = how(nperm=999)) # p = 0.072 .

vif.cca(rda.BDC.9)
# ave.air_temp   BarrenLand
# 1.350152     1.350152

head(BDC[,c(5:7,10,33:42)])
# we can use model selection instead of picking variables we think are important -- based on p values
rda.BDC.h1 = ordistep(rda(b.clr_BDC ~ 1, data = BDC[,c(5,33)]),
                      scope=formula(rda.BDC.9),
                      direction = "forward",
                      permutations = how(nperm=999))
# + BarrenLand    1 62.539 1.7548  0.063 .
# + ave.air_temp  1 62.680 1.6207  0.084 .

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.BDC.h2 = ordiR2step(rda(b.clr_BDC ~ 1, data = BDC[,c(5,33)]),
                        scope=formula(rda.BDC.9),
                        permutations = how(nperm=999))
# + BarrenLand  1 62.539 1.7548  0.062 .

# Final RDA for BDC: rda.BDC.9

#### Final RDAs ####
# RDA by sampling timepoint
head(meta.all.scaled)
head(b.clr)
rownames(b.clr) %in% rownames(meta.all.scaled) # sanity check 1

# all data
#rda.all2$call # best model for all data

rda.all<-rda(b.clr ~ Temp_DegC + Dissolved_OrganicMatter_RFU + DO_Percent_Local,data=meta.all.scaled)
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
rda.all.part<-varpart(b.clr, meta.all.scaled$Temp_DegC, meta.all.scaled$Dissolved_OrganicMatter_RFU,meta.all.scaled$DO_Percent_Local)
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
rownames(rda.sum.all$sites) %in% rownames(meta.all.scaled)
rda.axes.all<-data.frame(RDA1=rda.sum.all$sites[,1], RDA2=rda.sum.all$sites[,2], SampleID=rownames(rda.sum.all$sites), Depth_m=meta.all.scaled$Depth_m, SampDate=meta.all.scaled$SampDate)

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
