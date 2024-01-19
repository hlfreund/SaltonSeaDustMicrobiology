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
load("data/SSeawater_Data_Ready.Rdata") # save global env to Rdata file
#load("data/SSW_Amplicon_EnvDriver.Rdata")
load("data/SSW_Amplicon_EnvDriver_RDAsOnly.Rdata")

bac.dat.all[1:4,1:4]
bac.ASV_table[,1:4]
bac.ASV_table[(nrow(bac.ASV_table)-4):(nrow(bac.ASV_table)),(ncol(bac.ASV_table)-4):(ncol(bac.ASV_table))] # last 4 rows & cols

head(metadata)
head(meta_scaled)

# drop salinity from meta_scaled --> excluding this env variable
meta_scaled<-subset(meta_scaled, select=-c(Salinity_ppt))
head(meta_scaled)

# create column for Depth that is a numeric version of this variable, rather than a factor
meta_scaled$Depth.num<-as.numeric(as.character(meta_scaled$Depth_m))

#### Create Centered Log-Ratio Table from ASV table ####
bac.ASV_table[1:4,1:4]
b.clr<-decostand(bac.ASV_table[,-1],method = "clr", pseudocount = 1) #CLR transformation
b.clr[1:4,1:4]
# df must have rownames are SampleIDs, columns are ASV IDs for vegan functions below\

#### Separate All Data by Timepoints ####
# create metadata df that will contain scaled chemical data
head(metadata)
head(meta_scaled)

site_list<-unique(meta_scaled$SampDate) #define an array of string values
# go through metadata & create a list of data frames
## when metadata$Variable == element in site_list (aka x in this case), subset metadata by said element into elements of a list

# here the function(x) is using site_list aka x to subset metadata, when $Variable column == site_list
# Run the function so it's stored in Global Env
site_subsets<-lapply(site_list, function(x) {subset(meta_scaled, SampDate==x)})

site_subsets # sanity check1 (should see all elements in list)
site_subsets[[1]] # sanity check2 (see 1st element in list)
#rename the list elements

# name each element in list
names(site_subsets)<-site_list # * only do this if the order of names in site_list match order of the elements in site_subsets!
site_subsets$April.2022 # sanity check3 - should be able to pull dataframes by names rather than index now

# example of subsetting
site_subsets[[2]][1:3]
site_subsets$August.2021[1:3] # should produce same ouptut as line above

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

head(August.2021) # sanity check
August.2021[1:5,] # double check that our new Variable (here SampDate) data frames still have scaled chemical data
rownames(August.2021)

# matching data with user defined function -- here is the function, must run to store function in Global env
match_dat<-function(compdata, subset_metadata){
  subset_comp_data = pullrow<-(is.element(row.names(compdata), row.names(subset_metadata)))
  ### * comp data and metadata need to have row names - rownames should be Sample IDs
  subset_comp_data=compdata[pullrow,]
  return(subset_comp_data)
}

# double check that our data frames are ready for this function, aka that they both have the same rownames
## row #s do not have to be the same, but their row names should be in the same format and be able to match up
rownames(b.clr)
rownames(August.2021)

# run the function
b.clr_AUG21<-match_dat(b.clr,August.2021)
b.clr_DEC21<-match_dat(b.clr,December.2021)
b.clr_APR22<-match_dat(b.clr,April.2022)

# did the function work the way we wanted it to?

b.clr_AUG21[1:4,1:4]
rownames(August.2021) %in% rownames(b.clr_AUG21) # hopefully all of the rownames match, aka will get output of TRUE

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
b.dca #DCA1 axis length = 0.32591; use RDA
## The length of first DCA axis:
## > 4 indicates heterogeneous dataset on which unimodal methods should be used (CCA),
##  < 3 indicates homogeneous dataset for which linear methods are suitable (RDA)
## between 3 and 4 both linear and unimodal methods are OK.

# BY MONTH

b.clr_A21.pseudo<-b.clr_AUG21+1
b.A21.dca = decorana(b.clr_A21.pseudo)
b.A21.dca #DCA1 axis length = 0.211195; use RDA

b.clr_D21.pseudo<-b.clr_DEC21+1
b.D21.dca = decorana(b.clr_D21.pseudo)
b.D21.dca #DCA1 axis length = 0.229916; use RDA

b.clr_A22.pseudo<-b.clr_APR22+1
b.A22.dca = decorana(b.clr_A22.pseudo)
b.A22.dca #DCA1 axis length = 0.200988; use RDA

#### RDA w/ All Data ####

rownames(meta_scaled) %in% rownames(b.clr) # check order of DFs
head(meta_scaled)

rda.all.0<-rda(b.clr ~ DO_Percent_Local+ORP_mV+Temp_DegC+Dissolved_OrganicMatter_RFU+Depth.num+Sulfate_milliM+Sulfide_microM,data=meta_scaled)

# check summary of RDA
rda.all.0
summary(rda.all.0)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.all.0) # 51.74%
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.all.0, permutations = how(nperm=999)) # p = 0.001, significant

## we can also do a permutation test by RDA axis
#anova(rda.all.0, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.all.0, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                           Df Variance      F Pr(>F)
#DO_Percent_Local             1   261.23 11.2987  0.001 ***
#ORP_mV                       1    45.22  1.9558  0.060 .
#Temp_DegC                    1   100.51  4.3471  0.001 ***
#Dissolved_OrganicMatter_RFU  1   230.47  9.9680  0.001 ***

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.all.0)
# DO_Percent_Local               ORP_mV                   Temp_DegC    Dissolved_OrganicMatter_RFU                   Depth.num
# 11.560558                   28.792805                    7.844239                    3.475087                    3.026465
# Sulfate_milliM              Sulfide_microM
# 2.755198                   32.189492

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/

head(meta_scaled)
## we can use model selection instead of picking variables we think are important (by p values)
# more info on ordistep & ordiR2step here: https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
rda.all.a = ordistep(rda(b.clr ~ 1, data = meta_scaled[,c(8,10:11,14:16,18)]),
                     scope=formula(rda.all.0),
                     direction = "forward",
                     permutations = how(nperm=999))
# b.clr ~ Temp_DegC + Dissolved_OrganicMatter_RFU + DO_Percent_Local + Depth.num + Sulfate_milliM  = best model
rda.all.a$anova # see significance of individual terms in model
#                               Df    AIC      F Pr(>F)
#+ Temp_DegC                    1 162.38 9.6317  0.001 ***
# + Dissolved_OrganicMatter_RFU  1 158.36 5.9794  0.001 ***
#   + DO_Percent_Local             1 155.31 4.6845  0.001 ***
#   + Depth.num                    1 155.37 1.5993  0.013 *
#   + Sulfate_milliM               1 155.49 1.4705  0.027 *

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.all.a2 = ordiR2step(rda(b.clr ~ 1, data = meta_scaled[,c(8,10:11,14:16,18)]),
                        scope=formula(rda.all.0),
                        permutations = how(nperm=999))
# b.clr ~ Temp_DegC + Dissolved_OrganicMatter_RFU + DO_Percent_Local + Depth.num  = best model
rda.all.a2$anova # see significance of individual terms in model
#                               R2.adj Df    AIC      F Pr(>F)
# + Temp_DegC                   0.27288  1 162.38 9.6317  0.001 ***
# + Dissolved_OrganicMatter_RFU 0.40708  1 158.36 5.9794  0.001 ***
# + DO_Percent_Local            0.49558  1 155.31 4.6845  0.001 ***
# + Depth.num                   0.51026  1 155.37 1.5993  0.019 *

# check best fit model based on above results
anova(rda.all.a, permutations = how(nperm=999)) # p =  0.001, significant

# Let's double check by removing the variables with high VIF
rda.all1<-rda(b.clr ~ Temp_DegC + Dissolved_OrganicMatter_RFU + DO_Percent_Local + Depth.num + Sulfate_milliM,data=meta_scaled)
summary(rda.all1)
RsquareAdj(rda.all1) # how much variation is explained by our model? 52.21%
anova(rda.all1, by = "terms", permutations = how(nperm=999)) ### by variables
# Temp_DegC                    1   335.51 14.6543  0.001 ***
#   Dissolved_OrganicMatter_RFU  1   169.84  7.4183  0.001 ***
#   DO_Percent_Local             1   113.20  4.9444  0.001 ***

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.all1)
# Temp_DegC Dissolved_OrganicMatter_RFU            DO_Percent_Local                   Depth.num              Sulfate_milliM
# 7.187248                    3.237747                    9.958522                    2.771739                    2.491430
head(meta_scaled)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.all.b1 = ordistep(rda(b.clr ~ 1, data = meta_scaled[,c(8,11,14:15,18)]),
                      scope=formula(rda.all1),
                      direction = "forward",
                      permutations = how(nperm=999))
# b.clr ~ Temp_DegC + Dissolved_OrganicMatter_RFU + DO_Percent_Local +      Depth.num + Sulfate_milliM  = best model
rda.all.b1$anova # see significance of individual terms in model
#                               Df    AIC      F Pr(>F)
#+ Temp_DegC                    1 162.38 9.6317  0.001 ***
#+ Dissolved_OrganicMatter_RFU  1 158.36 5.9794  0.001 ***
  # + DO_Percent_Local             1 155.31 4.6845  0.001 ***
  # + Depth.num                    1 155.37 1.5993  0.016 *
  # + Sulfate_milliM               1 155.49 1.4705  0.037 *

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.all.b2 = ordiR2step(rda(b.clr ~ 1, data = meta_scaled[,c(8,11,14:15,18)]),
                        scope=formula(rda.all1),
                        permutations = how(nperm=999))
# b.clr ~ Dissolved_OrganicMatter_RFU + Sulfate_milliM + Temp_DegC + Depth.num + DO_%Local = best model
rda.all.b2$anova # see significance of individual terms in model
#                               R2.adj Df    AIC      F Pr(>F)
# + Temp_DegC                   0.27288  1 162.38 9.6317  0.001 ***
#   + Dissolved_OrganicMatter_RFU 0.40708  1 158.36 5.9794  0.001 ***
#   + DO_Percent_Local            0.49558  1 155.31 4.6845  0.001 ***
#   + Depth.num                   0.51026  1 155.37 1.5993  0.021 *
#   + Sulfate_milliM              0.52209  1 155.49 1.4705  0.037 *

# check best fit model based on above results
anova(rda.all.b1, permutations = how(nperm=999)) # p =  0.001, significant

# compare model fits to each other
anova(rda.all.0, rda.all.b1)

rda.all2<-rda(b.clr ~ Temp_DegC + Dissolved_OrganicMatter_RFU + DO_Percent_Local,data=meta_scaled)
summary(rda.all2)
RsquareAdj(rda.all2) # how much variation is explained by our model? 49.55%
anova(rda.all2, by = "terms", permutations = how(nperm=999)) ### by variables
# Temp_DegC                    1   335.51 13.8841  0.001 ***
# Dissolved_OrganicMatter_RFU  1   169.84  7.0285  0.001 ***
# DO_Percent_Local             1   113.20  4.6845  0.001 ***

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.all2)
# Temp_DegC Dissolved_OrganicMatter_RFU            DO_Percent_Local
# 2.957894                    2.122611                    3.906739
head(meta_scaled)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.all.c1 = ordistep(rda(b.clr ~ 1, data = meta_scaled[,c(8,11,14)]),
                      scope=formula(rda.all2),
                      direction = "forward",
                      permutations = how(nperm=999))
# b.clr ~ Temp_DegC + Dissolved_OrganicMatter_RFU + DO_Percent_Local +      Depth.num + Sulfate_milliM  = best model
rda.all.c1$anova # see significance of individual terms in model
#                               Df    AIC      F Pr(>F)

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.all.c2 = ordiR2step(rda(b.clr ~ 1, data = meta_scaled[,c(8,11,14)]),
                        scope=formula(rda.all2),
                        permutations = how(nperm=999))
# b.clr ~ Dissolved_OrganicMatter_RFU + Temp_DegC  + DO_%Local = best model
rda.all.c2$anova # see significance of individual terms in model
#                               R2.adj Df    AIC      F Pr(>F)
# + Temp_DegC                   0.27288  1 162.38 9.6317  0.001 ***
#   + Dissolved_OrganicMatter_RFU 0.40708  1 158.36 5.9794  0.001 ***
#   + DO_Percent_Local            0.49558  1 155.31 4.6845  0.001 ***
#   + Depth.num                   0.51026  1 155.37 1.5993  0.021 *
#   + Sulfate_milliM              0.52209  1 155.49 1.4705  0.037 *

#### RDA - August 2021 ####

rownames(August.2021) %in% rownames(b.clr_AUG21) # check order of DFs
head(August.2021)

rda.aug2021.0<-rda(b.clr_AUG21 ~ DO_Percent_Local+ORP_mV+Dissolved_OrganicMatter_RFU+Sulfate_milliM+Sulfide_microM+Depth.num,data=August.2021)

# check summary of RDA
rda.aug2021.0
summary(rda.aug2021.0)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.aug2021.0) # 9.47%
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.aug2021.0, permutations = how(nperm=999))

## we can also do a permutation test by RDA axis
#anova(rda.aug2021.0, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.aug2021.0, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                               Df Variance      F Pr(>F)
# ORP_mV                       1  146.426 1.6944  0.025 *

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.aug2021.0)
#  DO_Percent_Local           ORP_mV        Dissolved_OrganicMatter_RFU              Sulfate_milliM              Sulfide_microM
# 51.647076                   39.906439                  152.528418                    3.064175                   24.074013
# Depth.num
# 62.700985

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(August.2021)
## we can use model selection instead of picking variables we think are important (by p values)
rda.aug2021.a = ordistep(rda(b.clr_AUG21 ~ 1, data = August.2021[,c(8,10,14:16,18)]),
                         scope=formula(rda.aug2021.0),
                         direction = "forward",
                         permutations = how(nperm=999))
# b.clr_AUG21 ~ Dissolved_OrganicMatter_RFU = best model
rda.aug2021.a$anova # see significance of individual terms in model

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.aug2021.a2 = ordiR2step(rda(b.clr_AUG21 ~ 1, data = August.2021[,c(8,10,14:16,18)]),
                            scope=formula(rda.aug2021.0),
                            permutations = how(nperm=999))
# none

# check best fit model based on above results
anova(rda.aug2021.a, permutations = how(nperm=999)) # p =  0.001, significant
#anova(rda.aug2021.a2, permutations = how(nperm=999)) # p =  0.001, significant

# Let's double check by removing the variables with high VIF, & picking significant variables from ordistep
# dropped Sulfate because had smallest R^2 contribution, also not significant
rda.aug2021.1<-rda(b.clr_AUG21 ~ ORP_mV+Dissolved_OrganicMatter_RFU+DO_Percent_Local+Sulfide_microM,data=August.2021)
summary(rda.aug2021.1)
RsquareAdj(rda.aug2021.1) # how much variation is explained by our model? 11.77%
anova(rda.aug2021.1, by = "terms", permutations = how(nperm=999)) ### by variables
#                             Df Variance      F Pr(>F)
# ORP_mV                       1  179.150 2.1730  0.006 **

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.aug2021.1)
# ORP_mV Dissolved_OrganicMatter_RFU            DO_Percent_Local              Sulfide_microM
# 36.43692                    95.20177                    41.34070                    23.76813

head(August.2021)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.aug2021.b1 = ordistep(rda(b.clr_AUG21 ~ 1, data = August.2021[,c(8,10,14,16)]),
                          scope=formula(rda.aug2021.1),
                          direction = "forward",
                          permutations = how(nperm=999))
rda.aug2021.b1$anova
# b.clr_AUG21 ~ Dissolved_OrganicMatter_RFU = best model
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.aug2021.b2 = ordiR2step(rda(b.clr_AUG21 ~ 1, data = August.2021[,c(8,10,14,16)]),
                            scope=formula(rda.aug2021.1),
                            permutations = how(nperm=999))
# nothing significant; ORP, Sulfide have highest R2
rda.aug2021.b2$anova
# check best fit model based on above results
anova(rda.aug2021.b1, permutations = how(nperm=999)) # p =  0.001, significant

anova(rda.aug2021.0, rda.aug2021.1) # p =  0.003, significant

# choosing sig variables from ordistep & variables with highest variation
rda.aug2021.2<-rda(b.clr_AUG21 ~ Dissolved_OrganicMatter_RFU+ORP_mV+Sulfide_microM,data=August.2021)
summary(rda.aug2021.2)
RsquareAdj(rda.aug2021.2) # how much variation is explained by our model? 13.47%
anova(rda.aug2021.2, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)
#Dissolved_OrganicMatter_RFU  1   157.93 1.8998  0.005 **
#ORP_mV                       1   102.42 1.2320  0.148
#Sulfide_microM               1    77.98 0.9380  0.531
#Residual                     4   332.53
anova(rda.aug2021.2, by=NULL,permutations = how(nperm=999)) # p =  0.019, significant

vif.cca(rda.aug2021.2)
# Dissolved_OrganicMatter_RFU     ORP_mV              Sulfide_microM
#3.835695                        18.095850                   23.552496

# check if ORP & Sulfide are significantly correlated in August, which they are [strong, sig negative corr]
cor.test(meta_scaled[metadata$SampDate=="August.2021",]$Sulfide_microM, meta_scaled[metadata$SampDate=="August.2021",]$ORP_mV, method="pearson") # ******

head(August.2021)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.aug2021.c1 = ordistep(rda(b.clr_AUG21 ~ 1, data = August.2021[,c(10,14,16)]),
                          scope=formula(rda.aug2021.2),
                          direction = "forward",
                          permutations = how(nperm=999))
# b.clr_AUG21 ~ Sulfide_microM = best model
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.aug2021.c2 = ordiR2step(rda(b.clr_AUG21 ~ 1, data = August.2021[,c(10,14,16)]),
                            scope=formula(rda.aug2021.2),
                            permutations = how(nperm=999))
# no sig variables, but ORP & Sulfide have highest R^2

# check best fit model based on above results
anova(rda.aug2021.c1, permutations = how(nperm=999)) # p =  0.001, significant

anova(rda.aug2021.0, rda.aug2021.2) # p =  0.001, significant

rda.aug2021.3<-rda(b.clr_AUG21 ~ ORP_mV+Sulfide_microM,data=August.2021)
summary(rda.aug2021.3)
RsquareAdj(rda.aug2021.3) # how much variation is explained by our model? 13.53%
anova(rda.aug2021.3, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)
#ORP_mV          1   179.15 2.1707   0.01 **
#Sulfide_microM  1    76.38 0.9255   0.51

anova(rda.aug2021.3, by=NULL,permutations = how(nperm=999)) # p =  0.012, significant

vif.cca(rda.aug2021.3)
#ORP_mV Sulfide_microM
#17.40405       17.40405

head(August.2021)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.aug2021.d1 = ordistep(rda(b.clr_AUG21 ~ 1, data = August.2021[,c(10,16)]),
                          scope=formula(rda.aug2021.3),
                          direction = "forward",
                          permutations = how(nperm=999))
# b.clr_AUG21 ~ Sulfide = best model
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.aug2021.d2 = ordiR2step(rda(b.clr_AUG21 ~ 1, data = August.2021[,c(10,16)]),
                            scope=formula(rda.aug2021.3),
                            permutations = how(nperm=999))
# nothing sig, ORP is marginally higher

anova(rda(b.clr_AUG21 ~ Dissolved_OrganicMatter_RFU,data=August.2021)) # p =  0.001, significant

rda.aug2021.4<-rda(b.clr_AUG21 ~ Dissolved_OrganicMatter_RFU+Sulfide_microM,data=August.2021)
summary(rda.aug2021.4)
RsquareAdj(rda.aug2021.4) # how much variation is explained by our model? 14.28%
anova(rda.aug2021.4, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)
#Dissolved_OrganicMatter_RFU  1   157.93 1.9302  0.002 **
#Sulfide_microM               1   101.16 1.2363  0.154
#Residual                     4   329.89

anova(rda.aug2021.4, by=NULL,permutations = how(nperm=999)) # p =  0.002, significant

vif.cca(rda.aug2021.4)
#Dissolved_OrganicMatter_RFU              Sulfide_microM
#3.689057                    3.689057

head(August.2021)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.aug2021.e1 = ordistep(rda(b.clr_AUG21 ~ 1, data = August.2021[,c(14,16)]),
                          scope=formula(rda.aug2021.4),
                          direction = "forward",
                          permutations = how(nperm=999))
# b.clr_AUG21 ~ Dissolved_OrganicMatter_RFU = best model
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.aug2021.e2 = ordiR2step(rda(b.clr_AUG21 ~ 1, data = August.2021[,c(14,16)]),
                            scope=formula(rda.aug2021.4),
                            permutations = how(nperm=999))
# b.clr_AUG21 ~ Sulfide has higher R2, not sig

#### RDA - December 2021 ####

rownames(December.2021) %in% rownames(b.clr_DEC21) # check order of DFs
head(December.2021)

rda.dec2021.0<-rda(b.clr_DEC21 ~ DO_Percent_Local+ORP_mV+Dissolved_OrganicMatter_RFU+Sulfate_milliM+Sulfide_microM+Depth.num,data=December.2021)

# check summary of RDA
rda.dec2021.0
summary(rda.dec2021.0)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.dec2021.0) # 1.5%
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.dec2021.0, permutations = how(nperm=999)) # not significant

## we can also do a permutation test by RDA axis
#anova(rda.dec2021.0, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.dec2021.0, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)
# ORP_mV                       1   75.737 1.2457  0.017 *

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.dec2021.0)
# DO_Percent_Local            RP_mV.    Dissolved_OrganicMatter_RFU              Sulfate_milliM              Sulfide_microM
#160.90484                    18.01647                    24.78209                    11.06779                    16.82469
#Depth.num
#302.42627

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(December.2021)
## we can use model selection instead of picking variables we think are important (by p values)
rda.dec2021.a = ordistep(rda(b.clr_DEC21 ~ 1, data = December.2021[,c(8,10,14:16,18)]),
                         scope=formula(rda.dec2021.0),
                         direction = "forward",
                         permutations = how(nperm=999))
# b.clr_DEC21 ~ ORP_mV  - best model
rda.dec2021.a$anova # see significance of individual terms in model

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.dec2021.a2 = ordiR2step(rda(b.clr_DEC21 ~ 1, data = December.2021[,c(8,10,14:16,18)]),
                            scope=formula(rda.dec2021.0),
                            permutations = how(nperm=999))
# nothing sig

# check best fit model based on above results
anova(rda.dec2021.a, permutations = how(nperm=999))
#anova(rda.dec2021.a2, permutations = how(nperm=999)) # not significant

# Let's get rid of depth and rerun
rda.dec2021.1<-rda(b.clr_DEC21 ~ DO_Percent_Local+ORP_mV+Dissolved_OrganicMatter_RFU+Sulfate_milliM+Sulfide_microM,data=December.2021)
summary(rda.dec2021.1)
RsquareAdj(rda.dec2021.1) # how much variation is explained by our model? 4.11%
anova(rda.dec2021.1, by = "terms", permutations = how(nperm=999)) ### by variables
#                            Df Variance      F Pr(>F)
#ORP_mV          1   75.737 1.2794  0.001 ***

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.dec2021.1)
#DO_Percent_Local            ORP_mV Dissolved_OrganicMatter_RFU              Sulfate_milliM              Sulfide_microM
# 2.261942                    3.423645                    2.539956                    1.153821                    1.844823

head(December.2021)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.dec2021.b1 = ordistep(rda(b.clr_DEC21 ~ 1, data = December.2021[,c(8,10,14:16)]),
                          scope=formula(rda.dec2021.1),
                          direction = "forward",
                          permutations = how(nperm=999))
#b.clr_DEC21 ~ ORP_mV

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.dec2021.b2 = ordiR2step(rda(b.clr_DEC21 ~ 1, data = December.2021[,c(8,10,14:16)]),
                            scope=formula(rda.dec2021.1),
                            permutations = how(nperm=999))
# ORP has highest R^2 but not sig

# check best fit model based on above results
anova(rda.dec2021.b1, permutations = how(nperm=999))

anova(rda.dec2021.0, rda.dec2021.1) # no significant difference

rda.dec2021.2<-rda(b.clr_DEC21 ~ DO_Percent_Local+ORP_mV+Dissolved_OrganicMatter_RFU+Sulfate_milliM,data=December.2021)
summary(rda.dec2021.2)
RsquareAdj(rda.dec2021.2) # how much variation is explained by our model? 1.17%
anova(rda.dec2021.2, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                           Df Variance      F Pr(>F)
#ORP_mV          1   75.737 1.2413  0.013 *

vif.cca(rda.dec2021.2)
#DO_Percent_Local             ORP_mV      Dissolved_OrganicMatter_RFU              Sulfate_milliM
# 2.231614                    2.802203                    2.498753                    1.053136

head(December.2021)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.dec2021.c1 = ordistep(rda(b.clr_DEC21 ~ 1, data = December.2021[,c(8,10,14:15)]),
                          scope=formula(rda.dec2021.2),
                          direction = "forward",
                          permutations = how(nperm=999))
# b.clr_DEC21 ~ ORP_mV

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.dec2021.c2 = ordiR2step(rda(b.clr_DEC21 ~ 1, data = December.2021[,c(8,10,14:15)]),
                            scope=formula(rda.dec2021.2),
                            permutations = how(nperm=999))
# ORP has highest R^2 but not sig

# check best fit model based on above result
anova(rda.dec2021.0, rda.dec2021.2) # no significant difference

rda.dec2021.3<-rda(b.clr_DEC21 ~ ORP_mV+DO_Percent_Local+Sulfate_milliM,data=December.2021)
summary(rda.dec2021.3)
RsquareAdj(rda.dec2021.3) # how much variation is explained by our model? 4.49%
anova(rda.dec2021.3, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                           Df Variance      F Pr(>F)
#ORP_mV            1   77.391 1.3125  0.006 **

vif.cca(rda.dec2021.3)
#ORP_mV DO_Percent_Local   Sulfate_milliM
#1.230812         1.232172         1.007516

head(December.2021)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.dec2021.d1 = ordistep(rda(b.clr_DEC21 ~ 1, data = December.2021[,c(8,10,15)]),
                          scope=formula(rda.dec2021.3),
                          direction = "forward",
                          permutations = how(nperm=999))
# ORP is sig

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.dec2021.d2 = ordiR2step(rda(b.clr_DEC21 ~ 1, data = December.2021[,c(8,10,15)]),
                            scope=formula(rda.dec2021.3),
                            permutations = how(nperm=999))
# ORP is sig

# check best fit model based on above result
anova(rda.dec2021.0, rda.dec2021.3) # no significant difference
anova(rda.dec2021.2, rda.dec2021.3) # no significant difference

rda.dec2021.4<-rda(b.clr_DEC21 ~ ORP_mV+Sulfate_milliM,data=December.2021)
summary(rda.dec2021.4)
RsquareAdj(rda.dec2021.4) # how much variation is explained by our model? 5.3%
anova(rda.dec2021.4, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                           Df Variance      F Pr(>F)
#ORP_mV          1   77.391 1.3240  0.002 **
#Sulfate_milliM  1   62.509 1.0694  0.184
#Residual        5  292.258

vif.cca(rda.dec2021.4)
#ORP_mV Sulfate_milliM
#1.001605       1.001605

head(December.2021)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.dec2021.e1 = ordistep(rda(b.clr_DEC21 ~ 1, data = December.2021[,c(10,15)]),
                          scope=formula(rda.dec2021.4),
                          direction = "forward",
                          permutations = how(nperm=999))
#  ORP is sig

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.dec2021.e2 = ordiR2step(rda(b.clr_DEC21 ~ 1, data = December.2021[,c(10,15)]),
                            scope=formula(rda.dec2021.4),
                            permutations = how(nperm=999))
# ORP is sig

#### RDA - April 2022 ####

rownames(April.2022) %in% rownames(b.clr_APR22) # check order of DFs
head(April.2022)

rda.apr2022.0<-rda(b.clr_APR22 ~ DO_Percent_Local+ORP_mV+Dissolved_OrganicMatter_RFU+Sulfate_milliM+Sulfide_microM+Depth.num,data=April.2022)

# check summary of RDA
rda.apr2022.0
summary(rda.apr2022.0)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.apr2022.0) # -10%
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.apr2022.0, permutations = how(nperm=999)) # not significant

## we can also do a permutation test by RDA axis
#anova(rda.apr2022.0, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.apr2022.0, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#   nothing sig

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.apr2022.0)
# DO_Percent_Local                      ORP_mV Dissolved_OrganicMatter_RFU              Sulfate_milliM              Sulfide_microM
#33.781602                   27.495126                   37.344732                    2.421525                    3.771850
#Depth.num
#8.486058

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(April.2022)
## we can use model selection instead of picking variables we think are important (by p values)
rda.apr2022.a = ordistep(rda(b.clr_APR22 ~ 1, data = April.2022[,c(8,10,14:16,18)]),
                         scope=formula(rda.apr2022.0),
                         direction = "forward",
                         permutations = how(nperm=999))
# b.clr_APR22 ~ DOM
rda.apr2022.a$anova # see significance of individual terms in model

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.apr2022.a2 = ordiR2step(rda(b.clr_APR22 ~ 1, data = April.2022[,c(8,10,14:16,18)]),
                            scope=formula(rda.apr2022.0),
                            permutations = how(nperm=999))
# nothing

# check best fit model based on above results
#anova(rda.apr2022.a, permutations = how(nperm=999)) # p =  0.036, significant

rda.apr2022.1<-rda(b.clr_APR22 ~ DO_Percent_Local+ORP_mV+Dissolved_OrganicMatter_RFU+Sulfate_milliM+Sulfide_microM,data=April.2022)
summary(rda.apr2022.1)
RsquareAdj(rda.apr2022.1) # how much variation is explained by our model? -2.5%
anova(rda.apr2022.1, by = "terms", permutations = how(nperm=999)) ### by variables
#  nothing significant

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.apr2022.1)
# DO_Percent_Local                      ORP_mV Dissolved_OrganicMatter_RFU              Sulfate_milliM              Sulfide_microM
# 17.704118                   20.480959                   36.649047                    1.941701                    2.946815

head(April.2022)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.apr2022.b1 = ordistep(rda(b.clr_APR22 ~ 1, data = April.2022[,c(8,10,14:16)]),
                          scope=formula(rda.apr2022.1),
                          direction = "forward",
                          permutations = how(nperm=999))
# b.clr_APR22 ~ DOM
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.apr2022.b2 = ordiR2step(rda(b.clr_APR22 ~ 1, data = April.2022[,c(8,10,14:16)]),
                            scope=formula(rda.apr2022.1),
                            permutations = how(nperm=999))
# b.clr_APR22 - nothing significant

# check best fit model based on above results
anova(rda.apr2022.b1, permutations = how(nperm=999))

anova(rda.apr2022.0, rda.apr2022.1) # no significant difference

rda.apr2022.2<-rda(b.clr_APR22 ~ DO_Percent_Local+Dissolved_OrganicMatter_RFU+Sulfate_milliM,data=April.2022)
summary(rda.apr2022.2)
RsquareAdj(rda.apr2022.2) # how much variation is explained by our model? 0.47%
anova(rda.apr2022.2, by = "terms", permutations = how(nperm=999)) ### by variables
# ^ nothing significant
anova(rda.apr2022.2, by = NULL, permutations = how(nperm=999)) ### model not sig

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.apr2022.2)
#DO_Percent_Local Dissolved_OrganicMatter_RFU              Sulfate_milliM
#12.670210                   10.431717                    1.770423

head(April.2022)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.apr2022.c1 = ordistep(rda(b.clr_APR22 ~ 1, data = April.2022[,c(8,14:15)]),
                          scope=formula(rda.apr2022.2),
                          direction = "forward",
                          permutations = how(nperm=999))
# b.clr_APR22 ~ DOM
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.apr2022.c2 = ordiR2step(rda(b.clr_APR22 ~ 1, data = April.2022[,c(8,14:15)]),
                            scope=formula(rda.apr2022.2),
                            permutations = how(nperm=999))

# nothing
anova(rda.apr2022.c1, permutations = how(nperm=999)) # 0.04

anova(rda.apr2022.0, rda.apr2022.2) # no significant difference

rda.apr2022.3<-rda(b.clr_APR22 ~ Dissolved_OrganicMatter_RFU+Sulfate_milliM,data=April.2022)
summary(rda.apr2022.3)
RsquareAdj(rda.apr2022.3) # how much variation is explained by our model? 2.61%
anova(rda.apr2022.3, by = "terms", permutations = how(nperm=999)) ### by variables
#               Df Variance      F Pr(>F)
# Dissolved_OrganicMatter_RFU  1   60.077 1.1195  0.011 *
# Sulfate_milliM               1   57.312 1.0680  0.100 .
# Residual                     5  268.313

anova(rda.apr2022.3, by = NULL, permutations = how(nperm=999)) #0.044

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.apr2022.3)
# Dissolved_OrganicMatter_RFU              Sulfate_milliM
# 1.231086                    1.231086

head(April.2022)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.apr2022.e1 = ordistep(rda(b.clr_APR22 ~ 1, data = April.2022[,c(14:15)]),
                          scope=formula(rda.apr2022.3),
                          direction = "forward",
                          permutations = how(nperm=999))
# b.clr_APR22 ~ DOM
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.apr2022.e2 = ordiR2step(rda(b.clr_APR22 ~ 1, data = April.2022[,c(14:15)]),
                            scope=formula(rda.apr2022.3),
                            permutations = how(nperm=999))
# b.clr_APR22 ~ DOM

#### Final RDAs ####
# RDA by sampling timepoint
head(meta_scaled)
head(b.clr)
rownames(b.clr) %in% rownames(meta_scaled) # sanity check 1

# all data
#rda.all2$call # best model for all data

rda.all<-rda(b.clr ~ Temp_DegC + Dissolved_OrganicMatter_RFU + DO_Percent_Local,data=meta_scaled)
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

# August 2021
#rda.aug2021.4$call # best model

rda.aug2021<-rda(b.clr_AUG21 ~ Dissolved_OrganicMatter_RFU+Sulfide_microM,data=August.2021)
summary(rda.aug2021)
RsquareAdj(rda.aug2021) # how much variation is explained by our model? 14.28%
anova(rda.aug2021, permutations = how(nperm=999)) # p-value = 0.005 **
anova(rda.aug2021, by = "terms", permutations = how(nperm=999))
#                           Df Variance      F Pr(>F)
#Dissolved_OrganicMatter_RFU  1   157.93 1.9302  0.001 ***
#Sulfide_microM               1   101.16 1.2363  0.148
#Residual                     5   409.10
aov.rda.aug<-anova(rda.aug2021, by = "terms", permutations = how(nperm=999))
p.adjust(aov.rda.aug$`Pr(>F)`,method="bonferroni",n=3) # adjusted pvalues

# December 2021
#rda.dec2021.2$call # best model from above

rda.dec2021<-rda(b.clr_DEC21 ~ ORP_mV + Sulfate_milliM,data=December.2021)
summary(rda.dec2021)
RsquareAdj(rda.dec2021) # how much variation is explained by our model? 5.3%
anova(rda.dec2021, permutations = how(nperm=999)) # p-value = 0.005
anova(rda.dec2021, by = "terms", permutations = how(nperm=999))
#                 Df Variance      F Pr(>F)
# ORP_mV          1   77.391 1.3240  0.003 **
# Sulfate_milliM  1   62.509 1.0694  0.192
# Residual        5  292.258
aov.rda.dec<-anova(rda.dec2021, by = "terms", permutations = how(nperm=999))
p.adjust(aov.rda.dec$`Pr(>F)`,method="bonferroni",n=3) # adjusted pvalues

# April 2022
#rda.apr2022.3$call  #best mode

rda.apr2022<-rda(b.clr_APR22 ~ Dissolved_OrganicMatter_RFU + Sulfate_milliM,data=April.2022)
summary(rda.apr2022)
RsquareAdj(rda.apr2022) # how much variation is explained by our model? 2.61%
anova(rda.apr2022, permutations = how(nperm=999)) # p-value = 0.039
anova(rda.apr2022, by = "terms", permutations = how(nperm=999))
#                           Df Variance      F Pr(>F)
# Dissolved_OrganicMatter_RFU  1   60.077 1.1195  0.011 *
# Sulfate_milliM               1   57.312 1.0680  0.092 .
# Residual                     5  268.313
aov.rda.apr<-anova(rda.apr2022, by = "terms", permutations = how(nperm=999))
p.adjust(aov.rda.apr$`Pr(>F)`,method="bonferroni",n=3) # adjusted pvalues

# save RDAs as R object
save.image("data/SSW_Amplicon_EnvDriver_RDAsOnly.Rdata")

#### Plot RDA - ALL data ####
#plot(rda.aug2021) # depending on how many species you have, this step may take a while
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
rda.all.part<-varpart(b.clr, meta_scaled$Temp_DegC, meta_scaled$Dissolved_OrganicMatter_RFU,meta_scaled$DO_Percent_Local)
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
rownames(rda.sum.all$sites) %in% rownames(meta_scaled)
rda.axes.all<-data.frame(RDA1=rda.sum.all$sites[,1], RDA2=rda.sum.all$sites[,2], SampleID=rownames(rda.sum.all$sites), Depth_m=meta_scaled$Depth_m, SampDate=meta_scaled$SampDate)

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
  scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
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
  scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Seawater",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [30.80%]") + ylab("RDA2 [23.65%]")

ggsave(rda.plot3,filename = "figures/EnvDrivers/SSW_16S_RDA_AllData_bigger.png", width=15, height=15, dpi=600)

#### Plot RDA - Aug 2021 ####
#plot(rda.aug2021) # depending on how many species you have, this step may take a while
plot(rda.aug2021, scaling = 1)
## scaling = 1 -> emphasizes relationships among sites
plot(rda.aug2021, scaling = 2)
## scaling = 2 -> emphasizes relationships among species

# check summary of RDA
summary(rda.aug2021)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.aug2021) # 14.28%
## ^^ use this b/c chance correlations can inflate R^2

png('figures/EnvDrivers/SSW_Aug21_autoplot_rda_example.png',width = 700, height = 600, res=100)
autoplot(rda.aug2021, arrows = TRUE,data = rda.aug2021 ,layers=c("biplot","sites"),label = FALSE, label.size = 3, shape = FALSE, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3, scale= 0)+theme_classic()
dev.off()
## FOR AUTOPLOT -> must load packagve ggvegan first

# variance partitioning of RDA
rda.aug21.part<-varpart(b.clr_AUG21, August.2021$Dissolved_OrganicMatter_RFU, August.2021$Sulfide_microM)
rda.aug21.part$part
# plot variance partitioning results
png('figures/EnvDrivers/SSW_Aug21_RDA_VariancePartitioning.png',width = 900, height = 900, res=100)
plot(rda.aug21.part,
     Xnames = c("DOM (RFU)", "Sulfide (microM)"), # name the partitions
     bg = c("#ffbe0b", "darkgreen"), alpha = 80, # colour the circles
     digits = 3, # only show 3 digits
     cex = 1.5)
dev.off()

rda.sum.a21<-summary(rda.aug2021)
rda.sum.a21$sites[,1:2]
rda.sum.a21$cont #cumulative proportion of variance per axis
# RDA1=26.71%, RDA2=12.06%

# create data frame w/ RDA axes for sites
rda.axes.a21<-data.frame(RDA1=rda.sum.a21$sites[,1], RDA2=rda.sum.a21$sites[,2], SampleID=rownames(rda.sum.a21$sites), Depth_m=August.2021$Depth_m)

# create data frame w/ RDA axes for variables
arrows.a21<-data.frame(RDA1=rda.sum.a21$biplot[,1], RDA2=rda.sum.a21$biplot[,2], Label=rownames(rda.sum.a21$biplot))
#arrows.a21$Label[(arrows.a21$Label) == "ORP_mV"] <- "ORP (mV)"
arrows.a21$Label[(arrows.a21$Label) == "Dissolved_OrganicMatter_RFU"] <- "DOM (RFU)"
#arrows.a21$Label[(arrows.a21$Label) == "Sulfate_milliM"] <- "Sulfate (milliM)"
arrows.a21$Label[(arrows.a21$Label) == "Sulfide_microM"] <- "Sulfide (microM)"

rda.plot5<-ggplot(rda.axes.a21, aes(x = RDA1, y = RDA2)) + geom_point(size=2) +
  geom_segment(data = arrows.a21,mapping = aes(x = 0, y = 0, xend = RDA1, yend = RDA2),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.5,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.a21,aes(label = Label, x = RDA1, y = RDA2, fontface="bold"))+
  coord_fixed() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))

rda.plot6<-ggplot(rda.axes.a21, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=as.numeric(as.character(Depth_m))),size=4) +
  geom_segment(data = arrows.a21,mapping = aes(x = 0, y = 0, xend = RDA1*8, yend = RDA2*8),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.8,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.a21,aes(label = Label, x = RDA1*9.85, y = RDA2*9.5, fontface="bold"), size=4)+
  coord_fixed(ratio = 1, xlim = c(-5,15), ylim = c(-10,10)) + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Seawater, August 2021",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [26.71%]") + ylab("RDA2 [12.06%]")

ggsave(rda.plot6,filename = "figures/EnvDrivers/SSW_16S_RDA_Aug2021.png", width=16, height=12, dpi=600)

rda.plot6b<-ggplot(rda.axes.a21, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=as.numeric(as.character(Depth_m))),size=5) +
  geom_segment(data = arrows.a21,mapping = aes(x = 0, y = 0, xend = RDA1*8, yend = RDA2*8),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 1,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.a21,aes(label = Label, x = RDA1*9, y = RDA2*9.5, fontface="bold"), size=5)+
  coord_fixed(ratio = 1, xlim = c(-10,10), ylim = c(-10,10)) + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Seawater, August 2021",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [26.71%]") + ylab("RDA2 [12.06%]")

ggsave(rda.plot6b,filename = "figures/EnvDrivers/SSW_16S_RDA_Aug2021_bigger.png", width=15, height=15, dpi=600)

#### Plot RDA - Dec 2021 ####
#plot(rda.dec2021) # depending on how many species you have, this step may take a while
plot(rda.dec2021, scaling = 1)
## scaling = 1 -> emphasizes relationships among sites
plot(rda.dec2021, scaling = 2)
## scaling = 2 -> emphasizes relationships among species

# check summary of RDA
summary(rda.dec2021)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.dec2021) # 0.0532124
## ^^ use this b/c chance correlations can inflate R^2

png('figures/EnvDrivers/SSW_Dec21_autoplot_rda_example.png',width = 700, height = 600, res=100)
autoplot(rda.dec2021, arrows = TRUE,data = rda.dec2021 ,layers=c("biplot","sites"),label = FALSE, label.size = 3, shape = FALSE, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3, scale= 0)+theme_classic()
dev.off()
## FOR AUTOPLOT -> must load packagve ggvegan first

# variance partitioning of RDA
rda.dec21.part<-varpart(b.clr_DEC21, December.2021$ORP_mV, December.2021$Sulfate_milliM)
rda.dec21.part$part
# plot variance partitioning results
png('figures/EnvDrivers/SSW_Dec21_RDA_VariancePartitioning.png',width = 900, height = 900, res=100)
plot(rda.dec21.part,
     Xnames = c("ORP (mV)", "Sulfate (milliM)"), # name the partitions
     bg = c("#3a0ca3", "#8ac926"), alpha = 80, # colour the circles
     digits = 3, # only show 3 digits
     cex = 1.5)
dev.off()

rda.sum.d21<-summary(rda.dec2021)
rda.sum.d21$sites[,1:2]
rda.sum.d21$cont # cumulative proportion of variation per axis
# RDA1 = 18.19, RDA2 = 14.18

# create data frame w/ RDA axes for sites
rda.axes.d21<-data.frame(RDA1=rda.sum.d21$sites[,1], RDA2=rda.sum.d21$sites[,2], SampleID=rownames(rda.sum.d21$sites), Depth_m=December.2021$Depth_m)

# create data frame w/ RDA axes for variables
arrows.d21<-data.frame(RDA1=rda.sum.d21$biplot[,1], RDA2=rda.sum.d21$biplot[,2], Label=rownames(rda.sum.d21$biplot))
arrows.d21$Label[(arrows.d21$Label) == "ORP_mV"] <- "ORP (mV)"
arrows.d21$Label[(arrows.d21$Label) == "Sulfate_milliM"] <- "Sulfate (milliM)"
#arrows.d21$Label[(arrows.d21$Label) == "DO_Percent_Local"] <- "DO%"

rda.plot7<-ggplot(rda.axes.d21, aes(x = RDA1, y = RDA2)) + geom_point(size=2) +
  geom_segment(data = arrows.d21,mapping = aes(x = 0, y = 0, xend = RDA1, yend = RDA2),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.5,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.d21,aes(label = Label, x = RDA1, y = RDA2, fontface="bold"))+
  coord_fixed() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))

rda.plot8<-ggplot(rda.axes.d21, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=as.numeric(as.character(Depth_m))),size=4) +
  geom_segment(data = arrows.d21,mapping = aes(x = 0, y = 0, xend = RDA1*9, yend = RDA2*9),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.8,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.d21,aes(label = Label, x = RDA1*10.5, y = RDA2*10, fontface="bold"), size=4)+
  coord_fixed(ratio = 1, xlim = c(-10,11), ylim = c(-10,10)) + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Seawater",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [18.19%]") + ylab("RDA2 [14.18%]")

ggsave(rda.plot8,filename = "figures/EnvDrivers/SSW_16S_RDA_Dec2021.png", width=15, height=12, dpi=600)

rda.plot8b<-ggplot(rda.axes.d21, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=as.numeric(as.character(Depth_m))),size=5) +
  geom_segment(data = arrows.d21,mapping = aes(x = 0, y = 0, xend = RDA1*9, yend = RDA2*9),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 1,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.d21,aes(label = Label, x = RDA1*10.5, y = RDA2*10.5, fontface="bold"), size=5)+
  coord_fixed(ratio = 1, xlim = c(-10,11), ylim = c(-10,10)) + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Seawater, December 2021",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [18.19%]") + ylab("RDA2 [14.18%]")

ggsave(rda.plot8b,filename = "figures/EnvDrivers/SSW_16S_RDA_Dec2021_bigger.png", width=15, height=15, dpi=600)

#### Plot RDA - Apr 2022 ####
#plot(rda.dec2021) # depending on how many species you have, this step may take a while
plot(rda.apr2022, scaling = 1)
## scaling = 1 -> emphasizes relationships among sites
plot(rda.apr2022, scaling = 2)
## scaling = 2 -> emphasizes relationships among species

# check summary of RDA
summary(rda.apr2022)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.apr2022) # 2.61%
## ^^ use this b/c chance correlations can inflate R^2

png('figures/EnvDrivers/SSW_Apr22_autoplot_rda_example.png',width = 700, height = 600, res=100)
autoplot(rda.apr2022, arrows = TRUE,data = rda.apr2022 ,layers=c("biplot","sites"),label = FALSE, label.size = 3, shape = FALSE, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3, scale= 0)+theme_classic()
dev.off()
## FOR AUTOPLOT -> must load packagve ggvegan first

# variance partitioning of RDA
rda.apr22.part<-varpart(b.clr_APR22, April.2022$Dissolved_OrganicMatter_RFU, April.2022$Sulfate_milliM)
rda.apr22.part$part
# plot variance partitioning results
png('figures/EnvDrivers/SSW_Apr22_RDA_VariancePartitioning.png',width = 900, height = 900, res=100)
plot(rda.apr22.part,
     Xnames = c("DOM (RFU)", "Sulfate (milliM)"), # name the partitions
     bg = c("#ffbe0b", "#8ac926"), alpha = 80, # colour the circles
     digits = 3, # only show 3 digits
     cex = 1.5)
dev.off()

rda.sum.a22<-summary(rda.apr2022)
rda.sum.a22$sites[,1:2]
rda.sum.a22$cont
# RDA1 = 16.04, RDA2 = 14.40

# create data frame w/ RDA axes for sites
rda.axes.a22<-data.frame(RDA1=rda.sum.a22$sites[,1], RDA2=rda.sum.a22$sites[,2], SampleID=rownames(rda.sum.a22$sites), Depth_m=April.2022$Depth_m)

# create data frame w/ RDA axes for variables
arrows.a22<-data.frame(RDA1=rda.sum.a22$biplot[,1], RDA2=rda.sum.a22$biplot[,2], Label=rownames(rda.sum.a22$biplot))
arrows.a22$Label[(arrows.a22$Label) == "Dissolved_OrganicMatter_RFU"] <- "DOM (RFU)"
arrows.a22$Label[(arrows.a22$Label) == "Sulfate_milliM"] <- "Sulfate (milliM)"
#arrows.a22$Label[(arrows.a22$Label) == "DO_Percent_Local"] <- "DO%"
#arrows.a22$Label[(arrows.a22$Label) == "Temp_DegC"] <- "Temp (C)"

rda.plot9<-ggplot(rda.axes.a22, aes(x = RDA1, y = RDA2)) + geom_point(size=2) +
  geom_segment(data = arrows.a22,mapping = aes(x = 0, y = 0, xend = RDA1, yend = RDA2),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.5,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.a22,aes(label = Label, x = RDA1, y = RDA2, fontface="bold"))+
  coord_fixed() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))

rda.plot10<-ggplot(rda.axes.a22, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=as.numeric(as.character(Depth_m))),size=4) +
  geom_segment(data = arrows.a22,mapping = aes(x = 0, y = 0, xend = RDA1*8, yend = RDA2*8),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.8,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.a22,aes(label = Label, x = RDA1*9, y = RDA2*9, fontface="bold"), size=4)+
  coord_fixed(ratio = 1, xlim = c(-10,10), ylim = c(-10,10)) + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Seawater, April 2022",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [16.04%]") + ylab("RDA2 [14.40%]")

ggsave(rda.plot10,filename = "figures/EnvDrivers/SSW_16S_RDA_April2022.png", width=15, height=12, dpi=600)

rda.plot10b<-ggplot(rda.axes.a22, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=as.numeric(as.character(Depth_m))),size=5) +
  geom_segment(data = arrows.a22,mapping = aes(x = 0, y = 0, xend = RDA1*9, yend = RDA2*9),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 1,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.a22,aes(label = Label, x = RDA1*10, y = RDA2*10, fontface="bold"), size=5)+
  coord_fixed(ratio = 1, xlim = c(-10,10), ylim = c(-10,10)) + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Seawater, December 2021",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [16.04%]") + ylab("RDA2 [14.40]")

ggsave(rda.plot10b,filename = "figures/EnvDrivers/SSW_16S_RDA_April2022_bigger.png", width=15, height=15, dpi=600)

#### Save Progress ####

save.image("data/SSW_Amplicon_EnvDriver.Rdata")
