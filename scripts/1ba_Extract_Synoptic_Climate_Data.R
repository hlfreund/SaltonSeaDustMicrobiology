#### Load Packages ####

suppressPackageStartupMessages({ # load packages quietly
  library(phyloseq)
  library(ggplot2)
  library(vegan)
  library(ggpubr)
  library(tidyverse)
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
  library(fst)
  library(mesowest)
  library(data.table)
  library(geosphere)
  library(lubridate)
  library(circular)
  #library(weatherData)

})

#### Log Into MesoWest w/ API key ####
# mesowest package GitHub: https://github.com/fickse/mesowest

requestToken(apikey = "KunTxSCUqjGlu9CBtKJJRmxjYHnFc6PK4uPoDgcUOg")

getparams()
#### Retrieve Data for Sites of Interest ####
# MesoWest variables: https://developers.synopticdata.com/mesonet/v2/api-variables/
# air_temp (ºC)
# wind_speed (m/s)
# wind_gust (m/s)
# wind_direction (degrees)
# wind_cardinal_direction (char)

# Set coordinates for research site; Salton Sea: 33.3286° N, 115.8434° W
seaCoords <- c(-115.8434, 33.3286)

# Get all MesoWest sites in Imperial and Riverside counties (get info for sites in your research area)
sites <- bind_rows(
  mw(service = 'metadata', complete=1, state='CA', county='Riverside' )[[1]],
  mw(service = 'metadata', complete=1, state='CA', county='Imperial' )[[1]]) %>%
  mutate(LATITUDE = as.numeric(LATITUDE)) %>%
  mutate(LONGITUDE = as.numeric(LONGITUDE)) %>%
  mutate(start = ymd_hms(PERIOD_OF_RECORD$start)) %>%
  mutate(end = ymd_hms(PERIOD_OF_RECORD$end))

## more on the variables here here: https://docs.synopticdata.com/services/json-output-format ; notes on less obvious labels below
## STID = station abbrevation/identifier
## MNET_ID = numeric ID of the network
## ELEV_DEM = elevation produced by the Digital Elevation Model
## STATUS = string defining if the stations is active or not
## NWSZONE = NWS Forecast zone
## UNITS = units for elevation and sensor height (position); observation units are within the UNITS element

# Filter for sites within 40 km of Salton Sea

## First we calculate the distance between the sites we just pulled info for, and the coordinates of our research site (Salton Sea)
## ^ use haversine method to calculate shortest distance between site and Salton Sea coordinates and record that in sites$dist
sites <- sites %>%
  mutate(dist = distm(sites[, c('LONGITUDE', 'LATITUDE')], seaCoords, fun = distHaversine))
# distm - distance matrix of a set of points or between 2 sets of points
# distHaversine - shortest distance between 2 points (i.e., great-circle-distance or as the crow flies), according to the 'haversine method'
# ^^ this method assumes a spherical Earth

# ** NEED TO FILTER AND FIND SITE THAT IS CLOSEST TO EACH SITE, MONTH/YEAR

## Then filter the sites by the following info:
### increase distance to find stations near all sample sites
### site has data between time points of dust collections that we are looking at
sites <- sites %>%
  filter(dist <= 2000000) %>%
  filter(start <= as.Date('2020-06-01')) %>%
  filter(end >= as.Date('2021-12-31'))

# pull out only sites that we are interested in
# BDC’s closest weather station - UCDE
# PD’s closest weather station - CI200
# DP’s closest weather station - DPMC1
# WI’s closest weather station - CQ125 or UP614

our.site.list<-c("UCDE","CI200","DPMC1","CQ125")

sites<-sites[sites$STID %in% our.site.list,] # only keep sites of interest

#write_rds(sites, 'data/MesoWest_SaltonSea_2020-2021_SiteData.rds')

#### Import Salton Sea Dust Microbiome Metadata ####

dust_meta<-as.data.frame(read_excel("data/Metadata_EnvMiSeqPlate_Winter23.xlsx", sheet="SSea_Dust_Metadata_Updated"), header=TRUE)
head(dust_meta)

# pull out station IDs and dust collector deployment/retrieval date/times
stids<-unique(data.frame(dust_meta[,14:17]))
stids # sanity check

precip.stids<-unique(data.frame(dust_meta[,c(14,16:17,20)]))
#stids$STID<-gsub("UP614","CQ125",stids$STID) # just to see what data we can get from CQ125 - is it more than UP614
#stids
#### Pull Out Climate TimeSeries Data for Sites Near Salton Sea ####
## full timeseries variable list: https://demos.synopticdata.com/variables/index.html

i.STID <- sites$STID[1] # checking the index (1)

time_series_data <- lapply(sites$STID, function(i.STID){
  dt <- mw(service = 'timeseries', stid =i.STID,
           vars = c('wind_speed','wind_gust','wind_direction','wind_cardinal_direction',
                    'air_temp','relative_humidity'),
           start = '202005130001', end = '202112090001', jsonsimplify= TRUE)

  obs <- dt[['STATION']][['OBSERVATIONS']]

  obs <- lapply(obs, unlist) %>%
    as.data.table

  obs[, STID:=i.STID]

  return(obs)
}) %>% rbindlist(fill=TRUE)

setnames(time_series_data, names(time_series_data) %>% str_remove_all('_set_1(d?)'))
# ^ remove '_set_1(d?)' from column names in time_series_data

time_series_data[, date_time:=ymd_hms(date_time)]
time_series_data[, date_time_hour:=ceiling_date(date_time, 'hour')]
# time_series_data[, DateLST:=ymd_hms(date_time, tz = 'America/Los_Angeles')]
# time_series_data[, DateLST:=ceiling_date(DateLST, 'hour')]

# the process of converting wind vector into a mathematical vector: http://colaweb.gmu.edu/dev/clim301/lectures/wind/wind-uv
# NOTES: Remember that a vector can be broken into component vectors, where the i unit vector runs parallel to the x asis, and the j vector runs parallel to the y axis.
# For winds, the u wind is parallel to the x axis. A positive u wind is from the west (to east). A negative u wind is from the east (to west).
## E is 90 degrees, W is 270
# The v wind runs parallel to the y axis. A positive v wind is from the south, and a negative v wind is from the north.
## S is 180, N is 360 or 0

# first, calculate mathematical direction
time_series_data$math_wind.dir<-270-time_series_data$wind_direction
#time_series_data[math_wind.dir < 0, math_wind.dir := math_wind.dir + 360]
time_series_data$math_wind.dir<-ifelse(time_series_data$math_wind.dir < 0, time_series_data$math_wind.dir + 360, time_series_data$math_wind.dir)
# now calculate component vector winds, u & v
#u = ws * cos(θ), where theta is wind direction (in radians, which you convert from degrees)
#v = ws * sin(θ)

time_series_data$u_E.W.wind<-time_series_data$wind_speed * cos(pi * time_series_data$math_wind.dir/180) #
time_series_data$v_N.S.wind<-time_series_data$wind_speed * sin(pi * time_series_data$math_wind.dir/180)

# we can use u and v to determine direciton and magnitude of the wind vector
# wind speed = sqrt(u2+v2)
# wind direciton = arctan(v/u)

# now for CIRCULAR wind direction
time_series_data[, wd:=circular(wind_direction, template='geographics', units='degrees')]
# ^ create column $wd which is the wind direction in degrees, not sure why we need to create this since wind_direction has the same info?

head(time_series_data)
# use time_series_data for your sanity checks to see if your u and v wind components and you wind direction are describing the same cardinal direction
# remember: a positive u is from the west, a negative u is from the east; a positive v is from the south, a negative v is from the north
# NorthEast winds should have negative u and v!

# note: [,.()] is data.table notation; just keeps data in data.table format
# https://stackoverflow.com/questions/33808705/dot-preceding-parentheses-in-data-table
# looks like we are taking the mean of air temp & wind speed to get averages per hour rather than values per 5 minutes
time_series_data_scalar <- time_series_data[, .(air_temp = mean(air_temp, na.rm=TRUE),
                              wind_speed=mean(wind_speed, na.rm=TRUE),
                              wind_gust=max(wind_gust, na.rm=TRUE),
                              relative_humidity=mean(relative_humidity,na.rm=TRUE),
                              u_E.W.wind=mean(u_E.W.wind,na.rm=TRUE),
                              v_N.S.wind=mean(v_N.S.wind,na.rm=TRUE)),
                          .(STID, date_time_hour)]
# ^ need to add accumulated precip & relative humidity in here, want to make sure I am treating the data the same way as Will


time_series_data_dir <- time_series_data[!is.na(wd) & !is.na(wind_speed) & wind_speed > 0,
                           .(wind_direction=weighted.mean.circular(wd,
                                                                   w=wind_speed,
                                                                   template='geographics',
                                                                   units='degrees')),
                           .(STID, date_time_hour)]
# ^ time_series_data_dir only contains STID, date_time_hour, & wind_direction in degrees

time_series_data_dir[wind_direction < 0, wind_direction := wind_direction + 360]

time_series_data_out <- full_join(time_series_data_scalar, time_series_data_dir)

head(time_series_data_out)

# #### Pull Out Precipitation TimeSeries Data for Sites Near Salton Sea ####
# ## full timeseries variable list: https://demos.synopticdata.com/variables/index.html
#
# i.STID <- sites$STID[1] # checking the index (1)
#
# precip_time_series_data <- lapply(sites$STID, function(i.STID){
#   dt <- mw(service = 'timeseries', stid =i.STID,
#            vars = c('precip_accum_24_hour'),
#            start = '202005130001', end = '202112090001', jsonsimplify= TRUE)
#
#   obs <- dt[['STATION']][['OBSERVATIONS']]
#
#   obs <- lapply(obs, unlist) %>%
#     as.data.table
#
#   obs[, STID:=i.STID]
#
#   return(obs)
# }) %>% rbindlist(fill=TRUE)
#
# setnames(precip_time_series_data, names(precip_time_series_data) %>% str_remove_all('_set_1(d?)'))
# # ^ remove '_set_1(d?)' from column names in precip_time_series_data
#
# precip_time_series_data[, date_time:=ymd_hms(date_time)]
# precip_time_series_data[, date_time_hour:=ceiling_date(date_time, 'hour')]
# # precip_time_series_data[, DateLST:=ymd_hms(date_time, tz = 'America/Los_Angeles')]
# # precip_time_series_data[, DateLST:=ceiling_date(DateLST, 'hour')]
#
# head(precip_time_series_data)
#
# # note: [,.()] is data.table notation; just keeps data in data.table format
# # https://stackoverflow.com/questions/33808705/dot-preceding-parentheses-in-data-table
# # looks like we are taking the mean of air temp & wind speed to get averages per hour rather than values per 5 minutes
# precip_time_series_data_scalar <- precip_time_series_data[, .(precip_24_accum=mean(precip_accum_24_hour,na.rm=TRUE)),
#                                             .(STID, date_time_hour)]
# # ^ need to add accumulated precip & relative humidity in here, want to make sure I am treating the data the same way as Will
#
# precip_time_series_data_dir <- precip_time_series_data[!is.na(wd) & !is.na(wind_speed) & wind_speed > 0,
#                                          .(wind_direction=weighted.mean.circular(wd,
#                                                                                  w=wind_speed,
#                                                                                  template='geographics',
#                                                                                  units='degrees')),
#                                          .(STID, date_precip_time_hour)]
# # ^ precip_time_series_data_dir only contains STID, date_precip_time_hour, & wind_direction in degrees
#
# precip_time_series_data_dir[wind_direction < 0, wind_direction := wind_direction + 360]
#
# precip_time_series_data_out <- full_join(precip_time_series_data_scalar, precip_time_series_data_dir)
#
#
#
#### Parse TimeSeries Data ####

# input for this function includes two data frames: the stations_df, and the timeseries_df
# the stations_df should be a dataframe containing the following information:
##STID column with station IDs, Deploy_dth which is a column of dust collector deployment date/times, and Collect_dth which is a column of dust collector retrieval date/times
head(stids) # we are going to use this df, stids, as our stations_df input
# the timeseries_df is the data frame containing the timeseries_df, and we are using the date_time_hour column to parse out the timeseries data
# if the timeseries_df does not contain the column date_time_hour, the function will go to the error message.

SubsetTimeSeries<-function(stations_df,timeseries_df){
  timeseries.aves.list<-list() # create empty list where we will store averages
  if ("STID" %in% colnames(stations_df) & "Deploy_dth" %in% colnames(stations_df) & "Collect_dth" %in% colnames(stations_df)){
    # ^ if specific column headers are found in stations_df, will proceed to loop
    for (i in seq_len(nrow(stations_df))){
      row<-stations_df[i,]
      STIDname<-row$STID # column with station ID info
      DeployDT<-row$Deploy_dth # dust collector deployment date-time info
      CollDT<-row$Collect_dth # dust collector collection date-time info
      ColNum<-row$CollectionNum # dust collection #

      if ("date_time_hour" %in% colnames(timeseries_df)){
        # below we subset the timeseries data frame by looking between our Deployment date/time & Collection date/time, and by station ID
        df<-subset(timeseries_df, date_time_hour >= ymd_hms(DeployDT) & date_time_hour <= ymd_hms(CollDT) & STID==STIDname)

        df.nums<-select_if(df, is.numeric) # from subsetted time series data, which columns are numeric are subsetted
        df.ave1<-data.frame(ave=(t(colMeans(df.nums,na.rm=TRUE)))) # take the column means of the time series data that is numeric, and add the "ave" tag before the column headers
        df.ave<-cbind(row,df.ave1) # combine the row data from our station df and our average time series data we just calculated

        assign(paste0(STIDname,".Coll",ColNum), df,envir = .GlobalEnv) # save each subsetted time series data in global env
        assign(paste0(STIDname,".Coll",ColNum,".Ave"), df.ave,envir = .GlobalEnv) # save averages of subsetted time series data in global env

        timeseries.aves.list[[i]]<-df.ave # add the average of subsetted time series data to the list we made above
        names(timeseries.aves.list)[[i]]<-paste0(STIDname,".Coll",ColNum,".Ave") # name each element in list with the station ID, collection #, and "Ave" so you know its time series averages

      } else {
        next
      }

    }

  } else {
    print("Please ensure your 1st input data frame has columns with the following names: STID, Deploy_dth, and Collect_dth")
    print("STID contains station IDs, Deploy_dth contains dust collector deployment date/time, and Collect_dth contains dust collection retrieval date/time")
    print("Also, please ensure your 2nd input data frame has column with the following name: date_time_hour")
  }
  assign(paste0("timeseries.aves.list"), timeseries.aves.list,envir = .GlobalEnv) # assign the list of subsetted time series averages to global Env
  timeseries.ave = do.call(rbind, timeseries.aves.list)
  assign(paste0("timeseries.aves"), timeseries.ave,envir = .GlobalEnv) # assign the list of subsetted time series averages to global Env
}

SubsetTimeSeries(stids,time_series_data_out) # this works
#SubsetTimeSeries(sites,sites) # this is a test to show what happens if you do not have the necessary column headers in your 1st input df

# collapse the list of time series averages into 1 df
#timeseries.ave = do.call(rbind, timeseries.aves.list)

### here is the loop I made first before building the function; this is for subsetting time_series data using input from metadata
# for (i in seq_len(nrow(stids))){
#   row<-stids[i,]
#   STIDname<-row$STID # column with station ID info
#   DeployDT<-row$Deploy_dth # dust collector deployment date-time info
#   CollDT<-row$Collect_dth # dust collector collection date-time info
#   ColNum<-row$CollectionNum # dust collection #
#
#   df<-subset(time_series_data_out, date_time_hour >= ymd_hms(DeployDT) & date_time_hour <= ymd_hms(CollDT) & STID==STIDname)
#   #print(ymd_hms(DeployDT))
#   #print(ymd_hms(CollDT))
#   #print(ColNum)
#
#   #assign(paste0(STIDname,".Coll",ColNum), df,envir = .GlobalEnv)
#   nums <- unlist(lapply(df, is.numeric), use.names = FALSE)
#
#   df
#   #df.ave1<-colMeans(df[df.nums],na.rm=TRUE)
#   #df.ave<-rbind(STIDname,DeployDT,CollDT,df.ave1)
#
#   #assign(paste0(STIDname,".Coll",ColNum), df,envir = .GlobalEnv)
#   #assign(paste0(STIDname,".Coll",ColNum,".Ave"), df.ave,envir = .GlobalEnv)
#
#
# }


### below is the code I was using to figure out how to subset before creating my function
# # UCDE - for BDC
# UCDE.Coll1<-subset(time_series_data_out, date_time_hour >= ymd_hms("2020-05-13 12:00:00") & date_time_hour <= ymd_hms("2020-07-09 17:00:00") & STID=="UCDE")
# UCDE.Coll2<-subset(time_series_data_out, date_time_hour >= ymd_hms("2020-07-09 12:00:00") & date_time_hour <= ymd_hms("2020-08-14 17:00:00") & STID=="UCDE")
# UCDE.Coll3<-subset(time_series_data_out, date_time_hour >= ymd_hms("2020-08-14 12:00:00") & date_time_hour <= ymd_hms("2020-10-08 17:00:00") & STID=="UCDE")
# UCDE.Coll4<-subset(time_series_data_out, date_time_hour >= ymd_hms("2020-10-08 12:00:00") & date_time_hour <= ymd_hms("2020-11-06 17:00:00") & STID=="UCDE")
# UCDE.Coll5<-subset(time_series_data_out, date_time_hour >= ymd_hms("2021-06-05 12:00:00") & date_time_hour <= ymd_hms("2021-07-27 17:00:00") & STID=="UCDE")
# UCDE.Coll6<-subset(time_series_data_out, date_time_hour >= ymd_hms("2021-07-27 12:00:00") & date_time_hour <= ymd_hms("2021-09-29 17:00:00") & STID=="UCDE")
# UCDE.Coll7<-subset(time_series_data_out, date_time_hour >= ymd_hms("2021-09-29 12:00:00") & date_time_hour <= ymd_hms("2021-12-08 17:00:00") & STID=="UCDE")
#
# #CI200 - for PD
# CI200.Coll1<-subset(time_series_data_out, date_time_hour >= ymd_hms("2020-05-13 12:00:00") & date_time_hour <= ymd_hms("2020-07-09 17:00:00") & STID=="CI200")
# CI200.Coll2<-subset(time_series_data_out, date_time_hour >= ymd_hms("2020-07-09 12:00:00") & date_time_hour <= ymd_hms("2020-08-14 17:00:00") & STID=="CI200")
# CI200.Coll3<-subset(time_series_data_out, date_time_hour >= ymd_hms("2020-08-14 12:00:00") & date_time_hour <= ymd_hms("2020-10-08 17:00:00") & STID=="CI200")
# CI200.Coll4<-subset(time_series_data_out, date_time_hour >= ymd_hms("2020-10-08 12:00:00") & date_time_hour <= ymd_hms("2020-11-06 17:00:00") & STID=="CI200")
# CI200.Coll5<-subset(time_series_data_out, date_time_hour >= ymd_hms("2021-06-05 12:00:00") & date_time_hour <= ymd_hms("2021-07-27 17:00:00") & STID=="CI200")
# CI200.Coll6<-subset(time_series_data_out, date_time_hour >= ymd_hms("2021-07-27 12:00:00") & date_time_hour <= ymd_hms("2021-09-18 17:00:00") & STID=="CI200")
# CI200.Coll7<-subset(time_series_data_out, date_time_hour >= ymd_hms("2021-09-18 12:00:00") & date_time_hour <= ymd_hms("2021-12-08 17:00:00") & STID=="CI200")
#
# # CQ125 - for WI
# CQ125.Coll1<-subset(time_series_data_out, date_time_hour >= ymd_hms("2020-06-01 12:00:00") & date_time_hour <= ymd_hms("2020-07-10 17:00:00") & STID=="CQ125")
# CQ125.Coll2<-subset(time_series_data_out, date_time_hour >= ymd_hms("2020-07-10 12:00:00") & date_time_hour <= ymd_hms("2020-08-30 17:00:00") & STID=="CQ125")
# CQ125.Coll3<-subset(time_series_data_out, date_time_hour >= ymd_hms("2020-08-30 12:00:00") & date_time_hour <= ymd_hms("2020-10-10 17:00:00") & STID=="CQ125")
# CQ125.Coll4<-subset(time_series_data_out, date_time_hour >= ymd_hms("2020-10-10 12:00:00") & date_time_hour <= ymd_hms("2020-11-05 17:00:00") & STID=="CQ125")
# CQ125.Coll5<-subset(time_series_data_out, date_time_hour >= ymd_hms("2021-06-08 12:00:00") & date_time_hour <= ymd_hms("2021-07-29 17:00:00") & STID=="CQ125")
# CQ125.Coll6<-subset(time_series_data_out, date_time_hour >= ymd_hms("2021-07-29 12:00:00") & date_time_hour <= ymd_hms("2021-09-18 17:00:00") & STID=="CQ125")
# CQ125.Coll7<-subset(time_series_data_out, date_time_hour >= ymd_hms("2021-09-18 12:00:00") & date_time_hour <= ymd_hms("2021-12-08 17:00:00") & STID=="CQ125")
#
# # DPMC1 - for DP
# DPMC1.Coll1<-subset(time_series_data_out, date_time_hour >= ymd_hms("2020-06-01 12:00:00") & date_time_hour <= ymd_hms("2020-07-10 17:00:00") & STID=="DPMC1")
# DPMC1.Coll2<-subset(time_series_data_out, date_time_hour >= ymd_hms("2020-07-10 12:00:00") & date_time_hour <= ymd_hms("2020-08-30 17:00:00") & STID=="DPMC1")
# DPMC1.Coll3<-subset(time_series_data_out, date_time_hour >= ymd_hms("2020-08-30 12:00:00") & date_time_hour <= ymd_hms("2020-10-10 17:00:00") & STID=="DPMC1")
# DPMC1.Coll4<-subset(time_series_data_out, date_time_hour >= ymd_hms("2020-10-10 12:00:00") & date_time_hour <= ymd_hms("2020-11-05 17:00:00") & STID=="DPMC1")
# DPMC1.Coll5<-subset(time_series_data_out, date_time_hour >= ymd_hms("2021-06-08 12:00:00") & date_time_hour <= ymd_hms("2021-08-19 17:00:00") & STID=="DPMC1")
# DPMC1.Coll6<-subset(time_series_data_out, date_time_hour >= ymd_hms("2021-08-19 12:00:00") & date_time_hour <= ymd_hms("2021-09-18 17:00:00") & STID=="DPMC1")
# DPMC1.Coll7<-subset(time_series_data_out, date_time_hour >= ymd_hms("2021-09-18 12:00:00") & date_time_hour <= ymd_hms("2021-12-08 17:00:00") & STID=="DPMC1")
#

#### Save Output ####
head(timeseries.aves)

#dust_meta_clim<-merge(dust_meta,timeseries.aves,by=c("CollectionNum","STID","Deploy_dth","Collect_dth"))
#head(dust_meta_clim)

saveRDS(timeseries.aves, file = "data/Climate/SaltonSea_SynopticClimateData_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)
save.image("data/Climate/SSD_Synoptic_ClimateData_All.Rdata")

write.table(timeseries.aves, "data/Climate/SSD_Synoptic_ClimateData_All.tsv",
            sep="\t", quote=F, col.names=T,row.names=F)
