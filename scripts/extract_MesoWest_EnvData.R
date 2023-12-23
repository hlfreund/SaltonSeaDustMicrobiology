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

})

#### Log Into MesoWest w/ API key ####
# mesowest package GitHub: https://github.com/fickse/mesowest

requestToken(apikey = "KunTxSCUqjGlu9CBtKJJRmxjYHnFc6PK4uPoDgcUOg")

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

## Then filter the sites by the following info:
### distance between site + Salton Sea coordinates is <= 40km
### site has data between time points of dust collections that we are looking at
sites <- sites %>%
  filter(dist <= 40000) %>%
  filter(start <= as.Date('2020-06-01')) %>%
  filter(end >= as.Date('2021-12-31'))

write_rds(sites, 'data/MesoWest_SaltonSea_2020-2021_SiteData.rds')

#### Pull Out TimeSeries Data for Sites Near Salton Sea ####
## full timeseries variable list: https://demos.synopticdata.com/variables/index.html

i.STID <- sites$STID[1]

metdata <- lapply(sites$STID, function(i.STID){
  dt <- mw(service = 'timeseries', stid =i.STID,
           vars = c('wind_speed','wind_gust','wind_direction',
                    'wind_cardinal_direction', 'air_temp','relative_humidity','precip_accum',
                    'PM_25_concentration','incoming_radiation_uv','PM_10_concentration'),
           start = '202006010001', end = '202201010001', jsonsimplify= TRUE)

  obs <- dt[['STATION']][['OBSERVATIONS']]

  obs <- lapply(obs, unlist) %>%
    as.data.table

  obs[, STID:=i.STID]

  return(obs)
}) %>% rbindlist(fill=TRUE)

setnames(metdata, names(metdata) %>% str_remove_all('_set_1(d?)'))

metdata[, date_time:=ymd_hms(date_time)]
metdata[, date_time_hour:=ceiling_date(date_time, 'hour')]
# metdata[, DateLST:=ymd_hms(date_time, tz = 'America/Los_Angeles')]
# metdata[, DateLST:=ceiling_date(DateLST, 'hour')]

metdata[, wd:=circular(wind_direction, template='geographics', units='degrees')]

metdata_scalar <- metdata[, .(air_temp = mean(air_temp, na.rm=TRUE),
                              wind_speed=mean(wind_speed, na.rm=TRUE),
                              wind_gust=max(wind_gust, na.rm=TRUE)),
                          .(STID, date_time_hour)]

metdata_dir <- metdata[!is.na(wd) & !is.na(wind_speed) & wind_speed > 0,
                           .(wind_direction=weighted.mean.circular(wd,
                                                                   w=wind_speed,
                                                                   template='geographics',
                                                                   units='degrees')),
                           .(STID, date_time_hour)]
metdata_dir[wind_direction < 0, wind_direction := wind_direction + 360]

metdata_out <- full_join(metdata_scalar, metdata_dir)

#### Pull Out Precipitation Data for Sites Near Salton Sea ####
## info on precipitation data can be found here: https://docs.synopticdata.com/services/precipitation
mw(service = 'precipitation', stid='UP612', start =  '202007010001', end ='202201010001', pmode='intervals', interval = 'hour', jsonsimplify = TRUE, returnURL =FALSE)


write_rds(metdata_out, 'data/MesoWest_SaltonSea_2013-2020.rds')
