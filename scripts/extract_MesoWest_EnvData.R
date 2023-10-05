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
# MesoWest variables: https://developers.synopticdata.com/mesonet/v2/api-variables/
# air_temp (ºC)
# wind_speed (m/s)
# wind_gust (m/s)
# wind_direction (degrees)
# wind_cardinal_direction (char)
# mesowest package GitHub: https://github.com/fickse/mesowest

requestToken(apikey = "KuRvN1RwhvIS6gHsnK5HunnVSIlvMFV7nW3I4AlmMp")

#### Retrieve Data for Sites of Interest ####

# Salton Sea: 33.3286° N, 115.8434° W
seaCoords <- c(-115.8434, 33.3286)

# Get all MesoWest sites in Imperial and Riverside counties
sites <- bind_rows(
  mw(service = 'metadata', complete=1, state='CA', county='Riverside' )[[1]],
  mw(service = 'metadata', complete=1, state='CA', county='Imperial' )[[1]]) %>%
  mutate(LATITUDE = as.numeric(LATITUDE)) %>%
  mutate(LONGITUDE = as.numeric(LONGITUDE)) %>%
  mutate(start = ymd_hms(PERIOD_OF_RECORD$start)) %>%
  mutate(end = ymd_hms(PERIOD_OF_RECORD$end))


# Filter for sites within 40 km of Salton Sea
sites <- sites %>%
  mutate(dist = distm(sites[, c('LONGITUDE', 'LATITUDE')], seaCoords, fun = distHaversine))

sites <- sites %>%
  filter(dist <= 40000) %>%
  filter(start <= as.Date('2019-01-01')) %>%
  filter(end >= as.Date('2020-12-31'))

write_rds(sites, 'RData/MesoWest_SaltonSea_2013-2020_SiteData.rds')


i.STID <- sites$STID[1]
metdata <- lapply(sites$STID, function(i.STID){
  dt <- mw(service = 'timeseries', stid =i.STID,
           vars = c('wind_speed','wind_gust','wind_direction',
                    'wind_cardinal_direction', 'air_temp'),
           start = '202101010001', end = '202101020001', jsonsimplify= TRUE)

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

write_rds(metdata_out, 'data/MesoWest_SaltonSea_2013-2020.rds')
