# Lane Scher
# clanescher@gmail.com
# 06/24/2020

# This script processes BBS dates to produce start and end
# dates making a four-week window around the survey date

library(dplyr)

#### read in data ####
BBSdates <- read.csv("X:/eBird-BBS/DATA/dataBBS/raw/weather.csv")

#### clean up data ####
BBSdates <- BBSdates[which(BBSdates$RunType == 1),]

BBSdates$geeID <- paste(BBSdates$CountryNum, BBSdates$StateNum, BBSdates$Route, sep = "_")
BBSdates$date <- paste(BBSdates$Year, BBSdates$Month, BBSdates$Day, sep = "-")
BBSdates$date <- as.Date(BBSdates$date)

BBSdates <- BBSdates[,c("geeID", "Year", "date")]

BBSdates <- BBSdates[which(BBSdates$Year > 2009),]

# if a route-year is recorded twice, use the median date
BBSdates1 <- aggregate(BBSdates$date,
                       by = list(BBSdates$geeID, BBSdates$Year),
                       FUN = "median")

colnames(BBSdates1) <- c("geeID", "year", "median")

# make start and end dates of window
BBSdates1$start <- BBSdates1$median - 14
BBSdates1$end <- BBSdates1$median + 14


#### save ####
BBSdates <- BBSdates1
save(BBSdates, file = "DATA/dataBBS/processed/dates.rdata")
