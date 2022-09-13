

# Lane Scher
# clanescher@gmail.com
# 06/24/2020

# This script combines ebird and BBS observations with
# covariates. The output is a dataframe with all observations
# and covariates. 



############################################
#                                          #
#              read in data                #
#                                          #
############################################

#### libraries ####
library(gjam)
library(dplyr)
library(lubridate)
library(data.table)
library(suncalc)
library(lutz)

if (region == "MW") {codes <- c(19)}
if (region == "PN") {codes <- c(5)}
if (region == "SE") {codes <- c(29)}



#### read in BBS data ####
# BBS counts
load(paste0("DATA/dataBBS/processed/ydata-", region, ".rdata"))

# BBS observer names
load(paste0("DATA/dataBBS/processed/observers-", region, ".rdata"))

# BBS lat lon - filter by BCR
startBBS <- read.csv("DATA/dataBBS/raw/routes.csv")
startBBS <- startBBS[which(startBBS$BCR %in% codes),]
startBBS$geeID <- paste0(startBBS$CountryNum, "_", startBBS$StateNum, "_", startBBS$Route)

# clip BCR to certain states
if (region == "SE") {startBBS <- startBBS[which(startBBS$StateNum %in% c(88, 63, 80, 27, 2)),]}
if (region == "PN") {startBBS <- startBBS[which(startBBS$StateNum %in% c(89, 69)),]}




# BBS observation times 
times <- read.csv("DATA/dataBBS/raw/weather.csv")
times$geeID <- paste0(times$CountryNum, "_", times$StateNum, "_", times$Route)
times <- times[which(times$geeID %in% startBBS$geeID),]



#### read in ebird data ####
readFile <- paste0("DATA/dataEbird/processed/ydata-", region, ".rdata")
load(readFile)
dataEBI_CL$StateNum <- substr(dataEBI_CL$geeID, 5, nchar(dataEBI_CL))
dataEBI_CL$StateNum <- gsub("_.*", '', dataEBI_CL$StateNum)


if (region == "SE") {dataEBI_CL <- dataEBI_CL[which(dataEBI_CL$StateNum %in% c(88, 63, 80, 27, 2)),]}
if (region == "PN") {dataEBI_CL <- dataEBI_CL[which(dataEBI_CL$StateNum %in% c(89, 69)),]}


#### read in species names and common names ####
traits <- read.csv("DATA/dataOther/birdTraits3.csv")


#### make list to store count into ####
countInfo <- list()


############################################
#                                          #
#         clean up eBird CL data           #
#                                          #
############################################

#### get ebird checklist covariates ####
# add source
dataEBI_CL$source <- "EBI"

# add geeIDyear
dataEBI_CL$geeIDyear <- paste0(dataEBI_CL$geeID, "_", dataEBI_CL$year)

#### add effort distance for stationary CL ####
dataEBI_CL$effort_distance_km[which(dataEBI_CL$protocol == "Stationary")] <- 0

#### add time since sunrise ####
# find sunrise and sunset time
ebdSunrise <- dataEBI_CL[,c("observation_date", "latitude", "longitude")]
colnames(ebdSunrise) <- c("date", "lat", "lon")
ebdSunrise$date <- as.Date(ebdSunrise$date)
tmp <- getSunlightTimes(data = ebdSunrise,
                        keep = c("sunrise", "sunset"))
dataEBI_CL$sunrise <- tmp$sunrise
dataEBI_CL$sunset <- tmp$sunset

# make time of observation in UTC
dataEBI_CL$timeString <- paste0(as.character(dataEBI_CL$observation_date),
                                              " ", dataEBI_CL$time)
dataEBI_CL$tz <- tz_lookup_coords(dataEBI_CL$latitude,
                                  dataEBI_CL$longitude,
                                   method = "accurate",
                                   warn = T)
dataEBI_CL$timeString <- as.POSIXct(dataEBI_CL$timeString,
                                     format = "%Y-%m-%d %H:%M:%S")
dataEBI_CL$timeUTC <- force_tzs(dataEBI_CL$timeString,
                                 tzones = dataEBI_CL$tz)

# calculate time since sunrise and sunset
dataEBI_CL$timeSinceSunrise <- as.numeric(difftime(dataEBI_CL$timeUTC, dataEBI_CL$sunrise,
                                         units = "mins"))
dataEBI_CL$timeSinceSunset <- as.numeric(difftime(dataEBI_CL$timeUTC, dataEBI_CL$sunset,
                                        units = "mins"))
dataEBI_CL$timeSinceSunrise2 <- dataEBI_CL$timeSinceSunrise + 0.5*(dataEBI_CL$duration_minutes)
dataEBI_CL$timeSinceSunset2 <- dataEBI_CL$timeSinceSunset + 0.5*(dataEBI_CL$duration_minutes)


# add time category
dataEBI_CL$timeCat <- "night"
dataEBI_CL$timeCat[which(dataEBI_CL$timeSinceSunrise2 > -60 &
                           dataEBI_CL$timeSinceSunrise2 <= 60)] <- "sunrise"
dataEBI_CL$timeCat[which(dataEBI_CL$timeSinceSunrise2 > 60 &
                           dataEBI_CL$timeSinceSunrise2 <= 360)] <- "aamorning"
dataEBI_CL$timeCat[which(dataEBI_CL$timeSinceSunrise2 > 360 &
                           dataEBI_CL$timeSinceSunset2 <= -180)] <- "afternoon"
dataEBI_CL$timeCat[which(dataEBI_CL$timeSinceSunset2 > -180 &
                            dataEBI_CL$timeSinceSunset2 <= 60)] <- "evening"



ebdCovar <- data.frame(obsID = dataEBI_CL$checklist_id,
                       observer = dataEBI_CL$observer_id,
                       effKm = dataEBI_CL$effort_distance_km,
                       effMin = dataEBI_CL$duration_minutes,
                       time = dataEBI_CL$timeSinceSunrise2,
                       timeCat = dataEBI_CL$timeCat,
                       source = dataEBI_CL$source,
                       geeIDyear = dataEBI_CL$geeIDyear)
ebdCovar <- ebdCovar[which(ebdCovar$effMin > 1),] # if effort = 1 min, log(eff) = 0 which is bad for gjam


#### count up effort before and after spatiotemporal filter ####

# load full checklist info for checklist counts
ebiCLcounts <- data.frame(year = years,
                          EBIall = NA,
                          EBIsubset = NA)

for (y in 1:length(years)) {
  load(paste0("DATA/dataEbird/processed/", years[y], "-", region, "-CL.rdata"))
  load(paste0("DATA/dataEbird/processed/", years[y], "-", region, ".rdata"))
  
  # clip
  if (region == "SE") {ebd_CL <- ebd_CL[which(ebd_CL$state %in% c("Virginia", "North Carolina", "South Carolina", "Georgia", "Alabama")),]}
  if (region == "PN") {ebd_CL <- ebd_CL[which(ebd_CL$state %in% c("Oregon", "Washington")),]}
  
  ebd_df <- ebd_df[which(ebd_df$checklist_id %in% ebd_CL$checklist_id),]
  
  #### apply consistency filters ####

  # take out checklists with Xs
  X_CL <- ebd_df[which(ebd_df$observation_count == "X"),]
  X_CL <- unique(X_CL[,c("checklist_id", "state")])
  ebd_CL <- ebd_CL[-which(ebd_CL$checklist_id %in% X_CL$checklist_id),]

  # take out checklists without effort
  ebd_CL <- ebd_CL[which(is.na(ebd_CL$duration_minutes) == F),]
  ebd_CL <- ebd_CL[which(ebd_CL$duration_minutes > 1),] # if effort is 1, log(effort) is 0 which doesn't work

  # take out checklists that are probably BBS routes
  bbsPC <- ebd_CL[which(ebd_CL$duration_minutes == 3),]
  bbsPC1 <- bbsPC[,c("checklist_id", "observer_id", "observation_date")]
  bbsPC1 <- unique(bbsPC1)

  rm <- as.data.frame(table(bbsPC1$observer_id, bbsPC1$observation_date))
  rm$obsDate <- paste0(rm$Var1, "-", rm$Var2)
  rm <- rm$obsDate[which(rm$Freq > 40 & rm$Freq < 60)]

  ebd_CL$obsDate <- paste0(ebd_CL$observer_id, "-", ebd_CL$observation_date)
  ebd_CL <- subset(ebd_CL, !(ebd_CL$obsDate %in% rm))

  # take out group checklists
  group_CL <- ebd_CL$checklist_id[which(substr(ebd_CL$checklist_id, 1, 1) == "G")]
  ebd_CL <- ebd_CL[-which(ebd_CL$checklist_id %in% group_CL),]


  ebiCLcounts[y,2] <- length(unique(ebd_CL$checklist_id))
}


############################################
#                                          #
#          clean up BBS PC data            #
#                                          #
############################################

#### get BBS point count covariates ####
# keep BBS routes with RunType = 1
keepBBS <- times$RouteDataID[which(times$RunType == 1 &
                                     times$Year > 2009)]

dataBBS$RouteDataID <- gsub("[.].*", "", dataBBS$geeID)
dataBBS <- dataBBS[which(dataBBS$RouteDataID %in% keepBBS),]

#### make BBS covariate df ####
dataBBS_PC <- data.frame(RouteDataID = dataBBS$RouteDataID,
                        eff_min = 3,
                        eff_km = 0,
                        source = "BBS")
dataBBS_PC <- unique(dataBBS_PC)
#dataBBS_PC$geeIDyes <- paste0(dataBBS_PC$geeID, "_", dataBBS_PC$)

# add geeIDyear to dataBBS_PC
times <- times[which(times$RouteDataID %in% keepBBS),]

times$geeIDyear <- paste(times$CountryNum,
                         times$StateNum,
                         times$Route,
                         times$Year, sep = "_")
dataBBS_PC <- merge(dataBBS_PC, times[,c("geeID", "RouteDataID", "RunType", "Year")], by = "RouteDataID")


#### get time from sunrise for each point along route ####
times$startTime <- paste0(substr(times$StartTime, 1, nchar(times$StartTime)-2),
                          ":", 
                          substr(times$StartTime, nchar(times$StartTime)-1, nchar(times$StartTime)))
times$endTime <- paste0(substr(times$EndTime, 1, nchar(times$EndTime)-2),
                        ":", 
                        substr(times$EndTime, nchar(times$EndTime)-1, nchar(times$EndTime)))
times$startTime[which(nchar(times$startTime) == 4)] <-
  paste0("0", times$startTime[which(nchar(times$startTime) == 4)])
times$endTime[which(nchar(times$endTime) == 4)] <-
  paste0("0", times$endTime[which(nchar(times$endTime) == 4)])
times$Month[which(nchar(times$Month) == 1)] <-
  paste0("0", times$Month[which(nchar(times$Month) == 1)])
times$Day[which(nchar(times$Day) == 1)] <-
  paste0("0", times$Day[which(nchar(times$Day) == 1)])

times$start <- paste0(times$Year, "-", times$Month, "-", times$Day,
                     " ", times$startTime)
times$end <- paste0(times$Year, "-", times$Month, "-", times$Day,
                     " ", times$endTime)
times$start1 <- as.POSIXct(times$start,
                          format = "%Y-%m-%d %H:%M")
times$end1 <- as.POSIXct(times$end,
                           format = "%Y-%m-%d %H:%M")

# add lat and lon
times$geeID <- paste0(times$CountryNum, "_", times$StateNum,
                      "_", times$Route)
times <- merge(times, startBBS[,c("Latitude", "Longitude",
                                  "geeID")], by = "geeID")

# put in UTC
times$tz <- tz_lookup_coords(times$Latitude,
                             times$Longitude,
                             method = "accurate",
                             warn = T)
times$startUTC <- force_tzs(times$start1,
                            tzones = times$tz)
times$endUTC <- force_tzs(times$end1,
                          tzones = times$tz)

# rearrange
times$total <- difftime(times$endUTC, times$startUTC, units = "mins")
times1 <- times[,c("geeIDyear", "total", "startUTC")]

# find time at each stop
times1$each <- times1$total/50

# time 1 at start time
times1$`1` <- times1$startUTC
# Each subsequent stop is (total time / 50) minutes later
for (i in 1:49) {
  times1[,5+i] <- times1[,4+i] + times1$each
}
colnames(times1)[5:54] <- paste0("stop", 1:50)

# reformat
timesM <- reshape2::melt(times1[,c(1,5:54)], id.vars = "geeIDyear")
timesM$value <- as_datetime(timesM$value)

timesM$obsID <- paste0(timesM$geeIDyear, ".",
                       gsub("stop", "", timesM$variable))
timesM <- timesM[,c("obsID", "value")]
colnames(timesM) <- c("obsID", "time")

timesM$geeIDyear <- gsub("[.].*", "", timesM$obsID)


# find sunrise and sunset
bbsSunrise <- times[,c("start1", "Latitude", "Longitude")]
colnames(bbsSunrise) <- c("date", "lat", "lon")
bbsSunrise$date <- (substr(bbsSunrise$date, 1, 10))
bbsSunrise$date <- as.Date(bbsSunrise$date)
tmp <- getSunlightTimes(data = bbsSunrise,
                        keep = c("sunrise", "sunset"))
times$sunrise <- tmp$sunrise
times$sunset <- tmp$sunset

timesM <- merge(timesM, times[,c("geeIDyear", "sunrise", "sunset")],
                by = "geeIDyear")

# calculate time since sunrise and sunset
timesM$timeSinceSunrise <- as.numeric(difftime(timesM$time, timesM$sunrise,
                                                    units = "mins"))
timesM$timeSinceSunset <- as.numeric(difftime(timesM$time, timesM$sunset,
                                                   units = "mins"))
timesM$timeSinceSunrise2 <- timesM$timeSinceSunrise + 1.5
timesM$timeSinceSunset2 <- timesM$timeSinceSunset + 1.5


# add time category
timesM$timeCat <- "night"
timesM$timeCat[which(timesM$timeSinceSunrise2 > -60 &
                       timesM$timeSinceSunrise2 <= 60)] <- "sunrise"
timesM$timeCat[which(timesM$timeSinceSunrise2 > 60 &
                       timesM$timeSinceSunrise2 <= 360)] <- "aamorning"
timesM$timeCat[which(timesM$timeSinceSunrise2 > 360 &
                       timesM$timeSinceSunset2 <= -180)] <- "afternoon"
timesM$timeCat[which(timesM$timeSinceSunset2 > -180 &
                       timesM$timeSinceSunset2 <= 60)] <- "evening"
#timesM$timeCat <- as.factor(timesM$timeCat)

# add time since sunrise to df
dataBBS_PC$geeIDyear <- paste0(dataBBS_PC$geeID, "_", dataBBS_PC$Year)

dataBBS_PC <- merge(dataBBS_PC, timesM[,c("geeIDyear", "obsID", "timeSinceSunrise2",
                                          "timeCat")], 
                   by = "geeIDyear")

#### add BBS observers ####
dataBBS_PC <- merge(dataBBS_PC, obs, by = "geeIDyear")

#### fix obsID ####
dataBBS_PC$obsID2 <- paste0(dataBBS_PC$RouteDataID, ".",
                           gsub(".*[.]", "", dataBBS_PC$obsID))

bbsCovar <- data.frame(obsID = dataBBS_PC$obsID2,
                       observer = dataBBS_PC$observer,
                       effKm = 0,
                       effMin = 3,
                       time = dataBBS_PC$timeSinceSunrise2,
                       timeCat = dataBBS_PC$timeCat,
                       source = "BBS",
                       geeIDyear = dataBBS_PC$geeIDyear)





#### count up effort before spatiotemporal filter ####

bbsPCcounts <- data.frame(year = years,
                          BBSall = NA,
                          BBSsubset = NA)

tmp <- substr(bbsCovar$geeIDyear, nchar(bbsCovar$geeIDyear)-3, nchar(bbsCovar$geeIDyear))
bbsPCcounts$BBSall <- as.vector(table(tmp))



############################################
#                                          #
#        subset by space and time          #
#                                          #
############################################


#### only keep obs that are in the covariate data ####
dataBBS <- dataBBS[which(dataBBS$geeID %in% bbsCovar$obsID),]
dataEBI <- dataEBI[which(dataEBI$checklist_id %in% ebdCovar$obsID),]

#### make sure to only keep observations that have matching data ####
bbsLocs <- as.character(unique(bbsCovar$geeIDyear))
ebiLocs <- as.character(unique(ebdCovar$geeIDyear))

sameLocs <- intersect(bbsLocs, ebiLocs)

bbsKeep <- as.character(bbsCovar$obsID[which(bbsCovar$geeIDyear %in% sameLocs)])
ebiKeep <- as.character(ebdCovar$obsID[which(ebdCovar$geeIDyear %in% sameLocs)])

dataBBS <- dataBBS[which(dataBBS$geeID %in% bbsKeep),]
dataEBI <- dataEBI[which(dataEBI$checklist_id %in% ebiKeep),]

bbsCovar <- bbsCovar[which(bbsCovar$obsID %in% bbsKeep),]
ebdCovar <- ebdCovar[which(ebdCovar$obsID %in% ebiKeep),]


#### now take out afternoon, evening, and night ####
# but first record counts during each time period
ebdTime <- as.data.frame(table(ebdCovar$timeCat))
bbsTime <- as.data.frame(table(bbsCovar$timeCat))

timePeriods <- merge(ebdTime, bbsTime, by = "Var1", all = T)
colnames(timePeriods) <- c("period", "ebd", "bbs")

countInfo[[2]] <- timePeriods

# now get rid of afternoon, evening, and night
bbsCovar <- bbsCovar[which(bbsCovar$timeCat %in% c("sunrise",
                                                   "aamorning")),]

ebdCovar <- ebdCovar[which(ebdCovar$timeCat %in% c("sunrise",
                                                   "aamorning")),]


# and get rid of ydata from afternoon, evening, and night
dataBBS <- dataBBS[which(dataBBS$geeID %in% bbsCovar$obsID),]
dataEBI <- dataEBI[which(dataEBI$checklist_id %in% ebdCovar$obsID),]


##### now make sure all observations still have matches after removing
##### afternoon, evening, and night
bbsLocs <- as.character(unique(bbsCovar$geeIDyear))
ebiLocs <- as.character(unique(ebdCovar$geeIDyear))

sameLocs <- intersect(bbsLocs, ebiLocs)

bbsKeep <- as.character(bbsCovar$obsID[which(bbsCovar$geeIDyear %in% sameLocs)])
ebiKeep <- as.character(ebdCovar$obsID[which(ebdCovar$geeIDyear %in% sameLocs)])

dataBBS <- dataBBS[which(dataBBS$geeID %in% bbsKeep),]
dataEBI <- dataEBI[which(dataEBI$checklist_id %in% ebiKeep),]

bbsCovar <- bbsCovar[which(bbsCovar$obsID %in% bbsKeep),]
ebdCovar <- ebdCovar[which(ebdCovar$obsID %in% ebiKeep),]



#### add to effort df ####
tmp <- substr(bbsCovar$geeIDyear, nchar(bbsCovar$geeIDyear)-3, nchar(bbsCovar$geeIDyear))
bbsPCcounts$BBSsubset <- as.vector(table(tmp))

tmp <- substr(ebdCovar$geeIDyear, nchar(ebdCovar$geeIDyear)-3, nchar(ebdCovar$geeIDyear))
ebiCLcounts$EBIsubset <- as.vector(table(tmp))


# combine dfs and save
annualEffort <- merge(bbsPCcounts, ebiCLcounts, by = "year")

countInfo[[1]] <- annualEffort






############################################
#                                          #
#              combine ydata               #
#                                          #
############################################

#### clean up ####
class(dataEBI) <- "data.frame"
dataBBS$RouteDataID <- NULL

# dataBBS1 <- dataBBS
# dataEBI1 <- dataEBI
#### fix column names ####

# change to obsID
colnames(dataBBS)[1] <- "obsID"
colnames(dataEBI)[1] <- "obsID"


#### keep only observations that are in bbsCovar and ebdCovar ####
dataBBS <- dataBBS[which(dataBBS$obsID %in% bbsCovar$obsID),]
dataEBI <- dataEBI[which(dataEBI$obsID %in% ebdCovar$obsID),]


# take out colnames with sp.
bbsrm <- grep("sp[.]", colnames(dataBBS))
if (length(bbsrm) > 0) {
  dataBBS <- dataBBS[,-bbsrm]
}

ebirm <- grep("sp[.]", colnames(dataEBI))
if (length(ebirm) > 0) {
  dataEBI <- dataEBI[,-ebirm]
}


# take out colnames with /
bbsrm <- grep("/", colnames(dataBBS))
if (length(bbsrm) > 0) {
  dataBBS <- dataBBS[,-bbsrm]
}

ebirm <- grep("/", colnames(dataEBI))
if (length(ebirm) > 0) {
  dataEBI <- dataEBI[,-ebirm]
}

# take out colnames with x
bbsrm <- grep(" x ", colnames(dataBBS))
if (length(bbsrm) > 0) {
  dataBBS <- dataBBS[,-bbsrm]
}

ebirm <- grep(" x ", colnames(dataEBI))
if (length(ebirm) > 0) {
  dataEBI <- dataEBI[,-ebirm]
}

#### change colnames to common names BBS ####

# first, combine great blue herons
gbh <- grep("Ardea herodias", colnames(dataBBS))
if (length(gbh) == 2) {
  dataBBS$`Ardea herodias` <- dataBBS$`Ardea herodias` + dataBBS$`Ardea herodias occidentalis`
  dataBBS$`Ardea herodias occidentalis` <- NULL
}

# then combine northern flickers
nf <- grep("Colaptes auratus", colnames(dataBBS))
if (length(nf) > 1) {
  dataBBS$`Colaptes auratus` <- rowSums(dataBBS[,nf], na.rm = T)
  
  dataBBS[,nf[2:length(nf)]] <- NULL
}

# then combine dark eyed juncos
dej <- grep("Junco hyemalis", colnames(dataBBS)) 
if (length(dej) > 1) {
  dataBBS$`Junco hyemalis` <- rowSums(dataBBS[,dej], na.rm = T)
  
  dataBBS[,dej[2:length(dej)]] <- NULL
}


# now change names
sciToCom <- function(df, 
                     sci = traits$altScientific, # the default refers to
                     com = traits$nameCom        # the bird trait table, but
                                                 # these can be replaced with
                                                 # any dataframe
) {
  spec <- colnames(df)
  comNew <- c()

  for (i in 1:length(spec)) {
    sp <- spec[i]

    com1 <- grep(sp, sci)

    if (length(com1) > 0) {
      comNew[i] <- as.character(com[com1])
    } else {             #### if the species isn't in the df,
      comNew[i] <- sp    #### keeps the scientific name
    }
  }

  return(comNew)
}
tmp <- sciToCom(dataBBS[,2:ncol(dataBBS)])
colnames(dataBBS)[2:ncol(dataBBS)] <- tmp


#### change colnames to common names eBird ####

tmp <- sciToCom(dataEBI[,2:ncol(dataEBI)])
colnames(dataEBI)[2:ncol(dataEBI)] <- tmp


#### look for duplicates ####
length(unique(colnames(dataBBS)))
ncol(dataBBS)

length(unique(colnames(dataEBI)))
ncol(dataEBI)

#### take out species that were never observed ####

ebdrm <- colSums(dataEBI[,2:ncol(dataEBI)], na.rm = T)
tmp <- colnames(dataEBI)[2:ncol(dataEBI)]
ebdrm <- tmp[which(is.na(ebdrm) == T |
                     ebdrm == 0)]

dataEBI <- dataEBI[,-which(colnames(dataEBI) %in% ebdrm)]


bbsrm <- colSums(dataBBS[,2:(ncol(dataBBS)-1)], na.rm = T)
tmp <- colnames(dataBBS)[2:(ncol(dataBBS)-1)]
bbsrm <- tmp[which(is.na(bbsrm) == T |
                     bbsrm == 0)]

dataBBS <- dataBBS[,-which(colnames(dataBBS) %in% bbsrm)]

spObs <- list(ebd = colnames(dataEBI[2:ncol(dataEBI)]),
              bbs = colnames(dataBBS[2:ncol(dataBBS)]))

countInfo[[3]] <- spObs


#### count number of observations in ebd and bbs where each species was seen ####
spAll <- unique(c(spObs[[1]], spObs[[2]]))

counts <- data.frame(species = spAll,
                     ebi = NA,
                     bbs = NA)

for (s in 1:nrow(counts)) {
  sp <- counts$species[s]
  
  tmp <- dataEBI[,grep(sp, colnames(dataEBI))]
  counts$ebi[grep(sp, counts$species)] <- length(tmp[which(tmp > 0)])
  
  tmp <- dataBBS[,grep(sp, colnames(dataBBS))]
  if (class(tmp) != "numeric") {
    tmp$total <- rowSums(tmp, na.rm = T)
    counts$bbs[grep(sp, counts$species)] <- length(tmp[which(tmp$total > 0),3])
    
  } else {
    counts$bbs[grep(sp, counts$species)] <- length(tmp[which(tmp > 0)])
    
  }
  
}

countInfo[[4]] <- counts

############################################
#                                          #
#             combine data                 #
#                                          #
############################################


#### only keep matching species ####
keepSp <- intersect(colnames(dataBBS), colnames(dataEBI))
dataBBS <- dataBBS[,which(colnames(dataBBS) %in% keepSp)]
dataEBI <- dataEBI[,which(colnames(dataEBI) %in% keepSp)]


#### combine ####
ydata <- bind_rows(dataEBI, dataBBS)
ydata[is.na(ydata)] <- 0

ebdCovar$observer <- as.character(ebdCovar$observer)
bbsCovar$observer <- as.character(bbsCovar$observer)
xdata <- bind_rows(ebdCovar, bbsCovar)


save(xdata, ydata, file = paste0("DATA/combined-", region, ".rdata"))
save(countInfo, file = paste0("DATA/stats-", region, ".rdata"))
