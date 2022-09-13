
# Lane Scher
# clanescher@gmail.com
# 06/24/2020

# Reads in BBS route shape file, adds buffer, and saves 
# bufferred routes as a shape file. These polygons are 
# used to download environmental covariate data



############################################
#                                          #
#               load data                  #
#                                          #
############################################

#### libraries ####
library(rgdal)
library(raster)
library(rgeos)
library(sf)

#### read in data ####

# route shapefiles
setwd("DATA/dataBBS/raw/routesShapefile")
routesBBS <- readOGR(dsn = ".", layer = "bbs_geog83")
# downloaded from: http://www.pasda.psu.edu/uci/DataSummary.aspx?dataset=599
setwd(wd)

routesBBS$ROUTE_ <- as.character(routesBBS$ROUTE_)
routesBBS$OBJECTID <- as.character(routesBBS$OBJECTID)

# start locations of BBS routes, this has route metadata and geeID
startBBS <- read.csv("DATA/dataBBS/raw/routes.csv")
startBBS$geeID <- paste0(startBBS$CountryNum, "_", startBBS$StateNum, "_", startBBS$Route)
startBBS <- startBBS[,c(2, 12, 4, 6:7)]

# BBS state codes
statesBBS <- read.csv("DATA/dataBBS/raw/BBSstateCodes.csv")

# state data (neighbors and UTM)
neighborsList <- read.csv("DATA/dataOther/stateInfo.csv")


############################################
#                                          #
#          clean up route paths            #
#                                          #
############################################


# take out routes longer than 30km and shorter than 50km
routesBBS <- routesBBS[which(routesBBS@data$LENGTH > 0.30 &
                               routesBBS@data$LENGTH < 0.50),]


dups <- as.data.frame(table(routesBBS$ROUTE_))
dups <- dups[which(dups$Freq > 1),]

### There are still 11 duplicates. 
### Four of these duplicates cover the same paths but with slightly
### different locations. Manually remove one of each of these duplicates
## 61101, 61120, 69147, 76136 can be removed

## 61101
loc <- grep("61101", routesBBS$ROUTE_)
plot(routesBBS@lines[[loc[1]]]@Lines[[1]]@coords)
plot(routesBBS@lines[[loc[2]]]@Lines[[1]]@coords)
# second one is much lower resolution, remove it
rm <- routesBBS$OBJECTID[loc[2]]
routesBBS <- routesBBS[-which(routesBBS$OBJECTID == rm),]


loc <- grep("61120", routesBBS$ROUTE_)
plot(routesBBS@lines[[loc[1]]]@Lines[[1]]@coords)
plot(routesBBS@lines[[loc[2]]]@Lines[[1]]@coords)
# second one is much lower resolution, remove it
rm <- routesBBS$OBJECTID[loc[2]]
routesBBS <- routesBBS[-which(routesBBS$OBJECTID == rm),]


loc <- grep("69147", routesBBS$ROUTE_)
plot(routesBBS@lines[[loc[1]]]@Lines[[1]]@coords)
plot(routesBBS@lines[[loc[2]]]@Lines[[1]]@coords)
# second one is much lower resolution, remove it
rm <- routesBBS$OBJECTID[loc[2]]
routesBBS <- routesBBS[-which(routesBBS$OBJECTID == rm),]


loc <- grep("76136", routesBBS$ROUTE_)
plot(routesBBS@lines[[loc[1]]]@Lines[[1]]@coords)
plot(routesBBS@lines[[loc[2]]]@Lines[[1]]@coords)
# second one is much lower resolution, remove it
rm <- routesBBS$OBJECTID[loc[2]]
routesBBS <- routesBBS[-which(routesBBS$OBJECTID == rm),]


dups <- as.data.frame(table(routesBBS$ROUTE_))
dups <- dups[which(dups$Freq > 1),]

### now there are 7 duplicates. No way to tell which of these is
### correct, so remove both of them!

rm <- as.character(dups$Var1)

routesBBS <- routesBBS[-which(routesBBS$ROUTE_ %in% rm),]


##### each row is a unique route now.
##### now need to give it its geeID. This is in startBBS, although
##### the two dataframes don't share any of the same columns to merge by


# add Lat and Long of start location to routesBBS@data. That will be
# used to merge routesBBS with startBBS
for (l in 1:length(routesBBS@lines)) {
  routesBBS$Latitude[l] <- routesBBS@lines[[l]]@Lines[[1]]@coords[[1,2]]
  routesBBS$Longitude[l] <- routesBBS@lines[[l]]@Lines[[1]]@coords[[1,1]]
}

# save routesBBS as routesBBS1
routesBBS1 <- routesBBS
routes <- routesBBS@data


############################################
#                                          #
#        merge path and metadeta           #
#                                          #
############################################


# merge routes and startBBS by route name
tmp1 <- merge(routes, startBBS, by.x = "ROUTENAME", by.y = "RouteName",
              all = T)

# multiple routes have the same name, so only keep rows where 
# lat and long also match
tmp2 <- tmp1[which(abs(tmp1$Latitude.x - tmp1$Latitude.y) < 1 &
                     abs(tmp1$Longitude.x - tmp1$Longitude.y) < 1),]

#tmp2$routeGEEID <- paste0(tmp2$ROUTENAME, "-", tmp2$geeID)

### there are still some duplicated geeIDs
dups <- as.data.frame(table(tmp2$geeID))
dups <- dups[which(dups$Freq > 1),]

### get rid of them
tmp3 <- tmp2[-which(tmp2$geeID %in% dups$Var1),]

### there are still some duplicated OBJECTIDs
dups <- as.data.frame(table(tmp3$OBJECTID))
dups <- dups[which(dups$Freq > 1),]

# don't know which is the correct geeID for these, so get rid of them
tmp4 <- tmp3[-which(tmp3$OBJECTID %in% dups$Var1),]

tmp4 <- tmp4[,c("ROUTE_", "StateNum", "geeID")]
tmp4 <- merge(tmp4, statesBBS, 
              by.x = "StateNum", by.y = "code")
tmp4$stateNoSpace <- gsub(" ", "", tmp4$state)

# merge tmp4 with routesBBS@data based on ROUTE_
routesBBS2 <- merge(routesBBS, tmp4, by = "ROUTE_", all = F)



# save routes -- these are used to aggregate eBird data
save(routesBBS2, file = "DATA/dataBBS/processed/routesFinal.rdata")



############################################
#                                          #
#               add buffers                #
#                                          #
############################################

#### These will be used to get environmental covariates from GEE

#### add buffers ####

# make state vector
states <- toupper(neighborsList$state)

# start loop with first state
state <- states[1]
stateCode <- statesBBS$code[grep(state, statesBBS$state)]

# keep routes in state 1
test <- routesBBS2[which(routesBBS2@data$StateNum == stateCode),]

# select correct UTM zone for state 1
zone <- neighborsList$utm[which(toupper(neighborsList$state) == state)]
crsUTM <- paste0("+proj=utm +zone=", zone, " +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

# give routes that UTM
test <- spTransform(test, CRS(crsUTM))

# add buffer
test2 <- gBuffer(test, byid = T, width = 5000)

# put back in WGS84 for GEE
test2 <- spTransform(test2, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))

BBSbuffer5k <- test2
for (s in 2:length(states)){
  # select state
  state <- states[s]
  stateCode <- statesBBS$code[grep(state, statesBBS$state)]
  
  if (state == "VIRGINIA") {stateCode <- stateCode[1]}
  if (state == "KANSAS") {stateCode <- stateCode[2]}

  if (state != "HAWAII"){
    # keep routes in that state
    test <- routesBBS2[which(routesBBS2@data$StateNum == stateCode),]
    
    # select UTM zone
    zone <- neighborsList$utm[which(toupper(neighborsList$state) == state)]
    crsUTM <- paste0("+proj=utm +zone=", zone, " +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
    
    # give routes that UTM
    test <- spTransform(test, CRS(crsUTM))
    
    # add buffer
    test2 <- gBuffer(test, byid = T, width = 5000)
    
    # put back in WGS84 for GEE
    test2 <- spTransform(test2, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
    
    # add this state to the previous states
    BBSbuffer5k <- rbind(BBSbuffer5k, test2)
  }
  
  
}

setwd("DATA/dataBBS/processed/routesBufferShapefile")
writeOGR(obj = BBSbuffer5k, dsn = ".", driver = "ESRI Shapefile", layer = "BBSbuffer5k")
setwd(wd)

