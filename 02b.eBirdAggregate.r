
# Lane Scher
# clanescher@gmail.com
# 06/08/19

# This script aggregates eBird checklists to match BBS routes

library(rgeos)
library(lubridate)
library(reshape2)
library(dplyr)
library(data.table)

bufferDist <- 5000

# #### read in files ####
stateInfo <- read.csv("DATA/dataOther/stateInfo.csv")
load("DATA/dataBBS/processed/dates.rdata")
load("DATA/dataBBS/processed/routesFinal.rdata")

bbsRoutes <- rgdal::readOGR(
  dsn = paste0("DATA/dataBBS/processed/routesBufferShapefile") , 
  layer = "BBSbuffer5k",
  verbose = T
)
bbsRoutes <- sf::st_as_sf(bbsRoutes)


if (region == "SE") {bcr <- c(29)}
if (region == "MW") {bcr <- c(19)}
if (region == "PN") {bcr <- c(5)}

bbsRoutes1 <- bbsRoutes[which(bbsRoutes$BCR_REG == bcr),]

#states <- statesUse

dataEBI <- list()
dataEBI_CL <- list()


for (y in 1:length(years)) {
  year <- years[y]
  print(year)
  
  #### read in ebird checklists ####
  load(paste0("DATA/dataEbird/processed/", year, "-", region, "-CL.rdata"))
  load(paste0("DATA/dataEbird/processed/", year, "-", region, ".rdata"))
  
  # set up dataframes that will be filled in
  keep_CL <- as.data.frame(matrix(nrow = 0, ncol = 1))
  colnames(keep_CL) <- c("geeID")
  keep_CL$geeID <- as.character(keep_CL$geeID)
  
  keep_count <- as.data.frame(matrix(nrow = 0, ncol = 1))
  colnames(keep_count) <- "checklist_id"
  keep_count$checklist_id <- as.character(keep_count$checklist_id) 
  
  
  #### get rid of checklists we don't want ####
  
  # take out checklists with Xs
  X_CL <- ebd_df[which(ebd_df$observation_count == "X"),]
  X_CL <- unique(X_CL$checklist_id)
  ebd_CL <- ebd_CL[-which(ebd_CL$checklist_id %in% X_CL),]
  
  
  # take out checklists without effort
  ebd_CL <- ebd_CL[which(is.na(ebd_CL$duration_minutes) == F),]
  ebd_CL <- ebd_CL[which(ebd_CL$duration_minutes > 0),]
  
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
  
  ############################################
  #                                          #
  #   assign each checklist to a BBS route   #
  #                                          #
  ############################################
  
  if(nrow(ebd_CL) > 0){
    
    
    
    #### make checklists spatial points class ####
    ebd_CL_sf <- ebd_CL[,c("checklist_id", "longitude", "latitude")]
    ebd_CL_sf <- unique(ebd_CL_sf) # should already be unique
    
    ebd_CL_sf <- sf::st_as_sf(ebd_CL_sf, coords = c("longitude", 'latitude'), crs = 4326)
    
    
    ebd_CL_sf <- sf::st_join(ebd_CL_sf, bbsRoutes1)
    
    ebd_CL_match <- ebd_CL_sf[,c("checklist_id", "geeID")]
    ebd_CL_match <- ebd_CL_match[which(is.na(ebd_CL_match$geeID) == F),]
    
    ebd_CL <- ebd_CL[which(ebd_CL$checklist_id %in% ebd_CL_match$checklist_id),]
    ebd_CL <- merge(ebd_CL, ebd_CL_match[,c("checklist_id", 'geeID')], by = "checklist_id")
    
    # some checklists are within the buffer for two routes, choose one
    ebd_CL <- ebd_CL[match(unique(ebd_CL$checklist_id), ebd_CL$checklist_id),]
    
    
    
    #### keep only the correct dates ####
    
    # add start and end date to stateDF based on geeID and year
    ebd_CL <- merge(ebd_CL, BBSdates, by = c("geeID", "year"))
    
    # only keep checklists that match the dates
    ebd_CL <- ebd_CL[which(ebd_CL$observation_date < ebd_CL$end &
                             ebd_CL$observation_date > ebd_CL$start),]
    
    ebd_CL <- ebd_CL[,c("geeID", "bcr_code", "checklist_id",
                        "observer_id", "year", 
                        "duration_minutes", "effort_area_ha", "effort_distance_km",
                        "latitude", "longitude", 
                        "protocol_type", "time_observations_started", 
                        "observation_date")]
    
    
    
    #### now get counts from those checklists ####
    
    ebd_counts <- ebd_df[which(ebd_df$checklist_id %in% ebd_CL$checklist_id),]
    ebd_counts$observation_count <- as.numeric(as.character(ebd_counts$observation_count))
    
    ebd_counts <- reshape2::dcast(ebd_counts,
                                  checklist_id ~ scientific_name,
                                  value.var = "observation_count",
                                  FUN = "sum")
  }
  
  ### add to lists
  
  dataEBI[[y]] <- ebd_counts
  dataEBI_CL[[y]] <- ebd_CL
  
}

dataEBI <- data.table::rbindlist(dataEBI, fill = T, use.names = T)
dataEBI_CL <- data.table::rbindlist(dataEBI_CL, fill = T, use.names = T)


save(dataEBI,
     dataEBI_CL,
     file = paste0("DATA/dataEbird/processed/ydata-", region, ".rdata"))
