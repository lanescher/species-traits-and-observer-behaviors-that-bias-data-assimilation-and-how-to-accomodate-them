# Lane Scher
# clanescher@gmail.com
# 05/18/2020

# This script processes BBS observer data

weather <- read.csv("DATA/dataBBS/raw/weather.csv")
weather <- weather[which(weather$RunType == 1),]
weather$geeID <- paste0(weather$CountryNum, "_", weather$StateNum, "_", weather$Route)


# get BCR
load("DATA/dataBBS/processed/routesFinal.rdata")


regions <- c("MW", "PN", "SE")

for (r in 1:length(regions)) {
  reg <- regions[r]
  
  
  if (reg == "MW") {codes <- 19}
  if (reg == "PN") {codes <- 5}
  if (reg == "SE") {codes <- 29}
  
  
  #### keep only correct BCR ####
  keep <- routesBBS2[which(routesBBS2$BCR_REGION == codes),]
  weather1 <- weather[which(weather$geeID %in% keep$geeID),]
  
  #### keep only US ####
  weather1 <- weather1[which(weather1$CountryNum == 840),]
  
  
  #### keep only correct year ####
  weather1 <- weather1[which(weather1$Year %in% years),]
  
  
  weather1$geeIDyear <- paste(weather1$CountryNum, weather1$StateNum, weather1$Route,
                              weather1$Year, sep = "_")
  weather1$geeIDyearObs <- paste(weather1$geeIDyear, weather1$ObsN, sep = "_")
  
  
  obs <- weather1[,c("geeIDyear", "ObsN")]
  
  colnames(obs) <- c("geeIDyear", "observer")
  
  save(obs, file = paste0("DATA/dataBBS/processed/observers-", reg, ".rdata"))
  
}

