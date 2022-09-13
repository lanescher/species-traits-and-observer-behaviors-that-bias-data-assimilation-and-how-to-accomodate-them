# Lane Scher
# clanescher@gmail.com
# 06/24/2020

# This script processes BBS data and returns a dataframe in
# the correct format for ydata in GJAM


library(reshape2)


#### read in raw data ####
dpath <- 'DATA/dataBBS/raw/counts'
files <- list.files(dpath)

data <- read.csv(paste0(dpath, "/", files[1]))
colnames(data) <- tolower(colnames(data))

for (j in 2:length(files)) {
  print(j)
  tmp <- read.csv(paste0(dpath, "/", files[j]))
  colnames(tmp) <- tolower(colnames(tmp))
  data <- rbind(data, tmp)
}

# load migrants
bbsMigrants <- read.csv("X:/data-BBS/Feb-2021/MigrantNonBreeder/MigrantNonBreeder/Migrants/Migrants.csv")

# combine migrants with the others
colnames(bbsMigrants) <- tolower(colnames(bbsMigrants))
data <- bind_rows(data, bbsMigrants)
data$geeID <- paste0(data$countrynum, "_", data$statenum, "_", data$route)



# get BCR
load("DATA/dataBBS/processed/routesFinal.rdata")



#### read in AOU codes ####
spList <- read.fwf(file = "DATA/dataBBS/raw/SpeciesList.txt",
                   widths = c(7, 5, 50, 50, 50, 
                              50, 50, 50, 50),
                   skip = 11,
                   col.names = c("Seq", "AOU", "English",
                                 "French", "Spanish",
                                 "Order", "Family", "Genus", "Species"),
                   strip.white = T)

spList$gs <- paste(spList$Genus, spList$Species, sep = " ")
spList <- spList[,c("AOU", "gs")]




regions <- c("MW", "PN", "SE")

for (r in 1:length(regions)) {
  reg <- regions[r]
  
  # if (reg == "MW") {codes <- c(38, 54, 67)}
  # if (reg == "PN") {codes <- c(69, 89)}
  # if (reg == "SE") {codes <- c(63, 80, 27)}
  
  
  if (reg == "MW") {codes <- 19}
  if (reg == "PN") {codes <- 5}
  if (reg == "SE") {codes <- 29}
  
  
  #### keep only correct BCR ####
  keep <- routesBBS2[which(routesBBS2$BCR_REGION == codes),]
  data1 <- data[which(data$geeID %in% keep$geeID),]
  
  #### keep only US ####
  data1 <- data1[which(data1$countrynum == 840),]
  
  #### keep only correct year ####
  data1 <- data1[which(data1$year %in% years),]

  #### format ####
  dataM <- reshape2::melt(data1, id.vars = c("routedataid", "countrynum", "statenum",
                                            "route", "geeID", "rpid", "year", "aou"))
  
  dataM$variable <- gsub("stop", "", dataM$variable)
  dataM$routedataidP <- paste0(dataM$routedataid, ".", dataM$variable)
  
  #### replace AOU with species name ####
  dataM <- dataM[,c("routedataidP", "aou", "value")]
  
  dataM <- merge(dataM, spList, by.x = "aou", by.y = "AOU")
  dataM <- dataM[,c("routedataidP", "gs", "value")]
  dataM$count <- as.numeric(dataM$value)
  
  # combine species that were reported with multiple AOU
  dataM <- aggregate(dataM[,4], by = list(dataM$routedataidP, dataM$gs), FUN = sum)
  
  #### format for gjam ####
  dataBBS <- reshape2::dcast(dataM, Group.1 ~ Group.2, value.var = "x")
  colnames(dataBBS)[1] <- "geeID"
  
  #### save ####
  save(dataBBS, file = paste0("DATA/dataBBS/processed/ydata-", reg, ".rdata"))
  
}


