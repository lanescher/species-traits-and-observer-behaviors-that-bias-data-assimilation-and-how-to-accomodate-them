
# Lane Scher
# clanescher@gmail.com
# 03/25/2021

library(rgdal)
library(sf)
library(raster)
library(reshape2)
library(tidyverse)


devtools::install_github("ropensci/FedData")


freqNLCD_stack <- function(points, raster, buffer, max = F, numPix = F,
                           type = "nlcd") {
  
  .processNLCD <- function(df) {
    num_pix <- nrow(df)
    
    .tabProp <- function(values) {
      n <- length(values)
      tab1 <- table(values)
      tab1 <- as.data.frame(tab1/n)
      
      rownames(tab1) <- tab1$values
      tab1$values <- NULL
      
      tab1 <- as.data.frame(t(tab1))
      
      if (type == "nlcd") {
        lcNums <- c("11", "12", "21", "22", "23", "24", "31", "41", "42", "43",
                    "51", "52", "71", "81", "82", "90", "95")
      }
      
      if (type == "foresce") {
        lcNums <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                    "11", "12", "13", "14", "15", "16")
      }
      
      lcNumsMissing <- lcNums[which(lcNums %in% as.character(colnames(tab1)) == F)]
      if (length(lcNumsMissing > 0)) {
        for (i in 1:length(lcNumsMissing)){
          tab1[lcNumsMissing[i]] <- 0
        }
      }
      
      tab1 <- tab1[,lcNums]
      
      return(tab1)
    }
    tab <- apply(df, 2, .tabProp)
    
    #tab_filled <- lapply(tab, .fillMissingLC)
    nms <- names(tab)
    
    tab <- bind_rows(tab)
    
    tab$lc <- colnames(tab)[max.col(tab, ties.method = "random")]
    tab$year <- substr(nms, nchar(nms)-3, nchar(nms))
    tab$numPix <- num_pix
    
    if (type == "nlcd") {
      allFreq <- tab[,c("year", "numPix", "lc",
                        "11", "12", "21", "22", "23", "24", "31", "41", "42", "43",
                        "51", "52", "71", "81", "82", "90", "95")]
      
    }
    if (type == "foresce") {
      allFreq <- tab[,c("year", "numPix", "lc",
                        "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                        "11", "12", "13", "14", "15", "16")]
    }
    
    
    return(allFreq)
  }
  
  
  start <- Sys.time()
  tmp1 <- raster::extract(raster, points)
  
  
  test1 <- lapply(tmp1, .processNLCD)
  
  nlcd_out <- bind_rows(test1)
  
  # add locality_ids to nlcd_out
  names <- points$geeID
  
  nlcd_out <- cbind(locality_id = rep(names, each = raster::nlayers(raster)),
                    nlcd_out)
  
  if (max == F) {nlcd_out$max <- NULL}
  if (numPix == F) {nlcd_out$numPix <- NULL}
  
  return(nlcd_out)
  
  
  print(Sys.time() - start)
  
}


bbsRoutes <- rgdal::readOGR(
  dsn = paste0("DATA/dataBBS/processed/routesBufferShapefile") , 
  layer = "BBSbuffer5k",
  verbose = T
)

region <- "all"



# need to get previous winter too
years1 <- c(years[1]-1, years)


# get locations
locs <- bbsRoutes[which(bbsRoutes$BCR_REG %in% c(29, 19, 5)),]

# some states we are excluding
locs <- locs[-which(locs$state %in% c('ALASKA', 'MARYLAND', 'NEW JERSEY', 'PENNSYLVANIA')),]


names(locs@polygons) <- locs$geeID

#### check if data has already been downloaded
envD <- list.files("DATA/dataDaymet1/")
envDname <- paste0(region, "-weatherData.rdata")

if (envDname %in% envD) {
  
  load(paste0("DATA/dataDaymet1/", region, "-weatherData.rdata"))
  
} else {
  
  startTime <- Sys.time()
  
  #### pull monthly value for each BBS route
  prcpAll <- data.frame()
  tminAll <- data.frame()
  tmaxAll <- data.frame()
  vpAll <- data.frame()
  
  for (i in 1:length(locs)) {
    print(Sys.time() - startTime)
    print(paste0("Location ", i))
    # download data
    locsRaw <- FedData::get_daymet(template = locs[i,], 
                                   label = "geeID", 
                                   elements = c("prcp", 
                                                "tmin", 
                                                "tmax",
                                                "vp"),
                                   years = years1,
                                   tempo = "mon",
                                   extraction.dir = paste0("DATA/dataDaymet1/", region, "processed"))
    
    
    rasters_proj <- list()
    for (r in 1:length(locsRaw)) {
      rast <- raster::projectRaster(locsRaw[[r]], crs = sp::CRS('+init=epsg:4326'))
      rasters_proj[[r]] <- rast 
    }
    
    
    for (y in 1:length(years1)) {

      
      yearPrcp <- rasters_proj[1][[1]]
      yearPrcp <- yearPrcp[[grep(years1[y], names(yearPrcp))]]
      
      yearTmin <- rasters_proj[[2]]
      yearTmin <- yearTmin[[grep(years1[y], names(yearTmin))]]
      
      yearTmax <- rasters_proj[[3]]
      yearTmax <- yearTmax[[grep(years1[y], names(yearTmax))]]
      
      yearVP <- rasters_proj[[4]]
      yearVP <- yearVP[[grep(years1[y], names(yearVP))]]
      
      prcpMons <- data.frame()
      tminMons <- data.frame()
      tmaxMons <- data.frame()
      vpMons <- data.frame()
      for (m in 1:(12)) {
        
        # precipitation
        rast <- yearPrcp[[m]]

        # get mean prcp per month
        prcpMon <- raster::extract(rast, locs[i,])
        prcpMon <- lapply(prcpMon, FUN = mean, na.rm = T)
        prcpMon <- unlist(prcpMon)
        
        prcp <- data.frame(geeID = locs$geeID[i],
                           year = years1[y],
                           month = m,
                           prcp = prcpMon)
        prcpMons <- dplyr::bind_rows(prcpMons, prcp)
        
        #### tmin ####
        rast <- yearTmin[[m]]

        # get mean prcp per month
        tminMon <- raster::extract(rast, locs[i,])
        tminMon <- lapply(tminMon, FUN = mean, na.rm = T)
        tminMon <- unlist(tminMon)
        
        tmin <- data.frame(geeID = locs$geeID[i],
                           year = years1[y],
                           month = m,
                           tmin = tminMon)
        tminMons <- dplyr::bind_rows(tminMons, tmin)
        
        #### tmax ####
        rast <- yearTmax[[m]]

        # get mean prcp per month
        tmaxMon <- raster::extract(rast, locs[i,])
        tmaxMon <- lapply(tmaxMon, FUN = mean, na.rm = T)
        tmaxMon <- unlist(tmaxMon)
        
        tmax <- data.frame(geeID = locs$geeID[i],
                           year = years1[y],
                           month = m,
                           tmax = tmaxMon)
        tmaxMons <- dplyr::bind_rows(tmaxMons, tmax)
        
        #### vp ####
        rast <- yearVP[[m]]

        # get mean prcp per month
        vpMon <- raster::extract(rast, locs[i,])
        vpMon <- lapply(vpMon, FUN = mean, na.rm = T)
        vpMon <- unlist(vpMon)
        
        vp <- data.frame(geeID = locs$geeID[i],
                         year = years1[y],
                         month = m,
                         vp = vpMon)
        vpMons <- dplyr::bind_rows(vpMons, vp)
      }
      prcpAll <- dplyr::bind_rows(prcpAll, prcpMons)
      tminAll <- dplyr::bind_rows(tminAll, tminMons)
      vpAll <- dplyr::bind_rows(vpAll, vpMons)
      tmaxAll <- dplyr::bind_rows(tmaxAll, tmaxMons)
      
    }
    
    save(prcpAll, tminAll, vpAll, tmaxAll, file = paste0("DATA/dataDaymet1/", region, "-weatherData.rdata"))
  }
  
}


#### now calculate covariates

### annual prcp
prcpAnn <- prcpAll

prcpAnn$yearAnn <- prcpAnn$year
prcpAnn$yearAnn[which(prcpAnn$year %in% c("6", "7", "8", "9", "10", "11", "12"))] <-
  prcpAnn$year[which(prcpAnn$year %in% c("6", "7", "8", "9", "10", "11", "12"))] + 1

prcpAnn <- aggregate(prcpAnn[,c("prcp")], 
                     by = list(prcpAnn$geeID,
                               prcpAnn$yearAnn),
                     FUN = "sum")
colnames(prcpAnn) <- c("geeID", "year", "prcpAnnual")

### spring prcp
prcpSpr <- prcpAll[which(prcpAll$month %in% c("3", "4", "5")),]
prcpSpr <- aggregate(prcpSpr[,c("prcp")],
                     by = list(prcpSpr$geeID,
                               prcpSpr$year),
                     FUN = "sum")
colnames(prcpSpr) <- c("geeID", "year", "prcpSpring")

### winter prcp
prcpWin <- prcpAll[which(prcpAll$month %in% c("12", "1", "2")),]


prcpWin$yearWin <- prcpWin$year
prcpWin$yearWin[which(prcpWin$month %in% c("12"))] <- prcpWin$year[which(prcpWin$month %in% c("12"))] + 1

prcpWin <- aggregate(prcpWin[,c("prcp")],
                     by = list(prcpWin$geeID,
                               prcpWin$yearWin),
                     FUN = "sum")
colnames(prcpWin) <- c("geeID", "year", "prcpWinter")

### min winter temp
tminWin <- tminAll[which(tminAll$month %in% c("12", "1", "2")),]

tminWin$yearWin <- tminWin$year
tminWin$yearWin[which(tminWin$month %in% c("12"))] <- tminWin$year[which(tminWin$month %in% c("12"))] + 1

tminWin <- aggregate(tminWin[,c("tmin")],
                     by = list(tminWin$geeID,
                               tminWin$yearWin),
                     FUN = "min")
colnames(tminWin) <- c("geeID", "year", "tminWinter")

### min spring temp
tminSpr <- tminAll[which(tminAll$month %in% c("3", "4", "5")),]

tminSpr <- aggregate(tminSpr[,c("tmin")],
                     by = list(tminSpr$geeID,
                               tminSpr$year),
                     FUN = "min")
colnames(tminSpr) <- c("geeID", "year", "tminSpring")

### mean annual temp
tmeanAnn <- merge(tminAll, tmaxAll, by = c("geeID", "year", "month"))
tmeanAnn$tmean <- rowMeans(tmeanAnn[,c("tmin", "tmax")])

tmeanAnn$yearAnn <- tmeanAnn$year
tmeanAnn$yearAnn[which(tmeanAnn$year %in% c("6", "7", "8", "9", "10", "11", "12"))] <-
  tmeanAnn$year[which(tmeanAnn$year %in% c("6", "7", "8", "9", "10", "11", "12"))] + 1


tmeanAnn <- aggregate(tmeanAnn[,c("tmean")],
                      by = list(tmeanAnn$geeID,
                                tmeanAnn$yearAnn),
                      FUN = "mean")
colnames(tmeanAnn) <- c("geeID", "year", "tmeanAnn")


#### Land cover

nlcd_years <- c("2011", "2013", "2016", "2019")

locs <- bbsRoutes[which(bbsRoutes$BCR_REG %in% c(5, 19, 29)),]

### keep only continental US
sf::sf_use_s2(FALSE)
us_sf <- sf::st_as_sf(maps::map("state", plot = FALSE, fill = TRUE)) %>%
  sf::st_transform(crs = sf::st_crs(bbsRoutes)) %>%
  sf::st_buffer(dist = 0)

locs_sf <- sf::st_as_sf(locs)

locs_cont <- sf::st_join(locs_sf, us_sf)
rm_locs <- locs_cont[which(is.na(locs_cont$ID)),]

locs <- locs[-which(locs$geeID %in% rm_locs$geeID),]


### now get nlcd

if (file.exists("DATA/nlcd/final/all.rdata")) {
  load("DATA/nlcd/final/all.rdata")
} else {
  all_nlcd <- c()
  startTime <- Sys.time()
  for (r in 1:length(locs)) {
    print(r)
    print(Sys.time() - startTime)
    
    locs1 <- locs[r,]
    
    
    nlcd_stack <- raster::stack()
    for (y in 1:length(nlcd_years)) {
      year1 <- nlcd_years[y]
      print(year1)
      
      
      
      # download nlcd for this bbox
      tmp1 <- FedData::get_nlcd(template = locs1,
                                label = "test",
                                year = year1,
                                dataset = "landcover",
                                force.redo = F,
                                extraction.dir = paste0("DATA/nlcd/raw/",
                                                        locs1$geeID, "-", year1))
      
      # started at 12:32
      nlcd_stack <- raster::addLayer(nlcd_stack, tmp1)
      
      
    } # year
    
    tmp2 <- raster::projectRaster(nlcd_stack, crs="+proj=longlat +datum=WGS84",
                                  method = "ngb")
    
    
    print("Calculating NLCD")

    nlcd1 <- freqNLCD_stack(locs1, tmp2, buffer = 0, max = T, numPix = T)
    
    all_nlcd <- bind_rows(all_nlcd, nlcd1)
    
  }
  
  save(all_nlcd, file = paste0("DATA/nlcd/final/all.rdata"))
  
}


all_nlcd$water <- all_nlcd$`11` + all_nlcd$`12`
all_nlcd$developed <- all_nlcd$`21` + all_nlcd$`22` + all_nlcd$`23`
all_nlcd$forest <- all_nlcd$`41` + all_nlcd$`42` + all_nlcd$`43`
all_nlcd$shrub <- all_nlcd$`51` + all_nlcd$`52`
all_nlcd$grass <- all_nlcd$`71`
all_nlcd$wetland <- all_nlcd$`81` + all_nlcd$`82`




all_nlcd$lcYear <- all_nlcd$year
all_nlcd$geeID <- all_nlcd$locality_id


### download elevation

if (file.exists("DATA/elevation/elevation.rdata")) {
  load("DATA/elevation/elevation.rdata")
} else {
  locs <- bbsRoutes[which(bbsRoutes$BCR_REG %in% c(5, 19, 29)),]
  elev <- elevatr::get_elev_point(locs)
  elev <- data.frame(geeID = elev$geeID,
                     elev = elev$elevation)
  save(elev, file = "DATA/elevation/elevation.rdata")
}

#### squish them all together

locsData <- locs@data
locsData <- locsData[,c("geeID", "Latitud", "Longitd")]
colnames(locsData)[2:3] <- c("LAT", "LON")

covar <- merge(locsData, elev, by = c("geeID"))
covar <- merge(covar, prcpAnn, by = c("geeID"))
covar <- merge(covar, prcpSpr, by = c("geeID", "year"))
covar <- merge(covar, prcpWin, by = c("geeID", "year"))
covar <- merge(covar, tminWin, by = c("geeID", "year"))
covar <- merge(covar, tminSpr, by = c("geeID", "year"))
covar <- merge(covar, tmeanAnn, by = c("geeID", "year"))

#### add LC
covar$lcYear <- covar$year
covar$lcYear[which(covar$year %in% c(2009, 2010, 2011, 2012))] <- 2011
covar$lcYear[which(covar$year %in% c(2013, 2014))] <- 2013
covar$lcYear[which(covar$year %in% c(2015, 2016, 2017))] <- 2016
covar$lcYear[which(covar$year %in% c(2018, 2019, 2020))] <- 2019

covar <- merge(covar, all_nlcd[,c("geeID", "lcYear", "water", "developed", "forest", "shrub", "grass", "wetland")], by = c("geeID", "lcYear"))
covar$lcYear <- NULL

covar$geeIDyear <- paste0(covar$geeID, "_", covar$year)

save(covar, file = paste0("DATA/dataCovariates/covar-", region, ".rdata"))


