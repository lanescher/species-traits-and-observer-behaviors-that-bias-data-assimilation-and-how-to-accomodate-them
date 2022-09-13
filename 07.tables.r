
# Lane Scher
# clanescher@gmail.com
# 8/24/2021

# This script makes tables for the manuscript

library(dplyr)


# format beta table from gjam output
makeBetas <- function(gjamOut) {
  
  ### first do standardized betas
  betaStand <- gjamOut$parameters$betaStandXTable
  
  # make species column
  betaStand$species <- as.vector(rownames(betaStand))
  betaStand$species <- sub("_.*", "", betaStand$species)
  
  # make predictor column
  betaStand$predictor <- as.vector(rownames(betaStand))
  betaStand$predictor <- sub(".*_", "", betaStand$predictor)
  
  # add sig column
  betaStand$sig95 <- ""
  betaStand$sig95[which(betaStand$CI_025 > 0 &
                          betaStand$CI_975 > 0)] <- "*"
  betaStand$sig95[which(betaStand$CI_025 < 0 &
                          betaStand$CI_975 < 0)] <- "*"
  
  ### now do unstandardized betas
  betaUnstand <- gjamOut$parameters$betaTable
  
  # make species column
  betaUnstand$species <- as.vector(rownames(betaUnstand))
  betaUnstand$species <- sub("_.*", "", betaUnstand$species)
  
  # make predictor column
  betaUnstand$predictor <- as.vector(rownames(betaUnstand))
  betaUnstand$predictor <- sub(".*_", "", betaUnstand$predictor)
  
  # add sig column
  betaUnstand$sig95 <- ""
  betaUnstand$sig95[which(betaUnstand$CI_025 > 0 &
                            betaUnstand$CI_975 > 0)] <- "*"
  betaUnstand$sig95[which(betaUnstand$CI_025 < 0 &
                            betaUnstand$CI_975 < 0)] <- "*"
  
  ### put them into a list
  bts <- list(stand = betaStand,
              unstand = betaUnstand)
  
  return(bts)
}

############################################
#                                          #
#               trait table                #
#                                          #
############################################
regs <- c("SE", "MW", "PN")

# get SE
load(paste0("OUT/", regs[1], "-traitdata.rdata"))
colnames(traits2)[grep("obsSeen", colnames(traits2))] <- paste0("obsSeen-", regs[1])

traits2$`obsSeen-SE`[which(traits2$`obsSeen-SE` == -Inf)] <- NA

traitsAll <- traits2

# add MW
load(paste0("OUT/", regs[2], "-traitdata.rdata"))
colnames(traits2)[grep("obsSeen", colnames(traits2))] <- paste0("obsSeen-", regs[2])

traitsAll <- merge(traitsAll, traits2, all = T)

# add PN
load(paste0("OUT/", regs[3], "-traitdata.rdata"))
colnames(traits2)[grep("obsSeen", colnames(traits2))] <- paste0("obsSeen-", regs[3])

traitsAll <- merge(traitsAll, traits2, all = T)


traitsAll <- traitsAll[,c("nameCom", "nameSci", "Order",
                          "Major.Breeding.Habitats", "forHeight", "popularity", 
                          "percentVisDet",
                          "max.color.contrast", "BodyMass.Value",
                          "obsSeen-SE",
                          "obsSeen-MW", 
                          "obsSeen-PN")]

colnames(traitsAll) <- c("Common Name", "Scientific Name", "Order",
                         "Breeding habitats", "Foraging height", "Popularity",
                         "Sqrt(Percent visual detections)",
                         "Color contrast", "Log(Body mass)",
                         "Log(No. observations seen), SE",
                         "Log(No. observations seens), MW", 
                         "Log(No. observations seen), PN")

traitsAll$`Common Name` <- str_to_title(traitsAll$`Common Name`)

traitsAll[traitsAll == "NaN*"] <- NA

write.csv(traitsAll, file = "OUT/tables/traits.csv",
          row.names = F)


mn <- mean(as.numeric(traitsAll$Popularity), na.rm = T)
sd <- sd(as.numeric(traitsAll$Popularity), na.rm = T)
turkeyPop <- as.numeric(traitsAll$Popularity[which(traitsAll$`Common Name` == "Wild Turkey")])

(turkeyPop - mn)/sd

############################################
#                                          #
#              data quantity               #
#                                          #
############################################
# (table 2)


rnames <- c("NoObs", "NoObser", 
            "MeanObsObser",
            "MeanMin",
            "TotalMin",
            "MeanKm",
            "TotalKm",
            "PerStat", "PerSun", "PerMorn")
ns <- 150

dataQuantity <- data.frame(tmp = rep(NA, length(rnames)))
rownames(dataQuantity) <- rnames

regs <- c("PN", "MW", "SE")
for (r in 1:length(regs)) {
  region <- regs[r]
  
  
  
  dataQuan <- data.frame(eBird = rep(NA, length(rnames)),
                         BBS = rep(NA, length(rnames)))
  rownames(dataQuan) <- rnames
  

  source("CODE/04a.gjamPrep.r")
  
  if (region == "PN") {
    tmp <- substr(x$geeID, 5, 6)
    
    x <- x[tmp == "69",]
    y <- y[tmp == "69",]
  }
  
  xbbs <- x[x$source == "BBS",]
  xebi <- x[x$source == "EBI",]
  
  
  # number observations
  dataQuan[1, 1] <- nrow(xebi)
  dataQuan[1, 2] <- nrow(xbbs)
  
  # number observers
  dataQuan[2, 1] <- length(unique(xebi$observer))
  dataQuan[2, 2] <- length(unique(xbbs$observer))
  
  # mean and median obs/obser
  tmp <- as.data.frame(table(xebi$observer))
  tmp <- tmp[tmp$Freq > 0,]
  mn <- mean(tmp$Freq)
  std <- sd(tmp$Freq)
  dataQuan[3, 1] <- paste0(round(mn, 2), " (", round(std, 2), ")")
  #dataQuan[3, 1] <- mean(tmp$Freq)
  #dataQuan[4, 1] <- median(tmp$Freq)
  
  tmp <- as.data.frame(table(xbbs$observer))
  tmp <- tmp[tmp$Freq > 0,]
  mn <- mean(tmp$Freq)
  std <- sd(tmp$Freq)
  dataQuan[3, 2] <- paste0(round(mn, 2), " (", round(std, 2), ")")
  #dataQuan[3, 2] <- mean(tmp$Freq)
  #dataQuan[4, 2] <- median(tmp$Freq)
  
  # mean and median min
  mn <- mean(xebi$effMin)
  std <- sd(xebi$effMin)
  dataQuan[4, 1] <- paste0(round(mn, 2), " (", round(std, 2), ")")
  #dataQuan[5, 1] <- mean(xebi$effMin)
  #dataQuan[6, 1] <- median(xebi$effMin)
  
  mn <- mean(xbbs$effMin)
  std <- sd(xbbs$effMin)
  dataQuan[4, 2] <- paste0(round(mn, 2), " (", round(std, 2), ")")
  #dataQuan[5, 2] <- mean(xbbs$effMin)
  #dataQuan[6, 2] <- median(xbbs$effMin)
  
  # total minutes
  tot <- sum(xebi$effMin)
  dataQuan[5, 1] <- paste0(round(tot, 2))
  
  tot <- sum(xbbs$effMin)
  dataQuan[5, 2] <- paste0(round(tot, 2))
  
  # mean and median dist
  mn <- mean(xebi$effKm)
  std <- sd(xebi$effKm)
  dataQuan[6, 1] <- paste0(round(mn, 2), " (", round(std, 2), ")")
  #dataQuan[7, 1] <- mean(xebi$effKm)
  #dataQuan[8, 1] <- median(xebi$effKm)
  
  mn <- mean(xbbs$effKm)
  std <- sd(xbbs$effKm)
  dataQuan[6, 2] <- paste0(round(mn, 2), " (", round(std, 2), ")")
  #dataQuan[7, 2] <- mean(xbbs$effKm)
  #dataQuan[8, 2] <- median(xbbs$effKm)
  
  
  # total distance
  tot <- sum(xebi$effKm)
  dataQuan[7, 1] <- paste0(round(tot, 2))
  
  tot <- sum(xbbs$effKm)
  dataQuan[7, 2] <- paste0(round(tot, 2))
  
  
  # percent stationary
  dataQuan[8, 1] <- round(nrow(xebi[which(xebi$effKm == 0),])/nrow(xebi), 2)
  dataQuan[8, 2] <- round(nrow(xbbs[which(xbbs$effKm == 0),])/nrow(xbbs), 2)
  
  # time of day
  dataQuan[9, 1] <- round(nrow(xebi[which(xebi$timeCat == "sunrise"),])/nrow(xebi), 2)
  dataQuan[10, 1] <- round(nrow(xebi[which(xebi$timeCat == "aamorning"),])/nrow(xebi), 2)
  
  dataQuan[9, 2] <- round(nrow(xbbs[which(xbbs$timeCat == "sunrise"),])/nrow(xbbs), 2)
  dataQuan[10, 2] <- round(nrow(xbbs[which(xbbs$timeCat == "aamorning"),])/nrow(xbbs), 2)
  
  # merge with big df
  colnames(dataQuan) <- paste0(colnames(dataQuan), "-", region)
  dataQuantity <- bind_cols(dataQuantity, dataQuan)
}

dataQuantity$tmp <- NULL

rownames(dataQuantity) <- c("Number of observations",
                            "Number of observers",
                            "Observations per observer",
                            #"Median observations per observer",
                            "Minutes", "Total minutes",
                            "Distance (km)", "Total distance",
                            "Proportion stationary", 
                            "Proportion sunrise",
                            "Proportion morning")
write.csv(dataQuantity, file = "OUT/tables/dataQuantity.csv")





############################################
#                                          #
#               trait cors                 #
#                                          #
############################################

# remove species that haven't converged
sp.remove.time <- read.csv("OUT/tables/gelman_time.csv")
sp.remove.plat <- read.csv("OUT/tables/gelman_plat.csv")

traitsAll <- data.frame()
orderObs <- data.frame(Order = NA)
orderSun <- data.frame(Order = NA)
habitatObs <- data.frame(Major.Breeding.Habitats = NA)
habitatSun <- data.frame(Major.Breeding.Habitats = NA)

regs <- c("SE", "MW", "PN")
for (r in 1:length(regs)) {
  
  region <- region1 <- regs[r]
  
  #load("OUT/combined-parameters.rdata")
  load(paste0("OUT/gjamOutput/", propObs, "-",  regs[r], "_obs_plat-1-30000-censor50.rdata"))
  
  
  # load traits
  load(paste0("OUT/", region, "-traitdata.rdata"))
  
  betas <- makeBetas(out)$stand %>%
    filter(predictor == "timeCatsunrise")
  
  bt <- merge(betas, traits, by.x = "species", by.y = "bird")
  
  # remove species that haven't converged
  sp.rem <- sp.remove.time[,c("species", paste0(region, ".obs_plat"))]
  sp.rem <- sp.rem$species[which(sp.rem[,2] > 1.2)]
  
  if (length(sp.rem) > 0) {
    bt <- bt[-which(bt$species %in% sp.rem),]
  }
  
  # take wild turkey out from popularity
  btPop <- bt[which(bt$sp != "wildturkey"),]
  btPop <- bt[which(bt$sports.team == "N"),]
  
  traitsAll[r+3,1] <- round(cor(btPop$Estimate, btPop$popularity, use = "complete"), 2)
  traitsAll[r+3,2] <- round(cor(bt$Estimate, bt$obsSeen, use = "complete"), 2)
  traitsAll[r+3,3] <- round(cor(bt$Estimate, bt$percentVisDet, use = "complete"), 2)
  traitsAll[r+3,4] <- round(cor(bt$Estimate, bt$max.color.contrast, use = "complete"), 2)
  traitsAll[r+3,5] <- round(cor(bt$Estimate, bt$BodyMass.Value, use = "complete"), 2)
  traitsAll[r+3,6] <- round(cor(bt$Estimate, bt$forHeight, use = "complete"), 2)
  
  # discrete traits
  orderSunrise <- aggregate(bt[,"Estimate"], by = list(bt$Order), FUN = "mean")
  colnames(orderSunrise) <- c("Order", "Mean")
  orderSunrise$count <- table(bt$Order)
  
  colnames(orderSunrise)[2:3] <- paste0(colnames(orderSunrise)[2:3], "-", region)

  orderSun <- merge(orderSun, orderSunrise, by = "Order", all = T)
  
  habitatSunrise1 <- separate_rows(bt, 16, sep = "; ")
  habitatSunrise <- aggregate(habitatSunrise1[,"Estimate"], by = list(habitatSunrise1$Major.Breeding.Habitats), FUN = "mean")
  colnames(habitatSunrise) <- c("Major.Breeding.Habitats", "Mean")
  habitatSunrise$count <- table(habitatSunrise1$Major.Breeding.Habitats)
  
  colnames(habitatSunrise)[2:3] <- paste0(colnames(habitatSunrise)[2:3], "-", region)
  
  habitatSun <- merge(habitatSun, habitatSunrise, by = "Major.Breeding.Habitats", all = T)
  
  
  
  ## platform
  
  obs <- makeBetas(out)$stand %>%
    filter(predictor == "platformEBI")
  ot <- merge(obs, traits, by.x = "species", by.y = "bird")
  
  # remove species that haven't converged
  sp.rem <- sp.remove.plat[,c("species", paste0(region, ".obs_plat"))]
  sp.rem <- sp.rem$species[which(sp.rem[,2] > 1.2)]
  
  if (length(sp.rem) > 0) {
    ot <- ot[-which(ot$species %in% sp.rem),]
  }
  
  # take wild turkey out from popularity
  otPop <- ot[which(ot$species != "wildturkey" &
                      ot$sports.team == "N"),]
  
  traitsAll[r,1] <- round(cor(otPop$Estimate, otPop$popularity, use = "complete"), 2)
  traitsAll[r,2] <- round(cor(ot$Estimate, ot$obsSeen, use = "complete"), 2)
  traitsAll[r,3] <- round(cor(ot$Estimate, ot$percentVisDet, use = "complete"), 2)
  traitsAll[r,4] <- round(cor(ot$Estimate, ot$max.color.contrast, use = "complete"), 2)
  traitsAll[r,5] <- round(cor(ot$Estimate, ot$BodyMass.Value, use = "complete"), 2)
  traitsAll[r,6] <- round(cor(ot$Estimate, ot$forHeight, use = "complete"), 2)
  
  # discrete traits
  orderObserver <- aggregate(ot[,"Estimate"], by = list(ot$Order), FUN = "mean")
  colnames(orderObserver) <- c("Order", "Mean")
  orderObserver$count <- table(ot$Order)
  
  colnames(orderObserver)[2:3] <- paste0(colnames(orderObserver)[2:3], "-", region)
  
  orderObs <- merge(orderObs, orderObserver, by = "Order", all = T)
  
  habitatObserver1 <- separate_rows(ot, 16, sep = "; ")
  habitatObserver <- aggregate(habitatObserver1[,"Estimate"], by = list(habitatObserver1$Major.Breeding.Habitats), FUN = "mean")
  colnames(habitatObserver) <- c("Major.Breeding.Habitats", "Mean")
  habitatObserver$count <- table(habitatObserver1$Major.Breeding.Habitats)
  
  colnames(habitatObserver)[2:3] <- paste0(colnames(habitatObserver)[2:3], "-", region)
  
  habitatObs <- merge(habitatObs, habitatObserver, by = "Major.Breeding.Habitats", all = T)
}

colnames(traitsAll) <- c("Popularity", "Log(Obs Seen)", "Sqrt(Percent visual detections)",
                         "Color contrast", "Log(Body mass)", "Foraging stratum")
rownames(traitsAll) <- c("observer-SE", "observer-MW", "observer-PN",
                         "sunrise-SE", "sunrise-MW", "sunrise-PN")

write.csv(traitsAll, file = "OUT/tables/correlations.csv")

orderObs <- orderObs[-nrow(orderObs),]
orderObs[,c(2,4,6)] <- round(orderObs[,c(2,4,6)], 2)
write.csv(orderObs, file = "OUT/tables/order-observer.csv",
          row.names = F)

orderSun <- orderSun[-nrow(orderSun),]
orderSun[,c(2,4,6)] <- round(orderSun[,c(2,4,6)], 2)
write.csv(orderSun, file = "OUT/tables/order-sunrise.csv",
          row.names = F)

habitatObs <- habitatObs[-nrow(habitatObs),]
habitatObs[,c(2,4,6)] <- round(habitatObs[,c(2,4,6)], 2)
write.csv(habitatObs, file = "OUT/tables/habitat-observer.csv",
          row.names = F)

habitatSun <- habitatSun[-nrow(habitatSun),]
habitatSun[,c(2,4,6)] <- round(habitatSun[,c(2,4,6)], 2)
write.csv(habitatSun, file = "OUT/tables/habitat-sunrise.csv",
          row.names = F)



############################################
#                                          #
#             missing species              #
#                                          #
############################################

regs <- c("SE", "MW", "PN")

orderbbs <- data.frame(Order = NA)
orderebd <- data.frame(Order = NA)
habitatbbs <- data.frame(Major.Breeding.Habitats = NA)
habitatebd <- data.frame(Major.Breeding.Habitats = NA)

conTrait <- c()

for (r in 1:length(regs)) {
  region <- regs[r]
  
  # load data
  load(paste0("DATA/stats-", region, ".rdata"))
  
  
  # merge
  spEBD <- data.frame(sp = countInfo[[3]][[1]],
                      ebdPres = "yes")
  spBBS <- data.frame(sp = countInfo[[3]][[2]],
                      bbsPres = "yes")
  
  spAll <- merge(spEBD, spBBS, by = "sp",
                 all = T)
  
  spAll[is.na(spAll)] <- "no"
  
  spAll$species <- gsub(" ",  "", spAll$sp)
  spAll$species <- gsub("'",  "", spAll$species)
  spAll$species <- gsub("-",  "", spAll$species)
  
  # merge with traits
  load(paste0("OUT/", region, "-traitdata.rdata"))
  traits$`obsSeen`[which(traits$`obsSeen` == -Inf)] <- NA
  
  spAll <- merge(spAll, traits, by.x = "species", by.y = "bird")
  
  
  # do t tests
  
  # pull out important columns
  trts <- c("popularity", "percentVisDet",
            "max.color.contrast", "BodyMass.Value", "forHeight", "Major.Breeding.Habitats")
  
  small <- spAll[,c("sp", "bbsPres", trts)]
  
  small1 <- small[-grep("Coast|Wetland", small$Major.Breeding.Habitats),]
  
  ts <- data.frame(trait = trts[1:(length(trts)-1)],
                   meanNo = NA,
                   meanYes = NA,
                   p = NA,
                   platform = "bbs")
  
  for (t in 1:(length(trts)-1)) {
    
    if (trts[t] %in% c("BodyMass.Value", "forHeight")) {
      use <- small1
    } else {
      use <- small
    }
    
    
    tmp <- use[,c("bbsPres", trts[t])]
    no <- tmp[which(tmp$bbsPres == "no"),]
    yes <- tmp[which(tmp$bbsPres == "yes"),]
    
    tmp1 <- Bolstad::bayes.t.test(yes[,trts[t]], no[,trts[t]],
                                  mu = 0)
    
    ts[t, 4] <- tmp1$p.value
    ts[t, 3] <- tmp1$estimate[1]
    ts[t, 2] <- tmp1$estimate[2]
  }
  
  ts$trait <- paste0(ts$trait, "-", region)
  
  conTrait <- bind_rows(conTrait, ts)
  
  
  # count missing species by order
  spAllOrderBBS <- data.frame(table(spAll$Order,
                                    spAll$bbsPres))
  spAllOrderBBS <- pivot_wider(spAllOrderBBS,
                               names_from = "Var2",
                               values_from = "Freq")
  spAllOrderBBS$freqMissing <- spAllOrderBBS$no/(spAllOrderBBS$no + spAllOrderBBS$yes)
  
  colnames(spAllOrderBBS)[2:4] <- paste0(colnames(spAllOrderBBS)[2:4], "-", region)
  colnames(spAllOrderBBS)[1] <- "Order"
  orderbbs <- merge(orderbbs, spAllOrderBBS, by = "Order", all = T)
  

  
  # count missing species by habitat
  spAllHabitatBBS <- separate_rows(spAll, 13, sep = "; ")
  spAllHabitatBBS <- data.frame(table(spAllHabitatBBS$Major.Breeding.Habitats,
                                      spAllHabitatBBS$bbsPres))
  
  spAllHabitatBBS <- pivot_wider(spAllHabitatBBS,
                                 names_from = "Var2",
                                 values_from = "Freq")
  spAllHabitatBBS$freqMissing <- spAllHabitatBBS$no/(spAllHabitatBBS$no + spAllHabitatBBS$yes)
  
  colnames(spAllHabitatBBS)[1] <- "Major.Breeding.Habitats"
  colnames(spAllHabitatBBS)[2:4] <- paste0(colnames(spAllHabitatBBS)[2:4], "-", region)
  habitatbbs <- merge(habitatbbs, spAllHabitatBBS, by = "Major.Breeding.Habitats", all = T)
  
  
  
  
  
  
  
  # Now do for ebd
  
  # do t tests
  
  # pull out important columns
  trts <- c("popularity", "percentVisDet",
            "max.color.contrast", "BodyMass.Value", "forHeight", "Major.Breeding.Habitats")
  
  small <- spAll[,c("sp", "ebdPres", trts)]
  
  small1 <- small[-grep("Coast|Wetland", small$Major.Breeding.Habitats),]
  
  ts <- data.frame(trait = trts[1:(length(trts)-1)],
                   meanNo = NA,
                   meanYes = NA,
                   p = NA,
                   platform = "ebd")
  
  for (t in 1:(length(trts)-1)) {
    
    if (trts[t] %in% c("BodyMass.Value", "forHeight")) {
      use <- small1
    } else {
      use <- small
    }
    
    
    tmp <- use[,c("ebdPres", trts[t])]
    no <- tmp[which(tmp$ebdPres == "no"),]
    yes <- tmp[which(tmp$ebdPres == "yes"),]
    
    tmp1 <- Bolstad::bayes.t.test(yes[,trts[t]], no[,trts[t]],
                                  mu = 0)
    
    ts[t, 4] <- tmp1$p.value
    ts[t, 3] <- tmp1$estimate[1]
    ts[t, 2] <- tmp1$estimate[2]
  }
  
  ts$trait <- paste0(ts$trait, "-", region)
  
  conTrait <- bind_rows(conTrait, ts)
  
  # count missing species by order
  spAllOrderEBD <- data.frame(table(spAll$Order,
                                    spAll$ebd))
  spAllOrderEBD <- pivot_wider(spAllOrderEBD,
                               names_from = "Var2",
                               values_from = "Freq")
  spAllOrderEBD$freqMissing <- spAllOrderEBD$no/(spAllOrderEBD$no + spAllOrderEBD$yes)
  
  colnames(spAllOrderEBD)[1] <- "Order"
  colnames(spAllOrderEBD)[2:4] <- paste0(colnames(spAllOrderEBD)[2:4], "-", region)
  orderebd <- merge(orderebd, spAllOrderEBD, by = "Order", all = T)
  
  
  # count missing species by habitat
  spAllHabitatEBD <- separate_rows(spAll, 13, sep = "; ")
  spAllHabitatEBD <- data.frame(table(spAllHabitatEBD$Major.Breeding.Habitats,
                                      spAllHabitatEBD$ebdPres))
  
  spAllHabitatEBD <- pivot_wider(spAllHabitatEBD,
                                 names_from = "Var2",
                                 values_from = "Freq")
  spAllHabitatEBD$freqMissing <- spAllHabitatEBD$no/(spAllHabitatEBD$no + spAllHabitatEBD$yes)
  
  colnames(spAllHabitatEBD)[1] <- "Major.Breeding.Habitats"
  colnames(spAllHabitatEBD)[2:4] <- paste0(colnames(spAllHabitatEBD)[2:4], "-", region)
  habitatebd <- merge(habitatebd, spAllHabitatEBD, by = "Major.Breeding.Habitats", all = T)
  
}

write.csv(conTrait, file = "OUT/tables/missingSpeciesTraits.csv",
          row.names = F)

orderbbs <- orderbbs[-nrow(orderbbs),]
write.csv(orderbbs, file = "OUT/tables/order-bbsMissing.csv",
          row.names = F)

orderebd <- orderebd[-nrow(orderebd),]
write.csv(orderebd, file = "OUT/tables/order-ebdMissing.csv",
          row.names = F)

habitatbbs <- habitatbbs[-nrow(habitatbbs),]
write.csv(habitatbbs, file = "OUT/tables/habitat-bbsMissing.csv",
          row.names = F)

habitatebd <- habitatebd[-nrow(habitatebd),]
write.csv(habitatebd, file = "OUT/tables/habitat-ebdMissing.csv",
          row.names = F)

