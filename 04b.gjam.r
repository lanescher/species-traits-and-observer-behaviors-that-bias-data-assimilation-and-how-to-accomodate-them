
# Lane Scher
# clanescher@gmail.com
# 06/08/2020

# This script fits the gjam using output from 04a.gjamPrep.r
# It runs the gjam model, which identifies the source effect


#### read libraries ####
library(gjam)
library(tidyverse)



############################################
#                                          #
#                run GJAM                  #
#                                          #
############################################


#### set up model ####
#x$lc <- as.factor(x$lc)
x$platform <- as.factor(x$source)
x$observer <- as.factor(x$observer)

effort <- list(columns=1:ncol(y), values = x$effMin)


#### make censor list
# censor
censor <- 50
y[y > censor] <- censor



if (obs == T & plat == F) {
  formulaI <- as.formula(~ prcpAnnual + tminWinter + elev +
                           water + developed + forest + shrub + grass + wetland +
                           timeCat)
  
  ml   <- list(ng = ng, burnin = burnin, typeNames = "DA", 
               effort = effort,
               random = "observer")
}


if (obs == F & plat == T) {
  formulaI <- as.formula(~ prcpAnnual + tminWinter + elev +
                           water + developed + forest + shrub + grass + wetland +
                           timeCat + platform)
  
  ml   <- list(ng = ng, burnin = burnin, typeNames = "DA", 
               effort = effort)
}

if (obs == T & plat == T) {
  formulaI <- as.formula(~ prcpAnnual + tminWinter + elev +
                           water + developed + forest + shrub + grass + wetland +
                           timeCat + platform)
  
  ml   <- list(ng = ng, burnin = burnin, typeNames = "DA", 
               effort = effort,
               random = "observer")
}



print(paste0("Number of species: ", ncol(y)))
print(paste0("Observer random effect: ", obs))
print(paste0("Platform fixed effect: ", plat))
print(paste0("Run: ", run))
startTime <- Sys.time()
print(paste0("Started at ", startTime))

out <- gjam(formulaI, xdata = x, ydata = y, modelList = ml)

endTime <- Sys.time()
print(paste0("Ended at ", endTime))

if (obs == T & plat == F) {
  save(out, file = paste0("OUT/gjamOutput/", propObs, '-', region, "_obs-", run, "-", ng, "-censor", censor, ".rdata"))
}

if (obs == F & plat == T) {
  save(out, file = paste0("OUT/gjamOutput/", propObs, '-', region, "_plat-", run, "-", ng, "-censor", censor, ".rdata"))
}

if (obs == T & plat == T) {
  save(out, file = paste0("OUT/gjamOutput/", propObs, '-', region, "_obs_plat-", run, "-", ng, "-censor", censor, ".rdata"))
}
