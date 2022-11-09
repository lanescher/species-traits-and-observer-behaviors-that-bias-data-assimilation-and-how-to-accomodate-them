
# Lane Scher
# clanescher@gmail.com
# 11/21/2020

# This script runs the entire eBird-BBS analysis

# set working directory
wd <- "X:/eBird-BBS/"
setwd(wd)

years <- c(2010:2019)




############################################
#                                          #
#           1. format BBS data             #
#                                          #
############################################

# This section processes all three regions at the same time

#### --------process BBS dates-------- ####

source("CODE/01a.BBSdates.r")

# This script reads in the dates that BBS routes were surveyed
# and creates a two-week window on either side that is used
# to filter eBird checklists.



#### --------process BBS route location data -------- ####

source("CODE/01b.BBSroutes.r")

# This script processes BBS route data and merges it with
# metadata required for eBird checklist aggregation



#### --------format BBS count data-------- ####

source("CODE/01c.BBScounts.r")

# This script processes BBS counts and saves them in 
# gjam format




#### --------process BBS observers-------- ####

source("CODE/01d.BBSobservers.r")

# This script processes BBS observer data



############################################
#                                          #
#           2: format eBird data           #
#                                          #
############################################


#### --------process eBird count data-------- ####

source("CODE/02a.eBirdFilter.r")

# This script filters the raw eBird data


#### --------aggregate checklists to BBS routes-------- ####

source("CODE/02b.eBirdAggregate.r")

# This script pulls raw eBird data, cleans it, and keeps the
# correct checklists within the buffer around BBS sites.
# It then reformats the eBird data for GJAM.


############################################
#                                          #
#               combine data               #
#                                          #
############################################



#### --------combine ebird and bbs data-------- ####
region <- c("SE")
years <- c(2010:2019)
source("CODE/03a.combineData.r")

# This script combines eBird and BBS data into xdata and ydata.
# It formats their count-specific covariates (like time, distance,
# observer, etc) and only keeps observations that match time and location


#### --------download covariates-------- ####

source("CODE/03b.downloadEnv.r")

# This script downloads data from Daymet and reformats
# it for GJAM


############################################
#                                          #
#              4: run gjam                 #
#                                          #
############################################

#### --------combine covariates and bird counts-------- ####
region <- c("MW")
propObs <- .05                          # proportion of observations species must be in to be included in the model
source("CODE/04a.gjamPrep.r")

# This script combines all data from eBird and BBS, and
# merges them with environmental covariates. It produces
# a dataframe with all the data required for 04b.gjam.r


#### --------run gjam-------- ####
plat <- T      # include platform as fixed effect?
obs <- T       # include observer as random effect? 
run <- 3       # which of the 3 chains is this?

ng <- 30000                          # iterations
burnin <- ng/2                         # burnin
source("CODE/04b.gjam.r")

# This script takes output from 04a.combineData.r and fits
# a gjam. It saves the output in OUT/gjamOutput/


#### --------compare model fits-------- ####
propObs <- .05                          # proportion of observations species must be in to be included in the model
source("CODE/04c.compareModelFit.r")

# This script compares the fit of the models produced by
# CODE/04b.gjam.r


############################################
#                                          #
#               5: output                  #
#                                          #
############################################




#### --------output-------- ####
propObs <- .05                          # proportion of observations species must be in to be included in the model
source("CODE/05.output.r")


#### --------figures-------- ####
propObs <- 0.05
source("CODE/06.figures.r")

# fixed map


#### --------tables-------- ####
propObs <- .05                          # proportion of observations species must be in to be included in the model
source("CODE/07.tables.r")



