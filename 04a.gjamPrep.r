# Lane Scher
# clanescher@gmail.com
# 03/03/2020

# Prep ydata and xdata for gjam

library(gjam)

############################################
#                                          #
#                load data                 #
#                                          #
############################################

#### load checklist data #####
load(paste0("DATA/combined-", region, ".rdata"))


#### load env covariates ####
load(paste0("DATA/dataCovariates/covar-all.rdata"))
covar[,5:11] <- apply(covar[,5:17],
                       MARGIN = 2,
                       FUN = as.numeric)

############################################
#                                          #
#              format data                 #
#                                          #
############################################

#### combine data####
# xdata and ydata
all <- merge(ydata, xdata, by = "obsID")

# add environmental covariates
covar$geeIDyear <- paste0(covar$geeID, "_", covar$year)
all <- merge(all, covar, by = "geeIDyear")


# add row names
rownames(all) <- all$obsID


#### clean data ####
all$observer <- as.factor(all$observer)

# time categories
all$timeCat <- as.factor(all$timeCat)
all$timeCat <- factor(all$timeCat, levels = c("aamorning",
                                               "sunrise"))


# separate x and y
startX <- grep("observer", colnames(all))
x <- all[,c(1, 2, startX:ncol(all))]
y <- all[,c(3:(startX - 1))]



#### clean ydata ####

# get rid of species we don't have trait data for
rm <- c("Artemisiospiza nevadensis",
        "Selasphorus calliope")
y <- y[,!colnames(y) %in% rm]


minObs <- round(nrow(y)*propObs)

y <- as.data.frame(gjam::gjamTrimY(y[,1:ncol(y)], minObs = minObs)$y)
y <- sapply(y, as.numeric)
y <- as.data.frame(y)
y$other <- NULL

y[is.na(y) == T] <- 0


# xdata names
xVarNames <- data.frame(variable = colnames(x),
                        description = c("geeID and year",
                                        "observation ID - this column is unique to each observation",
                                        "observer", "distance in km", "duration in minutes", "minutes since local sunrise",
                                        "time of day category", "source (ebird or BBS)",
                                        "geeID", "year",
                                        "latitude", "longitude",
                                        "elevation",
                                        "total annual precipitation (June_y-1 through May_y)",
                                        "total spring precipitation (March_y through May_y)",
                                        "total winter precipitation (Dec_y-1 through Feb_y)",
                                        "mean temperature of coldest winter month (Dec_y-1 through Feb_y)",
                                        "mean temperature of coldest spring month (March_y through May_y)",
                                        "mean annual temperature (June_y-1 through May_y)",
                                        "% water habitat",
                                        "% developed habitat",
                                        "% forest habitat",
                                        "% shrub habitat",
                                        "% grassland habitat",
                                        "% wetland habitat"))

