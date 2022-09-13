
# Lane Scher
# clanescher@gmail.com
# 11/21/2020

# This script analyzes output and makes preliminary figures

# update 2/3/22: modified significance for observer effect

library(tidyverse)


all.cors <- data.frame(trait = c("pop", "com", "vis", "col", "siz", "hei"))


regs <- c("SE", "MW", "PN")

for (r in 1:length(regs)) {
  
  #load("OUT/combined-parameters.rdata")
  
  region <- region1 <- regs[r]

  #### load gjam output ####
  load(paste0("OUT/gjamOutput/0.05-", region1, "_plat-1-20000-censor50.rdata"))
  out.plat <- out
  
  load(paste0("OUT/gjamOutput/0.05-", region1, "_obs_plat-1-30000-censor50.rdata"))
  out.obsplat <- out
  
  
  #### also get inputs ####
  source("CODE/04a.gjamPrep.r")
  
  
  ############################################
  #                                          #
  #           get traits together            #
  #                                          #
  ############################################
  
  #### calculate commonness ####
  
  load(paste0("X:/eBird-BBS/DATA/stats-", region1, ".rdata"))
  counts <- countInfo[[4]]
  counts$obsSeen <- counts$ebi + counts$bbs
  
  counts$bird <- gsub(" ", "", counts$species)
  counts$bird <- gsub("'", "", counts$bird)
  counts$bird <- gsub("-", "", counts$bird)
  
  
  
  #### read in trait data ####
  traitsOriginal <- read.csv(file = "DATA/dataOther/birdTraits3.csv")
  
  traitsOriginal$bird <- gsub(" ", "", traitsOriginal$nameCom)
  traitsOriginal$bird <- gsub("-", "", traitsOriginal$bird)
  traitsOriginal$bird <- gsub("'", "", traitsOriginal$bird)
  
  #### merge traits and counts ####
  traits <- merge(traitsOriginal, counts, by = "bird",
                  all.x = T)
  
  #### get mean foraging height ####
  for_height <- traits[,c(grep("nameCom", colnames(traits)), grep("ForStrat", colnames(traits)))]
  for_height$ForStrat.Source <- for_height$ForStrat.SpecLevel <- NULL
  colnames(for_height) <- c("nameCom", -1, 0, 0, 1, 2, 3, 4)
  
  for_height_mean <- pivot_longer(for_height, cols = !nameCom) %>%
    mutate(name = as.numeric(name)) %>%
    group_by(nameCom) %>%
    summarize(forHeight = weighted.mean(name, value))
  
  traits <- merge(traits, for_height_mean, by = "nameCom")
  
  #### subset traits ####
  # select important traits
  trts <- c("popularity", "obsSeen",
            "percentVisDet", "max.color.contrast",
            "BodyMass.Value", "forHeight", "Major.Breeding.Habitats",
            "Order", "sports.team")
  
  # subset
  traits <- traits[, c("nameCom", "bird", "nameSci", trts)]
  traits$obsSeen <- log(traits$obsSeen)
  traits$percentVisDet <- sqrt(traits$percentVisDet)
  traits$BodyMass.Value <- log(traits$BodyMass.Value)
  
  
  #### calculate genus mean (this will be used to fill in species with missing data) ####
  # calculate genus mean
  traits$genus <- gsub(' .*', "", traits$nameSci)
  traitsGen <- aggregate(traits[,trts[1:6]],
                         by = list(traits$genus),
                         FUN = "mean", na.rm = T)
  colnames(traitsGen) <- c("genus", "popMean", "comMean", "visMean",
                           "colMean", "bodMean", "forMean")
  traits <- merge(traits, traitsGen, by = "genus")
  
  #### now just keep species from analysis ####
  spKeep <- countInfo[[4]]$species
  
  #traits <- traits[which(traits$nameCom %in% colnames(y)),]
  traits <- traits[which(traits$nameCom %in% spKeep),]
  
  
  #### fill in means ####
  # make new df
  traits2 <- traits
  
  # fill in means
  traits2$popularity[which(is.na(traits2$popularity))] <- traits2$popMean[which(is.na(traits2$popularity))]
  traits2$obsSeen[which(is.na(traits2$obsSeen))] <- traits2$comMean[which(is.na(traits2$obsSeen))]
  traits2$percentVisDet[which(is.na(traits2$percentVisDet))] <- traits2$visMean[which(is.na(traits2$percentVisDet))]
  traits2$max.color.contrast[which(is.na(traits2$max.color.contrast))] <- traits2$colMean[which(is.na(traits2$max.color.contrast))]
  traits2$BodyMass.Value[which(is.na(traits2$BodyMass.Value))] <- traits2$bodMean[which(is.na(traits2$BodyMass.Value))]
  traits2$forHeight[which(is.na(traits2$forHeight))] <- traits2$forMean[which(is.na(traits2$forHeight))]
  
  # save df as traits3 to use for analysis
  traits3 <- traits2[, c("nameCom", "bird", "nameSci", trts)]
  
  # add asterisks and save for table
  traits2 <- traits2[,c("nameCom", "nameSci", trts)]
  traits2[,3:7] <- round(traits2[,3:7], digits = 3)
  
  traitsAst <- traits[,c("nameCom", "nameSci", trts[1:6])]
  traitsAst[is.na(traitsAst)==F] <- ""
  traitsAst[is.na(traitsAst)] <- "*"
  
  traits2$popularity <- paste0(traits2$popularity, traitsAst$popularity)
  traits2$obsSeen <- paste0(traits2$obsSeen, traitsAst$obsSeen)
  traits2$percentVisDet <- paste0(traits2$percentVisDet, traitsAst$percentVisDet)
  traits2$max.color.contrast <- paste0(traits2$max.color.contrast, traitsAst$max.color.contrast)
  traits2$BodyMass.Value <- paste0(traits2$BodyMass.Value, traitsAst$BodyMass.Value)
  
  
  #### use traits3 ####
  traits <- traits3
  
  #### save traits and traits2 (traits2 has asterisks showing genus means) ####
  save(traits, traits2, file = paste0("OUT/", region1, "-traitdata.rdata"))
  
  
  #### start writing to txt file ####
  sink(file = paste0("OUT/text/", region1, "-summary.txt"))
  
  
  ############################################
  #                                          #
  #             summarize data               #
  #                                          #
  ############################################
  
  xBBS <- x[which(x$source == "BBS"),]
  xEBI <- x[which(x$source == "EBI"),]
  
  xBBS$observer <- as.character(xBBS$observer)
  xEBI$observer <- as.character(xEBI$observer)
  
  
  #### BBS data ####
  writeLines("\n")
  print("BBS DATA")
  print(paste0("number BBS point counts: ", nrow(xBBS)))
  print(paste0("from ", length(unique(xBBS$geeIDyear)), " route-years"))
  
  print(paste0("number unique observers: ", length(unique(xBBS$observer))))
  tmp <- as.data.frame(table(xBBS$observer))
  print(paste0("mean point counts per observer: ", mean(tmp$Freq)))
  print(paste0("median point counts per observer: ", median(tmp$Freq)))
  
  tmp <- table(xBBS$timeCat)
  tmp <- tmp/nrow(xBBS)
  print(tmp)
  
  
  #### EBI data ####
  writeLines("\n")
  print("EBIRD DATA")
  print(paste0("number EBI checklists: ", nrow(xEBI)))
  print(paste0("from ", length(unique(xEBI$geeIDyear)), " route-years"))
  tmp <- as.data.frame(table(xEBI$geeIDyear))
  print(paste0("mean checklists per route-year: ", mean(tmp$Freq)))
  print(paste0("median checklists per route-year: ", median(tmp$Freq)))
  
  print(paste0("number unique observers: ", length(unique(xEBI$observer))))
  tmp <- as.data.frame(table(xEBI$observer))
  print(paste0("mean checklists per observer: ", mean(tmp$Freq)))
  print(paste0("median checklists per observer: ", median(tmp$Freq)))
  
  tmp <- table(xEBI$timeCat)
  tmp <- tmp/nrow(xEBI)
  print(tmp)
  
  print(paste0("shortest duration: ", min(xEBI$effMin)))
  print(paste0("longest duration: ", max(xEBI$effMin)))
  print(paste0("mean duration: ", mean(xEBI$effMin)))
  print(paste0("median duration: ", median(xEBI$effMin)))
  
  print(paste0("shortest distance traveled: ", min(xEBI$effKm)))
  print(paste0("longest distance traveled: ", max(xEBI$effKm)))
  print(paste0("mean distance traveled: ", mean(xEBI$effKm)))
  print(paste0("median distance traveled: ", median(xEBI$effKm)))
  tmp <- as.data.frame(table(xEBI$effKm))
  print(paste0("percent stationary checklists: ", (tmp$Freq[1]/nrow(xEBI))*100))
  
  
  ############################################
  #                                          #
  #               time of day                #
  #                                          #
  ############################################
  
  
  #### output ####
  betas.plat <- out.plat$parameters$betaStandXTable

  # make bird column
  betas.plat$species <- as.vector(rownames(betas.plat))
  betas.plat$species <- sub("_.*", "", betas.plat$species)

  # make predictor column
  betas.plat$predictor <- as.vector(rownames(betas.plat))
  betas.plat$predictor <- sub(".*_", "", betas.plat$predictor)

  betas.plat <- betas.plat %>%
    filter(predictor == "timeCatsunrise")
  
  
  #### output ####
  betas.obsplat <- out.obsplat$parameters$betaStandXTable
  
  # make bird column
  betas.obsplat$species <- as.vector(rownames(betas.obsplat))
  betas.obsplat$species <- sub("_.*", "", betas.obsplat$species)
  
  # make predictor column
  betas.obsplat$predictor <- as.vector(rownames(betas.obsplat))
  betas.obsplat$predictor <- sub(".*_", "", betas.obsplat$predictor)
  
  betas.obsplat <- betas.obsplat %>%
    filter(predictor == "timeCatsunrise")
  
 
  corbetas <- cor(betas.obsplat$Estimate, betas.plat$Estimate)
  print(paste0("Correlation between sunrise betas: ", round(corbetas, 3)))
  
  
  # pull out the species that didn't converge
  sp.remove <- read.csv("OUT/tables/gelman_time.csv")
  sp.remove.col.obsplat <- grep(paste0(region1, ".obs_plat"), colnames(sp.remove))
  sp.remove.col.plat <- grep(paste0(region1, ".plat"), colnames(sp.remove))
  sp.remove <- sp.remove[,c(1, sp.remove.col.obsplat, sp.remove.col.plat)]
  sp.remove.obsplat <- sp.remove$species[which(sp.remove[,2] > 1.2)]
  sp.remove.plat <- sp.remove$species[which(sp.remove[,3] > 1.2)]
  sp.remove.both <- sp.remove$species[which(sp.remove[,2] > 1.2 |
                                              sp.remove[,3] > 1.2)]
  
  if (length(sp.remove.plat) > 0) {
    betas.plat <- betas.plat[-which(betas.plat$species %in% sp.remove.plat),]
  }
  if (length(sp.remove.obsplat) > 0) {
    betas.obsplat <- betas.obsplat[-which(betas.obsplat$species %in% sp.remove.obsplat),]
  }
  
  if (length(sp.remove.plat) > 0 |
      length(sp.remove.obsplat) > 0) {
    betas.plat1 <- betas.plat[-which(betas.plat$species %in% sp.remove.both),]
    betas.obsplat1 <- betas.obsplat[-which(betas.obsplat$species %in% sp.remove.both),]
  } else {
    betas.plat1 <- betas.plat
    betas.obsplat1 <- betas.obsplat
  }
  
  betas <- list(betas.plat, betas.obsplat, betas.plat1, betas.obsplat1)
  for (b in 1:length(betas)) {
    betas[[b]] <- merge(betas[[b]], traits, by.x = "species", by.y = "bird")
    
  }
  
  
  
  
  #### time categories with  continuous traits
  
  
  #preds <- c("timeCatsunrise")
  trts <- c("popularity", "obsSeen", "percentVisDet",
            "max.color.contrast", "BodyMass.Value", 'forHeight')
  
  
  writeLines("\n")
  print("TIME OF DAY")
  
  
  #### effect at different times of day ####
  for (d in 1:length(betas)) {
    if (d == 1) {print("PLATFORM ONLY")
      lab <- "plat"}
    if (d == 2) {print("PLATFORM AND OBSERVER")
      lab <- "obsplat"}
    if (d == 3) {print("PLATFORM ONLY - matching species")
      lab <- "plat-matching"}
    if (d == 4) {print("PLATFORM AND OBSERVER - matching species")
      lab <- "obsplat-matching"}
    
    time1 <- betas[[d]]
    
    ns <- length(unique(time1$species))
    
    sig <- nrow(time1[which(time1$sig95 == "*"),])
    print(paste0("Sunrise has a significant effect on ", round((sig/ns)*100, 2), " % of species"))
    time1Sig <- time1[which(time1$sig95 == "*"),]
    
    pos <- nrow(time1Sig[which(time1Sig$Estimate > 0),])
    neg <- nrow(time1Sig[which(time1Sig$Estimate < 0),])
    print(paste0("Of those, ", round((pos/nrow(time1Sig))*100, 2), "% are positive effects and ", 
                 round((neg/nrow(time1Sig))*100,2) , "% are negative effects"))
    
    
    print(paste0("Min effect at sunrise: ", min(time1Sig$Estimate)))
    print(paste0("Max effect at sunrise: ", max(time1Sig$Estimate)))
    print(paste0("Mean effect at sunrise: ", mean(time1Sig$Estimate)))
    print(paste0("Median effect at sunrise: ", median(time1Sig$Estimate)))
    
    writeLines("\n")

    summary <- aggregate(time1[,"Estimate"], by = list(time1$Order), FUN = "mean")
    colnames(summary) <- c("Order", "Mean")
    summary$count <- table(time1$Order)
    
    print(summary)
    
    writeLines("\n")
    
    
    time2 <- separate_rows(time1, 16, sep = "; ")
    summary <- aggregate(time2[,"Estimate"], by = list(time2$Major.Breeding.Habitats), FUN = "mean")
    colnames(summary) <- c("Habitat", "Mean")
    summary$count <- table(time2$Major.Breeding.Habitats)
    
    print(summary)
    
    writeLines("\n")
    
    
    
    cors <- data.frame(trait = c("pop", "com", "vis", "col", "siz", "hei"),
                       cor = NA)
    for (t in 1:nrow(cors)) {
      if (t == 1) { # take turkey and sports teams out of popularity
        time3 <- time1 %>%
          filter(species != 'wildturkey',
                 sports.team == "N")
      } 
      
      if (t %in% c(2:6)) {
        time3 <- time1
      }
      cors[t, 2] <- cor(time3[,"Estimate"], time3[,9+t], use = "complete.obs")
    }
    
    print(paste0("Effect at sunrise is correlated with: "))
    print(cors)
    
    all.cors <- inner_join(all.cors, cors, by = "trait")
    colnames(all.cors)[ncol(all.cors)] <- paste0(region, '-', lab, "-sunrise")
   
    
    writeLines("\n")
  }
  
  
  
  
  ############################################
  #                                          #
  #               platform                   #
  #                                          #
  ############################################
  
  #### output ####
  betas.plat <- out.plat$parameters$betaStandXTable
  
  # make bird column
  betas.plat$species <- as.vector(rownames(betas.plat))
  betas.plat$species <- sub("_.*", "", betas.plat$species)
  
  # make predictor column
  betas.plat$predictor <- as.vector(rownames(betas.plat))
  betas.plat$predictor <- sub(".*_", "", betas.plat$predictor)
  
  betas.plat <- betas.plat %>%
    filter(predictor == "platformEBI")
  
  
  #### output ####
  betas.obsplat <- out.obsplat$parameters$betaStandXTable
  
  # make bird column
  betas.obsplat$species <- as.vector(rownames(betas.obsplat))
  betas.obsplat$species <- sub("_.*", "", betas.obsplat$species)
  
  # make predictor column
  betas.obsplat$predictor <- as.vector(rownames(betas.obsplat))
  betas.obsplat$predictor <- sub(".*_", "", betas.obsplat$predictor)
  
  betas.obsplat <- betas.obsplat %>%
    filter(predictor == "platformEBI")
  
  
  corbetas <- cor(betas.obsplat$Estimate, betas.plat$Estimate)
  print(paste0("Correlation between ebird betas: ", round(corbetas, 3)))
  
  
  # pull out the species that didn't converge
  sp.remove <- read.csv("OUT/tables/gelman_plat.csv")
  sp.remove.col.obsplat <- grep(paste0(region1, ".obs_plat"), colnames(sp.remove))
  sp.remove.col.plat <- grep(paste0(region1, ".plat"), colnames(sp.remove))
  sp.remove <- sp.remove[,c(1, sp.remove.col.obsplat, sp.remove.col.plat)]
  sp.remove.obsplat <- sp.remove$species[which(sp.remove[,2] > 1.2)]
  sp.remove.plat <- sp.remove$species[which(sp.remove[,3] > 1.2)]
  sp.remove.both <- sp.remove$species[which(sp.remove[,2] > 1.2 |
                                              sp.remove[,3] > 1.2)]
  
  if(length(sp.remove.plat) > 0) {
    betas.plat <- betas.plat[-which(betas.plat$species %in% sp.remove.plat),]
  }
  if (length(sp.remove.obsplat) > 0) {
    betas.obsplat <- betas.obsplat[-which(betas.obsplat$species %in% sp.remove.obsplat),]
  }
  
  if (length(sp.remove.plat) > 0 |
      length(sp.remove.obsplat) > 0) {
    betas.plat1 <- betas.plat[-which(betas.plat$species %in% sp.remove.both),]
    betas.obsplat1 <- betas.obsplat[-which(betas.obsplat$species %in% sp.remove.both),]
  } else {
    betas.plat1 <- betas.plat
    betas.obsplat1 <- betas.obsplat
  }
  
  betas <- list(betas.plat, betas.obsplat, betas.plat1, betas.obsplat1)
  
  for (b in 1:length(betas)) {
    betas[[b]] <- merge(betas[[b]], traits, by.x = "species", by.y = "bird")
    
  }
  
  
  
  
  #### eBird with  continuous traits
  
  trts <- c("popularity", "obsSeen", "percentVisDet",
            "max.color.contrast", "BodyMass.Value", 'forHeight')
  
  
  writeLines("\n")
  print("PLATFORM")
  
  
  #### effect of eBird ####
  for (d in 1:length(betas)) {
    if (d == 1) {print("PLATFORM ONLY")
      lab <- "plat"}
    if (d == 2) {print("PLATFORM AND OBSERVER")
      lab <- "obsplat"}
    if (d == 3) {print("PLATFORM ONLY - matching species")
      lab <- "plat-matching"}
    if (d == 4) {print("PLATFORM AND OBSERVER - matching species")
      lab <- "obsplat-matching"}
    
    time1 <- betas[[d]]
    
    ns <- length(unique(time1$species))
    
    sig <- nrow(time1[which(time1$sig95 == "*"),])
    print(paste0("eBird has a significant effect on ", round((sig/ns)*100, 2), " % of species"))
    time1Sig <- time1[which(time1$sig95 == "*"),]
    
    pos <- nrow(time1Sig[which(time1Sig$Estimate > 0),])
    neg <- nrow(time1Sig[which(time1Sig$Estimate < 0),])
    print(paste0("Of those, ", round((pos/nrow(time1Sig))*100, 2), "% are positive effects and ", 
                 round((neg/nrow(time1Sig))*100,2) , "% are negative effects"))
    
    
    print(paste0("Min effect of eBird: ", min(time1Sig$Estimate)))
    print(paste0("Max effect of eBird: ", max(time1Sig$Estimate)))
    print(paste0("Mean effect of eBird: ", mean(time1Sig$Estimate)))
    print(paste0("Median effect of eBird: ", median(time1Sig$Estimate)))
    
    writeLines("\n")
    
    summary <- aggregate(time1[,"Estimate"], by = list(time1$Order), FUN = "mean")
    colnames(summary) <- c("Order", "Mean")
    summary$count <- table(time1$Order)
    
    print(summary)
    
    writeLines("\n")
    
    
    time2 <- separate_rows(time1, 16, sep = "; ")
    summary <- aggregate(time2[,"Estimate"], by = list(time2$Major.Breeding.Habitats), FUN = "mean")
    colnames(summary) <- c("Habitat", "Mean")
    summary$count <- table(time2$Major.Breeding.Habitats)
    
    print(summary)
    
    writeLines("\n")
    
    
    
    cors <- data.frame(trait = c("pop", "com", "vis", "col", "siz", "hei"),
                       cor = NA)
    for (t in 1:nrow(cors)) {
      if (t == 1) { # take turkey and sports teams out of popularity
        time3 <- time1 %>%
          filter(species != 'wildturkey',
                 sports.team == "N")
      } 
      
      if (t %in% c(2:6)) {
        time3 <- time1
      }
      cors[t, 2] <- cor(time3[,"Estimate"], time3[,9+t], use = "complete.obs")
    }
    
    print(paste0("Effect of eBird is correlated with: "))
    print(cors)
    
    all.cors <- inner_join(all.cors, cors, by = "trait")
    colnames(all.cors)[ncol(all.cors)] <- paste0(region, '-', lab, "-platform")
    
    
    writeLines("\n")
  }
  
  
  
  
  
  
  ############################################
  #                                          #
  #               observers                  #
  #                                          #
  ############################################
  



  ### get difference
  obs <- out.obsplat$parameters$randByGroupMu

  wb <- which(substr(colnames(obs), 1, 1) != "o")
  we <- which(substr(colnames(obs), 1, 1) == "o")

  muB <- apply(obs[,wb], 1, mean)
  muE <- apply(obs[,we], 1, mean)

  mdiff <- muE - muB

  vrB <- apply(obs[,wb], 1, var)
  vrE <- apply(obs[,we], 1, var)

  vdiff <- sqrt((vrB + vrE)/(length(wb) + length(we)))

  spObs <- data.frame(species = names(mdiff),
                      Estimate = mdiff,
                      se = vdiff)

  spObs <- spObs %>%
    mutate(hi = Estimate + 2*se,
           lo = Estimate - 2*se,
           sig = case_when(hi * lo > 0 ~ "*",
                           hi * lo < 0 ~ ""))

  
  spObs <- merge(spObs, traits, by.x = "species", by.y = "bird")
  
  
  writeLines("\n")
  print("OBSERVER EFFECT")
  
  # ebird effect
  sig <- spObs[which(spObs$sig == "*"),]
  print(paste0("Platform effect is significant for ", round(nrow(sig)/ns*100, 2), " % of species"))
  
  posEff <- nrow(sig[which(sig$Estimate > 0),])
  print(paste0("Platform effect is positive and significant for eBirders for ", round((posEff/ns)*100, 2), " % of species"))
  
  posEff <- nrow(sig[which(sig$Estimate < 0),])
  print(paste0("Platform effect is negative and significant for eBirders for ", round((posEff/ns)*100, 2), " % of species"))
  
  print(paste0("Max effect of eBird observers: ", round(max(spObs$Estimate), 2)))
  print(paste0("Min effect of eBird observers: ", round(min(spObs$Estimate), 2)))
  print(paste0("Mean effect of eBird observers: ", round(mean(spObs$Estimate), 2)))
  print(paste0("Median effect of eBird observers: ", round(median(spObs$Estimate), 2)))
  
  
  #### in each habitat
  writeLines("\n")
  
  spObs1 <- separate_rows(spObs, 15, sep = "; ")
  summary <- aggregate(spObs1[,"Estimate"], by = list(spObs1$Major.Breeding.Habitats), FUN = "mean")
  colnames(summary) <- c("Habitat", "Mean")
  summary$count <- table(spObs1$Major.Breeding.Habitats)
  print(summary)
  
  #### in each Order 
  writeLines("\n")
  
  summary <- aggregate(spObs[,"Estimate"], by = list(spObs$Order), FUN = "mean")
  colnames(summary) <- c("Order", "Mean")
  traitsTMP <- traits[traits$bird %in% spObs$species,]
  summary$count <- table(traitsTMP$Order)
  print(summary)
  
  #### with continuous traits ####
  writeLines("\n")
  
  spObs <- spObs[, c("Estimate", "popularity", "obsSeen",
                     "percentVisDet", "max.color.contrast",
                     "BodyMass.Value", "forHeight", "sports.team", "Major.Breeding.Habitats", "species")]
  
  cors <- data.frame(trait = c("pop", "com", "vis", "col", "siz", "hei"),
                     cor = NA)
  for (t in 2:(ncol(spObs)-3)) {
    if (t == 1) {
      #spObs2 <- spObs[which(spObs$sp != "wildturkey"),]
      spObs2 <- spObs %>%
        filter(sp != 'wildturkey',
               sports.team == "Y")
    } 
    
    if (t %in% 2:6) {
      spObs2 <- spObs
    }
    cors[t-1, 2] <- cor(spObs2[,1], spObs2[,t], use = "complete.obs")
  }
  
  # now get just most common half
  
  print(paste0("Effect of eBird observers is correlated with: "))
  print(cors)
  
  

  
  ############################################
  #                                          #
  #           missing  species               #
  #                                          #
  ############################################
  
  writeLines("\n")
  print("SPECIES REPORTED")
  
  # load data
  load(paste0("DATA/stats-", region1, ".rdata"))
  
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
  spAll <- merge(spAll, traits, 
                 by.x = "sp", by.y = "nameCom",
                 all.x = T)
  
  spAll <- merge(spAll, counts,
                 by = "bird")
  

  
  
  #### find species missing from BBS ####
  
  # count missing species by order
  spAllOrderBBS <- data.frame(table(spAll$Order,
                                    spAll$bbsPres))
  spAllOrderBBS <- pivot_wider(spAllOrderBBS,
                               names_from = "Var2",
                               values_from = "Freq")
  spAllOrderBBS$freqMissing <- spAllOrderBBS$no/(spAllOrderBBS$no + spAllOrderBBS$yes)
  
  writeLines("\n")
  
  print("Present in BBS: ")
  print(as.data.frame(spAllOrderBBS))
  
  
  # count missing species by habitat
  spAllHabitatBBS <- separate_rows(spAll, 13, sep = "; ")
  spAllHabitatBBS <- data.frame(table(spAllHabitatBBS$Major.Breeding.Habitats,
                                      spAllHabitatBBS$bbsPres))
  
  spAllHabitatBBS <- pivot_wider(spAllHabitatBBS,
                                 names_from = "Var2",
                                 values_from = "Freq")
  spAllHabitatBBS$freqMissing <- spAllHabitatBBS$no/(spAllHabitatBBS$no + spAllHabitatBBS$yes)
  
  writeLines("\n")
  
  print("Present in BBS: ")
  print(as.data.frame(spAllHabitatBBS))
  
  tot <- sum(spAllOrderBBS$no)
  print(paste0("Total number of species missing from BBS: ", tot))
  
  
  
  #### find species missing from eBird ####
  
  # count missing species by order
  spAllOrderEBD <- data.frame(table(spAll$Order,
                                    spAll$ebd))
  spAllOrderEBD <- pivot_wider(spAllOrderEBD,
                               names_from = "Var2",
                               values_from = "Freq")
  spAllOrderEBD$freqMissing <- spAllOrderEBD$no/(spAllOrderEBD$no + spAllOrderEBD$yes)
  
  writeLines("\n")
  
  print("Present in EBD: ")
  print(as.data.frame(spAllOrderEBD))
  
  
  # count missing species by habitat
  spAllHabitatEBD <- separate_rows(spAll, 13, sep = "; ")
  spAllHabitatEBD <- data.frame(table(spAllHabitatEBD$Major.Breeding.Habitats,
                                      spAllHabitatEBD$ebdPres))
  
  spAllHabitatEBD <- pivot_wider(spAllHabitatEBD,
                                 names_from = "Var2",
                                 values_from = "Freq")
  spAllHabitatEBD$freqMissing <- spAllHabitatEBD$no/(spAllHabitatEBD$no + spAllHabitatEBD$yes)
  
  writeLines("\n")
  
  print("Present in EBD: ")
  print(as.data.frame(spAllHabitatEBD))
  
  tot <- sum(spAllOrderEBD$no)
  print(paste0("Total number of species missing from eBird: ", tot))
  
  
  
  #### look at missing species with continuous traits ####
  
  ### bbs
  preds <- c("bbsPres")
  trts <- c("popularity", "percentVisDet",
            "max.color.contrast", "BodyMass.Value", "forHeight")
  
  small <- spAll[,c("sp", "bbsPres", trts)]
  
  
  
  ts <- data.frame(trait = trts,
                   meanNo = NA,
                   meanYes = NA,
                   p = NA)
  for (t in 1:length(trts)) {
    tmp <- small[,c(preds, trts[t])]
    no <- tmp[which(tmp$bbsPres == "no"),]
    yes <- tmp[which(tmp$bbsPres == "yes"),]
    
    tmp1 <- Bolstad::bayes.t.test(yes[,trts[t]], no[,trts[t]],
                                  mu = 0)
    
    ts[t, 4] <- tmp1$p.value
    ts[t, 2] <- tmp1$estimate[2]
    ts[t, 3] <- tmp1$estimate[1]
  }
  
  ts$trait[1] <- "obsSeenEBI"
  
  writeLines("\n")
  print("t-test of species seen vs. not seen in BBS: ")
  
  print(ts)
  
  ### ebd
  preds <- c("ebdPres")
  trts <- c("popularity", "percentVisDet",
            "max.color.contrast", "BodyMass.Value", "forHeight")
  
  small <- spAll[,c("sp", "ebdPres", trts)]
  
  
  
  ts <- data.frame(trait = trts,
                   meanNo = NA,
                   meanYes = NA,
                   p = NA)
  for (t in 1:length(trts)) {
    tmp <- small[,c(preds, trts[t])]
    no <- tmp[which(tmp$ebdPres == "no"),]
    yes <- tmp[which(tmp$ebdPres == "yes"),]
    
    tmp1 <- Bolstad::bayes.t.test(yes[,trts[t]], no[,trts[t]],
                                  mu = 0)
    
    ts[t, 4] <- round(tmp1$p.value, 3)
    ts[t, 2] <- tmp1$estimate[2]
    ts[t, 3] <- tmp1$estimate[1]
  }
  
  ts$trait[1] <- "obsSeenEBI"
  
  writeLines("\n")
  print("t-test of species seen vs. not seen in EBD: ")
  
  print(ts)
  
  sink()
  
  
}



save(all.cors, file = "OUT/all.cors.rdata")


