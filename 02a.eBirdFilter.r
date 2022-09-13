
# Lane Scher
# clanescher@gmail.com
# 06/08/19

# This script reads in raw ebird data and reformats it as a
# dataframe. It then saves a .rdata file for each state.

library(auk)
library(lubridate)




############################################
#                                          #
#    filter ebird data and save by year    #
#                                          #
############################################

regions <- c("SE", "PN", "MW")

for (r in 1:length(regions)) {
  
  if (regions[r] == "SE") {bcrCode <- 29}
  if (regions[r] == "MW") {bcrCode <- 19}
  if (regions[r] == "PN") {bcrCode <- 5}
  
  print(regions[r])
  if (length(years) > 0) {
    
    for (y in 1:length(years)) {
      
      
      year <- years[y]
      print(year)
      
      if (file.exists(paste0("X:/eBird-BBS/DATA/dataEbird/processed/", year, "-", regions[r], ".rdata"))) next
      
      #### set folders
      # read in files from
      ebd <- auk_ebd(file = "X:/data-eBird/Apr-2022/ebd_relApr-2022.txt",
                     file_sampling = "X:/data-eBird/Apr-2022/ebd_sampling_relApr-2022.txt")
      
      # save files to
      f_ebd <- "DATA/dataEbird/raw/ebd.txt"
      f_sampling <- "DATA/dataEbird/raw/sampling.txt"
      
      
      #### set filters
      # complete
      ebd_filters <- auk_complete(ebd)
      
      # date
      start <- paste0(year, "-04-01")
      end <- paste0(year, "-08-15")
      ebd_filters <- auk_date(ebd_filters,
                              date = c(start, end))
      
      # protocol
      ebd_filters <- auk_protocol(ebd_filters, protocol = c("Stationary", "Traveling"))
      
      # duration
      ebd_filters <- auk_duration(ebd_filters, duration = c(0, 300))
      
      # distance
      ebd_filters <- auk_distance(ebd_filters, distance = c(0, 5))
      
      # states
      #ebd_filters <- auk_state(ebd_filters, state = c(stateCodes))
      
      # BCR
      ebd_filters <- auk_bcr(ebd_filters, bcr = bcrCode)
      
      # country
      ebd_filters <- auk_country(ebd_filters, country = "Us")
      
      #### filter data
      ebd_filtered <- auk_filter(ebd_filters, 
                                 file = f_ebd, 
                                 file_sampling = f_sampling,
                                 overwrite = T)
      
      #### read filtered ebd and sampling
      ebd_zf <- auk_zerofill(f_ebd, f_sampling, collapse = T)
      
      #### clean ####
      
      ebd_zf$year <- as.numeric(substr(ebd_zf$observation_date, 1, 4))
      
      ebd_df <- ebd_zf[,c("checklist_id", "observer_id", "bcr_code", "state",
                          "latitude", "longitude", "observation_date", "year",
                          "time_observations_started", "duration_minutes",
                          "effort_distance_km", "effort_area_ha", "protocol_type",
                          "number_observers", "all_species_reported", "group_identifier",
                          "scientific_name", "observation_count")]
      
      # save
      outFile <- paste0("X:/eBird-BBS/DATA/dataEbird/processed/", year, "-", regions[r], ".rdata")
      save(ebd_df, file = outFile)
      
      
      
      #### save unique checklist info ####
      
      ebd_CL <- read_sampling(f_sampling)
      
      ebd_CL$year <- year(ebd_CL$observation_date)
      ebd_CL <- ebd_CL[,c("checklist_id", "observer_id", "bcr_code", "state", 
                          "latitude", "longitude", "observation_date", "year",
                          "time_observations_started", "duration_minutes", 
                          "effort_distance_km", "effort_area_ha", "protocol_type",
                          "number_observers", "all_species_reported",
                          "group_identifier")]
      
      outFile <- paste0("X:/eBird-BBS/DATA/dataEbird/processed/", year, "-", regions[r], "-CL.rdata")
      save(ebd_CL, file = outFile)
      
    }
    
    
  }
  
}

