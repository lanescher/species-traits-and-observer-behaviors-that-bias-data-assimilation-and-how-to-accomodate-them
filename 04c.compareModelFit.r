
# Lane Scher
# clanescher@gmail.com
# 8/22/2022

# This script analyzes model fits

library(tidyverse)


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


##### chain convergence ----

library(tidyverse)
regs <- c("MW", "SE", "PN")
types <- c("obs_plat", "plat")
#types <- "plat"

gelman_all_plat <- c() # save gelman indices
gelman_all_time <- c()
plat <- c() # save platform effect
sunrise <- c() # save sunrise effect
fits <- data.frame() # save fit metrics
df_by_spec <- c()

for (r in 1:length(regs)) {
  
  
  for (t in 1:length(types)) {
    if (types[t] == "obs_plat") {ng1 <- 30000}
    if (types[t] == "plat") {ng1 <- 20000}
    
    files <- list.files("OUT/gjamOutput/", pattern = regs[r])
    
    # get correct models
    files_plat <- files[grep(paste0(propObs, "-", regs[r], "_", types[t]), files)]
    files_plat <- files_plat[grep(paste0(ng1, "-censor50"), files_plat)]
    
    
    mods <- list()
    all_chains_plat <- list()
    for (run in 1:length(files_plat)) {
      load(paste0("OUT/gjamOutput/", files_plat[run]))
      
      if (exists("outRO")) {out <- outRO}
      
      #### save model output
      mods[[run]] <- out
      
      #### save chains
      all_chains_plat[[run]] <- out$chains$bgibbs
      
      #### model fit
      fits[nrow(fits)+1, 1] <- regs[r]
      fits[nrow(fits), 2] <- types[t]
      fits[nrow(fits), 3] <- run
      
      fits[nrow(fits), 4] <- out$fit$DIC
      fits[nrow(fits), 5] <- out$fit$rmspeAll

      
      if (exists("out")) {rm(out)}
      if (exists("outRO")) {rm(outRO)}
    }
    
    
    
    sp <- unique(colnames(all_chains_plat[[1]]))
    sp <- gsub("_.*", "", sp)
    sp <- unique(sp)
    
    #### pull out chains for platform effect
    
    chains_by_spec_plat <- list()
    chains_by_spec_time <- list()
    for (s in 1:length(sp)) {
      
      # get platform chains
      
      chain_sp_plat <- coda::mcmc.list()
      chain_sp_df_plat <- data.frame(tmp = rep(NA, ng1))
      
      for (run in 1:length(all_chains_plat)) {
        ind <- grep(paste0(sp[s], '_platformEBI'), colnames(all_chains_plat[[run]]))
        
        # for df
        chain <- all_chains_plat[[run]][,ind]
        chain_sp_df_plat[,run] <- chain
        
        # for mcmc.list
        chain <- coda::mcmc(all_chains_plat[[run]][,ind])
        chain_sp_plat[[run]] <- chain
        
      }
      
      colnames(chain_sp_df_plat) <- paste0("C", 1:length(all_chains_plat))
      chain_sp_df_plat <- chain_sp_df_plat %>%
        mutate(g = 1:ng1) %>%
        pivot_longer(!g) %>%
        mutate(region = regs[r],
               type = types[t],
               sp = sp[s],
               pred = "platformEBI")
      
      df_by_spec <- bind_rows(df_by_spec, chain_sp_df_plat)
      
      chains_by_spec_plat[[s]] <- chain_sp_plat
      
      
      # get sunrise chains
      
      chain_sp_time <- coda::mcmc.list()
      chain_sp_df_time <- data.frame(tmp = rep(NA, ng1))
      
      for (run in 1:length(all_chains_plat)) {
        ind <- grep(paste0(sp[s], '_timeCatsunrise'), colnames(all_chains_plat[[run]]))
        
        # for df
        chain <- all_chains_plat[[run]][,ind]
        chain_sp_df_time[,run] <- chain
        
        # for mcmc.list
        chain <- coda::mcmc(all_chains_plat[[run]][,ind])
        chain_sp_time[[run]] <- chain
        
      }
      
      colnames(chain_sp_df_time) <- paste0("C", 1:length(all_chains_plat))
      chain_sp_df_time <- chain_sp_df_time %>%
        mutate(g = 1:ng1) %>%
        pivot_longer(!g) %>%
        mutate(region = regs[r],
               type = types[t],
               sp = sp[s],
               pred = "timeCatsunrise")
      
      df_by_spec <- bind_rows(df_by_spec, chain_sp_df_time)
      
      chains_by_spec_time[[s]] <- chain_sp_time
      
      
    }
    
    names(chains_by_spec_plat) <- sp
    names(chains_by_spec_time) <- sp
    
    
    #### gelman index for each species
    for (s in 1:length(sp)) {
      
      # platform
      tmp <- coda::gelman.diag(chains_by_spec_plat[[s]])[[1]][,1]
      tmp1 <- data.frame(sp = sp[s],
                         #num = s,
                         index = tmp,
                         type = types[t],
                         region = regs[r])
      
      gelman_all_plat <- bind_rows(gelman_all_plat, tmp1)
      
      
      # time
      tmp <- coda::gelman.diag(chains_by_spec_time[[s]])[[1]][,1]
      tmp1 <- data.frame(sp = sp[s],
                         #num = s,
                         index = tmp,
                         type = types[t],
                         region = regs[r])
      
      gelman_all_time <- bind_rows(gelman_all_time, tmp1)

      
    }
    
    
    
  } # end type 
  
} # end regs

beepr::beep()

# format gelman indices
gelman_plat <- gelman_all_plat %>%
  mutate(regtype = paste0(region, "-", type)) %>%
  select(sp, index, regtype) %>%
  pivot_wider(values_from = index,
              names_from = regtype) %>%
  arrange(sp)
colnames(gelman_plat)[1] <- "species"

low <- c()
all <- c()
for (c in 2:7) {
  low <- c(low, nrow(gelman_plat[which(gelman_plat[,c] < 1.2),]))
  all <- c(all, nrow(gelman_plat[which(is.na(gelman_plat[,c]) == F),]))
}

gelman_plat <- as.data.frame(gelman_plat)
gelman_plat[nrow(gelman_plat)+1,] <- c("< 1.2", low/all)
gelman_plat[,2:7] <- apply(gelman_plat[,2:7], 2, as.numeric)
gelman_plat[,2:7] <- round(gelman_plat[,2:7], 3)
write.csv(gelman_plat, file = "OUT/tables/gelman_plat.csv", row.names = F)


gelman_time <- gelman_all_time %>%
  mutate(regtype = paste0(region, "-", type)) %>%
  select(sp, index, regtype) %>%
  pivot_wider(values_from = index,
              names_from = regtype) %>%
  arrange(sp)
colnames(gelman_time)[1] <- "species"

low <- c()
all <- c()
for (c in 2:7) {
  low <- c(low, nrow(gelman_time[which(gelman_time[,c] < 1.2),]))
  all <- c(all, nrow(gelman_time[which(is.na(gelman_time[,c]) == F),]))
}

gelman_time <- as.data.frame(gelman_time)
gelman_time[nrow(gelman_time)+1,] <- c("< 1.2", low/all)
gelman_time[,2:7] <- apply(gelman_time[,2:7], 2, as.numeric)
gelman_time[,2:7] <- round(gelman_time[,2:7], 3)
write.csv(gelman_time, file = "OUT/tables/gelman_time.csv", row.names = F)



# format fits
colnames(fits) <- c("region", "type", "run", "DIC", "RMSPE")
write.csv(fits, file = "OUT/tables/fits.csv", row.names = F)



#### Figures for appendix

# chain convergence

sp.plot <- c("cliffswallow", "commonpheasant", "easternmeadowlark", "americangoldfinch", "uplandsandpiper")

#sp.plot <- c("cassinssparrow")
subs <- df_by_spec %>%
  filter(sp %in% sp.plot)

gel_plat <- gelman_plat %>%
  filter(species %in% sp.plot) %>%
  select("MW-obs_plat", "MW-plat")
colnames(gel_plat) <- paste0("platformEBI:", colnames(gel_plat))

gel_time <- gelman_time %>%
  filter(species %in% sp.plot) %>%
  select("MW-obs_plat", "MW-plat")
colnames(gel_time) <- paste0("timeCatsunrise:", colnames(gel_time))

gel <- bind_cols(gel_plat, gel_time) %>%
  mutate(sp = sp.plot) %>%
  pivot_longer(!sp) %>%
  mutate(pred = gsub(":.*", "", name),
         type = gsub(".*-", "", name),
         lab = round(value, 2))


burnin <- gel[,c("sp", "pred", "type")] %>%
  mutate(burnin = case_when(type == "obs_plat" ~ 15000,
                            type == "plat" ~ 10000))

subs %>%
  filter(region == "MW") %>%
  ggplot() +
  geom_hline(yintercept = 0, color = "black") +
  geom_vline(data = burnin, aes(xintercept = burnin), linetype = "dashed") +
  geom_line(aes(x = g, y = value, color = name)) +
  geom_text(data = gel, aes(x = 3000, y = 1.8, label = lab)) +
  facet_grid(cols = vars(sp), rows = vars(type, pred)) +
  coord_cartesian(ylim = c(-0.5, 2)) +
  labs(x = 'ng', y = "Estimate", color = "Run")
# library(patchwork)
# a/b


ggsave("OUT/figures/chainconvergence.jpg", height = 6, width = 10)



# parameters with and without random observer

ggplot(platformEff1) +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  geom_point(aes(x = obs_plat, y = plat)) +
  labs(x = "Fixed platform effect and observer random effect", y = "Only fixed platform effect",
       title = "Platform effect")

ggplot(sunriseEff1) +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  geom_point(aes(x = obs_plat, y = plat)) +
  labs(x = "Fixed platform effect and observer random effect", y = "Only fixed platform effect",
       title = "Sunrise Effect")



# autocorrelation plots 
# plat
a1 <- subs %>%
  filter(region == "MW",
         type == "plat",
         sp == sp.plot[1],
         name == "C1",
         pred == "platformEBI") %>%
  select(value) %>%
  forecast::ggAcf() +
  labs(title = "Cliff Swallow")

a2 <- subs %>%
  filter(region == "MW",
         type == "plat",
         sp == sp.plot[2],
         name == "C1",
         pred == "platformEBI") %>%
  select(value) %>%
  forecast::ggAcf() +
  labs(title = "Common Pheasant")

a3 <- subs %>%
  filter(region == "MW",
         type == "plat",
         sp == sp.plot[3],
         name == "C1",
         pred == "platformEBI") %>%
  select(value) %>%
  forecast::ggAcf() +
  labs(title = "Eastern Meadowlark")

a4 <- subs %>%
  filter(region == "MW",
         type == "plat",
         sp == sp.plot[4],
         name == "C1",
         pred == "platformEBI") %>%
  select(value) %>%
  forecast::ggAcf() +
  labs(title = "American Goldfinch")

a5 <- subs %>%
  filter(region == "MW",
         type == "plat",
         sp == sp.plot[5],
         name == "C1",
         pred == "platformEBI") %>%
  select(value) %>%
  forecast::ggAcf() +
  labs(title = "Upland Sandpiper")

# obs_plat
b1 <- subs %>%
  filter(region == "MW",
         type == "obs_plat",
         sp == sp.plot[1],
         name == "C1",
         pred == "platformEBI") %>%
  select(value) %>%
  forecast::ggAcf() +
  labs(title = "Cliff Swallow")

b2 <- subs %>%
  filter(region == "MW",
         type == "obs_plat",
         sp == sp.plot[2],
         name == "C1",
         pred == "platformEBI") %>%
  select(value) %>%
  forecast::ggAcf() +
  labs(title = "Common Pheasant")

b3 <- subs %>%
  filter(region == "MW",
         type == "obs_plat",
         sp == sp.plot[3],
         name == "C1",
         pred == "platformEBI") %>%
  select(value) %>%
  forecast::ggAcf() +
  labs(title = "Eastern Meadowlark")

b4 <- subs %>%
  filter(region == "MW",
         type == "obs_plat",
         sp == sp.plot[4],
         name == "C1",
         pred == "platformEBI") %>%
  select(value) %>%
  forecast::ggAcf() +
  labs(title = "American Goldfinch")

b5 <- subs %>%
  filter(region == "MW",
         type == "obs_plat",
         sp == sp.plot[5],
         name == "C1",
         pred == "platformEBI") %>%
  select(value) %>%
  forecast::ggAcf() +
  labs(title = "Upland Sandpiper")

library(patchwork)
(a4 | a1 | a2 | a3 | a5) / (b4 | b1 | b2 | b3 | b5)
ggsave("OUT/figures/autocorrelation.jpg", height = 6, width = 10)
