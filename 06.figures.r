

# Lane Scher
# clanescher@gmail.com
# 11/21/2020

# This script produces final figures

library(ggpubr)


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
#                 figure 1                 #
#                                          #
############################################

### map and data quantity over time


#### A
colors <- c("darkolivegreen3", "tan3", "deepskyblue4")



# get BCRs
BCR <- sf::st_read("X:/eBird-BBS/DATA/BCR/bcr_shp/BCR.shp")

BCR$color <- "gray95"
BCR$color[which(BCR$BCRNumber == 29)] <- colors[1]
BCR$color[which(BCR$BCRNumber == 19)] <- colors[2]
BCR$color[which(BCR$BCRNumber == 5)] <- colors[3]


# get states
states <- map_data("usa")

sf::sf_use_s2(FALSE)

us_sf <- sf::st_as_sf(maps::map("state", plot = FALSE, fill = TRUE)) %>%
  sf::st_transform(crs = sf::st_crs(BCR)) %>%
  sf::st_buffer(dist = 0)


BCR <- sf::st_intersection(BCR, us_sf)


bcr29 <- BCR %>%
  filter(BCRNumber == 29,
         ID %in% c("virginia", 'north carolina', 'south carolina', 'georgia', 'alabama')) 

bcr19 <- BCR %>%
  filter(BCRNumber == 19)

bcr5 <- BCR %>%
  filter(BCRNumber == 5,
         ID %in% c('washington', 'oregon'))

a <- ggplot() +
  geom_polygon(data = states, aes(x = long, y = lat, group = group), fill = "gray97", color = "black") +
  geom_sf(data = us_sf, color = "gray20", fill = NA) +
  scale_fill_identity(guide = "legend",
                      labels = c("Southeast", "Pacific Northwest", "Central Plains"),
                      name = "Regions") +
  geom_sf(data = bcr5, #map = states,
          aes(fill = color), 
          color = NA) +
  geom_sf(data = bcr19, #map = states,
          aes(fill = color), 
          color = NA) +
  geom_sf(data = bcr29, #map = states,
          aes(fill = color), 
          color = NA) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.background = element_blank(),
        legend.position = "none") +
  labs(x = "", y = "", title = "") +
  coord_sf()





#### B

regions <- c("SE", "PN", "MW")
counts <- list()
for (r in 1:length(regions)) {
  load(paste0("DATA/stats-", regions[r], ".rdata"))
  
  countInfo[[1]]$region <- regions[r]

  
  counts[[r]] <- countInfo[[1]]
}
counts <- dplyr::bind_rows(counts)

countsM <- reshape2::melt(counts, id.vars = c("year", "region"))

countsAll <- countsM[which(countsM$variable %in% c("BBSall", "EBIall")),]

fig1base <- ggplot() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank()) +
  scale_x_continuous(breaks = c(2010, 2012, 2014, 2016, 2018)) +
  scale_color_manual(values = c(colors[2], colors[3], colors[1]),
                     name = "Region", 
                     labels = c("Plains", "Pacific Northwest", "Southeast")) +
  scale_linetype_manual(name = "Source", 
                        values = c("dotted", "solid"),
                        labels = c("BBS", "eBird")) +
  labs(x = "Year", y = "Count")




countsSubset <- countsM[which(countsM$variable %in% c("BBSsubset", "EBIsubset")),]

countsSubset$value[which(countsSubset$variable == "EBIsubset")] <- 
  countsSubset$value[which(countsSubset$variable == "EBIsubset")]*3

b <- fig1base +
  geom_path(data = countsSubset,
            aes(x = year,
                y = value,
                color = region,
                linetype = variable),
            size = 1.5) +
  scale_y_continuous(name = "BBS point counts",
                     sec.axis = sec_axis(trans = ~./3, name = "eBird checklists")) +
  coord_cartesian(ylim = c(0, 2000)) +
  theme(legend.position = "none")


fig1 <- ggarrange(a, b, ncol = 2,
                  labels = c("a", "b"))

ggsave(fig1, file = "OUT/figures/fig1.jpg",
       height = 4, width = 12)

############################################
#                                          #
#                 figure 2                 #
#                                          #
############################################


# also get inputs
region <- region1 <- "SE"
propObs <- 0.05
source("CODE/04a.gjamPrep.r")

load(paste0("OUT/", region, "-traitdata.rdata"))

load(paste0("OUT/gjamOutput/0.05-", region, "_obs_plat-1-30000-censor50.rdata"))

# #### get observer effect ####
spObs <- makeBetas(out)$stand %>%
  filter(predictor == "platformEBI")


# merge with traits
spObs <- merge(spObs, traits, by.x = "species", by.y = "bird")

### betas for observer (a) and sunrise (b)

pointSize <- 1.3
errorBarSize <- 0.4
errorBarWidth <- 0.5
x.text <- 3.5
y.lab <- 3.2

fig2Base <- ggplot() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(t = 0.5,
                             r = 0.5,
                             b = 1,
                             l = 0.5, 'cm'))
  

#### A: observer

spObs <- spObs %>%
  mutate(group = case_when(sig95 == "*" &
                             Estimate > 0 ~ "pos",
                           sig95 == "*" &
                             Estimate < 0 ~ "neg",
                           sig95 == "" ~ "zero"))

a <- fig2Base +
  geom_vline(xintercept = 0, color = "gray30") +
  geom_point(data = spObs, aes(y = reorder(nameCom, Estimate),
                               x = Estimate,
                               color = group),
             size = pointSize) +
  geom_errorbar(data = spObs, aes(y = reorder(nameCom, Estimate),
                                  xmin = CI_975, xmax = CI_025,
                                  color = group),
                size = errorBarSize, width = errorBarWidth) +
  scale_color_manual(values = c("pos" = "darkblue",
                                "neg" = "brown4",
                                "zero" = "gray70")) +
  labs(y = "", x = "Platform effect") +
  geom_text(data = spObs, aes(y = reorder(nameCom, Estimate),
                              x = CI_025-0.02,
                              label = nameCom,
                              angle = 0,
                              hjust = 1,
                              vjust = 0.2,
                              color = group),
            size = x.text) + 
  geom_segment(aes(y = -2.5, yend = -2.5, 
                   x = 0.1, xend = 0.6),
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text(aes(y = -3.5, x = 0.38), label = "reported more \nin eBird",
            size = y.lab) + 
  geom_segment(aes(y = -2.5, yend = -2.5, 
                  x = -0.1, xend = -0.8),
              arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text(aes(y = -3.5, x = -0.4), label = "reported more \nin BBS",
            size = y.lab) +
  coord_cartesian(xlim = c(-0.85, 0.65), 
                  ylim = c(0.5, nrow(spObs) + 0.5),
                  clip = "off") +
  theme(plot.margin = margin(t = 0.5,
                             r = 1,
                             b = 2,
                             l = 0.5, 'cm'))



#### B: sunrise

betas <- makeBetas(out)$stand %>%
  filter(predictor == "timeCatsunrise")


betas <- merge(betas, traits, by.x = "species", by.y = "bird")

betas <- betas %>%
  mutate(group = case_when(sig95 == "*" &
                             Estimate > 0 ~ "pos",
                           sig95 == "*" &
                             Estimate < 0 ~ "neg",
                           sig95 == "" ~ "zero"))

b <- fig2Base +
  geom_vline(xintercept = 0, color = "gray30") +
  geom_point(data = betas, aes(y = reorder(nameCom, Estimate),
                                x = Estimate,
                                color = group),
             size = pointSize) +
  geom_errorbar(data = betas, aes(y = reorder(nameCom, Estimate),
                                xmin = CI_025, xmax = CI_975,
                                color = group),
                size = errorBarSize, width = errorBarWidth) +
  scale_color_manual(values = c("pos" = "darkblue",
                                "neg" = "brown4",
                                "zero" = "gray70")) +
  labs(y = "", x = "Sunrise effect") +
  geom_text(data = betas, aes(y = reorder(nameCom, Estimate),
                               x = CI_025-0.03,
                               label = nameCom,
                               angle = 0,
                               hjust = 1,
                               vjust = 0.2,
                               color = group),
            size = x.text) + 
  geom_segment(aes(y = -2.5, yend = -2.5, 
                   x = 0.1, xend = 0.3),
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text(aes(y = -3.5, x = 0.2), label = "reported more \nat sunrise",
            size = y.lab) + 
  geom_segment(aes(y = -2.5, yend = -2.5, 
                   x = -0.1, xend = -1.3),
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text(aes(y = -3.5, x = -0.65), label = "reported more \nafter sunrise",
            size = y.lab) + 
  coord_cartesian(xlim = c(-1.3, 0.3), 
                  ylim = c(0.5, nrow(spObs)+0.5),
                  clip = "off") +
  theme(plot.margin = margin(t = 0.5,
                             r = 1,
                             b = 2,
                             l = 0.5, 'cm'))


#### TOGETHER

fig2 <- ggarrange(a, b, 
          nrow = 1, ncol = 2,
          labels = c("a", "b"))

ggsave(fig2, filename = "OUT/figures/fig2.jpg",
       height = 10, width = 7.5)




############################################
#                                          #
#                 figure 3                 #
#                                          #
############################################


load("OUT/all.cors.rdata")
all.cors1 <- all.cors %>%
  pivot_longer(!trait) %>%
  mutate(region = substr(name, 1, 2),
         pred = gsub(".*-", "", name),
         model = substr(name, 4, nchar(name)),
         model = gsub("-platform|-sunrise", "", model),
         Model = case_when(model == "obsplat" ~ "Platform + observer\n ",
                           model == "plat" ~ "Platform\n ",
                           model == "obsplat-matching" ~ "Platform + observer\nMatching",
                           model == "plat-matching" ~ "Platform\nMatching"),
         Pred = case_when(pred == "platform" ~ "Platform",
                          pred == "sunrise" ~ "Sunrise"),
         Region = case_when(region == "SE" ~ "Southeast",
                            region == "MW" ~ "Plains",
                            region == "PN" ~ "Pacific"),
         Trait = case_when(trait == "vis" ~ "Percent visual detections",
                           trait == "siz" ~ "Body size",
                           trait == 'pop' ~ "Popularity",
                           trait == "hei" ~ "Foraging height",
                           trait == "com" ~ "Commonness",
                           trait == "col" ~ "Color contrast"),
         Value = round(value, 2),
         text.color = case_when(value > 0.45 ~ "light",
                                value < -0.6 ~ "light",
                                T ~ "dark"))


# full figure goes in the appendix
ggplot(all.cors1) +
  geom_point(aes(x = Region, y = Trait, fill = value), size = 20, shape = 22) +
  geom_text(aes(x = Region, y = Trait, label = Value, color = text.color), hjust = 0.5, vjust = 0.5, show.legend = F) +
  facet_grid(rows = vars(Pred), cols = vars(Model)) +
  scale_fill_gradient2(high = "darkblue", low = "red4", mid = "white", name = "Correlation") +
  scale_color_manual(values = c('gray10', 'gray90')) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file = "OUT/figures/allcorrelations.jpg", height = 8, width = 10)



# this figure goes in the main text

labs <- data.frame(lab = c("a", "b"),
                   x = 2,
                   y = 5,
                   Pred = c("Platform", "Sunrise"))


plat <- all.cors1 %>%
  filter(model %in% c("obsplat"),
         Pred == "Platform") %>%
  ggplot() +
  geom_point(aes(x = Region, y = Trait, fill = value), size = 20, shape = 22) +
  geom_text(aes(x = Region, y = Trait, label = Value, color = text.color), hjust = 0.5, vjust = 0.5, show.legend = F) +
  facet_wrap(~Pred, nrow = 1) +
  scale_fill_gradient2(high = "darkblue", low = "red4", mid = "white", name = "Correlation") +
  scale_color_manual(values = c('gray10', 'gray90')) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = '', y = "") +
  coord_cartesian(clip = "off")

sunrise <- all.cors1 %>%
  filter(model %in% c("obsplat"),
         Pred == "Sunrise") %>%
  ggplot() +
  geom_point(aes(x = Region, y = Trait, fill = value), size = 20, shape = 22) +
  geom_text(aes(x = Region, y = Trait, label = Value, color = text.color), hjust = 0.5, vjust = 0.5, show.legend = F) +
  facet_wrap(~Pred, nrow = 1) +
  scale_fill_gradient2(high = "darkblue", low = "red4", mid = "white", name = "Correlation") +
  scale_color_manual(values = c('gray10', 'gray90')) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_blank()) +
  labs(x = '', y = "") + 
  coord_cartesian(clip = "off")

ggarrange(egg::tag_facet_outside(plat, tag_pool = "a"), 
          nrow = 1, ncol = 2,
          \legend = "right", common.legend = T)

a <- egg::tag_facet(plat, open = "", close = "",
                    hjust = 2, vjust = -0.4,
                    tag_pool = "a")
b <- egg::tag_facet(sunrise, open = "", close = "",
                    hjust = 2, vjust = -0.4,
                    tag_pool = "b")

ggarrange(a, b, nrow = 1, 
          legend = "right", common.legend = T,
          widths = c(1.68,1))+
  theme(plot.margin = margin(1,0,0,0, "cm")) 


ggsave(file = "OUT/figures/fig3.jpg", height = 4.6, width = 6.2)





############################################
#                                          #
#                 figure 4                 #
#                                          #
############################################

# 3 panels, observer effect vs. commonness


sp.remove <- read.csv("OUT/tables/gelman_plat.csv")

regs <- c("SE", "MW", "PN")

obsEff <- list()
for (r in 1:length(regs)) {
  load(paste0("OUT/gjamOutput/", propObs, "-",  regs[r], "_obs_plat-1-30000-censor50.rdata"))
  plat <- makeBetas(out)$stand %>%
    filter(predictor == "platformEBI")
  
  sp.rem <- sp.remove %>%
    select(species, paste0(regs[r], ".obs_plat"))
  sp.rem <- sp.rem$species[which(sp.rem[,2] > 1.2)]
  
  if (length(sp.rem > 0)) {
    plat <- plat[-which(plat$species %in% sp.rem),]
  }
  
  plat$region <- regs[r]
  
  load(paste0("OUT/", regs[r], "-traitdata.rdata"))
  
  plat <- merge(plat, traits, by.x = "species", by.y = "bird")
  
  obsEff[[r]] <- plat
}


obsEff1 <- bind_rows(obsEff)


obsEff1$bins <- cut(obsEff1$obsSeen,breaks = 12)

obsEff1 <- obsEff1 %>%
  mutate(Region = case_when(region == "PN" ~ "Pacific",
                            region == "MW" ~ "Plains",
                            region == "SE" ~ "Southeast"))

arr.df <- data.frame(x = c(0.9, 0.9),
                     y = c(0.1, -0.1),
                     xend = c(0.9, 0.9),
                     yend = c(0.55, -0.55),
                     Region = "Pacific")
lab.df <- data.frame(x = c(1.8, 1.8),
                     y = c(0.325, -0.325),
                     label = c("reported more \nin eBird",
                               "reported more \nin BBS"),
                     Region = "Pacific")

fig3a <- ggplot(obsEff1) +
  geom_hline(yintercept = 0, color = "gray50") +
  geom_boxplot(aes(x = bins,
                   y = Estimate),
               fill = "gray90",
               outlier.shape = NA) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 11)) +
  labs(x = "Commonness", y = "Platform effect",
       color = "Species") +
  geom_text(aes(x = 2.5, y = -0.6),
            label = "rare",
            size = 3.5) +
  geom_text(aes(x = 10, y = -0.6),
            label = "common", 
            size = 3.5) +
  facet_wrap(~factor(Region,
                     levels = c("Pacific",
                                "Plains",
                                "Southeast"))) +
  coord_cartesian(ylim = c(-0.6, 0.6), 
                  xlim = c(1, 12)) +
  geom_segment(data = arr.df, 
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text(data = lab.df,
            aes(x = x, y = y, label = label),
            angle = 90,
            size = y.lab)

obsEff1$bins1 <- cut(obsEff1$percentVisDet, breaks = 12)

fig3b <- ggplot(obsEff1) +
  geom_hline(yintercept = 0, color = "gray50") +
  geom_boxplot(aes(x = bins1,
                   y = Estimate),
               fill = "gray90",
               outlier.shape = NA) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 11)) +
  labs(x = "Typical detection method", y = "Platform effect",
       color = "Species") +
  geom_text(aes(x = 2.5, y = -0.58),
            label = "detected\nby sound",
            size = 3.5) +
  geom_text(aes(x = 10, y = -0.58),
            label = "detected\nby sight", 
            size = 3.5) +
  facet_wrap(~factor(Region,
                     levels = c("Pacific",
                                "Plains",
                                "Southeast"))) +
  coord_cartesian(ylim = c(-0.6, 0.6), 
                  xlim = c(1, 12)) +
  geom_segment(data = arr.df, 
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text(data = lab.df,
            aes(x = x, y = y, label = label),
            angle = 90,
            size = y.lab)

fig3 <- ggarrange(fig3a, fig3b,
                  nrow = 2, ncol = 1,
                  labels = c("a", "b"))


ggsave(fig3, filename = "OUT/figures/fig4.jpg",
       height = 8, width = 8.5)



