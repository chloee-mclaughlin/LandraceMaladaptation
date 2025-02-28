## Substitution analysis for Barley ##

##Adapted from Caproni et al. 2022 DOI: 10.1111/gcb.16560

# Set working directory 
setwd("~/OneDrive - The Pennsylvania State University/Research/Food Resilience/NatComm_R1/data/")

# Load libraries 
library("geosphere")
library("dbscan")
library("fields")
library("maps")
library("shape")
library("conformal")
library("extendedForest")
library("gradientForest")
library("ggplot2")
library("tidyr")
library("dplyr")
library("cowplot")

## Includes GF offset values, which we will use to identify what locations are the most vulnerable 
b.150offset <- read.csv("../output/Offset/BarleyOffset_ur_150_FINAL.csv")
b.150offset$PixelID <- paste(b.150offset$Longitude, b.150offset$Latitude,sep="_")

b.accession <- b.150offset[,c(2:4,6,21)] ##for later in the script 

## Aspects of climate in the GF model
climate.aspect <- c("maturity_days", 
                    "water_cum_adj", "water_veg_adj", "water_repro_adj",
                    "tmpavg_cum_adj", "tmpavg_veg_adj", "tmpavg_repro_adj",
                    "cold_cum_adj", "cold_veg_adj", "cold_repro_adj", 
                    "solar_cum_adj", "solar_veg_adj", "solar_repro_adj") 

## Climate data for the 150 Tg scenario
target <- read.csv("../data/GF_data/Barley_targets150_FINAL.csv")

i=7 # Offset for 2 years post strike (soot injected in year 5, so year 7 is 2 years post-strike)

# Format the target climate data for the loop of finding substitutions
store.target <- cbind(target[,3:5], target[paste(climate.aspect,'_', i, sep='')])
store.target$PixelID <- paste(store.target$Longitude, store.target$Latitude,sep="_")
target.sub <- store.target
colnames(target.sub) <- sub(i , 'avg', colnames(target.sub)) ##Change "_7" for to "_avg" for GF later
colnames(target.sub) <- gsub("veg_adj", "veg", colnames(target.sub))

## Control data for identifying substitutions
control <- read.csv("../data/GF_data/Barley_control_FINAL.csv")

for (k in 1:length(climate.aspect)){
  cols <- control[grep(climate.aspect[k], colnames(control))]
  control$mean <- rowMeans(cols, na.rm = TRUE, dims = 1)
  
  names(control)[names(control) == 'mean'] <-
    paste(climate.aspect[k],'_avg', sep='')
}

# Format the control climate data for the loop of finding substitutions
store.control <- control[paste(climate.aspect,'_avg', sep='')]
store.control <- cbind(control[,3:5], store.control)
store.control$PixelID <- paste(store.control$Longitude, store.control$Latitude,sep="_")
control.sub <- store.control
colnames(control.sub) <- gsub("veg_adj", "veg", colnames(control.sub))

## GF control model 
gf.Barley <- readRDS("../output/GF_RDS/GF_Barley.rds")

# Variables in the model
imp.vars <- names(gf.Barley$overall.imp)

# Predicted GF model outcomes for control conditions
pred.gf <- predict(gf.Barley, control.sub[, imp.vars])
pred.gf$PixelID <- control.sub$PixelID 

# Projected GF model outcomes for target (150 Tg) conditions
proj.gf <- predict(gf.Barley, target.sub[, imp.vars]) 
proj.gf$PixelID <- b.150offset$PixelID 

# For looping through conditions of the worst year (Year 2 post-strike) across pixels
listModL <- colnames(b.150offset[c(7)]) #"offset_7" worst year 

#############################################
# Infer and plot migration for climate model#
#############################################

RES_Migr <- NULL
modL <- "offset_7"  ##for now just doing the worst year

for (modL in listModL){ 
  
  # Masking points that do not reach maturity
  # We only want to look for substitutions where agriculture is possible 
  datModL_prep <- cbind(b.150offset, target.sub[,c("maturity_days_avg")])
  colnames(datModL_prep)[22] <- "maturity_days_avg"
  datModL_prep_1 <- datModL_prep %>% filter(maturity_days_avg < 364)
  datModL <- datModL_prep_1[, c("Longitude","Latitude", "PixelID", modL)]
  
  # Use all pixels, many are dropped in the masking step (barley Cycles model rarely reaches maturity from extreme conditions under the 150 Tg scenario)
  get.pc <-  quantile(datModL[,4], .00)
  dat_maxVuln <- datModL[which(datModL[,4] >= get.pc),] 
  
  # Calculate the distance between the most vulnerable pixels for clustering
  dist_pts_maxVuln <- distm(dat_maxVuln[,c(1,2)], fun=distGeo)/1000
  
  # dbscan clustering, identifies geographic clusters of most vulnerable pixels 
  db <- dbscan(dist_pts_maxVuln, eps = 1000, minPts = 1) # Low min points b/c so many accessions masked in previous step
  dat_maxVuln$groups <- db$cluster
  
  # Extract the highest offset value by cluster
  Extract <- do.call(rbind, lapply(split(dat_maxVuln, dat_maxVuln$groups), 
                                   function(x) {return(x[which.max(x[,modL]),])}))
  
  ## Identifies the point for each cluster that will be run through the next step
  selPixel <- Extract$PixelID

  colnames(Extract)[4] <- "GenVulnMax"
  Extract$model <- modL
  
  # Loops through process for each of the most vulnerable pixels 
  for (Pix in selPixel){
    
    # Extract the projected genomic composition (PGC) at the selected pixel in the target climate 
    ProjmaxVuln.gf <- proj.gf[proj.gf$PixelID==Pix,] # projected conditions at the selected pixel under the 150 Tg scenario
    
    ProjmaxVuln.gf <- as.data.frame(ProjmaxVuln.gf)
    
    temp <- vector("numeric", length = nrow(pred.gf))
    
    # Calculate the Euclidean distance between the PGC at the selected pixel in the target climate  
    # and the PGC over all pixels in the current climate (best matched current genotype for the target climate)
    for (i in imp.vars) {
      # projected genomic composition of selPixel - predicted gf of all points
      temp <- temp + (ProjmaxVuln.gf[,i] - pred.gf[,i])^2 
      }
    
    Dist <- cbind(pred.gf, sqrt(temp))
    
    colnames(Dist)[colnames(Dist) == "sqrt(temp)"] = "EDist_to_maxVuln"
    Dist <- cbind(control.sub[,c(1:3)], Dist)
    
    ## Make a country column for within country migration 
    Dist$country <- map.where("world", Dist$Longitude, Dist$Latitude)
    Dist$country[is.na(Dist$country)] <- "Unknown"
    
    ## For all ED, compute best w/in country substitution 
    # Extract the smallest ED within the same country 
    sel_country <- Dist[Dist$PixelID==Pix,]$country[1] 
    sel_country_df <- Dist %>% filter(country == sel_country)
    
    Dist2Pix <- NULL
    for (i in 1:nrow(sel_country_df)){
      Dist2Pix[i] <- distGeo(c(Dist[Dist$PixelID==Pix,]$Longitude[1],
                               Dist[Dist$PixelID==Pix,]$Latitude[1]), 
                             c(sel_country_df$Longitude[i], sel_country_df$Latitude[i]))/1000
    }
    
    dat_EDmin_country <- cbind(sel_country_df, Dist2Pix)

    # Save information on the EDist of best overall and within countrysubstitutions
    Dist_sub <- Dist[,c("PixelID","Latitude", "Longitude", "EDist_to_maxVuln")]
    dat_EDmin_country_sub <- dat_EDmin_country[,c("PixelID","Latitude","Longitude","EDist_to_maxVuln", "Dist2Pix", "country")]

    ResUniPixel <- as.data.frame(c(Extract[Extract$PixelID==Pix,],
                                   Dist_sub[which.min(Dist_sub$EDist_to_maxVuln),],
                                   distGeo(c(Extract[Extract$PixelID==Pix,]$Longitude, 
                                             Extract[Extract$PixelID==Pix,]$Latitude),
                                           c(Dist_sub[which.min(Dist_sub$EDist_to_maxVuln),]$Longitude,
                                             Dist_sub[which.min(Dist_sub$EDist_to_maxVuln),]$Latitude))/1000,
                                   dat_EDmin_country_sub[which.min(dat_EDmin_country_sub$EDist_to_maxVuln),]
                                                         )) 
    
    colnames(ResUniPixel) <- c("lon.VulnPix", "lat.VulnPix", "VulnPixID", "GenVulnMax", "groups", "model",
                               "MinVulnPixID", "lat.MinVul", "lon.MinVul", "EDist_to_maxVuln", "geoDist2Vuln", 
                               "MinVulnPixID_country", "lat.MinVul_country", "lon.MinVul_country", "EDist_to_maxVuln_country", "geoDist2Vuln_country", "country")
    
    RES_Migr <- rbind(RES_Migr, ResUniPixel)
  }
}

colramp = colorRampPalette(c("#838bc5", "#a0a6d2", "#ba9bc9", "#d0947e","#c58658", "#c9403b"))(30)

temp <- NULL
temp <- RES_Migr[RES_Migr$model == modL,]

#Identify quantile of maximum distances for optimal migration and assign each "bin" a color on the ramp
qtn <- quantile(temp$EDist_to_maxVuln)
qtn <- round(qtn, 5)
#qtn_vals <- c(temp$EDist_to_maxVuln, temp$EDist_to_maxVuln_country)
#qtn <- quantile(qtn_vals)
#qtn <- round(qtn, 5)

## Assign color for best substitution 
temp$Class <- colramp[30] 
temp$Class[temp$EDist_to_maxVuln < qtn[4]] <- colramp[23] 
temp$Class[temp$EDist_to_maxVuln < qtn[3]] <- colramp[15] 
temp$Class[temp$EDist_to_maxVuln < qtn[2]] <- colramp[8]  
temp$Class[temp$EDist_to_maxVuln < qtn[1]] <- colramp[1] 

## Assing color for best within country substitution 
temp$Class_country <- colramp[30] 
temp$Class_country[temp$EDist_to_maxVuln_country < qtn[4]] <- colramp[23] 
temp$Class_country[temp$EDist_to_maxVuln_country < qtn[3]] <- colramp[15] 
temp$Class_country[temp$EDist_to_maxVuln_country < qtn[2]] <- colramp[8] 
temp$Class_country[temp$EDist_to_maxVuln_country < qtn[1]] <- colramp[1]  

#### Migration plot
## Large black dot is the vulnerable pixel and end point of line is the ideal migration, colored by the Euclidean distance of migration (how good of a sub)

# Base map data
world_map <- map_data("world")

# Best substitution
best_sub_b <- 
  ggplot() + 
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "white", color = "grey15", size = 0.3) +
  # Add vulnerable points
  geom_point(data = temp, aes(x = lon.VulnPix, y = lat.VulnPix), 
             shape = 16, color = "black", size=3) +
  # Add arrows for outline
  geom_segment(data = temp, aes(x = lon.MinVul, y = lat.MinVul, 
                                xend = lon.VulnPix, yend = lat.VulnPix), 
               arrow = arrow(length = unit(0.1, "inches")), 
               color = "black", size = 1.4) +
  # Add arrows colored by Class
  geom_segment(data = temp, aes(x = lon.MinVul, y = lat.MinVul, 
                                xend = lon.VulnPix, yend = lat.VulnPix, 
                                color = Class), 
               arrow = arrow(length = unit(0.1, "inches")), 
               size = 1) +
  coord_fixed(xlim = c(-20, 70), ylim = c(0, 65), expand = FALSE) + #Barley
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "none") +
  xlab("") +
  ylab("") +
  theme(axis.text = element_blank(),                             
        axis.ticks = element_blank()) +
  # Map the Class colors
  scale_color_identity()  +
  labs(title = "Barley",
       subtitle = "Optimal substitution") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, face = "bold", hjust = 0.5))  

# Best within country substitution
country_sub_b <-
ggplot() + 
  # Map 
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "white", color = "grey15", size = 0.3) +
  # Add vulnerable points
  geom_point(data = temp, aes(x = lon.VulnPix, y = lat.VulnPix), 
             shape = 16, color = "black", size=3) +
  # Add arrows for outline
  geom_segment(data = temp, aes(x = lon.MinVul_country, y = lat.MinVul_country, 
                                xend = lon.VulnPix, yend = lat.VulnPix), 
               arrow = arrow(length = unit(0.1, "inches")), 
               color = "black", size = 1.4) +
  # Add arrows colored by Class
  geom_segment(data = temp, aes(x = lon.MinVul_country, y = lat.MinVul_country, 
                                xend = lon.VulnPix, yend = lat.VulnPix, 
                                color = Class_country), 
               arrow = arrow(length = unit(0.1, "inches")), 
               size = 1) +
  coord_fixed(xlim = c(-20, 70), ylim = c(0, 65), expand = FALSE) + #Barley
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "none") +
  xlab("") +
  ylab("") +
  theme(axis.text = element_blank(),                             
        axis.ticks = element_blank()) +
  # Map the Class colors
  scale_color_identity() +
  labs(title = "Barley",
       subtitle = "Within country best substitution") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, face = "bold", hjust = 0.5))

# Histogram for how common substitutions fall within different levels of remaining offset 
best_hist_b <-
  ggplot(temp, aes(x=geoDist2Vuln, fill=Class)) + #best sub
  geom_histogram(binwidth = 1000) + 
  xlab("Substitution distance (km)") +
  theme_bw() +
  scale_fill_identity()
  
country_hist_b <-
ggplot(temp, aes(x=geoDist2Vuln_country, fill=Class_country)) + #best within country sub
  geom_histogram(binwidth = 200) + 
  xlab("Substitution distance (km)") +
  theme_bw() +
  scale_fill_identity()

best_grob_b <- ggplotGrob((best_hist_b) +
                            theme(plot.background = element_rect(colour = "black")))

country_grob_b <- ggplotGrob((country_hist_b) +
                               theme(plot.background = element_rect(colour = "black")))


# Add the inset to the main plot
b_best_plot <-
  best_sub_b +
  annotation_custom(
    grob = best_grob_b,
    xmin = -20, xmax = 10,  # Adjust position and size of inset
    ymin = 0, ymax = 25
  )

b_country_plot <-
  country_sub_b +
  annotation_custom(
    grob = country_grob_b,
    xmin = -20, xmax = 10,  # Adjust position and size of inset
    ymin = 0, ymax = 25
  )

b_fig4 <- plot_grid(b_best_plot, b_country_plot, labels = c("(a)", "(b)"))

Remaining Offset after substitution 
"#838BC5" = "< 0.007", 
"#A5A3D0" = "0.007 - 0.010" , 
"#C398A9" = "0.010 - 0.012", 
"#C7885F" = "0.012 - 0.014",
"#C9403B" = "> 0.014                           

