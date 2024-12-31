############################################
##For predicted adaptation vs CG phenotype##
############################################

library(lme4)
library(cowplot)
library(tidyr)
library(dplyr)
library(rworldmap)
library(ggplot2)
library(ggpmisc)
library(conformal)
library(extendedForest)
library(gradientForest)

setwd("~/OneDrive - The Pennsylvania State University/Research/Food Resilience/NatComm_R1/data/")

#store.control_m object is built in the 2.GradientForest_AllCrops.Rmd script and can be loaded with the .RData object below
load('./store.control_all.RData') 

control.climate_m <- read.csv("../data/GF_data/Maize_control_FINAL.csv")

gf.m <- readRDS("../output/GF_RDS/GF_Maize.rds")

imp.vars <- names(gf.m$overall.imp)

pred.gf <- predict(gf.m, store.control_m[, imp.vars])

pred.gf$PixelID <- paste(control.climate_m$Latitude, control.climate_m$Longitude, sep="_")

# List of pixels Cycles env predicted for that fall in the same coordinate cell as each of the CGs 
selPixel_list = c("20.46666667_-97.33333333", "19.4_-99.55" ,"20.59568_-100.777428","27.83333333_-109.3333333",
                  "20.381461_-101.143579" , "19.521431_-98.857468", "16.458148_-93.370653", "18.35_-100.65", 
                  "20.290769_-102.319883", "21.221061_-104.739704", "19.8_-101.7833333", "18.71666667_-99.36666667",
                  "25.624139_-103.562492")

# List of the Common Garden 'locale names'
CG <- c("Agua Fria", "Amoloya de Juarez", "Celaya", "Ciudad Obregon", "Cortazar",
        "El Batan", "Guadalupe-Victoria", "Iguala", "Numaran", 
        "San Pedro Lagunillas", "Tarimbaro", "Tlaltizapan", "Torreon")

Dist_CG_all <- NULL

#For each CG, this function will predict how maladapted landraces accessions are to the CG they are grown in 

for (j in 1:length(selPixel_list)){
  selPixel = selPixel_list[j]
  
  for (Pix in selPixel){
    
    # Extract the projected genomic composition (PGC) at the selected pixel in the current climate 
    Pred.gf_Pix <- pred.gf[pred.gf$PixelID==Pix,][1,] #conditions at the selected pixel 
    Pred.gf_Pix <- as.data.frame(Pred.gf_Pix)
    
    temp_1 <- vector("numeric", length = nrow(pred.gf))
    
    # Across all env variables
    # calculate the Euclidean distance (ED) between the genomic composition (GC) at each common garden pixel
    # and the genomic composition across all pixels
    for (i in imp.vars) {
      
      #Projected GC of the selected CG Pixel - predicted GC of all locations 
      temp_1 <- temp_1 + (Pred.gf_Pix[,i]-pred.gf[,i])^2 
    }
    
    Dist_bind <- cbind(pred.gf, sqrt(temp_1))
    colnames(Dist_bind)[colnames(Dist_bind) == "sqrt(temp_1)"] = "EDist_to_CommonGarden"
    
    Dist_CG <- cbind(control.climate_m[,3:5], Dist_bind)
    Dist_CG$CG <- CG[j]
    
    Dist_CG_all <-   rbind(Dist_CG_all, Dist_CG)
    
  }}

write.csv(Dist_CG_all, "../output/CommonGarden/EDist_to_CommGarden.csv")
Dist_CG_all <- read.csv("../output/CommonGarden/EDist_to_CommGarden.csv")

## BLUPs from Gate Preprint (https://data.cimmyt.org/dataset.xhtml?persistentId=hdl:11529/10548233)
# Since these files are from a preprint, I will not include them in this Github directory. The files can be downloaded at the above link
# These are the same values as what is from CYMMIT's website, filtered for reported disease score and outlier phenotypes 

#Uncomment the CG.dat object and corresponding phenotpye for each phenotype/plot 

#CG.dat <- read.csv("./CommonGarden/Curated GWAS trial data/Anthesis silking interval_altitude_GWASResiduals.Rimage.csv", skip = 2)
#CG.dat <- read.csv("./CommonGarden/Curated GWAS trial data/Bare cob weight_altitude_GWASResiduals.Rimage.csv", skip = 2)
#CG.dat <- read.csv("./CommonGarden/Curated GWAS trial data/Days to anthesis_altitude_GWASResiduals.Rimage.csv", skip = 2)
#CG.dat <- read.csv("./CommonGarden/Curated GWAS trial data/Field weight_altitude_GWASResiduals.Rimage.csv", skip = 2)
#CG.dat <- read.csv("./CommonGarden/Curated GWAS trial data/Grain weight per hectare_altitude_GWASResiduals.Rimage.csv", skip = 2)
CG.dat <- read.csv("./CommonGarden/Curated GWAS trial data/Plant height_altitude_GWASResiduals.Rimage.csv", skip = 2)

# Added extra random effect of tester (not accounted for as individual random effect in the source data)
mod_CG <- lmer(value ~ (1|tester) +  (1|year), data = CG.dat)
resid.CG <- cbind(CG.dat, resid=residuals(mod_CG))

## Match predicted ED for each individual grown in a CG to the correct phenotype when grown in a CG 
CG.pred <- left_join(Dist_CG_all, resid.CG, by = c("Accession" = "SampleID", "CG" = "locale"), multiple = "all") %>% 
  drop_na(value) #this drops all ED predictions for individuals that weren't phenotyped in a given CG 

#ASI <- 
#BCW <- 
#DA <- 
#FW <- 
#GW <- 
PH <- 
CG.pred %>% 
  ggplot(aes(log(EDist_to_CommonGarden), resid)) + #,colour=CG, fill=CG
  geom_point(alpha=0.2, pch=1) +
  #geom_point(alpha=0.5, pch=1, aes(colour=CG, fill=CG)) +
  theme(legend.position = "none") +
  xlab("log(GF Offset) to common garden") +
  theme_bw() +
  geom_smooth(method = "lm", se = FALSE, color = "#c9403b") +  # Regression line
  stat_poly_eq(
    aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~~")),
    formula = y ~ x,  # Regression formula
    parse = TRUE) +
#  ylab("Anthesis silking interval (days)") 
#  ylab("Bare cob weight (grams)")
#  ylab("Days to anthesis (days)") 
#  ylab("Field weight (kg)") 
#  ylab("Grain weight per hectare (tons)")
  ylab("Plant height (cm)") 

plot_grid(ASI, BCW, DA, FW, GW, PH, labels= c("B","C","D","E","F","G"))



## Fig S6 Map
map.polygon <- getMap(resolution = "low") 

CG_locality <- data.frame(locale = c("Agua Fria", "Amoloya de Juarez", "Celaya", "Ciudad Obregon", "Cortazar", "El Batan", 
                                     "Guadalupe-Victoria", "Iguala", "Numaran", "San Pedro Lagunillas", "Tarimbaro", 
                                     "Tlaltizapan", "Torreon"),
                          Latitude = c(20.46, 19.41, 20.58, 27.37, 20.45, 19.53, 16.44, 18.35, 20.26, 21.22, 19.78, 18.69, 25.56), 
                          Longitude = c(-97.64, -99.74, -100.82, -109.93, -101.14, -98.85, -93.13, -99.51, -101.93, 
                                        -104.73, -101.18, -99.13, -103.37) )

mexi_map <- ggplot() +
  geom_path(data = map.polygon, aes(x = long, y = lat, group = group)) + # world map
  coord_sf(xlim = c(-120, -82), ylim = c(12, 34), expand = FALSE) + #focused area
  geom_point(data = CG.pred, aes(x = Longitude, y = Latitude), cex = 0.1) +
  geom_point(data = CG_locality, aes(x = Longitude, y = Latitude), cex = 3, color = "#c9403b") +
  ggtitle("Common garden locations") + 
  #force_panelsizes(rows = unit(4, "cm"),cols = unit(6, "cm")) + # set the panel size
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        #axis.title.y = element_text(size =4, vjust = 0.5,color = 'black',face = 'bold', family = 'sans'),
        axis.title = element_blank(),
        #axis.title.x = element_text(size =4, vjust = 0.5,color = 'black',face = 'bold', family = 'sans'),
        plot.margin = unit(c(0,0,0,0), "cm"),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = 'white'))

mexi_map


