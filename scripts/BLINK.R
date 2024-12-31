source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/FarmCPU/FarmCPU_functions.txt")
library(BLINK)
library(stringr)

# Sorghum genetic data
geno.s <- read.table("./Genotypes/LDpruned_Genotypes/sorghum_LD_pruned_recoded.txt", header=T)
geno.s$IID <- gsub("\\..*","", geno.s$IID)

# Contains control averaged conditions for each crop species
load('../data/store.control_all.RData') 

## Restrict and sort geno file by accessions in store.control_s and drop extra INFO columns
geno.s.sub <- subset(geno.s[-c(1,3:6)], (IID %in% rownames(store.control_s)))

store.control_s <- store.control_s %>%
  tibble::rownames_to_column(var = "Accession")

AllCrops <- read.csv("../data/FoodResilience_Accessions_Waypoint.csv")
SorghumMaster <- AllCrops %>% filter(Species == "Sorghum")

# Perform the join
env.coord.s <- left_join(store.control_s, SorghumMaster, by = "Accession")

geno.s.sub.t <- t(geno.s.sub[-c(1)])
colnames(geno.s.sub.t) <- geno.s.sub$IID

geno.s.sub.map <- data.frame(
  rs = rownames(geno.sb.sub.t),
  chr   = str_extract(rownames(geno.sb.sub.t), "\\d+"),
  pos   = str_extract(rownames(geno.sb.sub.t), "(?<=_)[^_]+(?=_)")
)

env.coord <- env.coord.s


# genotype information data
#myGM = read.table("myData.map", head=T)
myGM <- geno.s.sub.map[-1, ]

# genotype data
#myGD=read.big.matrix("myData.dat",head=F,sep="\t",type="char") 
myGD <- geno.s.sub.t

cov = prcomp(as.matrix(myGD))
myCV = cov$x[,1:3]
# if you only have one CV, the format should be matrix format. 
#myCV = as.matrix(cov$x[,1])

# phenotype data
#myY = read.table("myData.txt",head = T) 
myY <- env.coord.s

# association analysis
myBlink_md <- Blink(Y=myY[,c(1,2)], GD=myGD, GM=myGM, CV = myCV, maxLoop=10, time.cal=T)

myBlink_wc <- Blink(Y=myY[,c(1,3)], GD=myGD, GM=myGM, CV = myCV, maxLoop=10, time.cal=T)
myBlink_wv <- Blink(Y=myY[,c(1,4)], GD=myGD, GM=myGM, CV = myCV, maxLoop=10, time.cal=T)
myBlink_wr <- Blink(Y=myY[,c(1,5)], GD=myGD, GM=myGM, CV = myCV, maxLoop=10, time.cal=T)

myBlink_tc <- Blink(Y=myY[,c(1,6)], GD=myGD, GM=myGM, CV = myCV, maxLoop=10, time.cal=T)
myBlink_tv <- Blink(Y=myY[,c(1,7)], GD=myGD, GM=myGM, CV = myCV, maxLoop=10, time.cal=T)
myBlink_tr <- Blink(Y=myY[,c(1,8)], GD=myGD, GM=myGM, CV = myCV, maxLoop=10, time.cal=T)

myBlink_cc <- Blink(Y=myY[,c(1,9)], GD=myGD, GM=myGM, CV = myCV, maxLoop=10, time.cal=T)
myBlink_cv <- Blink(Y=myY[,c(1,10)], GD=myGD, GM=myGM, CV = myCV, maxLoop=10, time.cal=T)
myBlink_cr <- Blink(Y=myY[,c(1,11)], GD=myGD, GM=myGM, CV = myCV, maxLoop=10, time.cal=T)

myBlink_sc <- Blink(Y=myY[,c(1,12)], GD=myGD, GM=myGM, CV = myCV, maxLoop=10, time.cal=T)
myBlink_sv <- Blink(Y=myY[,c(1,13)], GD=myGD, GM=myGM, CV = myCV, maxLoop=10, time.cal=T)
myBlink_sr <- Blink(Y=myY[,c(1,14)], GD=myGD, GM=myGM, CV = myCV, maxLoop=10, time.cal=T)

outliers_blink_md <- myBlink_md$GWAS %>% 
  filter(P.value < 0.05)

outliers_blink_sc <- myBlink_sc$GWAS %>% 
  filter(P.value < 0.05)

outliers_blink_sv <- myBlink_sv$GWAS %>% 
  filter(P.value < 0.05)

outliers_blink_sr <- myBlink_sr$GWAS %>% 
  filter(P.value < 0.05)

outliers_blink_wc <- myBlink_wc$GWAS %>% 
  filter(P.value < 0.05)

outliers_blink_wv <- myBlink_wv$GWAS %>% 
  filter(P.value < 0.05)

outliers_blink_wr <- myBlink_wr$GWAS %>% 
  filter(P.value < 0.05)

outliers_blink_tc <- myBlink_tc$GWAS %>% 
  filter(P.value < 0.05)

outliers_blink_tv <- myBlink_tv$GWAS %>% 
  filter(P.value < 0.05)

outliers_blink_tr <- myBlink_tr$GWAS %>% 
  filter(P.value < 0.05)

outliers_blink_cc <- myBlink_cc$GWAS %>% 
  filter(P.value < 0.05)

outliers_blink_cv <- myBlink_cv$GWAS %>% 
  filter(P.value < 0.05)

outliers_blink_cr <- myBlink_cr$GWAS %>% 
  filter(P.value < 0.05)

outliers_blink <- rbind(outliers_blink_md, 
                        outliers_blink_sc,
                        outliers_blink_sv, 
                        outliers_blink_sr,
                        outliers_blink_wc,
                        outliers_blink_wv, 
                        outliers_blink_wr,
                        outliers_blink_tc,
                        outliers_blink_tv, 
                        outliers_blink_tr,
                        outliers_blink_cc,
                        outliers_blink_cv, 
                        outliers_blink_cr)

out_blink_dup <- outliers_blink[!duplicated(outliers_blink$rs), ]

outliers_rdadapt_env.s <- read.csv("./Genotypes/pRDA_Outliers/S_outliers_rdadapt_env.csv")

outliers_rda <- outliers_rdadapt_env.s
colnames(outliers_rda)[colnames(outliers_rda) == "x"] <- "rs"

overlap_loci <- inner_join(outliers_rda,  out_blink_dup)
