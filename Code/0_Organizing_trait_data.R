# M. Paola Barajas Barbosa et al (2023) Assembly of functional diversity in an oceanic island flora
# by M. Paola Barajas

library(tidyverse)
library(plyr)
library(tidyr)
library(dplyr)

setwd("C:/PUB_2/Code_Data_sub3")

# Organizing trait-data: from individual level trait to species level trait
# MASTER ID & TABLE
master  <- read.csv("Barajasetal_traits_status_20230323.csv")
master  <- subset(master, select= c(IDmaster,  Species1) )

Barajasetal  <- read.csv("Barajasetal_traits_status_20230323.csv")

# MAX HEIGHT
max_H <- read.csv("Trait_MaxPlantHeitght.csv")
max_H <- subset(max_H, select= c(IDmaster, Plant_height_max_meter) ) 

## LEAF TRAITS
# leaf area
leaf_area <- read.csv("Trait_leaf_area.csv")
grep("  $",leaf_area$Species1,value = TRUE)
grep(" $",leaf_area$Species1,value = TRUE)

leaf_area$Species1 <- gsub("  $","",leaf_area$Species1)
leaf_area$Species1 <- gsub(" $","",leaf_area$Species1)

leaf_area <- plyr::join(master, leaf_area, by="Species1", type="right") 
leaf_area <- subset(leaf_area, select= c(IDmaster, LeafArea.cm.2..of.a.single.leaf., Description)) 
leaf_area <- leaf_area[ !(leaf_area$Description %in% c("leaf (ephimeral)","leaf(small)")), ] # Delete ephimeral and small leaves
leaf_area$Description <- NULL

agg_leaf_area <- aggregate(. ~ IDmaster, data = leaf_area, mean) 
colnames(agg_leaf_area)[2] <- "leaf_area_cm2"

# LDMC 
LDMC <- read.csv("Trait_LDMC.csv")
LDMC <- subset(LDMC, select= c(IDmaster, Leaf.dry.mass..g...Formula..G2.E2, LDMC..mg.g.)) 
colnames(LDMC)[2] <- "Leaf.dry.mass.g"
colnames(LDMC)[3] <- "LDMC_mgg"

agg_leaf_dmc <- aggregate(. ~ IDmaster, data = LDMC, mean) 
agg_leaf_dmc <- subset(agg_leaf_dmc, select = c(IDmaster, LDMC_mgg, Leaf.dry.mass.g))

# Leaf LMA
LMA          <- plyr::join(agg_leaf_dmc, agg_leaf_area, by="IDmaster", type="right")# Leaf drymass + area
LMA$LMA_gm2  <- LMA$Leaf.dry.mass.g/(LMA$leaf_area_cm2*0.0001) # g/m2
LMA          <- subset(LMA, select = c(IDmaster, LMA_gm2))

# Leaf thickness
leaf_th       <- read.csv("Trait_Leaf_thickness.csv")
leaf_th       <- subset(leaf_th, select= c( IDmaster, indiv_Leaf.thickness.mm.)) 

agg_leaf_th   <- aggregate(.~ IDmaster, data=leaf_th, mean)                                                 
colnames(agg_leaf_th)[2] <-"leaf_th_mm"

# Leaf Nitrogen
N <- read.csv("Trait_Leaf_N.csv")
N <- subset(N, select= c(Species,IDmaster, Ntotal.mg.g.)) 
N$Species <-NULL 

agg_N <- aggregate(. ~ IDmaster, data = N, mean) 

# STEM SPECIFIC DENSITY
ssd     <- read.csv("Trait_Stem_specific_density.csv")
ssd     <- subset(ssd, select= c(IDmaster,  Stem.specific.density..SSD...mg.mm.3..))
colnames(ssd)[2] <- "ssd_mgmm3"

agg_ssd <- aggregate(. ~ IDmaster, data = ssd, mean)

# SEED MASS
seed    <- read.csv("Trait_Seed_mass.csv")
seed$seedmass_mg <- (seed$Seeds.dry.mass..g.*1000)/seed$Amount.of.seeds

seed     <- subset(seed, select= c( IDmaster, seedmass_mg)) 
agg_seed <- aggregate(. ~ IDmaster, data = seed, mean)

# Trait data all checked!
# FIN