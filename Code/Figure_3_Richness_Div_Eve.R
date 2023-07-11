# M. Paola Barajas Barbosa et al (2023) Assembly of functional diversity in an oceanic island flora
# by M. Paola Barajas & Pierre Denelle

# # Figure 3a

library(dplyr) 
library(tidyr) 
library(hypervolume) 
library(BAT) 
library(FactoMineR) 
library(tidyverse)

setwd("C:/PUB_2/Code_Data_sub3") 

# Tenerife data
Barajasetal <-  read.csv("Barajasetal_traits_status_20230323.csv")

# Nearly complete cases (5 % data missing)
data_Tenerife <- Barajasetal[, 1:16]
data_Tenerife$missing_trait <- rowSums(is.na(data_Tenerife))
data_Tenerife <- filter(data_Tenerife, missing_trait < 4)   
data_Tenerife_imp <- read.csv("Trait_imputed_phylo8.csv")
data_Tenerife_imp$Species1 <- NULL

data_Tenerife <- left_join (data_Tenerife, data_Tenerife_imp, by = "IDmaster")

data_Tenerife <- data_Tenerife %>%   # replace the few NAs with the imputed data
  mutate(LDMC.x = coalesce(LDMC.x, LDMC.y)) %>%
  mutate(LMA.x = coalesce(LMA.x, LMA.y)) %>%

  mutate(Lth.x = coalesce(Lth.x, Lth.y)) %>%
  mutate(Nmass.x = coalesce(Nmass.x, Nmass.y)) %>%

  mutate(SM.x = coalesce(SM.x, SM.y)) %>%
  mutate(SSD.x = coalesce(SSD.x, SSD.y)) %>%
  mutate(H.x = coalesce(H.x, H.y))

anyNA(data_Tenerife)

data_Tenerife$source <- 'Tenerife_data'
names(data_Tenerife)

data_Tenerife <- data.frame(
  Species1      = data_Tenerife$Species1,

  Leaf_area     = scale(log10(data_Tenerife$LA.x)),
  LMA           = scale(log10(data_Tenerife$LMA.x)),
  Leaf_N        = scale(log10(data_Tenerife$Nmass.x)),
  Seed_mass     = scale(log10(data_Tenerife$SM.x)),
  Stem_density  = scale(log10(data_Tenerife$SSD.x)),
  Height        = scale(log10(data_Tenerife$H.x)),
  LDMC          = scale(log10(data_Tenerife$LDMC.x)),
  Leaf_th       = scale(log10(data_Tenerife$Lth.x)),

  Biogeo_status= data_Tenerife$Biogeo_status,
  Ende_status  = data_Tenerife$Ende_status,
  Lineage      = data_Tenerife$Lineage,
  source       = data_Tenerife$source)

# Null model ----
# PCA analysis
pca <- FactoMineR::PCA(data_Tenerife[, c("Leaf_area", "LMA", "Leaf_N", "Seed_mass", "Stem_density", "Height", "LDMC", "Leaf_th")],
                       scale.unit = FALSE, graph = FALSE)
pca$eig
sp_coord <- as.data.frame(pca$ind$coord[, 1:3])
sp_coord <- cbind(sp_coord, data_Tenerife$Species1)
sp_coord$Species1 <- sp_coord$`data_Tenerife$Species1`
sp_coord[,4] <- NULL

data_Tenerife <- dplyr::left_join(data_Tenerife, sp_coord, by = "Species1")
rownames(data_Tenerife) <- data_Tenerife$Species1

# We calculate the three following indices: kernel.alpha, kernel.dispersion and kernel.evennes. 
# Indices based on hypervolume
all_bw <- estimate_bandwidth(data_Tenerife[, c("Dim.1", "Dim.2", "Dim.3")], # Estimate bandwidth for the island
                             method = "cross-validation")

# ####
min(table(data_Tenerife$Biogeo_status))
min_rar_bio <- 30
nbperm <- 999 # Number of permutations
####

ae  <- data_Tenerife[which(data_Tenerife$Biogeo_status == "CE"),  c("Dim.1", "Dim.2", "Dim.3")]
mac <- data_Tenerife[which(data_Tenerife$Biogeo_status == "MAC"), c("Dim.1", "Dim.2", "Dim.3")]
ns  <- data_Tenerife[which(data_Tenerife$Biogeo_status == "NNE"), c("Dim.1", "Dim.2", "Dim.3")]
sie <- data_Tenerife[which(data_Tenerife$Biogeo_status == "TE"),  c("Dim.1", "Dim.2", "Dim.3")]
all_sp <- data_Tenerife[, c("Dim.1", "Dim.2", "Dim.3")]

rar_bio_comp <- c()
for (i in 1:nbperm){
  # All species - null community
  all_sp_i <- all_sp[sample(1:nrow(all_sp), size = min_rar_bio),
                     c("Dim.1", "Dim.2", "Dim.3")]

  all_sp_i <-  hypervolume_gaussian(
    all_sp_i, kde.bandwidth = all_bw,
    quantile.requested = 0.95, quantile.requested.type = "probability",
    verbose = FALSE)

  all_sp_i <- data.frame(perm = i,
                         status = "All_species",
                         rich = kernel.alpha(all_sp_i),
                         eve = kernel.evenness(all_sp_i),
                         div = kernel.dispersion(all_sp_i))

  # Tenerife endemics
  SIE_i <- sie[sample(1:nrow(sie), size = min_rar_bio),
               c("Dim.1", "Dim.2", "Dim.3")]

  SIE_i <-  hypervolume_gaussian(
    SIE_i, kde.bandwidth = all_bw,
    quantile.requested = 0.95, quantile.requested.type = "probability",
    verbose = FALSE)

  SIE_i <- data.frame(perm = i,
                      status = "TE",
                      rich = kernel.alpha(SIE_i),
                      eve = kernel.evenness(SIE_i),
                      div = kernel.dispersion(SIE_i))

  # Canary endemics
  AE_i <- ae[sample(1:nrow(ae), size = min_rar_bio),
             c("Dim.1", "Dim.2", "Dim.3")]

  AE_i <-  hypervolume_gaussian(
    AE_i, kde.bandwidth = all_bw,
    quantile.requested = 0.95, quantile.requested.type = "probability",
    verbose = FALSE)

  AE_i <- data.frame(perm = i,
                     status = "CE",
                     rich = kernel.alpha(AE_i),
                     eve = kernel.evenness(AE_i),
                     div = kernel.dispersion(AE_i))

  # Macaronesia endemics
  MAC_i <- mac[sample(1:nrow(mac), size = min_rar_bio),
               c("Dim.1", "Dim.2", "Dim.3")]

  MAC_i <-  hypervolume_gaussian(
    MAC_i, kde.bandwidth = all_bw,
    quantile.requested = 0.95, quantile.requested.type = "probability",
    verbose = FALSE)

  MAC_i <- data.frame(perm = i,
                      status = "MAC",
                      rich = kernel.alpha(MAC_i),
                      eve = kernel.evenness(MAC_i),
                      div = kernel.dispersion(MAC_i))

  # Non endemic natives
  NS_i <- ns[sample(1:nrow(ns), size = min_rar_bio),
             c("Dim.1", "Dim.2", "Dim.3")]

  NS_i <-  hypervolume_gaussian(
    NS_i, kde.bandwidth = all_bw,
    quantile.requested = 0.95, quantile.requested.type = "probability",
    verbose = FALSE)

  NS_i <- data.frame(perm = i,
                     status = "NNE",
                     rich = kernel.alpha(NS_i),
                     eve = kernel.evenness(NS_i),
                     div = kernel.dispersion(NS_i))

  # Bind results
  rar_bio_comp <- rbind(rar_bio_comp, all_sp_i, SIE_i, AE_i, MAC_i, NS_i)

  cat(paste0(round(100*i/nbperm, 0), " %; "))
}


clado <- data_Tenerife[which(data_Tenerife$Ende_status == "Cla"), c("Dim.1", "Dim.2", "Dim.3")]
rar_end_comp <- c()
for (i in 1:nbperm){  

  # Cladogenetic species
  clado_i <- clado[sample(1:nrow(clado), size = min_rar_bio),
                   c("Dim.1", "Dim.2", "Dim.3")]

  clado_i <-  hypervolume_gaussian(
    clado_i, kde.bandwidth = all_bw,
    quantile.requested = 0.95, quantile.requested.type = "probability",
    verbose = FALSE)

  clado_i <- data.frame(perm = i,
                        status = "Cla",
                        rich = kernel.alpha(clado_i),
                        eve = kernel.evenness(clado_i),
                        div = kernel.dispersion(clado_i))

  # Bind results
  rar_end_comp <- rbind(rar_end_comp, clado_i)

  cat(paste0(round(100*i/nbperm, 0), " %; "))
}

# save(rar_bio_comp, rar_end_comp, file = "results_Figure_3a_nullmodel_20230323_2.RData") 
# # # 

load("H:/Code_Data_sub3/results_Figure_3a_nullmodel_20230323_2.RData")

rar_bio_tidy <- gather(rar_bio_comp, metric, val, c("rich", "eve", "div"))
rar_end_tidy <- gather(rar_end_comp, metric, val, c("rich", "eve", "div"))

rarefied     <- rbind(rar_bio_tidy, rar_end_tidy)
rarefied$cat <- rarefied$status

# Calculating mean, sd, se, and ci from the data
FR <- filter(rarefied, metric == "rich") 
FE <- filter(rarefied, metric == "eve") 
FD <- filter(rarefied, metric == "div") 

sum_fr <- FR %>%       
  group_by(cat) %>%
  summarise( n=n(),mean_pao =mean(val),sd_pao   =sd(val) ) %>%
  mutate( se = sd_pao/sqrt(n))  %>%
  mutate( ic = se * qt((1-0.05)/2 + .5, n-1))  %>%
  mutate( mean_pao, metric= "Functional Richness") 
sum_fe <- FE %>%
  group_by(cat) %>%
  summarise( n=n(),mean_pao =mean(val),sd_pao   =sd(val) ) %>%
  mutate( se = sd_pao/sqrt(n))  %>%
  mutate( ic = se * qt((1-0.05)/2 + .5, n-1)) %>%
  mutate( mean_pao, metric= "Functional xEvenness") 
sum_fd <- FD %>%
  group_by(cat) %>%
  summarise( n=n(),mean_pao =mean(val),sd_pao   =sd(val) ) %>%
  mutate( se = sd_pao/sqrt(n))  %>%
  mutate( ic = se * qt((1-0.05)/2 + .5, n-1)) %>%
  mutate( mean_pao, metric= "Functional zDispersion") 
sum_metrics <-  rbind(sum_fr, sum_fe, sum_fd)

# Figure 3 Null model / Standardized values 
p_my_theme2 <-  theme(axis.text.x=element_text(colour=c("black"),face="bold",size=7),
                      axis.text.y=element_text(colour=c("black"),face="bold",size=7),
                      panel.background =element_rect(fill="transparent",colour="black"),
                      panel.grid.minor=element_blank(),
                      panel.border=element_rect(fill=NA,colour="grey"))
         
coco = c("grey","gold2", "mediumorchid2","green4  ","black","dodgerblue3")

Figure_3a <- ggplot(sum_metrics) +
  geom_linerange(aes(factor(cat, level = c ("NNE","MAC","CE","TE","Cla", "All_species")), 
                     ymin= mean_pao-ic, ymax= mean_pao+ic, color = cat),  size = 0.3,  
                 show.legend = FALSE) +
  geom_point    (aes(x = cat,  
                     y=mean_pao, color = cat),  size = 1,  
                 show.legend = FALSE) +
  scale_colour_manual(values = coco) +
  facet_wrap(~ metric, scales = "free") +
  labs(x = "", y = "") +
  p_my_theme2 

Figure_3a

# FIN