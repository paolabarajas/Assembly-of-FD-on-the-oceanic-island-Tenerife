# M. Paola Barajas Barbosa et al (2023) Assembly of functional diversity in an oceanic island flora
# by M. Paola Barajas

# # Figure 2

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggExtra)

setwd("C:/PUB_2/Code_Data_sub3")

# Tenerife data
Barajasetal <-  read.csv("Barajasetal_traits_status_20230323.csv")

p_my_theme1 <-  theme( axis.title.x=element_text(colour="black",face="bold",size=7),
                       axis.title.y=element_text(colour="black",face="bold",size=7),
                       axis.text.x=element_text(colour=c("black"),face="bold",size=7),
                       axis.text.y=element_text(colour=c("black"),face="bold",size=7),
                       legend.position  = c(0.11,0.15), legend.direction="vertical",    # position
                       legend.key       = element_rect(fill="transparent"),
                       legend.key.size  = unit(.3,"line"),
                       legend.title     = element_blank(), 
                       legend.text      = element_text(size= 7, color="black"),
                       panel.background = element_rect(fill="transparent",colour="black"),
                       panel.grid.minor = element_blank(),
                       panel.border     = element_rect(fill=NA,colour="grey"))

p_my_theme2 <-  theme( axis.text.x=element_text(colour=c("white"),face="bold",size=7),
                       axis.text.y=element_text(colour=c("white"),face="bold",size=7),
                       legend.position  = c(0.11,0.15), legend.direction="vertical",    
                       legend.key       = element_rect(fill="transparent"),
                       legend.key.size  = unit(.3,"line"),
                       legend.title     = element_blank(), 
                       legend.text      = element_text(size= 7, color="black"),
                       panel.background = element_rect(fill="transparent",colour="black"),
                       panel.grid.minor = element_blank(),
                       panel.border     = element_rect(fill=NA,colour="grey"))

# Nearly complete cases: Species with min 5 traits ----
data_Tenerife <- Barajasetal[, 1:16]
data_Tenerife$missing_trait <- rowSums(is.na(data_Tenerife))
data_Tenerife <- filter(data_Tenerife, missing_trait < 4)   # 348 species

data_Tenerife_imp <- read.csv("Trait_imputed_phylo8.csv" ) 
data_Tenerife_imp$Species1 <- NULL
data_Tenerife <- left_join (data_Tenerife, data_Tenerife_imp, by = "IDmaster")

data_Tenerife <- data_Tenerife %>%   # replace NAs with the imputed data
  mutate(LDMC.x = coalesce(LDMC.x, LDMC.y)) %>%
  mutate(LMA.x = coalesce(LMA.x, LMA.y)) %>%
  mutate(Lth.x = coalesce(Lth.x, Lth.y)) %>%
  mutate(Nmass.x = coalesce(Nmass.x, Nmass.y)) %>%
  mutate(SM.x = coalesce(SM.x, SM.y)) %>%
  mutate(SSD.x = coalesce(SSD.x, SSD.y)) %>%
  mutate(H.x = coalesce(H.x, H.y))

data_Tenerife$source <- 'Tenerife_data'

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

data_Tenerife %>% group_by(Biogeo_status) %>% summarise(n = n())
data_Tenerife %>% group_by(Ende_status) %>% summarise(n = n())

data_Tenerife %>% filter(Lineage != "") %>% count(Lineage) # 63 lineages

# # PCA 
PCA             <- prcomp(data_Tenerife[, c("Leaf_area","LMA","Leaf_N","Leaf_th","LDMC","Seed_mass","Stem_density","Height")] )
summary(PCA)

PCAvalues       <- data.frame(Species = data_Tenerife$Species1, Biogeo_status = data_Tenerife$Biogeo_status, 
                              Ende_status = data_Tenerife$Ende_status, PCA$x)
PCAvalues$PC1   <- PCAvalues$PC1*-1 ; PCAvalues$PC2 <- PCAvalues$PC2*-1 # for visualization purposes

PCAloadings     <- data.frame(Variables = rownames(PCA$rotation), PCA$rotation)  
PCAloadings$PC1 <- PCAloadings$PC1*-1 ; PCAloadings$PC2 <- PCAloadings$PC2*-1

# # 
# 1. Plots Figure 1  
TENERIFE <-
  PCAvalues %>% ggplot(aes(PC1,  PC2), size = 1)+ # Plot Tenerife PCA 
  stat_density_2d(geom = "polygon", contour = TRUE, aes(fill = after_stat(level)),  
                  colour = "gray", bins = 34) +
  
  scale_fill_distiller(palette = "BuGn", direction = 1) +
  geom_jitter(alpha=0.6,  size = .7  , colour = "turquoise4") +     #  Display the points 
  geom_text(data = PCAloadings, aes(x = PC1*4.7, y = PC2*4.7, label = Variables), size = 2.3) +
  geom_segment(data = PCAloadings, size = 0.2,    # Plots the loadings, i.e., traits 
               aes(x = 0, xend = PC1*4.2, y = 0, yend = PC2*4.2),
               arrow = arrow(length = unit(0.1, "cm")),colour = "black")   + 
  xlab(" n = 348") + ylab("PC2 (25%)")   +
  xlim(-5 , 5) + ylim(-5, 4) +
  p_my_theme1
TENERIFE

# Plot each group separately
NNE <- dplyr::filter(PCAvalues, Biogeo_status == "NNE")
NNE<- NNE  %>% ggplot(aes(PC1,  PC2))+ 
  stat_density_2d(geom = "polygon", contour = TRUE, aes(fill = after_stat(level)), colour = "gray", bins = 10) +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  geom_jitter(alpha=0.5,  size = 0.7  , colour = "black") +    
  xlim(-5 , 5) + ylim(-5, 4) +
  xlab("n = 54") + ylab("")   + p_my_theme1
NNE

MAC <- dplyr::filter(PCAvalues, Biogeo_status == "MAC")
MAC <- MAC %>% ggplot(aes(PC1,  PC2))+ 
  stat_density_2d(geom = "polygon", contour = TRUE, aes(fill = after_stat(level)), colour = "gray", bins = 5) +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  geom_jitter(alpha=0.5,  size = 0.7  , colour = "green4") +  
  xlim(-5 , 5) + ylim(-5, 4) +
  xlab("n = 41") + ylab("PC2 (25%)")   +
  p_my_theme1
MAC

CE <- dplyr::filter(PCAvalues, Biogeo_status == "CE")
CE <- CE %>% ggplot(aes(PC1,  PC2))+ 
  stat_density_2d(geom = "polygon", contour = TRUE, aes(fill = after_stat(level)), colour = "gray", bins = 15) +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  geom_jitter(alpha=0.5,  size = 0.7  , colour = "gold2") +                      
  xlim(-5 , 5) + ylim(-5, 4) +
  xlab("n = 168") + ylab("")   +
  p_my_theme1
CE

TE <- dplyr::filter(PCAvalues, Biogeo_status == "TE")
TE <-  TE  %>% ggplot(aes(PC1,  PC2))+  
  stat_density_2d(geom = "polygon", contour = TRUE, aes(fill = after_stat(level)), colour = "gray", bins = 15) +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  geom_jitter(alpha=0.5,  size = 0.7   , colour = "dodgerblue3") +     
  xlim(-5 , 5) + ylim(-5, 4) +
  xlab("PC1 (30%) n = 85") + ylab("PC2 (25%)")   +
  p_my_theme1
TE

CLADO <- dplyr::filter(PCAvalues, Ende_status == "Cla")
CLADO <- CLADO %>% ggplot(aes(PC1,  PC2))+ 
  stat_density_2d(geom = "polygon", contour = TRUE, aes(fill = after_stat(level)), colour = "gray", bins = 20) +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  geom_jitter(alpha=0.5, size = 0.7  , colour = "mediumorchid2") +  
  xlim(-5 , 5) + ylim(-5, 4) +
  xlab("PC1 (30%) n = 195") + ylab("")   +
  p_my_theme1
CLADO

Ten<- ggExtra:: ggMarginal(TENERIFE, type = "density", fill="transparent", size = 15)
ns <- ggExtra:: ggMarginal(NNE, type = "density", fill="transparent", size = 15) 
mac<- ggExtra:: ggMarginal(MAC, type = "density", fill="transparent", size = 15)
ae <- ggExtra:: ggMarginal(CE, type = "density", fill="transparent", size = 15)
sie<- ggExtra:: ggMarginal(TE, type = "density", fill="transparent", size = 15)
cla<- ggExtra:: ggMarginal(CLADO, type = "density", fill="transparent", size = 15)

Figure_2 <- ggpubr::ggarrange(Ten, ns, mac, ae,  sie, cla, ncol =2, nrow = 3) 
Figure_2

pdf(file = "_Figure_2_traitspaces.pdf",
width =  3.8, # 89 mm (single column) is 3.5 inches.
height=  5,
paper = "a4r")
Figure_2
dev.off()

# FIN