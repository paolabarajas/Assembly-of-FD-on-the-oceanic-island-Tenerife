# M. Paola Barajas Barbosa et al (2023) Assembly of functional diversity in an oceanic island flora
# by M. Paola Barajas

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggExtra)
library(corrplot)
library(tidyr)
library(scales)
library(agricolae)
library(tibble)

# Extended Data Figure 7. Sensitivity analysis using species (n= 237) with complete empirical data for all 8 traits. 
# Plots set for complete cases
# # complete cases ----
# a. trait spaces
# b. functional diversity metrics: richness, eveness, dispersion
# c. functional diversity metrics: contribution, originality

setwd("C:/PUB_2/Code_Data_sub3")

Barajasetal <- read.csv("Barajasetal_traits_status_20230323.csv")

data_Tenerife <- Barajasetal[, 1:16]
data_Tenerife$complete <- complete.cases(data_Tenerife)
data_Tenerife %>% group_by(complete) %>% summarise(n = n())  # 237 species with complete cases
data_Tenerife <- filter (data_Tenerife, complete  == TRUE)
data_Tenerife$source <- 'Tenerife_data'
names(data_Tenerife)

data_Tenerife <- data.frame(
  Species1      = data_Tenerife$Species1,
  Leaf_area     = scale(log10(data_Tenerife$LA)),
  LMA           = scale(log10(data_Tenerife$LMA)),
  Leaf_N        = scale(log10(data_Tenerife$Nmass)),
  Seed_mass     = scale(log10(data_Tenerife$SM)),
  Stem_density  = scale(log10(data_Tenerife$SSD)),
  Height        = scale(log10(data_Tenerife$H)),
  LDMC          = scale(log10(data_Tenerife$LDMC)),
  Leaf_th       = scale(log10(data_Tenerife$Lth)),
  Biogeo_status= data_Tenerife$Biogeo_status,
  Ende_status  = data_Tenerife$Ende_status,
  Lineage      = data_Tenerife$Lineage,
  source       = data_Tenerife$source)

# # Extended Figure 7a ----
PCA         <- prcomp(data_Tenerife[, c("Leaf_area","LMA","Leaf_N","Leaf_th","LDMC","Seed_mass","Stem_density","Height")] )
summary(PCA)
PCAvalues   <- data.frame(Species = data_Tenerife$Species1, Biogeo_status = data_Tenerife$Biogeo_status,
                          Ende_status = data_Tenerife$Ende_status, PCA$x)
PCAloadings <- data.frame(Variables = rownames(PCA$rotation), PCA$rotation)  

PCAloadings$PC1 = PCAloadings$PC1*-1
PCAloadings$PC2 = PCAloadings$PC2*-1

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

TENERIFE <-  PCAvalues %>% ggplot(aes(PC1,  PC2), size = 1)+ 
  stat_density_2d(geom = "polygon", contour = TRUE, aes(fill = after_stat(level)),
                  colour = "gray", bins = 25) +
  
  scale_fill_distiller(palette = "BuGn", direction = 1) +
  geom_jitter(alpha=0.6,  size = 0.8  , colour = "turquoise4") +     
  geom_text(data = PCAloadings, aes(x = PC1*4.7, y = PC2*4.7, label = Variables), size = 2.3) +
  geom_segment(data = PCAloadings, size = 0.2,    
               aes(x = 0, xend = PC1*4.2, y = 0, yend = PC2*4.2),
               arrow = arrow(length = unit(0.1, "cm")),colour = "black")   +
  xlab("n = 237") + ylab("PC2 (26%)")   +
  xlim(-5.4 , 5.4) + ylim(-5, 5) +
  p_my_theme1
TENERIFE

# Plot each group separately
NNE <- dplyr::filter(PCAvalues, Biogeo_status == "NNE")
NNE<- NNE  %>% ggplot(aes(PC1,  PC2))+
  stat_density_2d(geom = "polygon", contour = TRUE, aes(fill = after_stat(level)), colour = "gray", bins = 15) +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  geom_jitter(alpha=0.5,  size = 0.7  , colour = "black") +
  xlim(-5.4 , 5.4) + ylim(-5, 5) +
  xlab("n = 37") + ylab("")   + p_my_theme1
NNE

MAC <- dplyr::filter(PCAvalues, Biogeo_status == "MAC")
MAC <- MAC %>% ggplot(aes(PC1,  PC2))+
  stat_density_2d(geom = "polygon", contour = TRUE, aes(fill = after_stat(level)), colour = "gray", bins = 3.7) +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  geom_jitter(alpha=0.5,  size = 0.8  , colour = "green4") +
  xlim(-5.4 , 5.4) + ylim(-5, 5) +
  xlab("n = 30") + ylab("PC2 (26%)")   +
  p_my_theme1
MAC

CE <- dplyr::filter(PCAvalues, Biogeo_status == "CE")
CE <- CE %>% ggplot(aes(PC1,  PC2))+
  stat_density_2d(geom = "polygon", contour = TRUE, aes(fill = after_stat(level)), colour = "gray", bins = 15) +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  geom_jitter(alpha=0.5,  size = 0.7   , colour = "gold2") +
  xlim(-5.4 , 5.4) + ylim(-5, 5) +
  xlab("n = 115") + ylab("")   +
  p_my_theme1
CE

TE <- dplyr::filter(PCAvalues, Biogeo_status == "TE")
TE <-  TE  %>% ggplot(aes(PC1,  PC2))+
  stat_density_2d(geom = "polygon", contour = TRUE, aes(fill = after_stat(level)), colour = "gray", bins = 12) +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  geom_jitter(alpha=0.5,  size = 0.7   , colour = "dodgerblue3") +
  xlim(-5.4 , 5.4) + ylim(-5, 5) +
  xlab("PC1 (30%) n = 55") + ylab("PC2 (26%)")   +
  p_my_theme1
TE

CLADO <- dplyr::filter(PCAvalues, Ende_status == "Cla")
CLADO <- CLADO %>% ggplot(aes(PC1,  PC2))+
  stat_density_2d(geom = "polygon", contour = TRUE, aes(fill = after_stat(level)), colour = "gray", bins = 13) +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  geom_jitter(alpha=0.5, size = 0.7  , colour = "mediumorchid2") +
  xlim(-5.4 , 5.4) + ylim(-5, 5) +
  xlab("PC1 (30%) n = 138") + ylab("")   +
  p_my_theme1
CLADO

Ten<- ggExtra:: ggMarginal(TENERIFE, type = "density", fill="transparent", size = 15)
ns <- ggExtra:: ggMarginal(NNE, type = "density", fill="transparent", size = 15)
mac<- ggExtra:: ggMarginal(MAC, type = "density", fill="transparent", size = 15)
ae <- ggExtra:: ggMarginal(CE, type = "density", fill="transparent", size = 15)
sie<- ggExtra:: ggMarginal(TE, type = "density", fill="transparent", size = 15)
cla<- ggExtra:: ggMarginal(CLADO, type = "density", fill="transparent", size = 15)

tmp <- ggpubr::ggarrange(Ten, ns, mac, ae,  sie, cla, ncol =2, nrow = 3)
tmp

pdf(file = "Ext_Figure_7a_traitspaces.pdf",
    width =  3.3, 
    height=  4.7,
    paper = "a4r")
tmp
dev.off()

# # Extended Figure 7b ----

library(hypervolume)
library(BAT)

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
all_bw <- hypervolume::estimate_bandwidth(data_Tenerife[, c("Dim.1", "Dim.2", "Dim.3")], # Estimate bandwidth for the island
                             method = "cross-validation")

# Rarefaction/Null models Functional Diversity metrics ----
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

# save(rar_bio_comp, rar_end_comp, file = "results_Ext_Figure_7_20230326.RData") 

# # #
load("results_Ext_Figure_7_20230326.RData")

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

Figure_7 <- ggplot(sum_metrics) +
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

Figure_7

## Functional contribution and originality ----

set.seed(8)
hv_island <- hypervolume_gaussian(
  data_Tenerife[, c("Dim.1", "Dim.2", "Dim.3")],
  name = "island_volume",
  kde.bandwidth = all_bw,
  quantile.requested = 0.95,
  quantile.requested.type = "probability", verbose = FALSE)

sp_contrib <- BAT::kernel.contribution(hv_island)
sp_orig    <- BAT::kernel.originality(hv_island)

data_Tenerife$contrib <- sp_contrib 
data_Tenerife$orig    <- sp_orig

data_Tenerife$Lineage[data_Tenerife$Lineage==""] <- "NS"

save(data_Tenerife, file = "sp_tra_20230327_completecas.Rdata")

tema <- theme(panel.background =element_rect(fill="transparent",colour="black"),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(fill=NA,colour="grey"),
              axis.text.y =element_text(colour=c("black"), size=7),
              axis.text.x =element_text(colour=c("black"), size=7))

# How much does each group contribute to the total island trait space?
load("sp_tra_20230327_completecas.Rdata")

sp_tra <- data_Tenerife 

a<-dplyr::filter(sp_tra, Biogeo_status == "MAC") #  MAC only
b<-dplyr::filter(sp_tra, !Biogeo_status == "MAC") # ALL except MAC
b$Biogeo_status <- "NonMac"

plot_MAC <- ggplot(a, aes(x=Biogeo_status, y=contrib))+
  geom_boxplot(color = "green4")+
  geom_boxplot(data = b, aes(x=Biogeo_status, y=contrib), color = "gray50") +
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  ylab("") + xlab("") + tema
plot_MAC

a<-dplyr::filter(sp_tra, Biogeo_status == "CE")
b<-dplyr::filter(sp_tra, !Biogeo_status == "CE")
b$Biogeo_status <- "NonCE"
plot_CE <- ggplot(a, aes(x=Biogeo_status, y=contrib))+
  geom_boxplot(color = "gold2")+
  geom_boxplot(data = b, aes(x=Biogeo_status, y=contrib), color = "gray50") +
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  ylab("") + xlab("") + tema

a<-dplyr::filter(sp_tra, Biogeo_status == "TE")
b<-dplyr::filter(sp_tra, !Biogeo_status == "TE")
b$Biogeo_status <- "zNonTE"
plot_TE <- ggplot(a, aes(x=Biogeo_status, y=contrib))+
  geom_boxplot(color = "dodgerblue3")+
  geom_boxplot(data = b, aes(x=Biogeo_status, y=contrib), color = "gray50") +
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  ylab("") + xlab("") + tema

a<-dplyr::filter(sp_tra, Biogeo_status == "NNE")
b<-dplyr::filter(sp_tra, !Biogeo_status == "NNE")
b$Biogeo_status <- "NonNNE"
plot_NNE <- ggplot(a, aes(x=Biogeo_status, y=contrib))+
  geom_boxplot(color = "black")+
  geom_boxplot(data = b, aes(x=Biogeo_status, y=contrib), color = "gray50") +
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  ylab("Contribution") + xlab("") + tema

a<-dplyr::filter(sp_tra, Ende_status == "Cla")
b<-dplyr::filter(sp_tra, !Ende_status == "Cla")
b$Ende_status <- "NonCla"
plot_clado <- ggplot(a, aes(x=Ende_status, y=contrib))+
  geom_boxplot(color = "mediumorchid2")+
  geom_boxplot(data = b, aes(x=Ende_status, y=contrib), color = "gray50") +
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  ylab("") + xlab("") + tema
contribution_plot <- ggpubr::ggarrange(plot_NNE, plot_MAC,  plot_CE,  plot_TE, plot_clado, nrow  = 1)
contribution_plot

a<-dplyr::filter(sp_tra, Biogeo_status == "MAC")
b<-dplyr::filter(sp_tra, !Biogeo_status == "MAC") 
b$Biogeo_status <- "NonMac"
plot_MAC <- ggplot(a, aes(x=Biogeo_status, y=orig))+
  geom_boxplot(color = "green4")+
  geom_boxplot(data = b, aes(x=Biogeo_status, y=orig), color = "gray50") +
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  ylab("") + xlab("") + tema
plot_MAC

a<-dplyr::filter(sp_tra, Biogeo_status == "CE")
b<-dplyr::filter(sp_tra, !Biogeo_status == "CE")
b$Biogeo_status <- "NonCE"
plot_CE <- ggplot(a, aes(x=Biogeo_status, y=orig))+
  geom_boxplot(color = "gold2")+
  geom_boxplot(data = b, aes(x=Biogeo_status, y=orig), color = "gray50") +
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  ylab("") + xlab("") + tema

a<-dplyr::filter(sp_tra, Biogeo_status == "TE")
b<-dplyr::filter(sp_tra, !Biogeo_status == "TE")
b$Biogeo_status <- "zNonTE"
plot_TE <- ggplot(a, aes(x=Biogeo_status, y=orig))+
  geom_boxplot(color = "dodgerblue3")+
  geom_boxplot(data = b, aes(x=Biogeo_status, y=orig), color = "gray50") +
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  ylab("") + xlab("") + tema

a<-dplyr::filter(sp_tra, Biogeo_status == "NNE")
b<-dplyr::filter(sp_tra, !Biogeo_status == "NNE")
b$Biogeo_status <- "NonNNE"
plot_NNE <- ggplot(a, aes(x=Biogeo_status, y=orig))+
  geom_boxplot(color = "black")+
  geom_boxplot(data = b, aes(x=Biogeo_status, y=orig), color = "gray50") +
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  ylab("Originality") + xlab("") + tema

a<-dplyr::filter(sp_tra, Ende_status == "Cla")
b<-dplyr::filter(sp_tra, !Ende_status == "Cla")
b$Ende_status <- "NonCla"
plot_clado <- ggplot(a, aes(x=Ende_status, y=orig))+
  geom_boxplot(color = "mediumorchid2")+
  geom_boxplot(data = b, aes(x=Ende_status, y=orig), color = "gray50") +
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  ylab("") + xlab("") + tema

originality_plot <- ggpubr::ggarrange(plot_NNE, plot_MAC, plot_CE,  plot_TE, plot_clado, nrow  = 1)
originality_plot

Figure_7_1 <- ggpubr::ggarrange(contribution_plot, originality_plot, nrow = 2)
Figure_7_1

tempi = ggpubr::ggarrange(Figure_7, Figure_7_1, ncol=1, nrow=2 ,  widths = c(3.5, 1))
tempi

pdf(file = "Ext_Figure_7b.pdf", 
    bg = "transparent",
    width=12/2.54, height= 10/2.54, #  1 inch = 2.54 cm
    paper = "a4")
tempi
dev.off()

# Kruskal-Wallis-Test contribution and originality of groups to the rest of the island.
biogeo_contri_kruski <- c()
for (category in unique(sp_tra$Biogeo_status)) {
  tmp <- sp_tra
  tmp[tmp$Biogeo_status != category , ]$Biogeo_status  <- "Other" # takes everything that is not current category e.g., MAC
  a_kruski <- kruskal(tmp$contrib, trt = as.factor(tmp$Biogeo_status))
  biogeo_contri_kruski[[category]]<-a_kruski
}
biogeo_orig_kruski <- c()
for (category in unique(sp_tra$Biogeo_status)) {
  tmp <- sp_tra
  tmp[tmp$Biogeo_status != category , ]$Biogeo_status  <- "Other"
  a_kruski <- kruskal(tmp$orig, trt = as.factor(tmp$Biogeo_status))
  biogeo_orig_kruski[[category]]<-a_kruski
}
Ende_contri_kruski <- c()
for (category in unique(sp_tra$Ende_status)) {
  tmp <- sp_tra
  tmp[tmp$Ende_status != category,]$Ende_status  <- "Other"
  a_kruski <- kruskal(tmp$contrib, trt = as.factor(tmp$Ende_status))
  Ende_contri_kruski[[category]]<-a_kruski
}
Ende_orig_kruski <- c()
for (category in unique(sp_tra$Ende_status)) {
  tmp <- sp_tra
  tmp[tmp$Ende_status != category,]$Ende_status  <- "Other"
  a_kruski <- kruskal(tmp$orig, trt = as.factor(tmp$Ende_status))
  Ende_orig_kruski[[category]]<-a_kruski
}

contri_NNE <- biogeo_contri_kruski$NNE$groups %>% rownames_to_column() # 
contri_MAC <- biogeo_contri_kruski$MAC$groups %>% rownames_to_column()
contri_CE  <- biogeo_contri_kruski$CE$groups  %>% rownames_to_column()
contri_TE  <- biogeo_contri_kruski$TE$groups  %>% rownames_to_column()
contri_cla   <- Ende_contri_kruski$Cla$groups  %>% rownames_to_column()
kruskal_contri <- rbind(contri_NNE, contri_MAC, contri_CE, contri_TE,
                        contri_cla)

orig_NNE <- biogeo_orig_kruski$NNE$groups %>% rownames_to_column()
orig_MAC <- biogeo_orig_kruski$MAC$groups %>% rownames_to_column()
orig_CE  <- biogeo_orig_kruski$CE$groups  %>% rownames_to_column()
orig_TE  <- biogeo_orig_kruski$TE$groups  %>% rownames_to_column()
orig_cla     <- Ende_orig_kruski$Cla$groups  %>% rownames_to_column()
kruskal_orig <- rbind(orig_NNE, orig_MAC, orig_CE, orig_TE,
                      orig_cla)

kw_grou <- cbind(kruskal_contri, kruskal_orig); kw_grou

# FIN