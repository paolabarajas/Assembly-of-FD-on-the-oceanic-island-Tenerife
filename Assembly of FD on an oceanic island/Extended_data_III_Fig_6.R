# M. Paola Barajas Barbosa et al (2023) Assembly of functional diversity in an oceanic island flora
# by M. Paola Barajas

# Extended Data Figure 6

library(hypervolume)
library(BAT) 
library(agricolae) 
library(colortools)
library(colortools)
library(tibble)
library(ggpubr)
library(ggplot2)
library(dplyr)

setwd("C:/PUB_2/Code_Data_sub3")

# Lineages functional originality and contribution 
Barajasetal <-  read.csv("Barajasetal_traits_status_20230323.csv")
status        <- Barajasetal[, c("IDmaster", "Biogeo_status","Ende_status","Lineage")]
data_i        <- Barajasetal[, c ("IDmaster","Species1","LA","LDMC","LMA","Lth","Nmass","SM","SSD","H",
                                    "Biogeo_status","Ende_status","Lineage")]

data_i$missing_trait <- rowSums(is.na(data_i))
data_i <- filter(data_i, missing_trait < 4)   

data_Tenerife_imp <- read.csv("Trait_imputed_phylo8.csv")
data_Tenerife_imp$Species1 <- NULL
data_i <- left_join (data_i, data_Tenerife_imp, by = "IDmaster")
data_i <- data_i %>%   # replace the very few NAs with the imputed data
  mutate(LDMC.x = coalesce(LDMC.x, LDMC.y)) %>%
  mutate(LMA.x = coalesce(LMA.x, LMA.y)) %>%
  mutate(Lth.x = coalesce(Lth.x, Lth.y)) %>%
  mutate(Nmass.x = coalesce(Nmass.x, Nmass.y)) %>%
  mutate(SM.x = coalesce(SM.x, SM.y)) %>%
  mutate(SSD.x = coalesce(SSD.x, SSD.y)) %>%
  mutate(H.x = coalesce(H.x, H.y))

anyNA(data_i)

lin <- filter(data_i, ! Biogeo_status == "NS" & ! Lineage == "") 
unique(lin$Lineage)

data_i_trans <- data.frame(
  IDmaster      = data_i$IDmaster,
  Species1      = data_i$Species1,
  Leaf_area     = scale(log10(data_i$LA.x)),
  LDMC          = scale(log10(data_i$LDMC.x)),
  LMA           = scale(log10(data_i$LMA.x)),
  Leaf_th       = scale(log10(data_i$Lth.x)),
  Leaf_N        = scale(log10(data_i$Nmass.x)),
  Seed_mass     = scale(log10(data_i$SM.x)),
  Stem_density  = scale(log10(data_i$SSD.x)),
  Height        = scale(log10(data_i$H.x)),
  Biogeo_status = data_i$Biogeo_status,
  Ende_status   = data_i$Ende_status,
  Lineage       = data_i$Lineage )

rownames(data_i_trans) <- data_i_trans$Species1

pca <- FactoMineR::PCA( data_i_trans[, c("Leaf_area","LMA","Leaf_N","Leaf_th","LDMC","Seed_mass","Stem_density","Height")]
                        , scale.unit = FALSE, graph = FALSE)
sp_coord <- as.data.frame(pca$ind$coord[, 1:3])
sp_coord$Species1 <- rownames(sp_coord)
sp_tra <- dplyr::left_join(data_i_trans, sp_coord, by = "Species1")
rownames(sp_tra) <- sp_tra$Species1

PCA         <- prcomp(sp_tra[, c("Leaf_area","LMA","Leaf_N","Leaf_th","LDMC","Seed_mass","Stem_density","Height")])
PCAloadings <- data.frame(Variables = rownames(PCA$rotation), PCA$rotation)

load("sp_tra_20230324_semi.RData")

sp_tra$Lineage = gsub("Spartocytisus","Cytisus", sp_tra$Lineage)

#_______________________________________________________________________________
# Extended FIGURE 5 ----
# Functional contribution and originality for radiated plant lineages present in Tenerife 

tema <- theme(panel.background =element_rect(fill="transparent",colour="black"),
               panel.grid.minor=element_blank(),
               panel.border=element_rect(fill=NA,colour="grey"),
               axis.text.y =element_text(colour=c("black"), size=7), 
               axis.text.x =element_text(colour=c("black"), size=7))

sp_tra$Lineage <- sub("NS", "zzzNS", sp_tra$Lineage)
num_lineage <- as.data.frame(table(factor(sp_tra$Lineage))) 

# plotting each lineage contribution to island trait space
plot_contrib <- list()
for (i in unique(sort(sp_tra$Lineage))) {
  
  a<-dplyr::filter(sp_tra,   Lineage == i)  
  b<-dplyr::filter(sp_tra, ! Lineage == i) 
  
  b$Lineage <- "zzNon" 
  
  plot_contrib[[i]]  = ggplot( a , aes( x = Lineage, y = contrib))+
    geom_boxplot(color = "deeppink3",  lwd = 0.2  ,  outlier.size = 0.2) +
    geom_boxplot(data  = b, aes( x = Lineage, y = contrib), color = "gray50", lwd = 0.2  ,  outlier.size = 0.2) +
    ylab("") + xlab("") + tema
}

plot_orig <- list()
for (i in unique(sort(sp_tra$Lineage))) {
  
  a<-dplyr::filter(sp_tra,   Lineage == i)  
  b<-dplyr::filter(sp_tra, ! Lineage == i) # ALL except a Lineage
  
  b$Lineage <- "zzNon" 
  
  plot_orig[[i]]  <- ggplot( a , aes( x = Lineage, y = orig))+
    geom_boxplot(color = "deeppink3" , lwd = 0.2 ,  outlier.size = 0.2)+
    geom_boxplot(data  = b, aes( x = Lineage, y = orig), color = "gray50", lwd = 0.2 , outlier.size = 0.2) +
    ylab("") + xlab("") + tema
}

plot_contrib$zzzNS <- NULL
plot_orig$zzzNS <- NULL

li_contri  <- ggarrange(plotlist = plot_contrib, ncol = 7, nrow = 10)
li_orig    <- ggarrange(plotlist = plot_orig, ncol = 7, nrow = 10)

pdf(file = "Ext_Figure_5_lin_contri.pdf", bg = "transparent",
    width=6, 
    height=10,
    paper = "a4")
li_contri
dev.off()

pdf(file = "Ext_Figure_5_lin_orig.pdf", bg = "transparent",
    width=6, 
    height=10,
    paper = "a4")
li_orig
dev.off()

# Kruskal Wallis Test contribution and originality: 
# how much does each Lineage contributes to the total island trait space

Lin_contri_kruski <- c()
for (category in unique(sp_tra$Lineage)) {
  tmp <- sp_tra
  tmp[tmp$Lineage != category , ]$Lineage  <- "zOther" # takes everything that is not current category (i.e., Lineage)
  a_kruski <- kruskal(tmp$contrib, trt = as.factor(tmp$Lineage))
  Lin_contri_kruski[[category]] <- a_kruski
}

Lin_orig_kruski <- c()
for (category in unique(sp_tra$Lineage)) {
  tmp <- sp_tra
  tmp[tmp$Lineage != category , ]$Lineage  <- "zOther" 
  a_kruski <- kruskal(tmp$orig, trt = as.factor(tmp$Lineage))
  Lin_orig_kruski[[category]] <- a_kruski
}

df<- data.frame()
for (i in 2:length(Lin_contri_kruski)) { 
  z   <- as.data.frame(Lin_contri_kruski[[i]][5])
  df  <- rbind(df, z)
}

df1<- data.frame()
for (i in 2:length(Lin_orig_kruski)) { 
  z     <- as.data.frame(Lin_orig_kruski[[i]][5])
  df1   <- rbind(df1, z)
}

df$name <- rownames(df)
df1$name <- rownames(df1)

kw_lin <- cbind(df, df1)
View(kw_lin)

write.csv(kw_lin, "Extended_Fig5_kruskal.csv")

# Plotting lineages that increase island trait space: 

linlin1<- dplyr::filter(sp_tra, Lineage == "Aeonium") 
linlin2<- dplyr::filter(sp_tra, Lineage == "Polycarpaea") 

linlin1['new_col1'] <- "Aeonium_25"
linlin2['new_col1'] <- "Polycarpaea_5"

linlin_z <- dplyr::filter(sp_tra, ! Lineage == "Aeonium" & ! Lineage == "Polycarpaea") 
linlin_z['new_col1'] <- "zTenerife"

tmp <- rbind(linlin_z, linlin1, linlin2)
unique(tmp$new_col1)

# Plot
pdf(file = "Ext_Figure_5c_traitspace.pdf", bg = "transparent",
    width= 3.8,
    height= 2,
    paper = "a4")

ggplot(tmp) +   
  geom_point(size = 1, alpha=0.9, 
             aes(x = Dim.1, y = Dim.2, group = new_col1, colour = new_col1, size = new_col1),
             show.legend = TRUE) +
   scale_colour_manual(values=c( "mediumvioletred", "gold1", "turquoise4")) +
  geom_text(data = PCAloadings, aes(x = (PC1*4.7)*-1, y = (PC2*4.7)*-1, label = Variables), size = 2.3) +
  geom_segment(data = PCAloadings, size = 0.2,    # Plots the loading, i.e., traits 
               aes(x = 0, xend = (PC1*4.2)*-1, y = 0, yend = (PC2*4.2)*-1),
               arrow = arrow(length = unit(0.1, "cm")),colour = "black")   +
  xlab("PC1 (30%)") + ylab("PC2 (25%)")   +
  tema

dev.off()

# Proportion of lineages belonging to the endemic species groups 
# sp_tra object required 

num_lineage <- as.data.frame(table(factor(sp_tra$Lineage))) # counts number of spp in each lineage. There are 66 lineages. NS are species that do not belong to a lineage

LinTE <-dplyr::filter(sp_tra, Biogeo_status == "TE")
LinCE <-dplyr::filter(sp_tra, Biogeo_status == "CE")   
LinMAC<-dplyr::filter(sp_tra, Biogeo_status == "MAC")

l_TE <-unique(LinTE$Lineage)  # 31
l_CE <-unique(LinCE$Lineage)  # 54
l_MAC<-unique(LinMAC$Lineage) # 14

species_lin = dplyr::filter(sp_tra, Lineage != "zzzNS") # 195

data <- data.frame(
  category=c("Tenerife", "Canary Islands", "Macaronesia"),
  count=c(31, 54, 14))

pdf(file = "Ext_Figure_5_proport_lin.pdf", bg = "transparent",
    width= 2.5,
    height= 2,
    paper = "a4")
ggplot(data, aes(x=category, y=count, color=category)) +
  geom_bar(stat="identity", fill="white", show.legend = FALSE) +
  scale_colour_manual(values=c("gold2", "green4", "dodgerblue3" )) +
  xlab("") + ylab("Number of lineages")+
  tema
dev.off()

# FIN