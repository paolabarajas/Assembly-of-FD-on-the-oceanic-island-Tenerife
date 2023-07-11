# Barajas - Barbosa et al (2023) Assembly of functional diversity in an oceanic island flora
# by M. Paola Barajas

# Figure 3 b and c

library(hypervolume)
library(BAT) 
library(agricolae) 
library(colortools)
library(tibble)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(scales)

# Calculation of functional originality and contribution 

setwd("C:/PUB_2/Code_Data_sub3") 

# Tenerife data
Barajasetal <-  read.csv("Barajasetal_traits_status_20230323.csv")
data_i        <- Barajasetal[, c ("IDmaster","Species1","LA","LDMC","LMA","Lth","Nmass","SM","SSD","H",
                                    "Biogeo_status","Ende_status","Lineage")]

# Nearly complete cases (5 % data missing)
data_i$missing_trait <- rowSums(is.na(data_i))
data_i <- filter(data_i, missing_trait < 4)   

data_Tenerife_imp <- read.csv("Trait_imputed_phylo8.csv")
data_Tenerife_imp$Species1 <- NULL

data_i <- left_join (data_i, data_Tenerife_imp, by = "IDmaster")
data_i <- data_i %>%   # replace the few NAs with the imputed data
  mutate(LDMC.x = coalesce(LDMC.x, LDMC.y)) %>%
  mutate(LMA.x = coalesce(LMA.x, LMA.y)) %>%
  mutate(Lth.x = coalesce(Lth.x, Lth.y)) %>%
  mutate(Nmass.x = coalesce(Nmass.x, Nmass.y)) %>%
  mutate(SM.x = coalesce(SM.x, SM.y)) %>%
  mutate(SSD.x = coalesce(SSD.x, SSD.y)) %>%
  mutate(H.x = coalesce(H.x, H.y))

anyNA(data_i)

lin <- filter(data_i, ! Biogeo_status == "NS" & ! Lineage == "") # Lineages diversified at the Macaronesian level
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

# #
pca <- FactoMineR::PCA( data_i_trans[, c("Leaf_area","LMA","Leaf_N","Leaf_th","LDMC","Seed_mass","Stem_density","Height")]
                        , scale.unit = FALSE, graph = FALSE)

sp_coord <- as.data.frame(pca$ind$coord[, 1:3])

sp_coord$Species1 <- rownames(sp_coord)
sp_tra <- dplyr::left_join(data_i_trans, sp_coord, by = "Species1")
rownames(sp_tra) <- sp_tra$Species1

PCA         <- prcomp(sp_tra[, c("Leaf_area","LMA","Leaf_N","Leaf_th","LDMC","Seed_mass","Stem_density","Height")])
PCAloadings <- data.frame(Variables = rownames(PCA$rotation), PCA$rotation)

# Hypervolume
all_bw <- estimate_bandwidth(sp_tra[, c("Dim.1", "Dim.2", "Dim.3")],
                             method = "cross-validation")
set.seed(8)
hv_island <- hypervolume_gaussian(
  sp_tra[, c("Dim.1", "Dim.2", "Dim.3")],
  name = "island_volume",
  kde.bandwidth = all_bw,
  quantile.requested = 0.95,
  quantile.requested.type = "probability", verbose = FALSE)

sp_contrib <- BAT::kernel.contribution(hv_island)
sp_orig    <- BAT::kernel.originality(hv_island)

sp_tra$contrib <- sp_contrib 
sp_tra$orig    <- sp_orig

sp_tra$Lineage[sp_tra$Lineage==""] <- "NS"

save(sp_tra, file = "sp_tra_20230324_semi.Rdata")
# 
#______________________________________________________________________________
#    load("sp_tra_20230322_semi.Rdata")
# ggplot2 theme
tema <- theme(panel.background =element_rect(fill="transparent",colour="black"),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(fill=NA,colour="grey"),
              axis.text.y =element_text(colour=c("black"), size=7),
              axis.text.x =element_text(colour=c("black"), size=7))
#_______________________________________________________________________________
# FIGURE 3 ----
# Functional contribution and originality. Groups analysis. 

# c
# How much does each group contribute to the total island trait space?
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

Figure_3c <- ggpubr::ggarrange(contribution_plot, originality_plot, nrow = 2)
Figure_3c

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

write.csv (kw_grou, file = "results_Figure_3c_kw_con_orig.csv")

# Is the functional contribution and originality different across groups? 
sp_tra_cla <- dplyr::filter(sp_tra, Ende_status == "Cla")
sp_tra_cla$Biogeo_status <- "CLA"  # 264 species
sp_tra_tmp <- rbind(sp_tra, sp_tra_cla) # add species from cladogenetic class

unique(sp_tra_tmp$Biogeo_status)

coco = c("gold2","mediumorchid2","green4", "black","dodgerblue3") 

con <- ggplot(sp_tra_tmp) +
  geom_boxplot(aes(x= factor(Biogeo_status, level = c ("NNE","MAC","CE","TE", "CLA")),
                   y=contrib, color = Biogeo_status), show.legend = FALSE) +
  scale_colour_manual(values = coco) +
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  ylab("Contribution") + xlab("") + tema

ori <- ggplot(sp_tra_tmp) +
  geom_boxplot(aes(x= factor(Biogeo_status, level = c ("NNE","MAC","CE","TE", "CLA")),
                   y=orig, color = Biogeo_status), show.legend = FALSE) +
  scale_colour_manual(values = coco) +
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  ylab("Originality") + xlab("") + tema

Figure_3b <- ggarrange(con, ori, ncol = 1, nrow = 2)
Figure_3b

# Kruskal-Wallis-Test contribution and originality among groups
tmp_biog       <- kruskal(sp_tra_tmp$contrib, trt = as.factor(sp_tra_tmp$Biogeo_status))
a <- tmp_biog$groups %>% rownames_to_column()
tmp_biog       <- kruskal(sp_tra_tmp$orig, trt = as.factor(sp_tra_tmp$Biogeo_status))
b <- tmp_biog$groups %>% rownames_to_column()

ori_cont_groups <- cbind(a, b); ori_cont_groups
ori_cont_groups

# write.csv (ori_cont_groups, file = "results_Figure_3b_con_orig_2.csv", row.names = FALSE)
####

tempi = ggarrange(Figure_3a, Figure_3b, ncol=2 , nrow=1 ,  widths = c(3.5, 1))
tempi
Figure3 = ggarrange(tempi, Figure_3c, ncol=1, nrow=2, heights = c(2.5, 3))
Figure3

pdf(file = "_Figure_3.pdf", 
    bg = "transparent",
    width=15/2.54, height= 12/2.54, #  1 inch = 2.54 cm
    paper = "a4")
Figure3
dev.off()

# FIN