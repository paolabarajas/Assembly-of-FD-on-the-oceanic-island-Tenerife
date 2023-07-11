# M. Paola Barajas Barbosa et al (2023) Assembly of functional diversity in an oceanic island flora
# by M. Paola Barajas

# Extended Data Figures 1, 2 and 3.

library(dplyr)
library(GGally)
library(ggplot2)
library(reshape2)
library(tibble)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(agricolae)
library(visdat)
library(ggpubr)

setwd("C:/PUB_2/Code_Data_sub3")

traits_status <- read.csv("Barajasetal_traits_status_20230323.csv")
status        <- traits_status[, c("IDmaster", "Biogeo_status","Ende_status","Lineage")]
data_i        <- traits_status[, c ("IDmaster","Species1","LA","LDMC","LMA","Lth","Nmass","SM","SSD","H",
                                    "Biogeo_status","Ende_status","Lineage")]

my_theme <- theme(panel.background = element_blank(), legend.position = "none",
                  axis.title.y =element_blank(),
                  axis.text.x  =element_text(colour="black",size=7),
                  axis.text.y  =element_text(colour="black",size=7),
                  panel.border =element_rect(fill=NA,colour="grey"))

p_my_theme1 <-  theme( axis.title.x=element_text(colour="black",face="bold",size=7),
                       axis.title.y=element_text(colour="black",face="bold",size=7),
                       axis.text.x=element_text(colour=c("black"),face="bold",size=7),
                       axis.text.y=element_text(colour=c("black"),face="bold",size=7),
                       legend.position  = c(0.11,0.15), legend.direction="vertical",   
                       legend.key       = element_rect(fill="transparent"),
                       legend.key.size  = unit(.3,"line"),
                       legend.title     = element_blank(), 
                       legend.text      = element_text(size= 7, color="black"),
                       panel.background = element_rect(fill="transparent",colour="black"),
                       panel.grid.minor = element_blank(),
                       panel.border     = element_rect(fill=NA,colour="grey"))

theme1 <- theme( axis.title = element_text(size=7),
                 axis.text.x= element_text(size=7), axis.text.y = element_text(size=7),
                 legend.text= element_text(size= 7)) 

theme2 <- theme( axis.title= element_text(size=7),
                 panel.background =element_rect(fill="transparent",colour="black"),panel.grid.minor=element_blank(),
                 axis.text.x= element_text(colour= "black", size=7), axis.text.y= element_text(colour= "black", size=7),
                 legend.position = "none") 

# Near-complete cases'. ---
data_i$missing_trait <- rowSums(is.na(data_i))
data_i <- filter(data_i, missing_trait < 4)   # 348 species

data_Tenerife_imp <- read.csv("Trait_imputed_phylo8.csv")
data_Tenerife_imp$Species1 <- NULL
data_i <- left_join (data_i, data_Tenerife_imp, by = "IDmaster")
data_i <- data_i %>%   # replace NAs with the imputed data
  mutate(LDMC.x = coalesce(LDMC.x, LDMC.y)) %>%
  mutate(LMA.x = coalesce(LMA.x, LMA.y)) %>%
  mutate(Lth.x = coalesce(Lth.x, Lth.y)) %>%
  mutate(Nmass.x = coalesce(Nmass.x, Nmass.y)) %>%
  mutate(SM.x = coalesce(SM.x, SM.y)) %>%
  mutate(SSD.x = coalesce(SSD.x, SSD.y)) %>%
  mutate(H.x = coalesce(H.x, H.y))

anyNA(data_i)

lin <- filter(data_i, ! Biogeo_status == "NS" & ! Lineage == "") # Lineages diversified at the Canary islands & Macaronesian level
unique(lin$Lineage)

data_i_trans <- data.frame(
  IDmaster      = data_i$IDmaster,
  Species1      = data_i$Species1,
  Leaf_area     = scale(log10(data_i$LA.x)),
  LDMC          = scale(log10(data_i$LDMC.x)),
  LMA           = scale(log10(data_i$LMA.x)),
  Lth           = scale(log10(data_i$Lth.x)),
  Leaf_N        = scale(log10(data_i$Nmass.x)),
  Seed_mass     = scale(log10(data_i$SM.x)),
  Stem_density  = scale(log10(data_i$SSD.x)),
  Height        = scale(log10(data_i$H.x)),
  Biogeo_status = data_i$Biogeo_status,
  Ende_status   = data_i$Ende_status,
  Lineage       = data_i$Lineage )
#

# EXTENDED DATA FIGURE 1 ----
data_i_0 = data_i
#  a. 
pdf(file = "Ext_Figure_1a_1.pdf",   
    width =2.2,
    height=4, paper = "a4")
vis_miss(data_i_0[,3:10]) + theme1
dev.off()

plot <- list()
for (i in unique(sort(data_i_0$Biogeo_status))) {
  tmp = dplyr::filter(data_i_0, Biogeo_status == i )
  plot[[i]]  = vis_miss(tmp[,c("LA","LDMC","LMA","Lth","Nmass","SM","SSD","H")]) + labs(title = i) + theme1 
} 

plot1 <- list()
for (i in unique(data_i_0$Ende_status)) {
  tmp = dplyr::filter(data_i_0, Ende_status == i )
  plot1[[i]]  = vis_miss(tmp[,c("LA","LDMC","LMA","Lth","Nmass","SM","SSD","H")])  + labs(title = i) + theme1
} 

plot$Cla = plot1$Cla

pdf(file = "Ext_Figure_1a_2.pdf",
    width =5.5,
    height=5, paper = "a4")
ggarrange(plotlist = plot)
dev.off()

#  b. 
pdf(file = "Ext_Figure_1b_1.pdf",   
    width =2.2,
    height=4, paper = "a4")
vis_miss(data_i[,3:10]) + theme1
dev.off()

plot <- list()
for (i in unique(sort(data_i$Biogeo_status))) {
  tmp = dplyr::filter(data_i, Biogeo_status == i )
  plot[[i]]  = vis_miss(tmp[,c("LA","LDMC","LMA","Lth","Nmass","SM","SSD","H")]) + labs(title = i) + theme1 
} 

plot1 <- list()
for (i in unique(data_i$Ende_status)) {
  tmp = dplyr::filter(data_i, Ende_status == i )
  plot1[[i]]  = vis_miss(tmp[,c("LA","LDMC","LMA","Lth","Nmass","SM","SSD","H")])  + labs(title = i) + theme1
} 

plot$Cla = plot1$Cla

pdf(file = "Ext_Figure_1b_2.pdf",
    width =5.5,
    height=5, paper = "a4")
ggarrange(plotlist = plot)
dev.off()

#  c.
trait<- data_i # use object from line 22
LA<- ggplot(trait, aes(LA)) + geom_density()+ scale_x_log10() + p_my_theme1
LDMC<- ggplot(trait, aes(LDMC)) + geom_density()+ scale_x_log10() + p_my_theme1
LMA<- ggplot(trait, aes(LMA)) + geom_density()+ scale_x_log10() + p_my_theme1
Lth<- ggplot(trait, aes(Lth)) + geom_density()+ scale_x_log10() + p_my_theme1
Nmass<- ggplot(trait, aes(Nmass)) + geom_density()+ scale_x_log10() + p_my_theme1
SSD<- ggplot(trait, aes(SSD)) + geom_density()+ scale_x_log10() + p_my_theme1
H<- ggplot(trait, aes(H)) + geom_density()+ scale_x_log10() + p_my_theme1
SM<- ggplot(trait, aes(SM)) + geom_density()+ scale_x_log10() + p_my_theme1
trait_original<- ggpubr::ggarrange(LA, LDMC, LMA, Lth, Nmass, SSD, SM, H, 
                                   ncol = 8, nrow = 1,
                                   labels = c("", ""), hjust = -1.1)

trait_imputed_phylo<- data_Tenerife_imp
LA<- ggplot(trait_imputed_phylo, aes(LA)) + geom_density()+ scale_x_log10() + p_my_theme1
LDMC<- ggplot(trait_imputed_phylo, aes(LDMC)) + geom_density()+ scale_x_log10() + p_my_theme1
LMA<- ggplot(trait_imputed_phylo, aes(LMA)) + geom_density()+ scale_x_log10() + p_my_theme1
Lth<- ggplot(trait_imputed_phylo, aes(Lth)) + geom_density()+ scale_x_log10() + p_my_theme1
Nmass<- ggplot(trait_imputed_phylo, aes(Nmass)) + geom_density()+ scale_x_log10() + p_my_theme1
SSD<- ggplot(trait_imputed_phylo, aes(SSD)) + geom_density()+ scale_x_log10() + p_my_theme1
H<- ggplot(trait_imputed_phylo, aes(H)) + geom_density()+ scale_x_log10() + p_my_theme1
SM<- ggplot(trait_imputed_phylo, aes(SM)) + geom_density()+ scale_x_log10() + p_my_theme1
trait_imputed<- ggpubr::ggarrange(LA, LDMC, LMA, Lth, Nmass, SSD, SM, H, 
                                  ncol = 8, nrow = 1,
                                  labels = c("", ""), hjust = 2)

pdf(file = "Ext_Figure_1c.pdf",   
    width =7.5,
    height=2, paper = "a4")
ggpubr::ggarrange (trait_original, trait_imputed, ncol=1, nrow=2,
                   labels = c("trait original", "trait imputed"), hjust = -1)
dev.off()


# EXTENDED DATA FIGURE 2 ----
# Get 'comparison' object from sub_2_Figure_1_isl_glob_.R file

load("Comparison_20221202.RData")

#  a.
prcomp(comparison[,c("Leaf_area","LMA","Leaf_N","Seed_mass","Stem_density","Height" )])
res.pca <- PCA(comparison[,c("Leaf_area","LMA","Leaf_N","Seed_mass","Stem_density","Height" )], graph = FALSE)
eig.val <- get_eigenvalue(res.pca)
eig.val
a<- fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50),  barcolor = "gray50", barfill = "gray") + my_theme
a

pdf(file = "Ext_Figure_2a.pdf",
    width =3,     
    height=3, paper = "a4")
a
dev.off()

#  b.
var <- get_pca_var(res.pca)
corrplot(var$contrib, is.corr=FALSE, col= colorRampPalette(c("grey90","grey50", "grey10"))(30))

pdf(file = "Ext_Figure_2b.pdf",
    width =2.5,
    height=2.5, paper = "a4")
corrplot(var$contrib, is.corr=FALSE, col= colorRampPalette(c("grey90","grey50", "grey10"))(30))
dev.off()

#  c.
detach("package:plyr", unload = TRUE)

load("results_Figure_1_nullmodel_basic_2022117.RData") 
load("results_Figure_1_nullmodel_2022117.RData") 
coco = c("gold2","gray65","turquoise4") 

# Null mode controlling for growth form
null_mod_2 = as.data.frame(null_mod_2)
colnames(null_mod_2) <- c("jaccard","Similar_component_sorensen",
                          "Unique_component_Island","Unique_component_globe")

null_2_sorensen <- null_mod_2 %>%
  summarise( n=n(), mean = mean(Similar_component_sorensen), sd= sd(Similar_component_sorensen) ) %>%
  mutate( se = sd /sqrt(n))  %>%
  mutate( ic = se * qt((1-0.05)/2 + .5, n-1))  %>%
  mutate( mean, metric = "Similar_component_sorensen")
null_2_sorensen$ID <- "null model 2"
null_2_frac_unique_1 <- null_mod_2 %>%  # Unique_1 is Unique component of Island volume relative to global one
  summarise( n=n(), mean = mean(Unique_component_Island), sd= sd(Unique_component_Island) ) %>%
  mutate( se = sd /sqrt(n))  %>%
  mutate( ic = se * qt((1-0.05)/2 + .5, n-1))  %>%
  mutate( mean, metric = "Unique_component_Island")
null_2_frac_unique_1$ID <- "null_2_frac_unique_island"
null_2_frac_unique_2 <- null_mod_2 %>%  # Unique_2 is Unique component of global volume relative to the island one
  summarise( n=n(), mean = mean(Unique_component_globe), sd= sd(Unique_component_globe) ) %>%
  mutate( se = sd /sqrt(n))  %>%
  mutate( ic = se * qt((1-0.05)/2 + .5, n-1))  %>%
  mutate( mean, metric = "Unique_component_globe")
null_2_frac_unique_2$ID <- "null_2_frac_unique_globe"
null_models <- rbind( null_2_sorensen, null_2_frac_unique_1, null_2_frac_unique_2)

# Null model basic
null_mod_1 = as.data.frame(null_mod_1)
colnames(null_mod_1) <- c("jaccard","Similar_component_sorensen",
                          "Unique_component_Island","Unique_component_globe")

null_2_sorensen <- null_mod_1 %>%
  summarise( n=n(), mean = mean(Similar_component_sorensen), sd= sd(Similar_component_sorensen) ) %>%
  mutate( se = sd /sqrt(n))  %>%
  mutate( ic = se * qt((1-0.05)/2 + .5, n-1))  %>%
  mutate( mean, metric = "Similar_component_sorensen")
null_2_sorensen$ID <- "null model 2"
null_2_frac_unique_1 <- null_mod_1 %>%  
  summarise( n=n(), mean = mean(Unique_component_Island), sd= sd(Unique_component_Island) ) %>%
  mutate( se = sd /sqrt(n))  %>%
  mutate( ic = se * qt((1-0.05)/2 + .5, n-1))  %>%
  mutate( mean, metric = "Unique_component_Island")
null_2_frac_unique_1$ID <- "null_2_frac_unique_island"
null_2_frac_unique_2 <- null_mod_1 %>%  
  summarise( n=n(), mean = mean(Unique_component_globe), sd= sd(Unique_component_globe) ) %>%
  mutate( se = sd /sqrt(n))  %>%
  mutate( ic = se * qt((1-0.05)/2 + .5, n-1))  %>%
  mutate( mean, metric = "Unique_component_globe")
null_2_frac_unique_2$ID <- "null_2_frac_unique_globe"
null_models_basic <- rbind( null_2_sorensen, null_2_frac_unique_1, null_2_frac_unique_2)

overlap_stats <- ggplot(null_models, label= mean)+
  geom_linerange(aes(x = metric, y = mean*100, ymin= mean*100-ic*100, ymax= mean*100+ic*100, color = ID), show.legend = FALSE) +
  geom_point(aes(x = metric, y = mean*100, color = ID),  size = 2,show.legend = FALSE) +
  scale_y_continuous(breaks = seq(0, 100, by = 10))+
  scale_colour_manual(values = coco) +  
  labs(x = "", y = "Hypervolume overlap statistics (%) complete cases only") + my_theme

overlap_stats_basic <- ggplot(null_models_basic, label= mean)+
  geom_linerange(aes(x = metric, y = mean*100, ymin= mean*100-ic*100, ymax= mean*100+ic*100, color = ID), show.legend = FALSE) +
  geom_point(aes(x = metric, y = mean*100, color = ID),  size = 2,show.legend = FALSE) +
  scale_y_continuous(breaks = seq(0, 100, by = 10))+
  scale_colour_manual(values = coco) +  
  labs(x = "", y = "Hypervolume overlap statistics (%) complete cases only") + my_theme

ggpubr::ggarrange( overlap_stats, overlap_stats_basic)

pdf(file = "Ext_Figure_2d.pdf",
    bg = "transparent",
    width=2.5, 
    height=2.5,
    paper = "a4")
ggpubr::ggarrange( overlap_stats, overlap_stats_basic)
dev.off()

# EXTENDED DATA FIGURE 3 ----
#  a. Explained variances
res.pca <- PCA(data_i_trans[, 3:10], graph = FALSE)
eig.val <- get_eigenvalue(res.pca)
a <- fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50),  barcolor = "gray50", barfill = "gray") + my_theme
a

#  b. Contribution each trait
var <- get_pca_var(res.pca)
b   <- corrplot(var$contrib, is.corr=FALSE, col= colorRampPalette(c("grey90","grey50", "grey10"))(30)) 
b

#  c. Box plots + Kruskal wallis test 
PCA         <- prcomp(data_i_trans[, 3:10])
PCAvalues   <- data.frame(Species = data_i$Species1, 
                          Biogeo_status = data_i$Biogeo_status, 
                          Ende_status = data_i$Ende_status, PCA$x)

tmp_b <- PCAvalues
tmp_b[, c("Ende_status")] = NULL

tmp_c <- PCAvalues
tmp_c[,c("Biogeo_status")] <- NULL
names(tmp_c)[2]<-"Biogeo_status"  
tmp_c <- dplyr::filter (tmp_c, Biogeo_status == "Cla" ) 

tmp_a <- rbind(tmp_b, tmp_c)

unique (tmp_a$Biogeo_status)
tmp_a$Biogeo_status <- gsub("NNE"       , "aNNE"       , tmp_a$Biogeo_status)
tmp_a$Biogeo_status <- gsub("MAC"       , "bMAC"       , tmp_a$Biogeo_status)
tmp_a$Biogeo_status <- gsub("Cla"       , "zcla"       , tmp_a$Biogeo_status)
unique (tmp_a$Biogeo_status)

coco<- c( "black", "green4","gold2","dodgerblue3","mediumorchid2" )

boxi_PC1 <- ggplot(tmp_a, aes(x = Biogeo_status, y=PC1, color = Biogeo_status)) +
    scale_color_manual(values=coco)+
    geom_boxplot(show.legend = FALSE, lwd=0.3, outlier.size = 0.5, outlier.alpha = 0.5)+
    scale_x_discrete(breaks = NULL) + 
    labs(x = "", y = "PC1") + 
    p_my_theme1

boxi_PC2 <- ggplot(tmp_a, aes(x=Biogeo_status, y=PC2, color = Biogeo_status)) +
    scale_color_manual(values = coco)+
    geom_boxplot(show.legend = FALSE, lwd=0.3, outlier.size = 0.5, outlier.alpha = 0.5)+
    labs(x = "Status", y = "PC2") +
    p_my_theme1

c <- ggpubr::ggarrange(boxi_PC1, boxi_PC2, nrow = 2)
c

ks_tmp <- kruskal(tmp_a$PC1, trt = as.factor(tmp_a$Biogeo_status))
ks_PC1 <- ks_tmp$groups %>% rownames_to_column("status")
ks_PC1 <- ks_PC1[, c("status", "groups")]
ks_PC1$PC1 <- "pc1"

ks_tmp <- kruskal(tmp_a$PC2, trt = as.factor(tmp_a$Biogeo_status))
ks_PC2 <- ks_tmp$groups %>% rownames_to_column("status")
ks_PC2 <- ks_PC2[, c("status", "groups")]
ks_PC2$PC2 <- "pc2"

ks_pc1pc2<- cbind(ks_PC1, ks_PC2)
ks_pc1pc2

pdf(file = "Ext_Figure_3a.pdf",  
    width =5, 
    height=3, paper = "a4")
a
dev.off()

pdf(file = "Ext_Figure_3b.pdf", 
    width =6,
    height=3, paper = "a4")
corrplot(var$contrib, is.corr=FALSE, col= colorRampPalette(c("grey90","grey50", "grey10"))(30)) 
dev.off()

pdf(file = "Ext_Figure_3c.pdf",   
    width =2,
    height=2.5, paper = "a4")

c
dev.off()

# FIN