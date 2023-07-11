# M. Paola Barajas Barbosa et al (2023) Assembly of functional diversity in an oceanic island flora
# by M. Paola Barajas

# # Figure 1

library(plyr)
library(tidyr)
library(ggplot2)
library(ggfortify)
library(ggpubr)
library(devtools)
library(cowplot)
library(hypervolume)
library(dplyr)
library(visdat)
library(purrr)
library(tidyverse)

setwd("C:/PUB_2/Code_Data_sub3")

# # Global and Tenerife island data  
Diazetal <- read.csv("diaz_etal_new_20221202.csv")
Diazetal <- dplyr::select(Diazetal, Species1 = Species.name.standardized.TPL.2016, Accepted_name_TNRS2022 = Accepted_name_TNRS2022, growth_form,
                          LA               = Leaf.Area..mm2.,                     # mm²
                          LMA              = LMA..g.m2.,                          # g/m²
                          Nmass            = Nmass..mg.g. ,                       # mg/g
                          SM               = Seed.Mass..mg.,                      # mg
                          SSD              = SSD.combined..mg.mm3.,               # mg/mm³
                          H                = Plant.Height..m.  )                  # m

Diazetal$source <- 'Global_data'
Diazetal$dupes = duplicated(Diazetal$Species1)

data_Global <- data.frame(
  Species1      = Diazetal$Species1,
  Leaf_area     = scale(log10(Diazetal$LA)),
  LMA           = scale(log10(Diazetal$LMA)),
  Leaf_N        = scale(log10(Diazetal$Nmass)),
  Seed_mass     = scale(log10(Diazetal$SM)),
  Stem_density  = scale(log10(Diazetal$SSD)),
  Height        = scale(log10(Diazetal$H)),
  source        = Diazetal$source,
  growth_form   = Diazetal$growth_form)

row.names(data_Global) <- data_Global$Species1

# Tenerife data
Barajasetal <-  read.csv("Barajasetal_traits_status_20230323.csv")

# Tenerife island - with species having minimum 5 traits (i.e., 348 species)
data_Tenerife <- Barajasetal[, 1:16]
data_Tenerife$missing_trait <- rowSums(is.na(data_Tenerife))
data_Tenerife <- filter(data_Tenerife, missing_trait < 4)

data_Tenerife_imp <- read.csv("Trait_imputed_phylo8.csv")
data_Tenerife_imp$Species1 <- NULL
data_Tenerife <- left_join (data_Tenerife, data_Tenerife_imp, by = "IDmaster")

data_Tenerife <- data_Tenerife %>%
  mutate(LDMC.x = coalesce(LDMC.x, LDMC.y)) %>%
  mutate(LMA.x = coalesce(LMA.x, LMA.y)) %>%

  mutate(Lth.x = coalesce(Lth.x, Lth.y)) %>%
  mutate(Nmass.x = coalesce(Nmass.x, Nmass.y)) %>%

  mutate(SM.x = coalesce(SM.x, SM.y)) %>%
  mutate(SSD.x = coalesce(SSD.x, SSD.y)) %>%
  mutate(H.x = coalesce(H.x, H.y))
anyNA(data_Tenerife)
data_Tenerife$source <- 'Tenerife_data'; data_Tenerife$growth_form <- NA

data_Tenerife <- data.frame(
  Species1      = data_Tenerife$Species1,

  Leaf_area     = scale(log10(data_Tenerife$LA.x)),
  LMA           = scale(log10(data_Tenerife$LMA.x)),
  Leaf_N        = scale(log10(data_Tenerife$Nmass.x)),
  Seed_mass     = scale(log10(data_Tenerife$SM.x)),
  Stem_density  = scale(log10(data_Tenerife$SSD.x)),
  Height        = scale(log10(data_Tenerife$H.x)),

  source       = data_Tenerife$source,
  growth_form  = data_Tenerife$growth_form)
#

# Remove Tenerife species from Global Data
data_Tenerife_imp <- read.csv("Trait_imputed_phylo8.csv")
tmp0 <- data.frame(Species1      = data_Tenerife_imp$Species1,   Leaf_area     = data_Tenerife_imp$LMA,
  LMA           = data_Tenerife_imp$LMA,  Leaf_N        = data_Tenerife_imp$LMA,
  Seed_mass     = data_Tenerife_imp$LMA,  Stem_density  = data_Tenerife_imp$LMA,
  Height        = data_Tenerife_imp$LMA,  source        = as.character(data_Tenerife_imp$LMA),   
  growth_form   = as.character(data_Tenerife_imp$LMA))

tmp1 <- rbind(tmp0, data_Global)     # First full Tenerife data to be able to delete dupes from Diaz 
tmp1$dupes <- duplicated(tmp1$Species1) 
tmp1 <- filter(tmp1, dupes == FALSE)    
anyDuplicated(tmp1$Species1)

## 
gl <- dplyr::filter (tmp1, source == "Global_data" ) # 2278 nearly complete cases. 
gl$dupes <- NULL

te <- data_Tenerife

comparison <- rbind(gl, te)

save(comparison, file = "Comparison_20221202.RData")

# FIGURE 1 ----

PCA         <- prcomp(comparison[, c("Leaf_area","LMA","Leaf_N","Seed_mass","Stem_density","Height")])
summary(PCA) # variance explained
PCAvalues   <- data.frame(Species = comparison$Species1, source = comparison$source, growth_form= comparison$growth_form, PCA$x)# Extract PC axes for plotting
PCAloadings <- data.frame(Variables = rownames(PCA$rotation), PCA$rotation)          # Extract loading of the variables

PCAvalues$PC2 <- (PCAvalues$PC2*-1)
PCAloadings$PC2 <- (PCAloadings$PC2*-1)

p_my_theme <-  theme( legend.text = element_text(size=7), axis.title= element_text(size=7),
                      axis.title.x= element_text(colour="black", size=7),
                      axis.title.y= element_text(colour="black", size=7),
                      axis.text.x= element_text(colour= "black", size=7),
                      axis.text.y= element_text(colour= "black", size=7),
                      panel.background =element_rect(fill="transparent",colour="black"),panel.grid.minor=element_blank(),
                      panel.border=element_rect(fill=NA,colour="grey50"))

p_my_theme2 <- theme( legend.text = element_text(size=7), axis.title= element_text(size=7),
                      axis.text.x= element_text(colour= "black", size=7), axis.text.y= element_text(colour= "black", size=7))

# Island versus global trait space  
Figure1 <- ggplot(PCAvalues) +   
  geom_point(size = .7, alpha=0.7, 
             
             aes(x = PC1, y = PC2, group = source, colour = source, size = source, label = Species),  # label = T
             
             show.legend = FALSE) +
  
  # geom_text( aes(x = PC1, y = PC2, label = Species, colour = source, size = source), hjust = 0)+  # add species names
  
  coord_fixed() + 
  scale_colour_manual(values=c("gray65","turquoise4")) +
  geom_segment(data = PCAloadings, size = 0.25,
               aes(x = 0, xend = PC1*3.5, y = 0, yend = PC2*3.5),
               arrow = arrow(length = unit(0.1, "cm")),colour = "black")   +
  geom_text(data = PCAloadings, aes(x = PC1*3.6, y = PC2*3.6, label = Variables), size = 2.3,
            hjust=c(0, 0, 0, 0, 0, 0) , vjust=c(0, 0, 0, 0, 0, 1))    + 
  xlab("PC1 (47%)   island n = 237") + ylab("PC2 (25.3%) Globe n = 2274")   +
  p_my_theme 
Figure1

# Marginal density distribution island and global PCA  
xplot <- ggdensity(PCAvalues, "PC1", fill = "source", color = "source", palette = c("gray65","turquoise4")) + p_my_theme2
yplot <- ggdensity(PCAvalues, "PC2", fill = "source", color = "source", palette = c("gray65","turquoise4")) + rotate() + p_my_theme2

##
pdf(file = "_Figure_1.pdf",
    width =4,     
    height=5.3,
    paper = "a4")

cowplot::plot_grid(xplot, NULL, Figure1, yplot, ncol = 2, align = "hv", 
                   rel_widths = c(2, 1), rel_heights = c(1, 2))
dev.off()
##

# Figure 1 # Null model ---- 
# Controls for the difference of number of species between the global versus the island trait space.
# Controls for potential biased sampling intensity of growth form in the global data.
# It computes (999 times) the global hypervolume by randomly sampling the same number.
# of species of Tenerife, i.e., 300 from the the global data.     
# Using global and island hypervolume we compute overlap statistics.

### Null models with nearly complete cases (n = 348)
### Rarefy by 300 species
### number from global data: tree, herb, shrub: 64, 163, 73

overall_bandwidth <- estimate_bandwidth(PCAvalues[, c("PC1","PC2","PC3")])

nm_glob <- filter(PCAvalues, source == "Global_data")
nm_isla <- filter(PCAvalues, source == "Tenerife_data")

nbperm <- 999 # permutations

null_mod_2 <- c()
for (i in 1:nbperm){  
  
  tmp0 <- nm_isla %>% sample_n(300) # rarefying by minimum common species
  head(tmp0)
  
  set.seed(8)
  HV_Ten_tmp <-hypervolume_gaussian(tmp0[,c("PC1","PC2","PC3")], name = "Island volume",
                                    kde.bandwidth=overall_bandwidth, quantile.requested=0.95)
  
  # Samples 21% trees + 54% herbs + 24% shrubs  
  tmp = nm_glob %>%   
    group_by(growth_form) %>% 
    nest() %>%            
    ungroup() %>% 
    mutate(n = c(64, 163, 73)) %>% # tree, herb, shrub
    mutate(samp = map2(data, n, sample_n)) %>% 
    select(-data) %>%
    unnest(samp)
  
  HV_glob_tmp <- hypervolume_gaussian(tmp[, c("PC1","PC2","PC3")], name = "Global volume",
                                      kde.bandwidth=overall_bandwidth, quantile.requested=0.95)
  
  HVs         <- hypervolume_join(HV_Ten_tmp, HV_glob_tmp)
  HV_set      <- hypervolume_set(HV_Ten_tmp, HV_glob_tmp, num.points.max = NULL,
                                 verbose = TRUE, check.memory = F, distance.factor = 1)
  
  overlap     <- hypervolume_overlap_statistics(HV_set)
  
  null_mod_2  <- rbind(null_mod_2, overlap)
  
  cat(paste0(round(100*i/nbperm, 0), " %; "))
  
}

View(null_mod_2)

save(null_mod_2, file = "results_Figure_1_nullmodel_2022117.RData") 

# # ***
load("results_Figure_1_nullmodel_basic_2022117.RData") 
load("results_Figure_1_nullmodel_2022117.RData") 

coco = c("gold2","gray65","turquoise4") 

# Null mode controlling for growth form
null_mod_2 = as.data.frame(null_mod_2)
colnames(null_mod_2) <- c("jaccard","Similar_component_sorensen","Unique_component_Island","Unique_component_globe")

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

overlap_stats <- ggplot(null_models, label= mean)+
  geom_linerange(aes(x = metric, y = mean*100, ymin= mean*100-ic*100, ymax= mean*100+ic*100, color = ID), show.legend = FALSE) +
  geom_point(aes(x = metric, y = mean*100, color = ID),  size = 4,show.legend = FALSE) +
  scale_y_continuous(breaks = seq(0, 100, by = 10))+
  scale_colour_manual(values = coco) +  
  labs(x = "", y = "Hypervolume overlap statistics (%) complete cases only") + p_my_theme

# Null model basic
null_mod_1 = as.data.frame(null_mod_1)
colnames(null_mod_1) <- c("jaccard","Similar_component_sorensen","Unique_component_Island","Unique_component_globe")

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

overlap_stats_basic <- ggplot(null_models_basic, label= mean)+
  geom_linerange(aes(x = metric, y = mean*100, ymin= mean*100-ic*100, ymax= mean*100+ic*100, color = ID), show.legend = FALSE) +
  geom_point(aes(x = metric, y = mean*100, color = ID),  size = 4,show.legend = FALSE) +
  scale_y_continuous(breaks = seq(0, 100, by = 10))+
  scale_colour_manual(values = coco) +  
  labs(x = "", y = "Hypervolume overlap statistics (%) complete cases only") + p_my_theme

ggpubr::ggarrange( overlap_stats, overlap_stats_basic)

# FIN

# Brifing figure
PCAvalues$briefing <- ""
rownames(PCAvalues) <- 1:nrow(PCAvalues)

PCAvalues[2390, 10] <-"Echium w"  
PCAvalues[2376, 10] <-"DRACO"  

ggplot(PCAvalues) +   
  geom_point(size = 2, alpha=.5, 
             
             aes(x = PC1, y = PC2, group = source, colour = source, size = source, label = briefing),  # label = T
             
             show.legend = FALSE) +
  
  geom_text( aes(x = PC1, y = PC2, label = briefing, colour = source, size = source), hjust = 0)+  # add species names
  
  coord_fixed() + 
  scale_colour_manual(values=c("gray65","turquoise4")) +
  geom_segment(data = PCAloadings, size = 0.25,
               aes(x = 0, xend = PC1*3.5, y = 0, yend = PC2*3.5),
               arrow = arrow(length = unit(0.1, "cm")),colour = "black")   +
  geom_text(data = PCAloadings, aes(x = PC1*3.6, y = PC2*3.6, label = Variables), size = 2.3,
            hjust=c(0, 0, 0, 0, 0, 0) , vjust=c(0, 0, 0, 0, 0, 1))    + 
  xlab("PC1 (47%)   island n = 237") + ylab("PC2 (25.3%) Globe n = 2274")   +
  p_my_theme 





