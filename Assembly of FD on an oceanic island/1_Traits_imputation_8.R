# M. Paola Barajas Barbosa et al (2022) Assembly of functional diversity in an oceanic island flora
# by M. Paola Barajas

# NOTE: We used phylogenetic trait imputation to estimate missing trait values. 
# To run analysis and figures of this study, trait imputation has to be run first

library(tidyverse)
library(reshape2)
library(missForest)
library(adephylo)
library(phytools)
library(ape)
library(pez)

library(V.PhyloMaker)
library(TNRS)

setwd("C:/PUB_2/Code_Data_sub3") 

Barajasetal <-  read.csv("Barajasetal_traits_status_20230323.csv")

traits0    <- Barajasetal  # 8 traits
traits0$sp <- gsub(" ", "_", traits0$Accepted_name_TNRS2022)

# species by trait matrix
load("Tenerife_tree.RData")

tree <- ten_phyo
tree <- ten_phyo$scenario.2 
tree <- as(tree$run.1, "phylo")
class(tree)
length(tree$tip.label) # 436 species

sp_tree = as.data.frame(tree$tip.label)
duplicated(sp_tree)

traits0$dupes = duplicated(traits0$sp)
traits0 = filter( traits0, dupes == FALSE)

length( setdiff(sp_tree$`tree$tip.label`, traits0$sp) )
length( setdiff(traits0$sp, sp_tree$`tree$tip.label`) )

traits = traits0 %>% 
  filter(sp %in% sp_tree$`tree$tip.label`) %>% 
  select(sp, LA, LDMC, LMA, Lth, Nmass, SM, SSD, H) %>%
  column_to_rownames( var = "sp")

# Calculate Moran eigenvetors to use phylogenetic correlation structure to predict traits
# Eigenvectors eliminate features that have a strong correlation 
# between them and also help in reducing over-fitting.
# create phylogenetic proximity table. 
prox.Ab.all <- proxTips(tree, method = "Abouheif", normalize="none")
dim(prox.Ab.all) # 436

prox <- prop.table(prox.Ab.all, 1) # standardize by row
prox <- 0.5 * (prox + t(prox))     # make matrix symetric

ME <- me.phylo(prox = prox) # create Moran's eigenvectors based on phylogenetic distance matrix. closely related species will have similar ME values 
ME <- ME[rownames(traits),]

trait.imp <- cbind(traits, ME[,1:30]) # Morans Eigenvectors plus species by traits matrix (with column_to_rownames). 

#  ---- Phylogenetic imputation using missForest ----
set.seed(8) 
dfk <- data.frame(matrix(NA, nrow = 30, ncol = 9)) 
colnames(dfk) <- c("k", "OOB_LA","OOB_LDMC","OOB_LMA","OOB_Lth","OOB_Nmass","OOB_SM","OOB_SSD","OOB_H")  

for (n in 1:30) {
  dfimp <- trait.imp[, 1: ( 8 +n)] # per trait a loop is run for each phylogenetic eigenvector number (all number > 10 (i.e., traits number) until 30) 
  o <- missForest(dfimp, maxiter = 25, ntree = 100 , variablewise = TRUE) 
  dfk[n, 1] <- n
  dfk[n,2] <- o$OOBerror[1] # save OOBerror for target traits only
  dfk[n,3] <- o$OOBerror[2]    
  dfk[n,4] <- o$OOBerror[3]
  dfk[n,5] <- o$OOBerror[4]
  dfk[n,6] <- o$OOBerror[5]
  dfk[n,7] <- o$OOBerror[6]
  dfk[n,8] <- o$OOBerror[7]
  dfk[n,9] <- o$OOBerror[8]
  
}

dfk2<-dfk %>%                                                                          # imputation errors per trait
  summarize(min_LA = min(OOB_LA), k_min_LA = k[which.min(OOB_LA)], 
            min_LDMC = min(OOB_LDMC), k_min_LDMC = k[which.min(OOB_LDMC)],
            min_LMA = min(OOB_LMA), k_min_LMA = k[which.min(OOB_LMA)],
            min_Lth = min(OOB_Lth), k_min_Lth = k[which.min(OOB_Lth)],
            min_Nmass = min(OOB_Nmass), k_min_Nmass = k[which.min(OOB_Nmass)],
            min_SM = min(OOB_SM), k_min_SM = k[which.min(OOB_SM)],
            min_SSD = min(OOB_SSD), k_min_SSD = k[which.min(OOB_SSD)],
            min_H = min(OOB_H), k_min_H = k[which.min(OOB_H)]
  )

# phylogenetically-informed imputed data have generally similar/lower error rates
# Chose the number of eigenvectors that minimizes the imputation error (i.e., OOB_trait).

set.seed(8) 
LA_ideal   <-missForest(trait.imp[, 1: (10+dfk2$k_min_LA)], maxiter = 25, ntree = 100 ,   variablewise = TRUE)  
LDMC_ideal <-missForest(trait.imp[, 1: (10+dfk2$k_min_LDMC)], maxiter = 25, ntree = 100 ,   variablewise = TRUE) 
LMA_ideal  <-missForest(trait.imp[, 1: (10+dfk2$k_min_LMA)], maxiter = 25, ntree = 100 ,   variablewise = TRUE) 
Lth_ideal  <-missForest(trait.imp[, 1: (10+dfk2$k_min_Lth)], maxiter = 25, ntree = 100 ,   variablewise = TRUE) 
Nmass_ideal <-missForest(trait.imp[, 1: (10+dfk2$k_min_Nmass)], maxiter = 25, ntree = 100 ,   variablewise = TRUE) 
SM_ideal   <-missForest(trait.imp[, 1: (10+dfk2$k_min_SM)], maxiter = 25, ntree = 100 ,   variablewise = TRUE) 
SSD_ideal  <-missForest(trait.imp[, 1: (10+dfk2$k_min_SSD)], maxiter = 25, ntree = 100 ,   variablewise = TRUE) 
H_ideal    <-missForest(trait.imp[, 1: (10+ dfk2$k_min_H )], maxiter = 25, ntree = 100 ,   variablewise = TRUE) 

best_LA     <-tibble(LA=LA_ideal$ximp$LA, species=rownames(LA_ideal$ximp))
best_LDMC   <-tibble(LDMC=LDMC_ideal$ximp$LDMC, species=rownames(LDMC_ideal$ximp))
best_LMA    <-tibble(LMA=LMA_ideal$ximp$LMA, species=rownames(LMA_ideal$ximp))
best_Lth    <-tibble(Lth=Lth_ideal$ximp$Lth, species=rownames(Lth_ideal$ximp))
best_Nmass  <-tibble(Nmass=Nmass_ideal$ximp$Nmass, species=rownames(Nmass_ideal$ximp))
best_SM     <-tibble(SM = SM_ideal$ximp$SM, species=rownames(SM_ideal$ximp))
best_SSD    <-tibble(SSD = SSD_ideal$ximp$SSD, species=rownames(SSD_ideal$ximp))
best_H      <-tibble(H=H_ideal$ximp$H, species=rownames(H_ideal$ximp))

impute_out<-left_join(best_LA,best_LDMC,by="species")%>%   # traits imputed using phylogeny and OOBerror
  left_join(.,best_LMA, by="species")%>%
  left_join(.,best_Lth,by="species")%>%
  left_join(.,best_Nmass, by="species")%>%
  left_join(.,best_SM,by="species")%>%
  left_join(.,best_SSD,by="species")%>%
  left_join(.,best_H,by="species")%>%
  
  select(., species,LA,LDMC,LMA,Lth,Nmass,SM,SSD,H)

Trait_imputed_phylo <- impute_out 
Trait_imputed_phylo$Accepted_name_TNRS2022 <- gsub("_" , " "    , Trait_imputed_phylo$species)
Trait_imputed_phylo$species = NULL
z <- Barajasetal[, c("IDmaster" , "Accepted_name_TNRS2022", "Species1")]

Trait_imputed_phylo <- left_join(Trait_imputed_phylo, z,  by = "Accepted_name_TNRS2022") # add ID_master

Trait_imputed_phylo$dupes <- duplicated(Trait_imputed_phylo$Accepted_name_TNRS2022)

Trait_imputed_phylo <- filter(Trait_imputed_phylo, 
                                Species1 != "Hyparrhenia sinaica" & 
                                  Species1 != "Lolium edwardii" &
                                  Species1 != "Patellifolia webbiana"& 
                                  Species1 != "Pimpinella rupicola"&
                                  Species1 != "Scilla dasyantha") 

Trait_imputed_phylo$dupes <- NULL 
Trait_imputed_phylo$IDmaster <- as.character(Trait_imputed_phylo$IDmaster)

setdiff(sp_tree$`tree$tip.label`, traits0$sp) # 7 species missing

#  ---- Naive trait imputation ----
# here we do not used phylogenetic information. 

traits_1 <- Barajasetal[, c("IDmaster","LA","LDMC","LMA", "Lth", "Nmass","SM","SSD", "H")] # selects trait only
traits_1 <- column_to_rownames(traits_1, var ="IDmaster") 

naive <- missForest(traits_1, maxiter = 25, ntree = 100 , variablewise = TRUE) 

Trait_imputed_naive <- naive$ximp # ximp	contains the imputed data matrix. 
Trait_imputed_naive["IDmaster"] <- rownames(Trait_imputed_naive)

# # 

Trait_imputed_phylo_1 = left_join(Trait_imputed_naive, Trait_imputed_phylo, by = "IDmaster" )
summary(Trait_imputed_phylo_1)

Trait_imputed_phylo_2 <- Trait_imputed_phylo_1 %>%   # replace the very few NAs with the naive imputed data
  mutate(LDMC.y = coalesce(LDMC.y, LDMC.x)) %>%
  mutate(LMA.y = coalesce(LMA.y, LMA.x)) %>%
  mutate(Lth.y = coalesce(Lth.y, Lth.x)) %>%
  mutate(Nmass.y = coalesce(Nmass.y, Nmass.x)) %>%
  mutate(SM.y = coalesce(SM.y, SM.x)) %>%
  mutate(SSD.y = coalesce(SSD.y, SSD.x)) %>%
  mutate(H.y = coalesce(H.y, H.x))

Trait_imp_phylo <- Trait_imputed_phylo_2[, 9:19]
Trait_imp_phylo <- data.frame(LA = Trait_imp_phylo$LA.y, 
                              LDMC = Trait_imp_phylo$LDMC.y, 
                              LMA = Trait_imp_phylo$LMA.y, 
                              Lth = Trait_imp_phylo$Lth.y, 
                              Nmass = Trait_imp_phylo$Nmass.y, 
                              SM = Trait_imp_phylo$SM.y, 
                              SSD = Trait_imp_phylo$SSD.y, 
                              H = Trait_imp_phylo$H.y, 
                              IDmaster = as.integer(Trait_imp_phylo$IDmaster), 
                              IDmaster = Trait_imp_phylo$IDmaster)


Trait_imp_phylo$IDmaster.1 = NULL
write.csv(Trait_imp_phylo, "Trait_imputed_phylo8.csv",row.names=F) 

#  ---- Imputation errors ----
# Naive errors
OOB_multi<-data.frame(rbind(naive$OOBerror))[,1:8]
colnames(OOB_multi) <- c("OOB_LA","OOB_LDMC","OOB_LMA","OOB_Lth",
                         "OOB_Nmass","OOB_SM", "OOB_SSD", "OOB_H")
OOB_multi$trait<-"naive"   

errors_naive <-OOB_multi
errors_phylo_info <- dfk2
errors_phylo_info <- errors_phylo_info[c ("min_LA", "min_LDMC", "min_LMA", "min_Lth",    
                                          "min_Nmass", "min_SM", "min_SSD", "min_H" )]

# Phylogenetically informed errors
errors_phylo_info <- data.frame(OOB_LA  = errors_phylo_info$min_LA, 
                                OOB_LDMC= errors_phylo_info$min_LDMC, 
                                OOB_LMA =errors_phylo_info$min_LMA,
                                OOB_Lth =errors_phylo_info$min_Lth, 
                                OOB_Nmass =errors_phylo_info$min_Nmass, 
                                OOB_SM  =errors_phylo_info$min_SM, 
                                OOB_SSD =errors_phylo_info$min_SSD, 
                                OOB_H   =errors_phylo_info$min_H )

errors_phylo_info$trait <- "phylo"
OOB_errors = rbind(errors_phylo_info, errors_naive)
OOB_errors

# FIN