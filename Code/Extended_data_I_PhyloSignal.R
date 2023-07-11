# M. Paola Barajas Barbosa et al (2023) Assembly of functional diversity in an oceanic island flora

# By Dylan Craven
# Phylogenetic signal of each trait     #
# Phylogenetic signal of missing traits #

require(tidyverse)
require(phytools)
require(picante)
require(pez)
require(V.PhyloMaker)
require(TNRS)
require(caper)

#----------------------#
# Load & Organize data #
#----------------------#

trts_ten<- read.csv("Barajasetal_traits_status_20230323.csv")

# adjust species names
ten_spp<-trts_ten %>% 
         dplyr::filter(.,!Accepted_name_TNRS2022=="") %>%
         dplyr::select(.,Accepted_name_TNRS2022) %>% 
         distinct(.)
ten_spp<-ten_spp$Accepted_name_TNRS2022

ten_families<-TNRS(taxonomic_names=ten_spp,sources=c("tropicos","wfo"))
ten_familiess<-ten_families %>% 
               dplyr::select(.,Accepted_name_TNRS2022=Name_submitted, Name_matched_accepted_family) %>% 
               distinct(.)

# get family names (APGII)
trts_tenn<-left_join(trts_ten, ten_familiess, by="Accepted_name_TNRS2022")
ten_extra<-TNRS(taxonomic_names="Teline salsoloides",sources=c("tropicos","wfo"))

# adjustments
trts_tenn<-trts_tenn %>% 
           mutate(Name_matched_accepted_family=if_else(is.na(Name_matched_accepted_family)==TRUE,
                                                       "Fabaceae",Name_matched_accepted_family),
                  Accepted_name_TNRS2022=if_else(Accepted_name_TNRS2022=="", 
                                             "Teline salsoloides",Accepted_name_TNRS2022),
                  Accepted_name_TNRS2022=if_else(dupes_TNRS==TRUE,Species1,Accepted_name_TNRS2022)) %>% 
           dplyr::select(.,spp_TNRS=Accepted_name_TNRS2022,family_TNRS=Name_matched_accepted_family,
                           dupes_TNRS,
                           Biogeo_status,Ende_status,
                           LA, LDMC, LMA, Lth, Nmass, SM,SSD, H) %>% 
          mutate(dupes_spp=duplicated(spp_TNRS))

# adjust family names for V.PhyloMaker
trts_tenn<-trts_tenn %>% 
           mutate(spp_TNRS=ifelse(spp_TNRS=="Scilla haemorrhoidalis" & dupes_spp==FALSE,
                  "Scilla dasyantha",spp_TNRS),
                  spp_TNRS=ifelse(spp_TNRS=="Patellifolia procumbens" & dupes_spp==FALSE,
                  "Patellifolia patellaris",spp_TNRS),
                  genus_TRNS=word(spp_TNRS,1, sep=" "),
                  family_TNRS=ifelse(genus_TRNS=="Heliotropium","Boraginaceae",family_TNRS),
                  family_TNRS=ifelse(genus_TRNS=="Limonium","Plumbaginaceae",family_TNRS),
                  family_TNRS=ifelse(genus_TRNS=="Sambucus","Adoxaceae",family_TNRS),
                  family_TNRS=ifelse(genus_TRNS=="Semele","Asparagaceae",family_TNRS),
                  family_TNRS=ifelse(genus_TRNS=="Viburnum","Adoxaceae",family_TNRS))

#----------------------------#
# Build phylogenies          #
#----------------------------#

trts_spp<-trts_tenn %>% 
          dplyr::select(., spp_TNRS, genus_TRNS, family_TNRS) %>% 
          distinct(.) %>% 
          data.frame(.)

ten_phyo<-phylo.maker(trts_spp,tree = GBOTB.extended, scenarios="S2",r=100)

#---------------------------#
# Phylogenetic signal       #
# of each traits            #
#---------------------------#

ten_LA<-trts_tenn %>% 
        dplyr::select(spp_TNRS, LA) %>% 
        mutate(spp_TNRS=str_replace_all(spp_TNRS, " ","_")) %>% 
        drop_na(.) %>% 
        tibble::column_to_rownames(., var="spp_TNRS")

ten_LDMC<-trts_tenn %>% 
          dplyr::select(spp_TNRS, LDMC) %>% 
          mutate(spp_TNRS=str_replace_all(spp_TNRS, " ","_")) %>% 
          drop_na(.) %>% 
          tibble::column_to_rownames(., var="spp_TNRS")

ten_LMA<-trts_tenn %>% 
          dplyr::select(spp_TNRS, LMA) %>% 
          mutate(spp_TNRS=str_replace_all(spp_TNRS, " ","_")) %>% 
          drop_na(.) %>% 
          tibble::column_to_rownames(., var="spp_TNRS")

ten_Lth<-trts_tenn %>% 
          dplyr::select(spp_TNRS, Lth) %>% 
          mutate(spp_TNRS=str_replace_all(spp_TNRS, " ","_")) %>% 
          drop_na(.) %>% 
          tibble::column_to_rownames(., var="spp_TNRS")

ten_Nmass<-trts_tenn %>% 
            dplyr::select(spp_TNRS, Nmass) %>% 
            mutate(spp_TNRS=str_replace_all(spp_TNRS, " ","_")) %>% 
            drop_na(.) %>% 
            tibble::column_to_rownames(., var="spp_TNRS")

ten_SM<-trts_tenn %>% 
        dplyr::select(spp_TNRS, SM) %>% 
        mutate(spp_TNRS=str_replace_all(spp_TNRS, " ","_")) %>% 
        drop_na(.) %>% 
        tibble::column_to_rownames(., var="spp_TNRS")

ten_SSD<-trts_tenn %>% 
          dplyr::select(spp_TNRS, SSD) %>% 
          mutate(spp_TNRS=str_replace_all(spp_TNRS, " ","_")) %>% 
          drop_na(.) %>% 
          tibble::column_to_rownames(., var="spp_TNRS")

ten_H<-trts_tenn %>% 
        dplyr::select(spp_TNRS, H) %>% 
        mutate(spp_TNRS=str_replace_all(spp_TNRS, " ","_")) %>% 
        drop_na(.) %>% 
        tibble::column_to_rownames(., var="spp_TNRS")

PhySig_out<-list(); 

for (i in 1:100) {
  # LA
  phy_up_la <- drop.tip(ten_phyo$scenario.2[[i]], 
                        setdiff(ten_phyo$scenario.2[[i]]$tip.label, rownames(ten_LA)))
  
  la_match<-match.phylo.data(phy_up_la, ten_LA)
  la_phy<-phylosig(la_match$phy, la_match$data[,1], method="lambda", test=FALSE, nsim=1000)
  la_phy<-cbind.data.frame(Trait="LA",Lambda=la_phy$lambda)
  
  # LDMC
  phy_up_ldmc <- drop.tip(ten_phyo$scenario.2[[i]], 
                          setdiff(ten_phyo$scenario.2[[i]]$tip.label, rownames(ten_LDMC)))
  
  ldmc_match<-match.phylo.data(phy_up_ldmc, ten_LDMC)
  ldmc_phy<-phylosig(ldmc_match$phy, ldmc_match$data[,1], method="lambda", test=FALSE, nsim=1000)
  ldmc_phy<-cbind.data.frame(Trait="LDMC",Lambda=ldmc_phy$lambda)
  
  # LMA
  phy_up_lma <- drop.tip(ten_phyo$scenario.2[[i]], 
                         setdiff(ten_phyo$scenario.2[[i]]$tip.label, rownames(ten_LMA)))
  
  lma_match<-match.phylo.data(phy_up_lma, ten_LMA)
  lma_phy<-phylosig(lma_match$phy, lma_match$data[,1], method="lambda", test=FALSE, nsim=1000)
  lma_phy<-cbind.data.frame(Trait="LMA",Lambda=lma_phy$lambda)
  
  # Lth
  phy_up_lth <- drop.tip(ten_phyo$scenario.2[[i]], 
                         setdiff(ten_phyo$scenario.2[[i]]$tip.label, rownames(ten_Lth)))
  
  lth_match<-match.phylo.data(phy_up_lth, ten_Lth)
  lth_phy<-phylosig(lth_match$phy, lth_match$data[,1], method="lambda", test=FALSE, nsim=1000)
  lth_phy<-cbind.data.frame(Trait="Lth",Lambda=lth_phy$lambda)
  
  # Nmass
  phy_up_nmass <- drop.tip(ten_phyo$scenario.2[[i]], 
                           setdiff(ten_phyo$scenario.2[[i]]$tip.label, rownames(ten_Nmass)))
  
  nmass_match<-match.phylo.data(phy_up_nmass, ten_Nmass)
  nmass_phy<-phylosig(nmass_match$phy, nmass_match$data[,1], method="lambda", test=FALSE, nsim=1000)
  nmass_phy<-cbind.data.frame(Trait="Nmass",Lambda=nmass_phy$lambda)
  
  # SM
  phy_up_sm <- drop.tip(ten_phyo$scenario.2[[i]], 
                        setdiff(ten_phyo$scenario.2[[i]]$tip.label, rownames(ten_SM)))
  
  sm_match<-match.phylo.data(phy_up_sm, ten_SM)
  sm_phy<-phylosig(sm_match$phy, sm_match$data[,1], method="lambda", test=FALSE, nsim=1000)
  sm_phy<-cbind.data.frame(Trait="SM",Lambda=sm_phy$lambda)
  
  # SSD
  phy_up_ssd <- drop.tip(ten_phyo$scenario.2[[i]], 
                         setdiff(ten_phyo$scenario.2[[i]]$tip.label, rownames(ten_SSD)))
  
  ssd_match<-match.phylo.data(phy_up_ssd, ten_SSD)
  ssd_phy<-phylosig(ssd_match$phy, ssd_match$data[,1], method="lambda", test=FALSE, nsim=1000)
  ssd_phy<-cbind.data.frame(Trait="SSD",Lambda=ssd_phy$lambda)
  
  # H
  phy_up_H <- drop.tip(ten_phyo$scenario.2[[i]], 
                       setdiff(ten_phyo$scenario.2[[i]]$tip.label, rownames(ten_H)))
  
  H_match<-match.phylo.data(phy_up_H, ten_H)
  H_phy<-phylosig(H_match$phy, H_match$data[,1], method="lambda", test=FALSE, nsim=1000)
  H_phy<-cbind.data.frame(Trait="H",Lambda=H_phy$lambda)
  
  a<-bind_rows(la_phy, ldmc_phy) 
  a<-bind_rows(a, lma_phy) 
  a<-bind_rows(a,lth_phy )  
  a<-bind_rows(a,nmass_phy)  
  a<-bind_rows(a,sm_phy)  
  a<-bind_rows(a, ssd_phy)  
  out<-bind_rows(a, H_phy)
  out$Iteration<-i
  cat("progress", i, sep=' ','\n')
  PhySig_out[[i]]<-rbind.data.frame(out)
  }

PhySig_out<-bind_rows(PhySig_out)

#---------------------------#
# Phylogenetic signal       #
# of missing traits         #
#---------------------------#

binary_ten<-trts_tenn %>% 
            mutate(spp_TNRS=str_replace_all(spp_TNRS, " ","_"),
                   LA=if_else(is.na(LA)==FALSE,0,1),
                   LDMC=if_else(is.na(LDMC)==FALSE,0,1),
                   LMA=if_else(is.na(LMA)==FALSE,0,1),
                   Lth=if_else(is.na(Lth)==FALSE,0,1),
                   Nmass=if_else(is.na(Nmass)==FALSE,0,1),
                   SM=if_else(is.na(SM)==FALSE,0,1),
                   SSD=if_else(is.na(SSD)==FALSE,0,1),
                   H=if_else(is.na(H)==FALSE,0,1)) %>% 
             mutate(complete_cases= rowSums(.[6:13]),
                    complete_cases=if_else(complete_cases==0,1,0)) %>% 
            dplyr::select(., spp_TNRS, LA, LDMC, LMA, Lth, Nmass, SM, SSD, H,complete_cases) %>% 
            tibble::column_to_rownames(., var="spp_TNRS") %>% 
            tibble::rownames_to_column(., var="spp_TNRS")

binary_comm<-cbind.data.frame(Site=1, Spp=rownames(binary_ten),PA=1) %>% 
             pivot_wider(id_cols=Site, names_from = Spp, values_from=PA ) %>% 
             tibble::column_to_rownames(., var="Site")
binary_comm<-as.matrix(binary_comm)

PhySig_miss<-list(); 

for (i in 1:100) {

  phy_up_la <- drop.tip(ten_phyo$scenario.2[[i]], 
                        setdiff(ten_phyo$scenario.2[[i]]$tip.label, binary_ten$spp_TNRS))
  
  phy_up_la$node.label<-NULL
  
  bin_data <- comparative.data(data=binary_ten, phy=phy_up_la, names.col=spp_TNRS)
  
  la_D<-phylo.d(bin_data,binvar=LA,permut = 1000)
  LDMC_D<-phylo.d(bin_data,binvar=LDMC,permut = 1000)
  Lth_D<-phylo.d(bin_data,binvar=Lth,permut = 1000)
  LMA_D<-phylo.d(bin_data,binvar=LMA,permut = 1000)
  Nmass_D<-phylo.d(bin_data,binvar=Nmass,permut = 1000)
  SM_D<-phylo.d(bin_data,binvar=SM,permut = 1000)
  SSD_D<-phylo.d(bin_data,binvar=SSD,permut = 1000)
  H_D<-phylo.d(bin_data,binvar=H,permut = 1000)
  complete_D<-phylo.d(bin_data,binvar=complete_cases,permut = 1000)
  
  bin_physig<-cbind.data.frame(LA=la_D$DEstimate, LDMC=LDMC_D$DEstimate, LMA=LMA_D$DEstimate,
                             Nmass=Nmass_D$DEstimate,
                             SM=SM_D$DEstimate, SSD=SSD_D$DEstimate, H=H_D$DEstimate,
                             complete_cases=complete_D$DEstimate) %>%
              pivot_longer(.,cols=LA:complete_cases, names_to = "Trait",values_to = "Dstat") %>% 
              mutate(iteration=i)
            
    cat("progress", i, sep=' ','\n')
  PhySig_miss[[i]]<-rbind.data.frame(bin_physig)
}

PhySig_miss<-bind_rows(PhySig_miss)

save(trts_tenn,ten_phyo,PhySig_out,PhySig_miss, 
     file="Phylo_Sig_out.RData" )

# FIN