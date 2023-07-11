# M. Paola Barajas Barbosa et al (2023) Assembly of functional diversity in an oceanic island flora

# By Dylan Craven
#-----------------------------------#
#hot-node analysis                  #
#(sensu Molina-Venegas et al. 2022) #
#-----------------------------------#

require(tidyverse)
require(phytools)
require(picante)
require(pez)
require(V.PhyloMaker)
require(TNRS)
require(caper)
require(ggtree)
require(ape)
require(geiger)
require(ggsci)

setwd("C:/PUB_2")

source("hot_nodes.R")

load(file="Phylo_Sig_out.RData" ) # from results in code Extended_data_I_PhyloSignal

#----------------------#
# Load & Organize data #
#----------------------#

spp_nodes<-trts_tenn %>% 
           dplyr::select(.,spp_TNRS,Biogeo_status) %>% 
           mutate(PA=1,
                  spp_TNRS=str_replace_all(spp_TNRS, " ","_")) %>% 
           pivot_wider(., id_cols=spp_TNRS, 
                       names_from=Biogeo_status,
                       values_from=PA,
                       values_fill=0) %>% 
           drop_na(.) %>% 
           tibble::column_to_rownames(., var="spp_TNRS")

# create vector for each group
spp_CE<-spp_nodes$CE
names(spp_CE)<-rownames(spp_nodes)

spp_TE<-spp_nodes$TE
names(spp_TE)<-rownames(spp_nodes)

spp_NNE<-spp_nodes$NNE
names(spp_NNE)<-rownames(spp_nodes)

spp_MAC<-spp_nodes$MAC
names(spp_MAC)<-rownames(spp_nodes)

#----------------------------#
# Hot nodes                  #
#----------------------------#

hotnodes_out<-list(); 

for (i in 1:100) {
  
  # CE
  phy_CE <- drop.tip(ten_phyo$scenario.2[[i]], 
                        setdiff(ten_phyo$scenario.2[[i]]$tip.label, names(spp_CE)))
  
  phy_match_CE<-match.phylo.data(phy_CE, spp_CE)
  
  hotnodes_CE <- hot.nodes(phy_match_CE$data, phy_match_CE$phy, 
                          runs = 1000, signif = 1.96, min.size = 4, 
                          positive = TRUE, verbose=FALSE) 
  
  hotnodes_CE<-data.frame(hotnodes_CE$`Hot nodes`) %>% 
                dplyr::select(., Node, Size, Usable, SES) %>% 
                mutate(Biogeo_status="CE")
  
  # TE
  phy_TE <- drop.tip(ten_phyo$scenario.2[[i]], 
                     setdiff(ten_phyo$scenario.2[[i]]$tip.label, names(spp_TE)))
  
  phy_match_TE<-match.phylo.data(phy_TE, spp_TE)
  
  hotnodes_TE <- hot.nodes(phy_match_TE$data, phy_match_TE$phy, 
                           runs = 1000, signif = 1.96, min.size = 4, 
                           positive = TRUE, verbose=FALSE) 
  
  hotnodes_TE<-data.frame(hotnodes_TE$`Hot nodes`)%>% 
                dplyr::select(., Node, Size, Usable, SES) %>% 
                mutate(Biogeo_status="TE")
  # MAC
  
  # NNE
  phy_NNE <- drop.tip(ten_phyo$scenario.2[[i]], 
                     setdiff(ten_phyo$scenario.2[[i]]$tip.label, names(spp_NNE)))
  
  phy_match_NNE<-match.phylo.data(phy_NNE, spp_NNE)
  
  hotnodes_NNE <- hot.nodes(phy_match_NNE$data, phy_match_NNE$phy, 
                           runs = 1000, signif = 1.96, min.size = 4, 
                           positive = TRUE, verbose=FALSE) 
  
  hotnodes_NNE<-data.frame(hotnodes_NNE$`Hot nodes`) %>% 
                dplyr::select(., Node, Size, Usable, SES) %>% 
                mutate(Biogeo_status="NNE")
  # MAC
  phy_MAC <- drop.tip(ten_phyo$scenario.2[[i]], 
                     setdiff(ten_phyo$scenario.2[[i]]$tip.label, names(spp_MAC)))
  
  phy_match_MAC<-match.phylo.data(phy_MAC, spp_MAC)
  
  hotnodes_MAC <- hot.nodes(phy_match_MAC$data, phy_match_MAC$phy, 
                           runs = 1000, signif = 1.96, min.size = 4, 
                           positive = TRUE, verbose=FALSE) 
  
  hotnodes_MAC<-data.frame(hotnodes_MAC$`Hot nodes`) %>% 
               dplyr::select(., Node, Size, Usable, SES) %>% 
               mutate(Biogeo_status="MAC")
  
  # merge
  
  nodes<-bind_rows(hotnodes_CE, hotnodes_TE)
  nodes<-bind_rows(nodes, hotnodes_NNE)
  nodes<-bind_rows(nodes, hotnodes_MAC)
  nodes$iteration<-i
  
  cat("progress", i, sep=' ','\n')
  hotnodes_out[[i]]<-rbind.data.frame(nodes)
}

hotnodes_out<-bind_rows(hotnodes_out)

save(hotnodes_out, 
     file="Hotspot_out.RData" )

#-------------------#
# Visualisation     #
#-------------------#

load( file="Hotspot_out.RData" )

# remover Nodes <10 species
# remover Nodes que aparecen en menos de 5 filogenias

# get node information
spp_nodes<-trts_tenn %>% 
            dplyr::select(.,species=spp_TNRS,genus=genus_TRNS,family=family_TNRS) %>% 
            mutate(species=str_replace_all(species, " ","_")) 

ten_phyo$scenario.2$run.1$node.label<-paste("N",1:ten_phyo$scenario.2$run.1$Nnode,sep="")

nodes<-build.nodes.1(ten_phyo$scenario.2$run.1, spp_nodes)

spp_nodess<-spp_nodes %>% 
            mutate(spp_TNRSs=str_replace_all(species,"_"," "))

ten_phyo$scenario.2$run.1$tip.label <- with(spp_nodess, spp_TNRSs[match(ten_phyo$scenario.2$run.1$tip.label, species)])

nodess<-nodes %>% 
        dplyr::select(., Node=bn,family) %>% 
        distinct(.) 

hotnodes_outt<-hotnodes_out %>% 
              mutate(Node=paste("N",Node, sep="")) %>% 
              dplyr::filter(Size>9) %>% 
              group_by(Biogeo_status,Node) %>% 
              summarize(xSES=mean(SES),
                        n=length(unique(iteration))) %>% 
              dplyr::filter(.,Biogeo_status=="TE") %>% 
              left_join(.,nodess, by="Node") %>% 
              ungroup(.) %>% 
              mutate(dupe_Node=duplicated(Node),
                     family=if_else(Node=="N165" & Biogeo_status=="NNE","Aizoaceae/Portulacaceae",family),
                     family=if_else(Node=="N256"& Biogeo_status=="MAC","Cucurbitaceae/Myricaceae",family)) %>% 
              dplyr::select(.,-dupe_Node) %>%
              distinct(.)

hotnodes_outt<- hotnodes_outt%>%
                dplyr::select(., Node, Biogeo_status, family,xSES,n )

# Tree

biogeo_tree<-ggtree(ten_phyo$scenario.2$run.1, 
                    layout="circular",
                    size=0.25,
                    color="gray30")+
             geom_tiplab2(aes(angle=angle),size=1.5, hjust=-0.1, color="gray30")

biogeo_tree<-biogeo_tree %<+% 
              hotnodes_outt+
              geom_nodepoint(aes(color=Biogeo_status, size= xSES),
                             color="#5ab4ac",fill="#5ab4ac",alpha=0.8)+
              scale_color_npg(na.translate = F)+
              guides(colour=guide_legend(title="Biogeographical\nstatus"), 
                     size="none", alpha="none")+
              theme(legend.position=c(0.93,0.88),
                    legend.title = element_text(size=7, color="black",face="bold"),
                    legend.text = element_text(size=6, color="black"))

ggsave(biogeo_tree,filename="Hotspot_phylogeny.png", 
       units="cm", 
       width=21, 
       height=21, 
       #pointsize=2, 
       dpi=300)

# FIN