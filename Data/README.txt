DATA DESCRIPTION 
Martha Paola Barajas Barbosa et al (2023) Assembly of functional diversity in an oceanic island flora
paolabarajas@gmail.com

-- .csv Data description --

- Barajasetal_traits_status_20230323.csv
Contains species by 8 trait data for 436 Tenerife native seed plants
(including missing data), floristic status, lineage status (citation for lineage), family
and growth form information.  

- diaz_etal_new_20221202.csv
Contains species by 6 trait data for 2288 plant species,
family and growth form information.

- Trait_imputed_phylo8.csv
Contains species by 8 trait data for 436 Tenerife native seed plants
without missing data.

The following files contain trait information at the individual level. We
aggregated the following .csv at the species level to obtain:
Barajasetal_traits_status_20221202.csv
- Trait_LDMC.csv
- Trait_leaf_area.csv
- Trait_Leaf_N.csv
- Trait_Leaf_thickness.csv
- Trait_MaxPlantHeitght.csv
- Trait_Seed_mass.csv
- Trait_Stem_specific_density.csv

-- .Rdata description --

Comparison_20221202.RData
Contains combined data, from global (Díaz et al 2016) and Tenerife (Barajas et al)
for 6 traits.  

Tenerife_tree.RData
Contains phylogeny for Tenerife native seed plant species.

results_Figure_1_nullmodel_2022117.RData
results_Figure_1_nullmodel_basic_2022117.RData
Both .RData files contain results from the null model related to Figure 1
in the main text.

results_Figure_3a_nullmodel_20230323.RData
The.RData files contain results from the null model related to Figure 3
in the main text and for Extended Figure 5, respectively.

sp_tra_20230324_semi.Rdata
sp_tra_20230327_completecas.Rdata
Contains intermediate results from the calculation of functional contribution and
functional originality, for both nearly complete cases (e.g., Figures 1-3) and
complete cases for sensitity analysis (extended Data Figure 7)

hot_nodes.R
R function to perform hot node analysis

Tene_clip.shp
Shapefile of Tenerife area related to Extended Data Figure 4a.