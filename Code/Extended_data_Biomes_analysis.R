# M. Paola Barajas Barbosa et al (2022) Assembly of functional diversity in an oceanic island flora

# # Extended data Fig. 6
# Tenerife map and biome representativity on the island. 

# https://rawgit.com/valentinitnelav/plotbiomes/master/html/Whittaker_biomes_examples.html
# devtools::install_github("valentinitnelav/plotbiomes")

library(plotbiomes)
library(ggplot2)
library(raster)
library(maptools)
library(sf)
library(rgdal)

whittaker_base_plot() + theme_bw()

# temp_c, precp_cm,
summary(Whittaker_biomes)
names(Whittaker_biomes)
plot(density(Whittaker_biomes$temp_c))
plot(density(Whittaker_biomes$precp_cm))

# prec and temp from: https://datadryad.org/stash/dataset/doi:10.5061/dryad.qrfj6q5fs
prec = raster("Temp_Prec_Canaries_FernandezPalaciosetal_2020/P_ann.tif")
temp = raster("Temp_Prec_Canaries_FernandezPalaciosetal_2020/T_ann.tif")

tmp  = raster::stack(temp, prec)
names(tmp) <- c("temperature", "precipitation")
temp_pp = tmp

# ===== Generate random locations within Tenerife
wrld_simpl <- st_read("Tene_clip.shp") ; wrld_simpl <- as_Spatial(wrld_simpl)

plot(prec); plot(wrld_simpl, add = TRUE)

set.seed(66) # Create random locations within world's polygons.
points <- sp::spsample(x = wrld_simpl, n = 1000, type = "random")

setwd("C:/2022_Assemply_ms_PhDII/pub_R_code/submission_2")
writeOGR(obj=points, dsn="tempdir", layer="points", driver="ESRI Shapefile")

plot(prec); plot(wrld_simpl, add = TRUE); plot(points, add = TRUE)

# ===== Extract from raster stack: Extract temperature and precipitation values from the raster datasets
extractions <- raster::extract(temp_pp, points, df = TRUE)

# Adjust temperature values to normal scale because WorldClim temperature data
# has a scale factor of 10 (integer storage for saving space).
extractions$temperature <- extractions$temperature
extractions$precipitation <- extractions$precipitation/10 # Convert precipitation from mm to cm

my_theme <- theme( panel.background =element_rect(fill="grey97",colour="black"),
                   legend.text = element_text(size=7), axis.title= element_text(size=7),
                   axis.text.x= element_text(colour= "black", size=7), axis.text.y= element_text(colour= "black", size=7))

wp = whittaker_base_plot() +
  # add the temperature - precipitation data points
  geom_point(data = extractions, 
             aes(x = temperature, 
                 y = precipitation), 
             size   = .2,
             shape  = 21,
             colour = "black", 
             fill   = "black",
             stroke = .51,
             alpha  = 0.5) +
  my_theme

wp

pdf(file = "whittaker_plot_Tenerife.pdf", 
    bg = "transparent",
    width=5,   height=3,
    paper = "a4")
wp
dev.off()

# FIN