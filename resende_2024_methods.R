# This R script was developed by Ludmila S. Resende, Rejane Silva-Souza e Luísa
# Carvalheiro with the aim to estimate the contribution of pollinators to
# agricultural production. The procedures were applied here using a map of
# land use and cover and orange production values for the year 2022 for the
# municipality of Itaberaí (Brazil), a list of bees from the Cerrado biome
# that pollinate oranges (Citrus sinensis L.) and information on nesting
# suitability of different habitats around crop fields for these pollinators.

# For more details, see Resende, L. S. Integrating pollinator attributes into
# natural capital estimates for agricultural production. Master dissertation.
# Universidade Federal de Goiás, 2024. Link: <...>

# The procedures were divided into four steps:

# 1. Reclassify target crop plantings
# 2. Resize raster resolution
# 3. Spatially model pollinator density
# 4. Estimate the contribution of pollinators to orange production

# Warning: step 1.1 is not necessary if you have a raster with correct classification
# of target crop area, but the step 1.2 must be done.
# For detailed information on these steps, see Resende 2024.


# Required packages -------------------------------------------------------
library(terra)
library(tidyverse)
library(sf)

# 1. Reclassify target crop plantings -------------------------------------
# 1.1 Delimitation of the cultivation area of the agricultural species of interest

landscape <- terra::rast("lulc_map_itaberai_2022_utm.tif") # LULC map of the study area
orange_fields <- terra::vect("orange_fields_utm_files/area_laranja_utm.shp") # shapefile of the target crop fields
orange_fields_raster <- rasterize(orange_fields, landscape, touches = TRUE) # transform shapefile in raster
landscape_rec <- merge(orange_fields_raster, landscape, first = TRUE) # join the orange fields with the LULC map
landscape_rec[landscape_rec == 1] <- 47 # assigning the specific Mapbiomas classification number to the Citrus class
print(landscape_rec) # view raster information
plot(landscape_rec) # visualize the raster

# 1.2 Reclassify raster to simplify land use and land cover classes

freq_raster <- terra::freq(landscape_rec) # obtain frequency of each lulc class
lulc_classes_before <- freq_raster$value # obtain lulc classes code
print(lulc_classes_before) # see sequence to replace some values according to following:

# 0 == NA; 
# 11 (Wetland), 25 (Other non Vegetated Areas) and 33 (River, Lake and Ocean) == 0
# 19 (Temporary Crop), 21 (Mosaic of Uses) and 47 (Citrus) == 18 (Agriculture)

lulc_classes_after <- c(NA,  3,  4,  9, 0, 12, 15, 18, 18, 24, 0, 0, 18) # new class values
classes <- cbind(lulc_classes_before,
                 lulc_classes_after) # reclassification matrix
landscape_rec_2 <- terra::classify(landscape_rec,
                                   classes) # obtain new raster with eight class values

# 2. Resize raster resolution ------------------------------------------------
# 2.1 Create a place to store the suitability maps that will be generated 

spp_pollin <- readr::read_csv("spp_pollin.csv") # habitat suitability values for pollinator nesting
spp_pollin <- spp_pollin %>%
  mutate(spp_suitab = gsub(" ", "_", spp_name) %>%
           paste("_suitab", sep = "")) # create names for the list
spp_suitab <- unique(spp_pollin$spp_suitab) # identify unique names
spp_suitab_list <- vector(mode = "list",
                            length = length(spp_suitab)) # create the list
names(spp_suitab_list) <- spp_suitab # assign names to the list

# 2.2 Create nesting suitability map for each species
rm(i)
for(i in 1 : length(unique(spp_pollin$spp_suitab))){
  spp_suitab_i <- unique(spp_pollin$spp_suitab)[i]
  
  # filter suitability values for species i by excluding class values that do not exist in landscape_rec_2
  freq_raster_rec <- terra::freq(landscape_rec_2)
  tab_spp <- spp_pollin %>%
    filter(spp_suitab == spp_suitab_i & (class_num %in% freq_raster_rec$value)) %>%
    arrange(class_num)
  
  # create reclassification matrix
  lulc_classes_raster <- freq_raster_rec$value[-1] # select all classes except class 0
  lulc_classes_suit <- tab_spp$suitability
  classes_rec <- cbind(lulc_classes_raster,
                       lulc_classes_suit)
  
  # create suitability map and store in spp_suitab_list
  spp_suitab_list[[i]] <- terra::classify(landscape_rec_2,
                                            classes_rec)
  
  print(i)
  
};rm(i)

# 2.3 Change map resolution from 10m to 100m
spp_suitab_list_100m <- spp_suitab_list # create place to store maps with 100m resolution
rm(b)
for(b in 1 : length(spp_suitab_list)){
  spp_suitab_list_100m[[b]] <- terra::aggregate(spp_suitab_list[[b]],
                                                  fun = "mean",
                                                  fact = 10)

  names(spp_suitab_list_100m[[b]]) <- unique(spp_pollin$spp_suitab)[b]
  
  print(b)
  
};rm(b)


# 3. Spatially model pollinator density -----------------------------------

# 3.1 Estimate foraging distance
dist_forage <- read_csv("dist_forage.csv")
ITD <- dist_forage$ITD_mm
dist_forage$distances_km <- pollimetry::foragedist(ITD, "GrMfd")
dist_forage$distances_m <- dist_forage$distances_km * 1000

# set the maximum foraging distance at 5000 meters
for (i in 1:length(dist_forage$distances_m)) {
  if (dist_forage$distances_m[i] >= 5000) {
    dist_forage$distances_m[i] <- 5000
  } else {
    dist_forage$distances_m[i] <- dist_forage$distances_m[i]
  }
}


# 3.2 Obtain coordinate information of plantation pixels
coord <- read.csv("coord_pixels_100m.csv") #Insert the attribute table of cropland fields shapefile containing coordinates of each pixel
coord_pixels <- coord[ ,c(2,3)] # select left and top
points_list <- list()
for (i in 1:nrow(coord_pixels)) {
  x_coord <- as.numeric(coord_pixels[i, "left"])
  y_coord <- as.numeric(coord_pixels[i, "top"])
  point <- st_point(cbind(x_coord, y_coord))
  points_list[[i]] <- point
}

# 3.3 Create a table to store the densities
densities <- data.frame(coord_pixels)
names <- gsub("adequab", "densid", spp_suitab)
num_cols <- length(names)
for (i in 1:num_cols) {
  densities[[names[i]]] <- NA
}

# 3.4 Run the model
# setting a progress bar
pb <- txtProgressBar(min = 0, max = length(spp_suitab_list_100m) * length(points_list), style = 3)
cont <- 0

# selecting the species i
for (i in 1:length(spp_suitab_list_100m)){
  raster_spp <- spp_suitab_list_100m[[i]] # selecting suitability raster for species i
  dist_for_spp <- dist_forage$distances_m[i] # selecting forage distance for species i
  
  # selecting the cropland pixel j
  for(j in 1:length(points_list)){
    pixel_focal <- points_list[[j]]
    
    # update the progress bar
    cont <- cont + 1
    setTxtProgressBar(pb, cont)
    
    # create a buffer around agricultural pixel j with a radius equal to the foraging distance of species i
    buffer <- sf::st_buffer(pixel_focal, dist = dist_for_spp) # draw a buffer
    buffer_coords <- sf::st_coordinates(buffer) # get buffer boundary coordinates
    buffer_coords <- buffer_coords[ ,c(1,2)] # keep only columns x and y
    poligono <- terra::vect(buffer_coords, type = "polygon") # transform buffer to polygon
    terra::crs(poligono) <- "EPSG:32722" # define the polygon CRS
    raster_spp_masked <- mask(raster_spp, poligono) # crop the raster by the extent of the mask
    raster2points <- terra::as.points(raster_spp_masked, values = TRUE) # transform all pixels inside the buffer into point vector
    all_pixels_buffer <- terra::extract(raster_spp, raster2points, cells = TRUE) # get positions of points in the raster
    celulas <- all_pixels_buffer$cell # create vector with positions
    coord_all_pixels_buffer <- terra::xyFromCell(raster_spp, celulas) # get coordinates from positions
    dados_buffer <- data.frame(all_pixels_buffer, coord_all_pixels_buffer) # gather information on land use and land cover classes, position and coordinates of agricultural pixels

    # calculate the distances between pixel i and each of the nesting pixels
    dist_matrix <- raster::pointDistance(pixel_focal, coord_all_pixels_buffer, lonlat = FALSE)
    dados_buffer$distances <- dist_matrix
    
    # calculate the pollinator density in the pixel i
    Dens_spp <- dados_buffer[,2]*exp(-(dados_buffer$distances+1)/dist_for_spp)
    densities[j, i+2] <- sum(Dens_spp, na.rm = TRUE)
  }
}

close(pb) # close the progress bar

densities$densTotal <- rowSums(densities[ , c(3:length(densities))]) # calculate the total pollinator density for each pixel

# optional: add position information to open table as map in QGis
densities <- data.frame(coord[ ,1], densities)

# 4. Estimate the contribution of pollinators to orange production --------

# 4.1 Calculate the fraction of optimal pollinator density
opt_dens <- max(densities$densTotal)
densities$dens_fraction <- densities$densTotal/opt_dens

# set the production parameters
QPT = 27996 * 1000 # total quantity produced (kg) in 2022
VPT = 31780 * 1000 # total production value (reais) in 2022

#set the dependence rates
PDmin=0.06
PDavg=0.19
PDmax=0.31

# prepare table to save results
results <- data.frame(
  Parameters = numeric(3),
  Total_production = numeric(3),
  Pollin_contr_PDmin_0.06 = numeric(3),
  Pollin_contr_PDavg_0.19 = numeric(3),
  Pollin_contr_PD_max_0.31 = numeric(3)
)
results[c(1:3),1] <- c("Quant. prod. (kg)",
                          "Prod. value (reais)",
                          "Polin. contr. (%)")
results[c(1,2),2] <- c(QPT, VPT)
results[3,2] <- "-"

# set the function to calculate pollinator contribution
calc_contrib <- function (parameter, PD) {
  V <-  parameter/ (nrow(densities)+PD*sum(densities$dens_fraction))
  without_pollin <- V * (nrow(densities))
  pollin_contrib <- parameter - without_pollin
  return(pollin_contrib)
}

# calculate the pollinator contribution to crop production

# 1. absolute values
results[1,c(3:5)] <- calc_contrib(QPT, c(PDmin, PDavg, PDmax))
results[2,c(3:5)] <- calc_contrib(VPT, c(PDmin, PDavg, PDmax))

# 2. percentage values
calc_contrib_perc <- function(parameter, contrib) {
  resu <- round((contrib*100/parameter),2)
  return(resu)
}
# only one production parameter is necessary (here, QPT)
results[3, c(3:5)] <- calc_contrib_perc(QPT, c(results[1,3], 
                                                     results[1,4],
                                                     results[1,5]))

write.csv(results, file="pollin_contrib_13-02-24.csv")


# estimate the gain in production in case of maximum level of pollination
# and loss in the absence of pollination

# function to calculate the production without pollinator
Prod_without_pollin <- function (parameter, PD){
  V <-  parameter/ (nrow(densities)+PD*sum(densities$dens_fraction))
  resu <- V * nrow(densities)
  return(resu)
}

# create data.frame to store the results
results2 <- data.frame(
  Cenary = rep(c("Without", "Maximum"), each = 2),
  Parameter = rep(c("Quant. prod.", "Prod. value"), times = 2),
  Total_values_PDmin_0.06 = numeric(4),
  Total_values_PDavg_0.19 = numeric(4),
  Total_values_PD_max_0.31 = numeric(4),
  stringsAsFactors = FALSE
)

# calculate
results2[1,c(3:5)] <- Prod_without_pollin(QPT, c(PDmin, PDavg, PDmax))
results2[2,c(3:5)] <- Prod_without_pollin(VPT, c(PDmin, PDavg, PDmax))


# calculate the production with maximum level of pollination

# obtain the complement of dependency rates on pollinators
complement <- 1 - c(PDmin, PDavg, PDmax)

results2[3, c(3:5)] <- results2[1,c(3:5)] / complement
results2[4, c(3:5)] <- results2[2,c(3:5)] / complement

# add information about gains and losses in percentage values

diff_without_pollin <- round((results2[1,c(3:5)]*100/QPT),2) - 100
diff_max_pollin <- round((results2[3, c(3:5)]*100/QPT),2) - 100

results2[5,c(1:5)] <- c("Loss without pollin. (%)", "Quant. and Prod. value", diff_without_pollin)
results2[6,c(1:5)] <- c("Maximum gain (%)", "Quant. and Prod. value", diff_max_pollin)

write.csv(results2, file="loss_and_gain_13-02-24.csv")

###

