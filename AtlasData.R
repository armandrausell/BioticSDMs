
library(igraph)
#library(blockmodels)
library(NetIndices)
library(ggplot2)
library(RColorBrewer)
library(readr)
library(tidyr)
library(magrittr)
library(dplyr)
#library(deSolve)
library(geodata)
library(sdm)
library(usdm)
library(terra)
library(raster)
library(rasterVis)
library(sf)
#library(sbm)
#library(networkD3)

# Loading species data from ATLAS data
Amph <- read_delim("Atlas_data/DB_Amphibians_IP.txt",
                   delim = "\t", escape_double = FALSE, 
                   trim_ws = TRUE,
                   show_col_types = FALSE)
Bird <- read_delim("Atlas_data/DB_Birds_IP.txt", 
                   delim = "\t", escape_double = FALSE, 
                   trim_ws = TRUE,
                   show_col_types = FALSE)
Mamm <- read_delim("Atlas_data/DB_Mammals_IP.txt", 
                   delim = "\t", escape_double = FALSE, 
                   trim_ws = TRUE,
                   show_col_types = FALSE)
Rept <- read_delim("Atlas_data/DB_Reptiles_IP.txt", 
                   delim = "\t", escape_double = FALSE, 
                   trim_ws = TRUE,
                   show_col_types = FALSE)

# Loading ATLAS's shapefile
spain_vector <- vect("Atlas_data/UTMgrid_IberianPeninsula/UTM_GRID.shp")
spain_vector_latlong <- project(spain_vector, "+proj=longlat +datum=WGS84")

# Merging species data and the grid to obtain a vector with data inside
Amph_spain <- merge(spain_vector, Amph, 
                    by.x = "UTMCODE", by.y = "UTMCODE") # Amphibians
Rept_spain <- merge(spain_vector, Rept, 
                    by.x = "UTMCODE", by.y = "UTMCODE") # Reptiles
Bird_spain <- merge(spain_vector, Bird, 
                    by.x = "UTMCODE", by.y = "UTMCODE") # Birds
Mamm_spain <- merge(spain_vector, Mamm, 
                    by.x = "UTMCODE", by.y = "UTMCODE") # Mammals
Species_spain <- merge(Amph_spain, Rept_spain, 
                       by.x = "UTMCODE", by.y = "UTMCODE") %>%
  merge(Bird_spain, 
        by.x = "UTMCODE", by.y = "UTMCODE") %>%
  merge(Mamm_spain, 
        by.x = "UTMCODE", by.y = "UTMCODE")

Species_spain_df <- as.data.frame(Species_spain) # Transform the vector into a df

plot(Species_spain) # Chech the UTM-based grid

spain <- geodata::world(path = "countries.shp") 
spain <- spain[spain$GID_0 %in% c("ESP","PRT"), ]
spain <- project(spain, crs(Species_spain))
e1 <- ext(-199683.152327576, 1116052.05195627, 
          3798402.70715166, 4950541.20614095)
spain_cropped <- terra::crop(spain, e1)


#Transform the vector to a raster
# Define the raster template (resolution and extent should be appropriate)
spain_raster <- rast(spain_vector, resolution = 10000)  # Adjust resolution as needed

# Fill the raster with 0
values(spain_raster) <- 0

# Plot the output
plot(spain_raster)
plot(spain_cropped, col = "white", border = "black", axes = F, add=T)


# Create a raster filled with NA
values(spain_raster) <- NA

# Rasterize the vector with a buffer (sets cells touching the vector to 0)
buffered_raster <- rasterize(spain_vector, spain_raster, field=1)  # Temporary raster with 1s

# Assign 0 to cells where the vector is present, keep others as NA
spain_raster[!is.na(buffered_raster)] <- 0

# Plot the output to verify
plot(spain_raster, col=c("gray"), main="Spain Raster")
plot(spain_cropped, add=TRUE, border="black", lwd=2)  # Overlay vector

# Save the modified raster
writeRaster(spain_raster, "spain_rasterUTM.tif", overwrite=TRUE)

############ PLOT POINTS TO RASTER ##############
species_to_plot <- names(Species_spain[,2]) # Choose from 2 to 293
plot(spain_cropped, col = "white", border = "black", axes = F)
sp.points <- Species_spain[Species_spain[[species_to_plot]] == 1, ]
#plot(sp.points, col = "blue", add=T)

# Rasterize the points (assign 1 to cells where points exist)
points_raster <- rasterize(sp.points, spain_raster, field=1)

# Assign 1 where points are located
spain_raster[!is.na(points_raster)] <- 1

# Plot the output to verify
plot(spain_raster, col=c("gray", "red"), main="Vector Points Set to 1")
#plot(sp.points, add=TRUE, col="blue", pch=20)  # Overlay points


