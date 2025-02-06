
library(igraph)
#library(blockmodels)
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

library(data.table)

# Merging species data and the grid to obtain a vector with data inside
Amph_spain <- merge(spain_vector, Amph, by.x = "UTMCODE", by.y = "UTMCODE") # Amphibians
Rept_spain <- merge(spain_vector, Rept, by.x = "UTMCODE", by.y = "UTMCODE") # Reptiles
Bird_spain <- merge(spain_vector, Bird, by.x = "UTMCODE", by.y = "UTMCODE") # Birds
Mamm_spain <- merge(spain_vector, Mamm, by.x = "UTMCODE", by.y = "UTMCODE") # Mammals

# Merging all groups step-by-step using Reduce() to avoid memory issues
Species_spain <- Reduce(function(x, y) merge(x, y, by = "UTMCODE", all = TRUE),
                        list(Amph_spain, Rept_spain, Bird_spain, Mamm_spain))

Species_spain_df <- as.data.frame(Species_spain) # Transform the vector into a df


#Get shapefile for spain
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

plot_species <- function(species_name) {
  # Load necessary library
  library(terra)
  
  # Get all species names from the dataset
  existing_names <- names(Species_spain)
  
  # Generate possible name formats "." or " ". 
  species_to_plot_dot <- gsub(" ", ".", species_name)  # Replace space with dot
  species_to_plot_space <- gsub("\\.", " ", species_name)  # Replace dot with space
  
  # Find the correct species name
  if (species_to_plot_dot %in% existing_names) {
    species_to_plot <- species_to_plot_dot
  } else if (species_to_plot_space %in% existing_names) {
    species_to_plot <- species_to_plot_space
  } else {
    stop(paste("Species not found in the dataset:", species_name))
  }
  
  # Filter the points for the selected species
  sp.points <- Species_spain[Species_spain[[species_to_plot]] == 1, ]
  
  # Check if there are points for this species
  if (nrow(sp.points) == 0) {
    stop(paste("No occurrence data for", species_to_plot))
  }
  
  # Plot base map
  plot(spain_cropped, col = "white", border = "black", axes = FALSE, main = species_to_plot)
  
  # Rasterize the species' occurrence points
  points_raster <- rasterize(sp.points, spain_raster, field=1)
  
  # Assign 1 where points are located
  spain_raster[!is.na(points_raster)] <- 1
  
  # Plot raster with species occurrences highlighted
  plot(spain_raster, col=c("gray", "red"), main=paste("Species:", species_to_plot))
  
  # Add occurrence points to the plot
  #plot(sp.points, col = "blue", pch = 20, add = TRUE)
}
plot_species("Lissotriton boscai")
