library(terra)

# Function to compute species richness map
compute_species_richness <- function(species_list) {
  # Initialize an empty raster with the same extent/resolution as the first species raster
  richness_raster <- spain_raster
  values(richness_raster) <- 0  # Start with all values as 0
  
  # Loop through each species and add its presence to the richness raster
  for (species in species_list) {
    cat("Adding", species, "to species richness map...\n")
    
    # Get species raster
    species_raster <- tryCatch(
      get.sp.raster(species),
      error = function(e) {
        cat("Skipping", species, "- no valid raster found.\n")
        return(NULL)
      }
    )
    
    # If the raster exists, add it to the cumulative species richness map
    if (!is.null(species_raster)) {
      richness_raster <- richness_raster + species_raster
    }
  }
  
  # Plot the species richness map
  plot(richness_raster, main = "Species Richness Map", col = terrain.colors(100))
  
  return(richness_raster)
}

# Example: Compute species richness for all unique species in the dataset
species_richness_map <- compute_species_richness(names_sps)
