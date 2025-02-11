
VisOverlap <- function(species_A, species_B) {
  # Get rasters for both species
  raster_A <- plot_species(species_A) #plot_species is created in AtlasData.R
  raster_B <- plot_species(species_B)
  
  # Sum both rasters (cells with 2 indicate overlap)
  combined_raster <- raster_A + raster_B
  plot(combined_raster)
}

VisOverlap("Ursus arctos", "Canis lupus")

