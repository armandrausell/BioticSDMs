
#Finding the overlapping distributions

#Remove UTM column
SpeciesNames<-Species_spain_df[,-1]

names_sps<-names(SpeciesNames) #Extract names

#Remove "." from the names to create clean vector of names
names_sps <- gsub("\\.", " ", names_sps)  # Replace dot with space
print(names_sps)

library(dplyr)

# Filter only interactions where both predators exist in names_sps, in Spain.
interaction_df_spain <- interaction_df %>%
  filter(Predator_A %in% names_sps & Predator_B %in% names_sps)

library(terra)
library(dplyr)

#Define function to plot raster of each species when called
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
  
  # Rasterize the species' occurrence points
  points_raster <- rasterize(sp.points, spain_raster, field=1)
  
  # Assign 1 where points are located
  spain_raster[!is.na(points_raster)] <- 1
  return(spain_raster)
  }


# Function to check overlap
check_overlap <- function(species_A, species_B) {
  # Get rasters for both species
  raster_A <- plot_species(species_A) #plot_species is created in AtlasData.R
  raster_B <- plot_species(species_B)
  
  # Sum both rasters (cells with 2 indicate overlap)
  combined_raster <- raster_A + raster_B
  
  # Count overlapping cells (value = 2)
  overlap_cells <- sum(values(combined_raster) == 2, na.rm = TRUE)
  
  # Count total occupied cells (value = 1 or 2)
  total_cells <- sum(values(raster_A) >= 1, na.rm = TRUE)
  
  # Compute proportion of overlap
  overlap_proportion <- ifelse(total_cells > 0, overlap_cells / total_cells, 0)
  
  # Return TRUE if overlap exists, FALSE otherwise
  return(list(overlap = overlap_cells > 0, proportion = overlap_proportion))
}

interaction_df_spain <- as.data.frame(interaction_df_spain)

interaction_df_spain[c("Overlap", "Overlap_Proportion")] <- t(apply(
  interaction_df_spain, 1, function(row) {
    result <- check_overlap(row["Predator_A"], row["Predator_B"])
    return(c(result$overlap, result$proportion))
  }
))

# Save to CSV if needed
write.csv(interaction_df_spain, "Complete_overlap_spain.csv", row.names = FALSE)

