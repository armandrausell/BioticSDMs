
#Finding the overlapping distributions

#Remove UTM column
SpeciesNames<-Species_spain_df[,-1]

names_sps<-names(SpeciesNames) #Extract names
names_sps<-sort(names_sps)

#Remove "." from the names to create clean vector of names
names_sps <- gsub("\\.", " ", names_sps)  # Replace dot with space
print(names_sps)

library(dplyr)
library(NetIndices)

# Filter only interactions where both predators exist in names_sps, in Spain.
interaction_df_spain <- interaction_df %>%
  filter(Predator_A %in% names_sps & Predator_B %in% names_sps)

# Filter only interactions where both PREY exist in names_sps, in Spain.
prey_df_spain <- prey_interaction_df %>%
  filter(Prey_A %in% names_sps & Prey_B %in% names_sps)

# Filter from the total list the species that are in Spain
prey_count_spain<-prey_df %>% filter(Species %in% names_sps)
predator_count_spain<- total_prey_df %>% filter(Predator %in% names_sps)

proportion_prey_predator<-cbind(prey_count_spain,predator_count_spain)
proportion_prey_predator<-proportion_prey_predator[,-3]

#This creates a column with predatory ranking
proportion_prey_predator$TrophicPosition<-
  (proportion_prey_predator$Times_Preyed_On+1)/(proportion_prey_predator$TotalPrey+1)*-1

range01 <- function(x){(x-min(x))/(max(x)-min(x))} # Range standardization 0-1

proportion_prey_predator$TrophicPosition<- range01(proportion_prey_predator$TrophicPosition)

#This will give us an order of species according the trophic position they occupy
proportion_prey_predator$TrophicPosition<-rank(proportion_prey_predator$TrophicPosition)

library(terra)
library(FAMEFMR)
# Create a cache to store raster layers
species_raster_cache <- list()


get.sp.raster <- function(species_name) {
  # Load necessary library
  library(terra)
  
  # Check if species raster is already cached
  if (species_name %in% names(species_raster_cache)) {
    return(species_raster_cache[[species_name]])  # Return cached raster
  }
  
  # Get all species names from the dataset
  existing_names <- names(Species_spain)
  
  # Generate possible name formats "." or " " 
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
  
  # Store the computed raster in cache
  species_raster_cache[[species_name]] <<- spain_raster
  print(paste("Cached raster for:", species_name))
  
  return(spain_raster)
}

#saveSpatRasterList(species_raster_cache,
 #                  filePath = "Species_raster_cache.qs") #Saving cache

if (file.exists("Species_raster_cache.qs")) {
  species_raster_cache <- readSpatRasterList("Species_raster_cache.qs")  # Load if file exists
  print("Loaded species_raster_cache from file.")
} else {
  print("File does not exist. Skipping load.")
}
names_species_cache<-names(species_raster_cache)

# Function to check overlap between two species
check_overlap <- function(species_A, species_B) {
  # Get rasters for both species
  raster_A <- get.sp.raster(species_A) #plot_species is created in AtlasData.R
  raster_B <- get.sp.raster(species_B)
  
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

# Check if the CSV file exists
if (file.exists("Complete_overlap_spain.csv")) {
  
  cat("CSV file found. Loading data...\n")
  interaction_df_spain <- read.csv("Complete_overlap_spain.csv")  # Load the existing data
  
} else {
  
  cat("CSV file not found. Running overlap computation...\n")
  
  # Compute overlap if the file does not exist
  interaction_df_spain[c("Overlap", "Overlap_Proportion")] <- t(apply(
    interaction_df_spain, 1, function(row) {
      result <- check_overlap(row["Predator_A"], row["Predator_B"])
      return(c(result$overlap, result$proportion))
    }
  ))
  
  # Save to CSV
  write.csv(interaction_df_spain, "Complete_overlap_spain.csv", row.names = FALSE)
}

# Remove unlikely interactions (known herbivores, insectivores...)
# Filter out rows where either TotalPrey_A or TotalPrey_B is less than 3
interaction_df_spain_clean <- interaction_df_spain %>%
  filter(!(TotalPrey_A <= 3 | TotalPrey_B <= 3))

# Save if needed
write.csv(interaction_df_spain_clean, "interaction_df_spain_clean.csv", row.names = FALSE)

####### REALIZED COMPETITION ############

# Function to get the names of the preys from binary_matrix_df while filtering names_sps
get_prey_species <- function(predator) {
  prey <- colnames(binary_matrix_df)[binary_matrix_df[predator, ] == 1]
  
  # Keep only prey species that exist in names_sps
  prey <- prey[prey %in% names_sps]
  
  return(prey)
}


# Function to check if a predator and prey co-occur using get.sp.raster()
get_co_occurring_prey <- function(predator, prey_list) {
  co_occurring_prey <- c()
  
  predator_raster <- get.sp.raster(predator)  # Get predator's raster
  
  for (prey in prey_list) {
    prey_raster <- get.sp.raster(prey)  # Get prey's raster

    combined_raster <- predator_raster + prey_raster
    if (sum(values(combined_raster) == 2, na.rm = TRUE) > 0) {
      co_occurring_prey <- c(co_occurring_prey, prey)
    }
  }
  
  return(co_occurring_prey)
}

library(dplyr)
library(purrr)

# Check if the CSV file exists
if (file.exists("Competition_data_spain.csv")) {
  
  cat("CSV file found. Loading data...\n")
  Competition_spain_df <- read.csv("Competition_data_spain.csv")  # Load the existing data
  
} else {
  
  cat("CSV file not found. Running computation...\n")
  
  # Compute Competition_spain_df if the file does not exist
  Competition_spain_df <- interaction_df_spain_clean %>%
    mutate(
      # Get prey species for all Predator_A and Predator_B at once
      Prey_A = map(Predator_A, ~ get_prey_species(.x)),
      Prey_B = map(Predator_B, ~ get_prey_species(.x)),
      
      # Get co-occurring prey for all Predator_A and Predator_B at once
      Co_Occurring_A = map2(Predator_A, Prey_A, ~ get_co_occurring_prey(.x, .y)),
      Co_Occurring_B = map2(Predator_B, Prey_B, ~ get_co_occurring_prey(.x, .y)),
      
      # Count shared co-occurring prey
      Shared_Co_Occurring_Prey = map2_int(Co_Occurring_A, Co_Occurring_B, ~ length(intersect(.x, .y)))
    ) %>%
    dplyr::select(-Prey_A, -Prey_B, -Co_Occurring_A, -Co_Occurring_B)  # Remove temp columns
  
  # Convert back to dataframe
  Competition_spain_df <- as.data.frame(Competition_spain_df)
  
  # Save the computed data to CSV
  write.csv(Competition_spain_df, "Competition_data_spain.csv", row.names = FALSE)
}

Competition_spain_df$RealizedDietOverlapA<-
  Competition_spain_df$Shared_Co_Occurring_Prey/Competition_spain_df$TotalPrey_A

Competition_spain_df$RealizedDietOverlapB<-
  Competition_spain_df$Shared_Co_Occurring_Prey/Competition_spain_df$TotalPrey_B


#Calculate alpha based on species generalism, a parameter to relate Realized and potential overlaps
Competition_spain_df <- Competition_spain_df %>%
  mutate(alphaA = (TotalPrey_A - min(TotalPrey_A, na.rm = TRUE)) /
           (max(TotalPrey_A, na.rm = TRUE) - min(TotalPrey_A, na.rm = TRUE)))

Competition_spain_df <- Competition_spain_df %>%
  mutate(alphaB = (TotalPrey_B - min(TotalPrey_B, na.rm = TRUE)) /
           (max(TotalPrey_B, na.rm = TRUE) - min(TotalPrey_B, na.rm = TRUE)))

#Calculate the amount of pressure competitors exert over the others
Competition_spain_df$PressureBtoA<-Competition_spain_df$alphaA*Competition_spain_df$Proportion_A_Shares+
  (1-Competition_spain_df$alphaA)*Competition_spain_df$RealizedDietOverlapA

Competition_spain_df$PressureAtoB<-Competition_spain_df$alphaB*Competition_spain_df$Proportion_B_Shares+
  (1-Competition_spain_df$alphaB)*Competition_spain_df$RealizedDietOverlapB


