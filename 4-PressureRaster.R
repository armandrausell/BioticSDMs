
library(dplyr)

library(dplyr)

# Define function to get top competitors for any species
get_top_competitors <- function(species_name, n = 5)  { #amount of species to get
  top_pressure_species <- Competition_spain_df %>%
    filter(Predator_A == species_name | Predator_B == species_name) %>%
    mutate(
      # Select the correct pressure column based on position
      Pressure_on_Species = ifelse(Predator_A == species_name, PressureBtoA, PressureAtoB),
      
      # Identify the species exerting pressure on the target species
      Pressure_Exerting_Species = ifelse(Predator_A == species_name, Predator_B, Predator_A)
    ) %>%
    arrange(desc(Pressure_on_Species)) %>%  # Sort by highest pressure
    slice(1:n)  # Select the top `n` species
  
  # Extract unique competitor species
  top_competitors <- unique(top_pressure_species$Pressure_Exerting_Species)
  print(top_competitors)
  return(list(Top_Pressure_Table = top_pressure_species, Top_Competitors = top_competitors))
}

# Example usage for Lynx pardinus
result <- get_top_competitors("Lynx pardinus", n = 5)

library(terra)
library(dplyr)

# Function to compute the total competitive pressure raster for a given species
compute_competitor_pressure <- function(species_name, n = 10) {
  # Step 1: Get top competitors and their pressures
  competitors_data <- get_top_competitors(species_name, n)
  top_pressure_table <- competitors_data$Top_Pressure_Table
  
  # Step 2: Initialize an empty raster (same extent/resolution as Spain raster)
  total_pressure_raster <- spain_raster
  values(total_pressure_raster) <- 0  # Start with all values as 0
  
  # Step 3: Iterate over competitors, get their raster, and weight it
  for (i in 1:nrow(top_pressure_table)) {
    competitor <- top_pressure_table$Pressure_Exerting_Species[i]
    pressure_value <- top_pressure_table$Pressure_on_Species[i]
    
    # Get competitor's presence raster
    competitor_raster <- get.sp.raster(competitor)
    
    # Multiply raster by pressure value
    weighted_raster <- competitor_raster * pressure_value
    
    # Add weighted raster to the total pressure raster
    total_pressure_raster <- total_pressure_raster + weighted_raster
  }
  
  # Step 4: Plot the final raster
  plot(total_pressure_raster, main = paste("Total Competitive Pressure on", species_name),
       col = terrain.colors(100))
  
  return(total_pressure_raster)
}

# Example: Compute and plot the competitive pressure on Lynx pardinus
sp<-"Coronella girondica"
pressure_raster <- compute_competitor_pressure(paste0(sp), n = 5)
r<-get.sp.raster(sp)
plot(r,main = paste("Distribution of", sp))

#############################
library(terra)
library(dplyr)

# Function to assess Spearman correlation between competitive pressure and actual distribution
assess_correlation <- function(species_name, n = 5) {
  # Step 1: Compute the competitive pressure raster
  pressure_raster <- compute_competitor_pressure(species_name, n)
  
  # Step 2: Get the actual species distribution raster
  actual_raster <- get.sp.raster(species_name)
  
  # Step 3: Extract values from both rasters where species is present
  pressure_values <- values(pressure_raster)
  actual_values <- values(actual_raster)
  
  # Remove NA values from both datasets (only keep matched data)
  valid_indices <- !is.na(pressure_values) & !is.na(actual_values)
  pressure_values <- pressure_values[valid_indices]
  actual_values <- actual_values[valid_indices]
  
  # Step 4: Compute Spearman correlation
  correlation_result <- cor(pressure_values, actual_values, method = "spearman")
  
  # Print and return correlation
  print(paste("Spearman correlation for", species_name, ":", correlation_result))
  
  return(correlation_result)
}

# Example: Compute correlation for Lynx pardinus
Sp.correlation <- assess_correlation("Lynx pardinus", n = 5)

# Get all the unique names
unique_species <- unique(c(Competition_spain_df$Predator_A, Competition_spain_df$Predator_B))

# Initialize an empty dataframe to store results
Sp.correlation <- data.frame(Species = character(), Spearman_Correlation = numeric())

# Loop over selected species and compute correlation
for (species in unique_species) {
  
  # Check if species exists in the dataset
  if (!(species %in% Competition_spain_df$Predator_A | species %in% Competition_spain_df$Predator_B)) {
    cat("Skipping species (not found in dataset):", species, "\n")
    next  # Skip this species and move to the next iteration
  }
  
  cat("Computing correlation for:", species, "\n")  # Print progress
  
  # Try to compute correlation, handle errors
  correlation_value <- tryCatch(
    assess_correlation(species, n = 5),
    error = function(e) {
      cat("Error computing correlation for", species, ":", e$message, "\n")
      return(NA)  # Return NA if there is an error
    }
  )
  
  # Store results in dataframe only if successful
  if (!is.na(correlation_value)) {
    Sp.correlation <- rbind(Sp.correlation, data.frame(Species = species, Spearman_Correlation = correlation_value))
  }
}

# Print final correlation table
print(Sp.correlation)

# Compute R² and round to 2 decimal places
Sp.correlation$R2 <- round(Sp.correlation$Spearman_Correlation^2, 3)


