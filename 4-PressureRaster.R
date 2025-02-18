
library(terra)
library(dplyr)
library(psych)

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
       col = rev(terrain.colors(100)))
  
  return(total_pressure_raster)
}

# Example: Compute and plot the competitive pressure on any species
sp<-"Coronella girondica" #Species name
pressure_raster <- compute_competitor_pressure(paste0(sp), n = 5)
r<-get.sp.raster(sp)
plot(r,main = paste("Distribution of", sp))

#############################

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

library(dplyr)
library(terra)

# Function to assess correlation for a species while varying competitor numbers
evaluate_competitor_correlation <- function(species_name, max_competitors = 15) {
  
  # Get species richness map (same for all competitor numbers)
  richness_map <- compute_species_richness(unique_species)
  
  # Get actual species distribution
  actual_distribution <- get.sp.raster(species_name)
  
  # Initialize results dataframe
  correlation_results <- data.frame(
    Species = character(),
    Num_Competitors = integer(),
    Correlation_Actual = numeric(),
    Correlation_Richness = numeric(),
    Corr_Act_Biserial = numeric()
  )
  
  # Loop through 1 to max_competitors
  for (n in 1:max_competitors) {
    cat("Computing for", species_name, "with", n, "competitors...\n")
    
    # Compute competition pressure map with n competitors
    pressure_map <- compute_competitor_pressure(species_name, n)
    
    # Extract values while removing NA
    pressure_values <- values(pressure_map)
    actual_values <- values(actual_distribution)
    richness_values <- values(richness_map)
    
    # Remove NA values from all three datasets
    valid_indices <- !is.na(pressure_values) & !is.na(actual_values) & !is.na(richness_values)
    pressure_values <- pressure_values[valid_indices]
    actual_values <- actual_values[valid_indices]
    richness_values <- richness_values[valid_indices]
    
    # Compute Spearman correlation with actual species distribution
    correlation_actual <- cor(pressure_values, actual_values, method = "spearman")
    
    # Compute Spearman correlation with species richness map
    correlation_richness <- cor(pressure_values, richness_values, method = "spearman")
    
    # Compute biserial correlation, similar to Pearson but for binary/continuous comparisons
    r_pb <- biserial(pressure_values, actual_values)
    
    # Store results
    correlation_results <- rbind(correlation_results, data.frame(
      Species = species_name,
      Num_Competitors = n,
      Correlation_Actual = correlation_actual,
      Correlation_Richness = correlation_richness,
      Corr_Act_Biserial = r_pb
    ))
  }
  
  return(correlation_results)
}



# Check if the CSV file exists
if (file.exists("Correlation_distribution_w_competition.csv")) {
  cat("CSV file found. Loading data...\n")
  all_correlation_results <- read.csv("Correlation_distribution_w_competition.csv")  # Load the existing data
} else {
  cat("CSV file not found. Running evaluation...\n")
  # Compute correlation results if the file does not exist
  all_correlation_results <- data.frame()
  for (species in unique_species) {
    cat("Evaluating:", species, "\n")
    species_results <- evaluate_competitor_correlation(species, max_competitors = 15)
    all_correlation_results <- rbind(all_correlation_results, species_results)
  }
  
  # Save the computed data to CSV
  write.csv(all_correlation_results, "Corr_distribution_w_competition.csv", row.names = FALSE)
}

all_correlation_results$difference<-all_correlation_results$Correlation_Actual-all_correlation_results$Correlation_Richness

# Select the best number of competitor count for each species
Opt_comptetitor_per_species <- all_correlation_results %>%
  group_by(Species) %>%
  filter(abs(Corr_Act_Biserial) == max(abs(Corr_Act_Biserial))) %>%
  mutate(Corr_Act_Biserial = abs(Corr_Act_Biserial) * sign(Corr_Act_Biserial)) %>%  # Preserve sign
  dplyr::select(Species, Num_Competitors, Corr_Act_Biserial, Correlation_Richness) %>%
  arrange(desc(Corr_Act_Biserial))

library(ggplot2)
library(ggrepel)  # For better label placement

# Ensure both size and color have the same scale and breaks
ggplot(Opt_comptetitor_per_species, aes(x = Corr_Act_Biserial, y = Correlation_Richness, label = Species)) +
  geom_point(alpha = 0.8, aes(size = Num_Competitors, color = Num_Competitors)) +  # Maintain point visibility
  scale_color_continuous(limits=c(1, 15), breaks=seq(1, 15, by=2),low = "#f7b9b5", high = "darkred") +
  scale_size_continuous(name = "Number of Competitors", breaks = seq(1, 15, by = 2), limits = c(1, 15)) +  
  #scale_color_viridis_c(name = "Number of Competitors", option = "plasma") +  # Using Viridis for better color mapping
  geom_vline(xintercept = 0, linetype="dashed")+
  geom_text_repel(aes(label = Species), size = 2, color = "black") +  # Keep text readable
  labs(
    title = "Relationship Between competitive pressure with Species Richness and actual dist",
    x = "Correlation with Actual Distribution",
    y = "Correlation with Species Richness"
  ) +
  guides(
    size = guide_legend(title = "Number of Competitors"),
    color = guide_legend(title = "Number of Competitors")  # Force merging
  ) +
  theme_minimal() +  
  theme(legend.position = "right")
