
library(terra)
library(dplyr)
library(psych)
#Competition_spain_df<-Competition_weighted

# Define function to get top competitors for any species
get_top_competitors <- function(species_name, n = 5, sizeWeight = TRUE, potentialDiet = FALSE)  { 
  # Check if size weighting should be applied
  if(sizeWeight == TRUE) {
    
    # Join bodymass_df to assign mass to Predator_A
    Competition_weighted <- Competition_spain_df %>%
      left_join(bm_combined, by = c("Predator_A" = "Species")) %>%
      rename(BodyMass_A = BodyMass_g, Class_A = Class)
    
    # Join again for Predator_B
    Competition_weighted <- Competition_weighted %>%
      left_join(bm_combined, by = c("Predator_B" = "Species")) %>%
      rename(BodyMass_B = BodyMass_g, Class_B = Class)
    
    # Convert new columns to numeric (if not already)
    Competition_weighted <- Competition_weighted %>%
      mutate(BodyMass_A = as.numeric(BodyMass_A),
             BodyMass_B = as.numeric(BodyMass_B))
    
    # Define Gaussian function
    gaussian_similarity <- function(mass_A, mass_B, sigma = 1.8) {
      exp(-((log(mass_A) - log(mass_B))^2 / (2 * sigma^2)))  # Gaussian decay
    }
    
    if(potentialDiet == TRUE){ #Only use the potential diet
      # Define Gaussian function
      gaussian_similarity <- function(mass_A, mass_B, sigma = 1.8) {
        exp(-((log(mass_A) - log(mass_B))^2 / (2 * sigma^2)))  # Gaussian decay
      }
      
      # Apply Gaussian weighting to favor similar-sized competitors
      Competition_weighted <- Competition_weighted %>%
        mutate(
          MassSimilarity = gaussian_similarity(BodyMass_A, BodyMass_B, sigma = 1.8),  # Adjust sigma
          TPSimilarity = gaussian_similarity(TPA, TPB, sigma = 1), # Calculate Trophic similarity
          PressureBtoA = round(
            (Proportion_A_Shares * MassSimilarity)*TPSimilarity, 3
          ),
          PressureAtoB = round(
            (Proportion_B_Shares* MassSimilarity)*TPSimilarity, 3
          )
        )
    } else {
    
    # Apply Gaussian weighting to favor similar-sized competitors
    Competition_weighted <- Competition_weighted %>%
      mutate(
        MassSimilarity = gaussian_similarity(BodyMass_A, BodyMass_B, sigma = 1.8),  # Adjust sigma
        TPSimilarity = gaussian_similarity(TPA, TPB, sigma = 1), # Calculate Trophic similarity
        PressureBtoA = round(
          ((RealizedDietOverlapA) * MassSimilarity)*TPSimilarity, 3
        ),
        PressureAtoB = round(
          ((RealizedDietOverlapB) * MassSimilarity)*TPSimilarity, 3
        )
        #PressureBtoA = round(
         # ((alphaA * Proportion_A_Shares + (1 - alphaA) * RealizedDietOverlapA) * MassSimilarity)*TPSimilarity, 3
        #),
        #PressureAtoB = round(
        #  ((alphaB * Proportion_B_Shares + (1 - alphaB) * RealizedDietOverlapB) * MassSimilarity)*TPSimilarity, 3
        
      )
    }
  } else {
    
    if(potentialDiet == TRUE){
      # Calculate pressure for the standard method
      Competition_weighted <- Competition_spain_df %>%
        mutate(
          PressureBtoA = Proportion_A_Shares*TPSimilarity,
          PressureAtoB = Proportion_B_Shares*TPSimilarity)
    } else {
    # Calculate pressure for the standard method
    Competition_weighted <- Competition_spain_df %>%
      mutate(
        #PressureBtoA = alphaA * Proportion_A_Shares + (1 - alphaA) * RealizedDietOverlapA,
        #PressureAtoB = alphaB * Proportion_B_Shares + (1 - alphaB) * RealizedDietOverlapB
        PressureBtoA = RealizedDietOverlapA*TPSimilarity,
        PressureAtoB = RealizedDietOverlapB*TPSimilarity
      )
  }
}
  # Find the top competitors
  top_pressure_species <- Competition_weighted %>%
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
result <- get_top_competitors("Lynx pardinus", n = 5, potentialDiet = F)

# Function to compute the total competitive pressure raster for a given species
compute_competitor_pressure <- function(species_name, n = 10, sizeWeight = TRUE, potentialDiet = FALSE) {
  # Step 1: Get top competitors and their pressures
  
  competitors_data <- get_top_competitors(species_name, n = n, sizeWeight = sizeWeight, potentialDiet = potentialDiet)
  
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
  
  if(sizeWeight == TRUE) {
    plot(total_pressure_raster, main = paste("Weighted Total Competitive Pressure on", species_name),
         col = rev(terrain.colors(100)))  
    } else {
      plot(total_pressure_raster, main = paste("Non-weighted Total Competitive Pressure on", species_name),
           col = rev(terrain.colors(100)))
  }
  
  
  return(total_pressure_raster)
}

#############################


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
  plot(richness_raster, main = "Species Richness Map", col = rev(terrain.colors(100)))
  
  return(richness_raster)
}

library(dplyr)
library(terra)
.pa<-unique(Competition_weighted$Predator_A)
.pb<-unique(Competition_weighted$Predator_B)
unique_species<-unique(c(.pa,.pb))

#Load hfp
dir_vars<-"~/Desktop/proyecto_GANSO"
hfp<-rast(paste0(dir_vars,"/variables/hfp_esp_2.5km.tif")) #1 km

# Step 1: Reproject hfp to match spain_raster's CRS
hfp_projected <- project(hfp, spain_raster)

# Step 2: Resample hfp_projected to match the resolution and extent of spain_raster
hfp_resampled <- resample(hfp_projected, spain_raster, method = "bilinear")  # or "near" for categorical

# Optional: plot to check result
plot(hfp_resampled)

hfp_resampled <- mask(hfp_resampled, anyNA(spain_raster), maskvalue=TRUE)

# Function to assess correlation for a species while varying competitor numbers
evaluate_competitor_correlation <- function(species_name, max_competitors = 15, sizeWeight=TRUE, potentialDiet = FALSE) {
  
  # Get species richness map (same for all competitor numbers)
  richness_map <- compute_species_richness(names_sps)
  
  # Get actual species distribution
  actual_distribution <- get.sp.raster(species_name)
  
  # Initialize results dataframe
  correlation_results <- data.frame(
    Species = character(),
    Num_Competitors = integer(),
    Correlation_Actual = numeric(),
    Correlation_Richness = numeric(),
    Corr_Act_Biserial = numeric(),
    Corr_hfp = numeric()
  )
  
  # Loop through 1 to max_competitors
  for (n in 1:max_competitors) {
    cat("Computing for", species_name, "with", n, "competitors...\n")
    
    pressure_map <- compute_competitor_pressure(
      species_name = species_name,
      n = n,
      sizeWeight = sizeWeight,
      potentialDiet = potentialDiet
    )
    
    # Extract values while removing NA
    pressure_values <- values(pressure_map)
    actual_values <- values(actual_distribution)
    richness_values <- values(richness_map)
    hfp_values <- values(hfp_resampled)
    
    # Remove NA values from all datasets
    valid_indices <- !is.na(pressure_values) & !is.na(actual_values) & !is.na(richness_values) &!is.na(hfp_values)
    pressure_values <- pressure_values[valid_indices]
    actual_values <- actual_values[valid_indices]
    richness_values <- richness_values[valid_indices]
    hfp_values <- hfp_values[valid_indices]
    
    # Compute Spearman correlation with actual species distribution
    correlation_actual <- cor(pressure_values, actual_values, method = "spearman")
    
    # Compute Spearman correlation with species richness map
    correlation_richness <- cor(pressure_values, richness_values, method = "spearman")
    
    # Compute Spearman correlation with human footprint
    correlation_hfp <- biserial(hfp_values, actual_values)
    
    # Compute biserial correlation, similar to Pearson but for binary/continuous comparisons
    r_pb <- biserial(pressure_values, actual_values)
    
    # Store results
    correlation_results <- rbind(correlation_results, data.frame(
      Species = species_name,
      Num_Competitors = n,
      Correlation_Actual = correlation_actual,
      Correlation_Richness = correlation_richness,
      Corr_Act_Biserial = r_pb,
      Corr_hfp = correlation_hfp
    ))
  }
  
  return(correlation_results)
}

# Set parameters
potentialDiet <- TRUE
sizeWeight <- TRUE

# Choose a name based on the parameters
result_name <- if (potentialDiet) {
  "all_correlation_results_potential"
} else {
  "all_correlation_results_realized"
}

# Initialize the object
all_correlation_results <- data.frame()

# Run loop
for (species in unique_species) {
  cat("Evaluating:", species, "\n")
  species_results <- evaluate_competitor_correlation(
    species,
    max_competitors = 25,
    sizeWeight = sizeWeight,
    potentialDiet = potentialDiet
  )
  all_correlation_results <- rbind(all_correlation_results, species_results)
}

# Assign it to the dynamic name
assign(result_name, all_correlation_results)

  
# Save the computed data to CSV
write.csv(all_correlation_results_potential, "Corr_competition_potential.csv", row.names = FALSE)
write.csv(all_correlation_results_realized, "Corr_competition_realized.csv", row.names = FALSE)


##---------------------------------POTENTIAL------------------------------------------##
# Select the best number of competitor count for each species
Opt_competitor_per_species_potential <- all_correlation_results_potential %>%
  group_by(Species) %>%
  filter(abs(Corr_Act_Biserial) == max(abs(Corr_Act_Biserial))) %>%
  slice_min(order_by = Num_Competitors, with_ties = FALSE) %>%
  mutate(Corr_Act_Biserial = abs(Corr_Act_Biserial) * sign(Corr_Act_Biserial)) %>%  # Preserve sign
  dplyr::select(Species, Num_Competitors, Corr_Act_Biserial, Correlation_Richness,Corr_hfp) %>%
  arrange(desc(Corr_Act_Biserial))

library(dplyr)

Opt_competitor_per_species_potential <- Opt_competitor_per_species_potential %>%
  left_join(proportion_prey_predator %>%
              dplyr::select(Species, TrophicPosition),
            by = "Species")

Opt_competitor_per_species_potential <- Opt_competitor_per_species_potential %>%
  left_join(bm_combined %>%
              dplyr::select(Species, Class,BodyMass_g),
            by = "Species")

Opt_competitor_per_species_potential <- Opt_competitor_per_species_potential %>%
  mutate(BodyMass_g = as.numeric(BodyMass_g))

Fixed_competitor_per_species_potential <- all_correlation_results_potential %>%
  filter(Num_Competitors == 5) %>%
  mutate(Corr_Act_Biserial = abs(Corr_Act_Biserial) * sign(Corr_Act_Biserial)) %>%  # Keep sign
  dplyr::select(Species, Num_Competitors, Corr_Act_Biserial, Correlation_Richness) %>%
  arrange(desc(Corr_Act_Biserial))

##--------------------------------REALIZED-------------------------------------------##
# Select the best number of competitor count for each species
Opt_competitor_per_species_realized <- all_correlation_results_realized %>%
  group_by(Species) %>%
  filter(abs(Corr_Act_Biserial) == max(abs(Corr_Act_Biserial))) %>%
  slice_min(order_by = Num_Competitors, with_ties = FALSE) %>%
  mutate(Corr_Act_Biserial = abs(Corr_Act_Biserial) * sign(Corr_Act_Biserial)) %>%  # Preserve sign
  dplyr::select(Species, Num_Competitors, Corr_Act_Biserial, Correlation_Richness) %>%
  arrange(desc(Corr_Act_Biserial))

library(dplyr)

Opt_competitor_per_species_realized <- Opt_competitor_per_species_realized %>%
  left_join(proportion_prey_predator %>%
              dplyr::select(Species, TrophicPosition),
            by = "Species")

Opt_competitor_per_species_realized <- Opt_competitor_per_species_realized %>%
  left_join(bm_combined %>%
              dplyr::select(Species, Class,BodyMass_g),
            by = "Species")

Opt_competitor_per_species_realized <- Opt_competitor_per_species_realized %>%
  mutate(BodyMass_g = as.numeric(BodyMass_g))

Fixed_competitor_per_species_realized <- all_correlation_results_realized %>%
  filter(Num_Competitors == 5) %>%
  mutate(Corr_Act_Biserial = abs(Corr_Act_Biserial) * sign(Corr_Act_Biserial)) %>%  # Keep sign
  dplyr::select(Species, Num_Competitors, Corr_Act_Biserial, Correlation_Richness) %>%
  arrange(desc(Corr_Act_Biserial))
##---------------------------------------------------------------------------##

library(ggplot2)
library(ggrepel)  # For better label placement
library(paletteer)

# Ensure both size and color have the same scale and breaks
gg<-ggplot(Opt_competitor_per_species_potential, aes(x = Corr_Act_Biserial, y = Correlation_Richness, label = Species)) +
  geom_point(alpha = 0.8, aes(size = Num_Competitors, color = Num_Competitors)) +  # Maintain point visibility
  scale_color_continuous(limits=c(1, 25), breaks=seq(1, 25, by=2),low = "#f7b9b5", high = "darkred") +
  scale_size_continuous(name = "Number of Competitors", breaks = seq(1, 25, by = 2), limits = c(1, 25)) +  
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
print(gg)


ggsave("Figures/Correlation_tp_With_size_potential.jpg", plot = gg, width = 10, height = 6, dpi = 600)


library(ggplot2)
library(dplyr)

##### FIXED NUMBER OF COMPETITORS ########

ggf<-ggplot(Fixed_competitor_per_species_potential, aes(x = Corr_Act_Biserial, y = Correlation_Richness, label = Species)) +
  geom_point(alpha = 0.8, aes(size = Num_Competitors, color = Num_Competitors)) +  # Maintain point visibility
  scale_color_continuous(limits=c(1, 15), breaks=seq(1, 15, by=2),low = "#f7b9b5", high = "darkred") +
  scale_size_continuous(name = "Number of Competitors", breaks = seq(1, 15, by = 2), limits = c(1, 15)) +  
  #scale_color_viridis_c(name = "Number of Competitors", option = "plasma") +  # Using Viridis for better color mapping
  geom_vline(xintercept = 0, linetype="dashed")+
  geom_text_repel(aes(label = Species), size = 2, color = "black") +  # Keep text readable
  labs(
    title = "Relationship Potential competitive pressure with Species Richness and actual dist",
    x = "Correlation with Actual Distribution",
    y = "Correlation with Species Richness"
  ) +
  guides(
    size = guide_legend(title = "Number of Competitors"),
    color = guide_legend(title = "Number of Competitors")  # Force merging
  ) +
  theme_minimal() +  
  theme(legend.position = "right")
print(ggf)

ggsave("Figures/Correlation_tp_With_size_fixed.jpg", plot = ggf, width = 10, height = 6, dpi = 600)

#########################################



# Combine the data for both treatments into a single dataframe
df_combined <- bind_rows(
  Opt_competitor_per_species_potential %>%
    dplyr::select(Species, Corr_Act_Biserial) %>%
    mutate(Treatment = "Potential"),
  
  Opt_competitor_per_species_realized %>%
    dplyr::select(Species, Corr_Act_Biserial) %>%
    mutate(Treatment = "Realized")
)

# Order species based on ascending correlation in the Weighted treatment
species_order <- df_combined %>%
  filter(Treatment == "Realized") %>%
  arrange(Corr_Act_Biserial) %>%
  pull(Species) %>%
  unique()  # Ensure unique species names

df_combined <- df_combined %>%
  filter(Species != "Neophron percnopterus")

df_combined <- df_combined %>%
  filter(Species != "Parus major")

# Convert Species column into a factor with the specified order
df_combined <- df_combined %>%
  mutate(Species = factor(Species, levels = species_order, ordered = TRUE))


# Calculate mean correlation for each treatment
mean_values <- df_combined %>%
  group_by(Treatment) %>%
  summarise(mean_corr = mean(Corr_Act_Biserial))

# Plot lollipop chart with horizontal mean lines
l1<-ggplot(df_combined, aes(x = Species, y = Corr_Act_Biserial, color = Treatment)) +
  geom_segment(aes(xend = Species, y = 0, yend = Corr_Act_Biserial), size = 0.5, alpha = 0.7) +
  geom_point(size = 3) +
  geom_hline(data = mean_values, aes(yintercept = mean_corr, color = Treatment), 
             linetype = "dashed", size = 1, alpha = 0.8) +  # Add mean line for each treatment
  scale_color_manual(values = c("blue", "red")) +  # Different colors for treatments
  labs(
    title = "Correlation (Biserial) with actual distribution by Species and Treatment",
    x = "Species",
    y = "Correlation (Biserial) with actual distribution"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # Rotate species names for readability
    legend.position = "top"
  )
print(l1)

ggsave("Lollipop_potential_vs_realized_dynamic.jpg", plot = l1, width = 11, height = 6, dpi = 600)

############

# Combine the data for both treatments into a single dataframe
df_combined <- bind_rows(
  Fixed_competitor_per_species_potential %>%
    dplyr::select(Species, Corr_Act_Biserial) %>%
    mutate(Treatment = "Potential"),
  
  Fixed_competitor_per_species_realized %>%
    dplyr::select(Species, Corr_Act_Biserial) %>%
    mutate(Treatment = "Realized")
)

# Order species based on ascending correlation in the Weighted treatment
species_order <- df_combined %>%
  filter(Treatment == "Realized") %>%
  arrange(Corr_Act_Biserial) %>%
  pull(Species) %>%
  unique()  # Ensure unique species names

df_combined <- df_combined %>%
  filter(Species != "Neophron percnopterus")

df_combined <- df_combined %>%
  filter(Species != "Parus major")

# Convert Species column into a factor with the specified order
df_combined <- df_combined %>%
  mutate(Species = factor(Species, levels = species_order, ordered = TRUE))


# Calculate mean correlation for each treatment
mean_values <- df_combined %>%
  group_by(Treatment) %>%
  summarise(mean_corr = mean(Corr_Act_Biserial))

# Plot lollipop chart with horizontal mean lines
l2<-ggplot(df_combined, aes(x = Species, y = Corr_Act_Biserial, color = Treatment)) +
  geom_segment(aes(xend = Species, y = 0, yend = Corr_Act_Biserial), size = 0.5, alpha = 0.7) +
  geom_point(size = 3) +
  geom_hline(data = mean_values, aes(yintercept = mean_corr, color = Treatment), 
             linetype = "dashed", size = 1, alpha = 0.8) +  # Add mean line for each treatment
  scale_color_manual(values = c("blue", "red")) +  # Different colors for treatments
  labs(
    title = "Correlation with actual distribution by Species and Treatment",
    x = "Species",
    y = "Correlation (Biserial) with actual distribution"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # Rotate species names for readability
    legend.position = "top"
  )
print(l2)

ggsave("Lollipop_potential_vs_realized_fixed_comp.jpg", plot = l1, width = 11, height = 6, dpi = 600)
