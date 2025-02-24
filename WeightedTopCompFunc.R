

get_top_competitors <- function(species_name, n = 5, sizeWeight = TRUE)  { 
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
    
    # Apply Gaussian weighting to favor similar-sized competitors
    Competition_weighted <- Competition_weighted %>%
      mutate(
        MassSimilarity = gaussian_similarity(BodyMass_A, BodyMass_B, sigma = 1.8),  # Adjust sigma
        PressureBtoA = round(
          (alphaA * Proportion_A_Shares + (1 - alphaA) * RealizedDietOverlapA) * MassSimilarity, 3
        ),
        PressureAtoB = round(
          (alphaB * Proportion_B_Shares + (1 - alphaB) * RealizedDietOverlapB) * MassSimilarity, 3
        )
      )
    
  } else {
    
    # Calculate pressure for the standard method
    Competition_weighted <- Competition_spain_df %>%
      mutate(
        PressureBtoA = alphaA * Proportion_A_Shares + (1 - alphaA) * RealizedDietOverlapA,
        PressureAtoB = alphaB * Proportion_B_Shares + (1 - alphaB) * RealizedDietOverlapB
      )
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
get_top_competitors("Tetrax tetrax",n=12,sizeWeight = T)
get_top_competitors("Tetrax tetrax",n=12,sizeWeight = F)
