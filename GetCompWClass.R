

Competition_spain_df<-Competition_weighted

# Find competitors and prioritize same-class species
get_top_competitors<-function(species_name, n = 5)  { #amount of species to get
  # Get the class of the target species (Vipera aspis in this case)
  species_class <- bm_combined %>%
    filter(Species == species_name) %>%  # Filter by species name
    pull(Class)
  
  # Find competitors and determine their class
  top_pressure_species <- Competition_spain_df %>%
    filter(Predator_A == species_name | Predator_B == species_name) %>%
    mutate(
      # Select the correct pressure column based on position
      Pressure_on_Species = ifelse(Predator_A == species_name, PressureBtoA, PressureAtoB),
      
      # Identify the species exerting pressure on the target species
      Pressure_Exerting_Species = ifelse(Predator_A == species_name, Predator_B, Predator_A),
      
      # Assign class to the exerting species
      Pressure_Exerting_Species_Class = ifelse(Predator_A == species_name, Class_B, Class_A),
      
      # Create a flag: 1 if the competitor is from the same class, 0 otherwise
      SameClass = ifelse(Pressure_Exerting_Species_Class == species_class, 1, 0)
    ) %>%
    arrange(desc(SameClass), desc(Pressure_on_Species)) %>%  # Prioritize same class, then pressure
    slice(1:n)  # Select top n competitors
  
  # Print results
  print(top_pressure_species)
  top_competitors <- unique(top_pressure_species$Pressure_Exerting_Species)
  print(top_competitors)
  return(list(Top_Pressure_Table = top_pressure_species, Top_Competitors = top_competitors))
}
get.comp("Glis glis",n=15)
get_top_competitors("Glis glis",n=15)
