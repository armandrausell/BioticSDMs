

# Calculation of interaction strenghts prey-predator from binary matrix

# Filter the matrix to only include species present in names_sps (both rows and columns)
binary_matrix_filtered <- binary_matrix_df[
  rownames(binary_matrix_df) %in% names_sps,
  colnames(binary_matrix_df) %in% names_sps
]

# Trophic level with NetIndices
tl<-TrophInd(t(binary_matrix_filtered))
library(tibble)
tl <- tibble::rownames_to_column(tl, "Species")

# Step 1: Extract predator-prey pairs from binary matrix
get_predator_prey_pairs <- function(binary_matrix_filtered) {
  # Ensure species names are row and column names
  predator_species <- rownames(binary_matrix_filtered)
  prey_species <- colnames(binary_matrix_filtered)
  
  # Create an empty data frame to store pairs
  interaction_list <- data.frame(Predator = character(), Prey = character(), stringsAsFactors = FALSE)
  
  # Loop over the matrix and record where interaction == 1
  for (i in seq_along(predator_species)) {
    for (j in seq_along(prey_species)) {
      if (binary_matrix_filtered[i, j] == 1) {
        interaction_list <- rbind(
          interaction_list,
          data.frame(Predator = predator_species[i], Prey = prey_species[j])
        )
      }
    }
  }
  
  return(interaction_list)
}

# Apply the function to your binary matrix
predation_pairs <- get_predator_prey_pairs(binary_matrix_filtered)

# Add trophic positions from the proportion_prey_predator dataframe

# Rename for clarity and avoid duplicate column names during join
trophic_info_tl <- tl %>%
  dplyr::select(Species, TL)

trophic_info_prop <- proportion_prey_predator %>%
  dplyr::select(Species, TrophicPosition)

trophic_info$TrophicPosition<-trophic_info_prop$TrophicPosition*trophic_info_tl$TL
# Join TrophicPosition for predators
predation_pairs_tp <- predation_pairs %>%
  left_join(trophic_info, by = c("Predator" = "Species")) %>%
  rename(TrophicPositionPred = TrophicPosition)

# Join TrophicPosition for prey
predation_pairs_tp <- predation_pairs_tp %>%
  left_join(trophic_info, by = c("Prey" = "Species")) %>%
  rename(TrophicPositionPrey = TrophicPosition)

# Remove self-predation cases
predation_pairs_tp <- predation_pairs_tp %>%
  filter(Predator != Prey)
########################################################
# Create a short lookup table from proportion_prey_predator
pred_info <- proportion_prey_predator %>%
  dplyr::select(Species, TotalPrey)

# Add TotalPrey for each predator
predation_pairs_tp1 <- predation_pairs_tp %>%
  left_join(pred_info, by = c("Predator" = "Species")) %>%
  rename(TotalPrey_Pred = TotalPrey)

prey_info <- proportion_prey_predator %>%
  dplyr::select(Species, Times_Preyed_On)

# Add Times_Preyed_On for each prey
predation_pairs_tp1 <- predation_pairs_tp1 %>%
  left_join(prey_info, by = c("Prey" = "Species")) %>%
  rename(Times_Prey_Preyed = Times_Preyed_On)
#######################################################
library(dplyr)

# Compute the median TP for each predator
median_tp_by_predator <- predation_pairs_tp1 %>%
  group_by(Predator) %>%
  summarise(
    MedianTP = median(TrophicPositionPrey, na.rm = TRUE),
    .groups = "drop"
  )

# Join the median back into predation_pairs_tp1
predation_pairs_tp1 <- predation_pairs_tp1 %>%
  left_join(median_tp_by_predator, by = "Predator")

# Compute absolute deviation from the median
predation_pairs_tp1 <- predation_pairs_tp1 %>%
  mutate(AbsDeviation = abs(TrophicPositionPrey - MedianTP))

# Compute the max deviation per predator
max_deviation_by_predator <- predation_pairs_tp1 %>%
  group_by(Predator) %>%
  summarise(MaxDeviation = max(AbsDeviation, na.rm = TRUE), .groups = "drop")

# Join max deviation back
predation_pairs_tp1 <- predation_pairs_tp1 %>%
  left_join(max_deviation_by_predator, by = "Predator")

# Compute Trophich Core Proximity (TCP)
predation_pairs_tp1 <- predation_pairs_tp1 %>%
  mutate(
    TrophicCoreProximity = ifelse(MaxDeviation > 0,
                                  1 - (AbsDeviation / MaxDeviation),
                                  1)  # If no deviation (all TP same), assign 1
  )

#######################################################
# Calculate Relative Predator Pressure (rPP)
Relative_pred_pressure <- predation_pairs_tp1 %>%
  mutate(
    rPP = (abs(TrophicPositionPred - TrophicPositionPrey) * (1 / TotalPrey_Pred)
           * TrophicCoreProximity) / 
      (Times_Prey_Preyed + 1)
  )

Relative_test <- predation_pairs_tp1 %>%
  mutate(
    rPP = (log2(abs(TrophicPositionPred - TrophicPositionPrey)) * 
             (1 / (1+TotalPrey_Pred*(1-TrophicCoreProximity)))) 
    / 
      (Times_Prey_Preyed + 1)
  )


############################################################
compare_predators_for_prey <- function(predator1, predator2, prey, df = Relative_test) {
  # Check for column existence
  required_cols <- c("Predator", "Prey", "TrophicPositionPred", "TotalPrey_Pred", 
                     "Times_Prey_Preyed", "MedianTP", "AbsDeviation", 
                     "TrophicCoreProximity", "rPP")

  # Filter the data
  df_filtered <- Relative_test %>%
    dplyr::filter(Prey == prey, Predator %in% c(predator1, predator2)) %>%
    dplyr::select(all_of(required_cols))  # safer than listing columns manually
  
  # Output
  #print(df_filtered, row.names = FALSE)
  return(df_filtered)
}

compare_predators_for_prey("Dendrocopos medius", "Jynx torquilla", "Parus major")
