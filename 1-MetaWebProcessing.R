
setwd("~/Desktop/Github/BioticSDMs")

df<-read.csv2("metaweb/TetraEU_pairwise_interactions.csv", )

# Load required library
library(dplyr)

# Get unique species names
species <- unique(c(df$sourceTaxonName, df$targetTaxonName))

# Create an empty binary matrix
binary_matrix <- matrix(0, nrow = length(species), ncol = length(species), 
                        dimnames = list(species, species))

# Fill the matrix based on interactions
for (i in 1:nrow(df)) {
  predator <- df$sourceTaxonName[i]
  prey <- df$targetTaxonName[i]
  binary_matrix[predator, prey] <- 1
}

# Convert to data frame for better readability
binary_matrix_df <- as.data.frame(binary_matrix)

# View the binary matrix
print(binary_matrix_df)

# Write to CSV
write.csv(binary_matrix_df, "binary_interaction_matrix.csv", row.names = TRUE)

# Compute the shared predation matrix
shared_predation_matrix <- binary_matrix %*% t(binary_matrix)

# Convert to dataframe for better readability
shared_predation_df <- as.data.frame(shared_predation_matrix)

# Display the shared predation matrix
print(shared_predation_df)

# Get predator species names
predators <- rownames(shared_predation_df)

# Get the total number of prey each predator eats (diagonal values)
total_prey <- diag(as.matrix(shared_predation_df))
total_prey_df <- data.frame(Predator = predators, TotalPrey = total_prey)

# Create an empty dataframe to store interactions
interaction_df <- data.frame(Predator_A = character(), 
                             Predator_B = character(), 
                             SharedPreyCount = numeric(), 
                             Proportion_A_Shares = numeric(), 
                             Proportion_B_Shares = numeric(), 
                             TotalPrey_A = numeric(), 
                             TotalPrey_B = numeric(), 
                             stringsAsFactors = FALSE)

# Loop through unique predator pairs (takes around 1-2 minutes)
for (i in 1:(length(predators) - 1)) {
  for (j in (i + 1):length(predators)) {
    predator_A <- predators[i]
    predator_B <- predators[j]
    
    shared_prey_count <- shared_predation_df[predator_A, predator_B]
    
    if (shared_prey_count > 0) {  # Only keep pairs that share at least one prey
      total_prey_A <- total_prey[i]
      total_prey_B <- total_prey[j]
      
      proportion_A_shares <- shared_prey_count / total_prey_A
      proportion_B_shares <- shared_prey_count / total_prey_B
      
      interaction_df <- rbind(interaction_df, 
                              data.frame(Predator_A = predator_A, 
                                         Predator_B = predator_B, 
                                         SharedPreyCount = shared_prey_count, 
                                         Proportion_A_Shares = proportion_A_shares, 
                                         Proportion_B_Shares = proportion_B_shares, 
                                         TotalPrey_A = total_prey_A, 
                                         TotalPrey_B = total_prey_B))
    }
  }
}

# View the dataframe
print(interaction_df)
