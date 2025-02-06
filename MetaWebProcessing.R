
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

