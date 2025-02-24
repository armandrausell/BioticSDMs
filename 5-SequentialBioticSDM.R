
#SEQUENTIAL BioticSDMS ####

library(terra)
library(dplyr)
library(sdm)

# Function to run SDM for each species
run_sdm_for_species <- function(sp_name) {
  cat("\nProcessing:", sp_name, "\n")
  
  # Load species raster
  species_raster <- get.sp.raster(sp_name)
  
  # Reproject species raster to match vars.es (EPSG:4326)
  species_reprojected <- project(species_raster, vars.es, method = "near")
  
  # Sample presence points
  n1 <- freq(species_reprojected)[2,3]  # Count 1s
  n <- round(n1 * 0.3)  
  sp1 <- sampleRast(species_reprojected, n, adjArea = TRUE, replace = FALSE, prob = TRUE) %>%
    as.data.frame()
  
  # Sample absence points
  n0 <- freq(species_reprojected)[1,3]  # Count 0s
  n0 <- round(n0 * 0.05)  
  inverted_raster <- 1 - species_reprojected  
  sp0 <- sampleRast(inverted_raster, n0, adjArea = TRUE, replace = FALSE, prob = TRUE) %>%
    as.data.frame()
  
  # Merge presence/absence
  sp0$species <- 0
  sp1$species <- 1
  sp <- rbind(sp1, sp0)
  sp <- vect(sp, geom = c("x", "y"))
  
  # Get optimal number of competitors
  num_comp <- Opt_competitor_per_species %>%
    filter(Species == sp_name) %>%
    pull(Num_Competitors)
  
  # Compute biotic pressure and reproject it
  pressure <- compute_competitor_pressure(sp_name, n = num_comp)
  pressure_reprojected <- project(pressure, vars.es, method = "bilinear")
  names(pressure_reprojected) <- "pressure"
  
  # Define abiotic and biotic predictors
  biotic <- c(vars.es, pressure_reprojected)
  abiotic <- vars.es
  
  # Build models for Biotic SDM
  db <- sdmData(species ~ ., train = sp, predictors = biotic)
  mb <- sdm(species ~ ., data = db, methods = c('glm', 'gam'), replications = "sub", test.percent = 30, n = 10)
  enb <- ensemble(mb, newdata = biotic, setting = list(method = 'weighted', stat = 'AUC'))
  eb <- evaluates(db, enb)
  auc_b <- eb@statistics$AUC  
  
  # Build models for Abiotic SDM
  da <- sdmData(species ~ ., train = sp, predictors = abiotic)
  ma <- sdm(species ~ ., data = da, methods = c('glm', 'gam'), replications = "sub", test.percent = 30, n = 10)
  ena <- ensemble(ma, newdata = abiotic, setting = list(method = 'weighted', stat = 'AUC'))
  ea <- evaluates(da, ena)
  auc_a <- ea@statistics$AUC  
  
  # Store AUC values in a dataframe
  return(data.frame(Species = sp_name, AUC_Biotic = auc_b, AUC_Abiotic = auc_a))
}

# Run the function for all species
species_list <- Opt_competitor_per_species$Species  
auc_results <- bind_rows(lapply(species_list, run_sdm_for_species))

# Save AUC results
write.csv(auc_results, "SDM_AUC_Results.csv", row.names = FALSE)

# Print final results
print(auc_results)

# Convert to long format for ggplot2
df_long <- auc_results %>%
  pivot_longer(cols = c(AUC_Biotic, AUC_Abiotic), 
               names_to = "Model_Type", 
               values_to = "AUC")

# Create Boxplot
ggplot(df_long, aes(x = Model_Type, y = AUC, fill = Model_Type)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +  # Adds data points
  labs(title = "Comparison of Biotic and Abiotic SDMs",
       x = "Model Type",
       y = "AUC Value") +
  scale_fill_manual(values = c("AUC_Biotic" = "blue", "AUC_Abiotic" = "red")) +
  theme_minimal()

library(ggplot2)
library(tidyr)


# Convert to long format for ggplot2
df_long <- auc_results %>%
  pivot_longer(cols = c(AUC_Biotic, AUC_Abiotic), 
               names_to = "Model_Type", 
               values_to = "AUC")

# Convert species to a factor for ordered x-axis
df_long$Species <- factor(df_long$Species, levels = unique(df_long$Species))

# Create Scatter Plot
ggplot(df_long, aes(x = Species, y = AUC, color = Model_Type)) +
  geom_point(size = 4, alpha = 0.8, position = position_dodge(width = 0.5)) +  # Dodge to separate Biotic & Abiotic
  labs(title = "AUC Comparison Between Biotic and Abiotic SDMs",
       x = "Species",
       y = "AUC Value") +
  scale_color_manual(values = c("AUC_Biotic" = "blue", "AUC_Abiotic" = "red"), 
                     labels = c("AUC Biotic", "AUC Abiotic")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate species names if needed

