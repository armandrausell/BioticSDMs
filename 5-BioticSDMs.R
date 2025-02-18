
### biotic SDM models ###

library(enmSdmX)
library(sdm)
library(future.apply)

# Get worldclim vars
vars.es<-worldclim_global(var="bio", res=5, path="vars/worldclim_2.5.tif", version="2.1")
new_names <- paste0("bio", 1:19) #Rename them if not
names(vars.es) <- new_names

ext.es <- ext(-10, 4.3, 35.8, 43.8)  # Extent spain
vars.es <- vars.es[[c("bio5", "bio6","bio12")]]
vars.es<-crop(vars.es, ext.es) # Crop global extent to spain

sp_name<-"Tetrax tetrax"

# Function to run SDM for each species
run_sdm_for_species <- function(sp_name) {
  cat("\nProcessing:", sp_name, "\n")
  print("line 23")
  library(terra)
  library(dplyr)
  library(sdm)
  library(enmSdmX)
  #vars.es<-unwrap(vars.es)
  
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
  
  # Load species raster
  species_raster <-spain_raster
  print(species_raster)
  print("line 30")
  vars.es <- unwrap(vars.es_wrapped)
  print(vars.es)
  print("line 37")
  # Reproject species raster to match vars.es (EPSG:4326)
  species_reprojected <- project(species_raster, vars.es, method = "near")
  print("line 33")
  # Sample presence points
  n1 <- freq(species_reprojected)[2,3]  # Count 1s
  n <- round(n1 * 0.3)  
  sp1 <- sampleRast(species_reprojected, n, adjArea = TRUE, replace = FALSE, prob = TRUE) %>%
    as.data.frame()
  print("line 39")
  # Sample absence points
  n0 <- freq(species_reprojected)[1,3]  # Count 0s
  n0 <- round(n0 * 0.05)  
  inverted_raster <- 1 - species_reprojected  
  sp0 <- sampleRast(inverted_raster, n0, adjArea = TRUE, replace = FALSE, prob = TRUE) %>%
    as.data.frame()
  print("line 46")
  # Merge presence/absence
  sp0$species <- 0
  sp1$species <- 1
  sp <- rbind(sp1, sp0)
  sp <- vect(sp, geom = c("x", "y"))
  print("line 51")
  # Get optimal number of competitors
  num_comp <- Opt_competitor_per_species %>%
    filter(Species == sp_name) %>%
    pull(Num_Competitors)
  
  print("line 58")
  # Compute biotic pressure and reproject it
  pressure <- compute_competitor_pressure(sp_name, n = num_comp)
  pressure_reprojected <- project(pressure, vars.es, method = "bilinear")
  names(pressure_reprojected) <- "pressure"
  print("line 63")
  # Define abiotic and biotic predictors
  biotic <- c(vars.es, pressure_reprojected)
  abiotic <- vars.es
  print("line 67")
  # Build models for Biotic SDM
  db <- sdmData(species ~ ., train = sp, predictors = biotic)
  mb <- sdm(species ~ ., data = db, methods = c('glm', 'gam'), replications = "sub", test.percent = 30, n = 10)
  enb <- ensemble(mb, newdata = biotic, setting = list(method = 'weighted', stat = 'AUC'))
  eb <- evaluates(db, enb)
  auc_b <- eb@statistics$AUC  
  print("line 74")
  # Build models for Abiotic SDM
  da <- sdmData(species ~ ., train = sp, predictors = abiotic)
  ma <- sdm(species ~ ., data = da, methods = c('glm', 'gam'), replications = "sub", test.percent = 30, n = 10)
  ena <- ensemble(ma, newdata = abiotic, setting = list(method = 'weighted', stat = 'AUC'))
  ea <- evaluates(da, ena)
  auc_a <- ea@statistics$AUC  
  print("line 81")
  # Store AUC values in a dataframe
  return(data.frame(Species = sp_name, AUC_Biotic = auc_b, AUC_Abiotic = auc_a))
}

library(future.apply)

# Wrap vars.es before exporting
vars.es_wrapped <- wrap(vars.es)

# Set up parallel execution with available CPU cores
num_cores <- 3
cl <- makeCluster(num_cores)  # Create a cluster
plan(cluster, workers = cl)  # Assign cluster to future.apply

# Ensure all nodes load required packages
clusterEvalQ(cl, {
  library(terra)
  library(dplyr)
  library(sdm)
  library(enmSdmX)  # If needed
  print("Packages loaded in worker")
})

# Explicitly export `get.sp.raster()` and required objects
clusterExport(cl, varlist = c("names_species_cache", "get.sp.raster", "Species_spain", "spain_raster", "species_raster_cache",
                              "compute_competitor_pressure", "vars.es", "Opt_competitor_per_species"), 
              envir = .GlobalEnv)

# Select only first 3 species
species_list <- head(Opt_competitor_per_species$Species, 3)

# Run in parallel using future_lapply
auc_results <- bind_rows(future_lapply(species_list, run_sdm_for_species))

# Save results
write.csv(auc_results, "SDM_AUC_Results.csv", row.names = FALSE)

# Print results
print(auc_results)

