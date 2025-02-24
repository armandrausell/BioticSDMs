
### biotic SDM models ###

library(enmSdmX)
library(sdm)
library(future.apply)
library(parallel)

# Get worldclim vars
vars.es<-worldclim_global(var="bio", res=5, path="vars/worldclim_2.5.tif", version="2.1")
new_names <- paste0("bio", 1:19) #Rename them if not
names(vars.es) <- new_names

ext.es <- ext(-10, 4.3, 35.8, 43.8)  # Extent spain
vars.es <- vars.es[[c("bio5", "bio6","bio12")]]
vars.es<-crop(vars.es, ext.es) # Crop global extent to spain

sp_name<-"Zamenis longissimus"

# Function to run SDM for each species
run_sdm_for_species <- function(sp_name,n=15) {
  cat("\nProcessing:", sp_name, "\n")
  print("line 23")
  library(terra)
  library(dplyr)
  library(sdm)
  library(enmSdmX)
  #vars.es<-unwrap(vars.es)
  Species_spain<-unwrap(Species_spain_wrapped)
  vars.es <- unwrap(vars.es_wrapped)
  spain_raster<-unwrap(spain_raster_wrapped)
  spain_raster[spain_raster == 1] <- 0
  print("line 32")
  
  # Get all species names from the dataset
  existing_names <- names(Species_spain)
  
  # Generate possible name formats "." or " " 
  species_to_plot_dot <- gsub(" ", ".", sp_name)  # Replace space with dot
  species_to_plot_space <- gsub("\\.", " ", sp_name)  # Replace dot with space
  print("line 40")
  # Find the correct species name
  if (species_to_plot_dot %in% existing_names) {
    species_to_plot <- species_to_plot_dot
  } else if (species_to_plot_space %in% existing_names) {
    species_to_plot <- species_to_plot_space
  } else {
    stop(paste("Species not found in the dataset:", sp_name))
  }
  print("line 49")
  # Filter the points for the selected species
  sp.points <- Species_spain[Species_spain[[species_to_plot]] == 1, ]
  
  # Check if there are points for this species
  if (nrow(sp.points) == 0) {
    stop(paste("No occurrence data for", species_to_plot))
  }
  print("line 57")
  # Rasterize the species' occurrence points
  points_raster <- rasterize(sp.points, spain_raster, field=1)
  print("line 60")
  
  spain_raster_copy <- spain_raster  # Keep original raster
  spain_raster_copy[!is.na(points_raster)] <- 1  # Modify only the copy
  species_raster<-spain_raster_copy
  
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
  
  plot(species_reprojected)
  plot(vect(sp0, geom = c("x", "y")), add=T, cex=0.5, col="blue")
  plot(vect(sp1, geom = c("x", "y")), add=T, cex=0.5, col="red")
  print("line 51")
  # Get optimal number of competitors
  num_comp <- Opt_competitor_per_species %>%
    filter(Species == sp_name) %>%
    pull(Num_Competitors)
  
  print("line 58")
  # Compute biotic pressure and reproject it
  ######## compute competitor pressure function ####
  ##
  #
  top_pressure_species <- Competition_spain_df %>%
    filter(Predator_A == sp_name | Predator_B == sp_name) %>%
    mutate(
      # Select the correct pressure column based on position
      Pressure_on_Species = ifelse(Predator_A == sp_name, PressureBtoA, PressureAtoB),
      
      # Identify the species exerting pressure on the target species
      Pressure_Exerting_Species = ifelse(Predator_A == sp_name, Predator_B, Predator_A)
    ) %>%
    arrange(desc(Pressure_on_Species)) %>%  # Sort by highest pressure
    slice(1:num_comp)  # Select the top `n` species
  
  # Extract unique competitor species
  top_competitors <- unique(top_pressure_species$Pressure_Exerting_Species)
  competitors_data <- top_competitors
  top_pressure_table <- top_pressure_species
  
  ### COMPUTE get sp ###
  
  # Step 2: Initialize an empty raster (same extent/resolution as Spain raster)
  total_pressure_raster <- spain_raster
  values(total_pressure_raster) <- 0  # Start with all values as 0
  
  
  # Step 3: Iterate over competitors, get their raster, and weight it
  for (i in 1:num_comp) {
    competitor <- top_pressure_table$Pressure_Exerting_Species[i]
    pressure_value <- top_pressure_table$Pressure_on_Species[i]
    
    # Print competitor and corresponding pressure value
    print(paste("Competitor:", competitor, "- Pressure Value:", pressure_value))
    # Get all species names from the dataset
    existing_names <- names(Species_spain)
    
    # Generate possible name formats "." or " " 
    species_to_plot_dot <- gsub(" ", ".", competitor)  # Replace space with dot
    species_to_plot_space <- gsub("\\.", " ", competitor)  # Replace dot with space
    #print("line 40")
    # Find the correct species name
    if (species_to_plot_dot %in% existing_names) {
      species_to_plot <- species_to_plot_dot
    } else if (species_to_plot_space %in% existing_names) {
      species_to_plot <- species_to_plot_space
    } else {
      stop(paste("Species not found in the dataset:", sp_name))
    }
    #print("line 49")
    # Filter the points for the selected species
    sp.points <- Species_spain[Species_spain[[species_to_plot]] == 1, ]
    
    # Check if there are points for this species
    if (nrow(sp.points) == 0) {
      stop(paste("No occurrence data for", species_to_plot))
    }
    #print("line 57")
    # Rasterize the species' occurrence points
    points_raster <- rasterize(sp.points, total_pressure_raster, field=1)
   
      # Assign 1 where points are located
    spain_raster_copy <- spain_raster  # Keep original raster
    spain_raster_copy[!is.na(points_raster)] <- 1
      
      competitor_raster<-spain_raster_copy
    
    ##############################
    # Multiply raster by pressure value
    weighted_raster <- competitor_raster * pressure_value
    
    # Add weighted raster to the total pressure raster
    total_pressure_raster <- total_pressure_raster + weighted_raster
    pressure<-total_pressure_raster
    #plot(pressure, main=paste0("Pressure of: ", competitor))
  }
  ######################
  
  #
  ##
  ##################################################
  
  #print(pressure)
  print("line 101")
  pressure_reprojected <- project(pressure, vars.es, method = "near")
  names(pressure_reprojected) <- "pressure"
  
  #print(pressure_reprojected)
  print("line 63")
  # Define abiotic and biotic predictors
  biotic <- c(vars.es[[c(1,2)]], pressure_reprojected)
  #print(biotic)
  abiotic <- vars.es
  #print(abiotic)
  print("line 67")
  
  print(sp)
  # Build models for Biotic SDM
  db <- sdmData(species ~ ., train = sp, predictors = biotic)
  mb <- sdm(species ~ ., data = db, methods = c('glm', 'gam'), replications = "sub", test.percent = 30, n = 20)
  enb <- ensemble(mb, newdata = biotic, setting = list(method = 'weighted', stat = 'AUC'))
  eb <- evaluates(db, enb)
  auc_b <- eb@statistics$AUC 
  print(auc_b)
  print("line 74")
  # Build models for Abiotic SDM
  da <- sdmData(species ~ ., train = sp, predictors = abiotic)
  ma <- sdm(species ~ ., data = da, methods = c('glm', 'gam'), replications = "sub", test.percent = 30, n = 20)
  ena <- ensemble(ma, newdata = abiotic, setting = list(method = 'weighted', stat = 'AUC'))
  ea <- evaluates(da, ena)
  auc_a <- ea@statistics$AUC  
  print("line 81")
  # Store AUC values in a dataframe
  
  #return(list(DataFrame = auc_data, Plot = auc_plot))
  
  return(data.frame(Species = sp_name, AUC_Biotic = auc_b, AUC_Abiotic = auc_a))
}

library(future.apply)

# Wrap vars.es before exporting
vars.es_wrapped <- wrap(vars.es)
Species_spain_wrapped<-wrap(Species_spain)
spain_raster_wrapped<-wrap(spain_raster)

# Set up parallel execution with available CPU cores
num_cores <- 5
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
species_list <- head(Opt_competitor_per_species$Species, 9)

# Run in parallel using future_lapply
auc_results <- bind_rows(future_lapply(species_list, run_sdm_for_species, future.seed = T))

# Save results
#write.csv(auc_results, "SDM_AUC_Results.csv", row.names = FALSE)

# Print results
print(auc_results)

