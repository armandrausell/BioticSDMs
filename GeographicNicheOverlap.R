
library(terra)
library(raster)
#Geographic niche overlap

#load environmental variables
dir_vars<- "~/Desktop/vars/"
vars<-raster::brick(paste0(dir_vars,"variables_stack_current.tif"))
gpp<-rast(paste0(dir_vars, "GPP_1km_res.tif"))
elevation<-rast(paste0(dir_vars, "srtm_1km_res.tif"))
hfp<-rast(paste0(dir_vars,"HFP2009ddeg.tif")) #1 km
fractions<-rast("~/Desktop/proyecto_GANSO/variables/fractions_esp_2.5km.tif")
fractions<-fractions[[-c(2)]]

#Load vector and create raster of Spain
spain_raster <- rast("spain_rasterUTM.tif")
plot(spain_raster)

# Step 1: Reproject hfp to match spain_raster's CRS
vars_p <- project(rast(vars), spain_raster)

gpp_p <- project(gpp, spain_raster)
names(gpp_p)<-"gpp"

elevation_p <- project(elevation, spain_raster)
names(elevation_p)<-"elevation"

hfp_p <- project(hfp, spain_raster)
names(hfp_p)<-"hfp"

fractions_p <- project(fractions, spain_raster)
names(fractions_p)<-c("grass","crops", "shrub","water")

vars_p<-c(vars_p,gpp_p,elevation_p,hfp_p,fractions_p)

# Resample hfp_projected to match the resolution and extent of spain_raster
vars_res <- resample(vars_p, spain_raster, method = "bilinear")  # or "near" for categorical

# plot to check result
plot(vars_res)

vars_res <- terra::mask(vars_res, anyNA(spain_raster), maskvalue=T)
plot(vars_res)


###############################

shared_env_niche <- function(sp1_name, sp2_name) {
  # Step 1: Get species raster (binary presence/absence)
  sp1.rast <- get.sp.raster(sp1_name)
  sp2.rast <- get.sp.raster(sp2_name)
  
  # Step 2: Extract presence cells
  cells_sp1 <- which(!is.na(terra::values(sp1.rast)) & terra::values(sp1.rast) == 1)
  cells_sp2 <- which(!is.na(terra::values(sp2.rast)) & terra::values(sp2.rast) == 1)
  
  # Step 3: Extract environmental values at presence locations
  ex1 <- terra::extract(vars_res, cells_sp1, xy = TRUE)
  ex2 <- terra::extract(vars_res, cells_sp2, xy = TRUE)
  
  env1 <- ex1[, !(names(ex1) %in% c("ID", "x", "y"))]
  env2 <- ex2[, !(names(ex2) %in% c("ID", "x", "y"))]
  
  # Step 4: Calculate shared range for each variable
  shared_ranges <- do.call(rbind, lapply(colnames(env1), function(var) {
    data.frame(
      Variable = var,
      Min = max(min(env1[[var]], na.rm = TRUE), min(env2[[var]], na.rm = TRUE)),
      Max = min(max(env1[[var]], na.rm = TRUE), max(env2[[var]], na.rm = TRUE))
    )
  }))
  
  # Step 5: Create logical mask for all shared environmental conditions
  shared_env_mask <- rep(TRUE, ncell(vars_res))
  for (i in 1:nrow(shared_ranges)) {
    var <- shared_ranges$Variable[i]
    min_val <- shared_ranges$Min[i]
    max_val <- shared_ranges$Max[i]
    v <- terra::values(vars_res[[var]])
    shared_env_mask <- shared_env_mask & (v >= min_val & v <= max_val)
  }
  
  shared_env_raster <- vars_res[[1]]
  terra::values(shared_env_raster) <- as.numeric(shared_env_mask)
  
  # Step 6: Binary layers for each condition
  binary_layers <- list()
  for (i in 1:nrow(shared_ranges)) {
    var <- shared_ranges$Variable[i]
    min_val <- shared_ranges$Min[i]
    max_val <- shared_ranges$Max[i]
    layer <- vars_res[[var]]
    
    binary_layer <- classify(layer, rcl = matrix(c(
      -Inf, min_val-1e-6, 0,
      min_val-1e-6, max_val+1e-6, 1,
      max_val+1e-6, Inf, 0
    ), ncol = 3, byrow = TRUE))
    
    binary_layers[[var]] <- binary_layer
  }
  
  binary_stack <- rast(binary_layers)
  sum_binary <- sum(binary_stack)
  
  final_shared_raster <- classify(sum_binary, matrix(c(
    -Inf, nrow(shared_ranges) - 0.1, 0,
    nrow(shared_ranges) - 0.1, Inf, 1
  ), ncol = 3, byrow = TRUE))
  
  # Step 7: Plot results
  par(mfrow = c(1, 3))
  plot(sp1.rast, main = sp1_name)
  plot(sp2.rast, main = sp2_name)
  plot(final_shared_raster, main = "Shared Environmental Niche", col = c("grey80", "darkgreen"), legend = T)
  
  # Return outputs for further use
  return(list(
    shared_ranges = shared_ranges,
    shared_env_raster = shared_env_raster,
    final_shared_raster = final_shared_raster
  ))
}
shared_env_niche("Lanius collurio", "Lanius senator")

sp1_name<-"Lanius collurio"
sp2_name<-"Lanius senator"

#This function retrieves the raster for which two species have excluding competition
exclusion_zone<-function(sp1_name,sp2_name, plot=TRUE){
# Run shared environmental niche function
result <- shared_env_niche(sp1_name, sp2_name)

# Get species binary rasters
sp1.rast <- get.sp.raster(sp1_name)
sp2.rast <- get.sp.raster(sp2_name)

#  Conditional rasters
mask_sp1 <- as.numeric((sp1.rast == 1) & (sp2.rast != 1) & (result$final_shared_raster == 1))
mask_sp2 <- as.numeric((sp2.rast == 1) & (sp1.rast != 1) & (result$final_shared_raster == 1))

both_sp <- as.numeric((sp1.rast == 1) & (sp2.rast == 1) & (result$final_shared_raster == 1))
true_overlap <- as.numeric((sp1.rast == 1) & (sp2.rast == 1))

non_overlap <- as.numeric(result$final_shared_raster - mask_sp1 - mask_sp2)
pot_comp <- as.numeric(result$final_shared_raster - non_overlap)
# Create a clean raster with 1 for TRUE, 0 for FALSE, NA elsewhere
if(plot==TRUE){
par(mfrow = c(2, 2))
plot(sp1.rast, main = paste0(sp1_name), col = c("gray", "blue"), legend = FALSE)
plot(sp2.rast, main = paste0(sp2_name), col = c("gray", "red"), legend = FALSE)
plot(result$final_shared_raster, main = "Shared environmental space", col = c("gray", "darkgreen"), legend = F)
plot(pot_comp, main = "Where exclusion might occur (in shared env.)", col = c("gray", "darkgreen"), legend = F)

par(mfrow = c(1, 2))
plot(true_overlap, main = "Overlap in real space",col = c("gray", "darkgreen"))
plot(both_sp, main = "Overlap in shared environmental space",col = c("gray", "darkgreen"))

par(mfrow = c(1, 3))
plot(mask_sp1, main = paste0("Where ",sp1_name, " potentially excludes"), col = c("gray", "blue"), legend = FALSE)
plot(mask_sp2, main = paste0("Where ",sp2_name, " potentially excludes"), col = c("gray", "red"), legend = FALSE)
plot(pot_comp, main = "Where exclusion occurs (in shared env.)", col = c("gray", "darkgreen"), legend = F)

}
dev.off()
return(pot_comp)
}

exclusion_zone(sp1_name,sp2_name, plot=TRUE)

##########

#In this part, we will stack the exclusion areas of each competitor pair
# Initialize an empty list to store all exclusion rasters
exclusion_rasters <- list()

# Loop through species list
for (spA in unique_species) {
  # Get top 2 exclusive competitors
  competitors <- get_exclusive_comp(spA)
  message("Processing ", spA)
  
  # Skip if less than 1 competitors found
  if (length(competitors) < 1) next
  
  # Get top 2 competitors
  comp1 <- competitors[1]
  comp2 <- competitors[2]
  
  if (length(competitors) == 1) {
  comp1 <- competitors[1]
  rast1 <- try(exclusion_zone(spA, comp1, plot = FALSE), silent = TRUE)
  if (inherits(rast1, "SpatRaster")) exclusion_rasters[[length(exclusion_rasters) + 1]] <- rast1
  
  } else {
    comp1 <- competitors[1]
    comp2 <- competitors[2]
    
    rast1 <- try(exclusion_zone(spA, comp1, plot = FALSE), silent = TRUE)
    rast2 <- try(exclusion_zone(spA, comp2, plot = FALSE), silent = TRUE)
    
    if (inherits(rast1, "SpatRaster")) exclusion_rasters[[length(exclusion_rasters) + 1]] <- rast1
    if (inherits(rast2, "SpatRaster")) exclusion_rasters[[length(exclusion_rasters) + 1]] <- rast2
    
  }
  # If both runs succeeded, store results
}

# Stack and sum all exclusion rasters into one
if (length(exclusion_rasters) > 0) {
  exclusion_stack <- terra::rast(exclusion_rasters)
  global_exclusion_zone <- sum(exclusion_stack, na.rm = TRUE)
  
  # Plot final exclusion map
  plot(global_exclusion_zone, main = "Global Exclusion Zone", col = hcl.colors(100, "Inferno"), legend = TRUE)
} else {
  message("No exclusion zones could be computed.")
}

richness<-compute_species_richness(names_sps)
global_excl_corr<-global_exclusion_zone/richness
plot(global_excl_corr)
head(values(global_excl_corr, na.rm =T))
writeRaster(global_excl_corr, filename = "vars/GlobalExclusionZonesCorr.tif", overwrite = T)

# 1. Extract raster values as a vector
vals <- values(global_excl_corr)

# 2. Find the 10 highest non-NA values
top_10_vals <- sort(unique(vals[!is.na(vals)]), decreasing = TRUE)[1:10]

# 3. Set those values to NA in a copy of the raster
corrected_raster <- global_excl_corr
values(corrected_raster)[vals %in% top_10_vals] <- NA

plot(corrected_raster, main = "Potential exclusion areas", col = hcl.colors(100, "Inferno"))
