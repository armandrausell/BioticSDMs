
library(virtualspecies)
library(raster)
library(terra)

#This code creates copies of the real species in the real world, supposedly less
# biased by biotic interactions. It also computes  

#load environmental variables
dir_vars<- "~/Desktop/vars/"
vars<-raster::brick(paste0(dir_vars,"variables_stack_current.tif"))
gpp<-rast(paste0(dir_vars, "GPP_1km_res.tif"))
elevation<-rast(paste0(dir_vars, "srtm_1km_res.tif"))
hfp<-rast(paste0(dir_vars,"HFP2009ddeg.tif")) #1 km
fractions<-rast("~/Desktop/proyecto_GANSO/variables/fractions_esp_2.5km.tif")
# Step 1: Reproject hfp to match spain_raster's CRS
vars_p <- project(rast(vars), spain_raster)

gpp_p <- project(gpp, spain_raster)
names(gpp_p)<-"gpp"

elevation_p <- project(elevation, spain_raster)
names(elevation_p)<-"elevation"

hfp_p <- project(hfp, spain_raster)
names(hfp_p)<-"hfp"

fractions_p <- project(fractions, spain_raster)
names(fractions_p)<-c("grass", "urban","crops", "shrub","water")

vars_p<-c(vars_p,gpp_p,elevation_p,hfp_p,fractions_p)

# Resample hfp_projected to match the resolution and extent of spain_raster
vars_res <- resample(vars_p, spain_raster, method = "bilinear")  # or "near" for categorical

# plot to check result
plot(vars_res)

vars_res <- mask(vars_res, anyNA(spain_raster), maskvalue=TRUE)
plot(vars_res)

#Function to create virtual species based on the environmental data extracted from
#the actual species. Will create a virtual copy of the species. You can plug a single
#species name or a vector of names and plot both maps next to each other if plot = T
correlation_virtual_species <- function(species_names, plot = TRUE) {
  if (is.character(species_names)) {
    species_names <- as.vector(species_names)
  } else {
    stop("species_names must be a character vector.")
  }
  
  library(terra)
  library(virtualspecies)
  
  results <- data.frame(Species = character(), Correlation = numeric(), stringsAsFactors = FALSE)
  
  for (sp_name in species_names) {
    message("Processing: ", sp_name)
    
    # Step 1: Get species raster (binary presence/absence)
    sp.rast <- get.sp.raster(sp_name)
    
    # Step 2: Extract cells with presence (value == 1)
    cells_with_value_1 <- which(!is.na(values(sp.rast)) & values(sp.rast) == 1)
    
    if (length(cells_with_value_1) == 0) {
      warning("No presence cells for species: ", sp_name)
      results <- rbind(results, data.frame(Species = sp_name, Correlation = NA))
      next
    }
    
    # Step 3: Extract environmental variables at presence locations
    ex <- terra::extract(vars_res, cells_with_value_1, xy = TRUE)
    median_bio <- ex[, c(3:15)]  # Adjust this if your environmental variables are in different columns
    
    # Step 4: Calculate mean and sd for each variable
    medians <- sapply(median_bio, function(x) mean(x, na.rm = TRUE))
    sd <- sapply(median_bio, function(x) sd(x, na.rm = TRUE))
    
    # Step 5: Define parameters for virtual species generation
    my.parameters <- list(
      bio2 = list(fun = "dnorm", args = c(mean = unname(medians[1]), sd = unname(sd[1])^2)),
      bio4 = list(fun = "dnorm", args = c(mean = unname(medians[2]), sd = unname(sd[2])^2)),
      bio5 = list(fun = "dnorm", args = c(mean = unname(medians[3]), sd = unname(sd[3])^2)),
      bio6 = list(fun = "dnorm", args = c(mean = unname(medians[4]), sd = unname(sd[4])^2)),
      bio12 = list(fun = "dnorm", args = c(mean = unname(medians[5]), sd = unname(sd[5])^2)),
      gpp = list(fun = "dnorm", args = c(mean = unname(medians[6]), sd = unname(sd[6])^2)),
      elevation = list(fun = "dnorm", args = c(mean = unname(medians[7]), sd = unname(sd[7])^2)),
      hfp = list(fun = "dnorm", args = c(mean = unname(medians[8]), sd = unname(sd[8])^2)),
      grass = list(fun = "dnorm", args = c(mean = unname(medians[9]), sd = unname(sd[9])^2)),
      urban = list(fun = "dnorm", args = c(mean = unname(medians[10]), sd = unname(sd[10])^2)),
      crops = list(fun = "dnorm", args = c(mean = unname(medians[11]), sd = unname(sd[11])^2)),
      shrub = list(fun = "dnorm", args = c(mean = unname(medians[12]), sd = unname(sd[12])^2)),
      water = list(fun = "dnorm", args = c(mean = unname(medians[13]), sd = unname(sd[13])^2))
      
    )
    
    # Step 6: Generate virtual species suitability map
    vs <- generateSpFromFun(raster.stack = vars_res, parameters = my.parameters, plot = FALSE)
    suitability_map <- vs$suitab.raster
    
    if (plot == TRUE) {
      par(mfrow = c(1, 2))  # Arrange 2 plots side by side
      plot(suitability_map, main = paste("Suitability -", sp_name))
      plot(sp.rast, main = paste("Observed Presence -", sp_name))
      par(mfrow = c(1, 1))  # Reset plotting layout
    }
    
    
    # Step 7: Align rasters (resample to match)
    sp_rast_resampled <- resample(sp.rast, suitability_map, method = "near")
    # Step 8: Extract values and compute correlation
    common_mask <- !is.na(values(sp_rast_resampled)) & !is.na(values(suitability_map))
    pres_values <- values(sp_rast_resampled)[common_mask]
    suit_values <- values(suitability_map)[common_mask]
    
    cor_val <- if (length(pres_values) > 2) {
      cor(pres_values, suit_values, method = "pearson")
    } else {
      NA
    }
    
    # Step 9: Store results
    results <- rbind(results, data.frame(Species = sp_name, Correlation = cor_val))
  }
  
  return(results)
}

results<-correlation_virtual_species(names_sps, plot = F)
