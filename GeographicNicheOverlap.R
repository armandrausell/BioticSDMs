
#Geographic niche overlap

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

vars_res <- terra::mask(vars_res, anyNA(spain_raster), maskvalue=T)
plot(vars_res)


###############################
#message("Processing: ", sp_name)

# Step 1: Get species raster (binary presence/absence)
sp1.rast <- get.sp.raster(sp1_name)
sp2.rast<-  get.sp.raster(sp2_name)

# Step 2: Extract cells with presence (value == 1)
cells_sp1_value_1 <- which(!is.na(terra::values(sp1.rast)) & terra::values(sp1.rast) == 1)
cells_sp2_value_1 <- which(!is.na(terra::values(sp2.rast)) & terra::values(sp2.rast) == 1)

# Step 3: Extract environmental variables at presence locations
ex1 <- terra::extract(vars_res, cells_sp1_value_1, xy = TRUE)
ex2 <- terra::extract(vars_res, cells_sp2_value_1, xy = TRUE)

env1 <- ex1[, !(names(ex1) %in% c("ID", "x", "y"))]
env2 <- ex2[, !(names(ex2) %in% c("ID", "x", "y"))]

shared_ranges <- data.frame()

for (var in colnames(env1)) {
  min_shared <- max(min(env1[[var]], na.rm = TRUE), min(env2[[var]], na.rm = TRUE))
  max_shared <- min(max(env1[[var]], na.rm = TRUE), max(env2[[var]], na.rm = TRUE))
  shared_ranges <- rbind(shared_ranges, data.frame(Variable = var, Min = min_shared, Max = max_shared))
}

# Mask of shared environmental space (initially all TRUE)
shared_env_mask <- rep(TRUE, ncell(vars_res))

# Apply each environmental constraint
for (i in 1:nrow(shared_ranges)) {
  var <- shared_ranges$Variable[i]
  min_val <- shared_ranges$Min[i]
  max_val <- shared_ranges$Max[i]
  
  v <- terra::values(vars_res[[var]])
  shared_env_mask <- shared_env_mask & (v >= min_val & v <= max_val)
}

shared_env_raster <- vars_res[[1]]
terra::values(shared_env_raster) <- as.numeric(shared_env_mask)

plot(shared_env_raster, main = "Shared Environmental Space")
plot(sp1.rast)
plot(sp2.rast)

####

library(terra)

# Initialize an empty list to hold binary rasters for each variable
binary_layers <- list()

# Loop through each variable to create a binary layer
for (i in 1:nrow(shared_ranges)) {
  var <- shared_ranges$Variable[i]
  min_val <- shared_ranges$Min[i]
  max_val <- shared_ranges$Max[i]
  
  # Extract raster layer
  layer <- vars_res[[var]]
  
  # Create binary raster: 1 if within range, 0 otherwise
  binary_layer <- classify(layer, rcl = matrix(c(-Inf, min_val, 0,
                                                 min_val, max_val, 1,
                                                 max_val, Inf, 0), ncol = 3, byrow = TRUE))
  
  binary_layers[[var]] <- binary_layer
}

# Stack the binary rasters
binary_stack <- rast(binary_layers)

# Sum across layers to count how many conditions are met at each pixel
sum_binary <- sum(binary_stack)

# Create final raster: 1 if all variables = 1, otherwise 0
final_shared_raster <- classify(sum_binary, matrix(c(-Inf, nrow(shared_ranges) - 0.1, 0,
                                                     nrow(shared_ranges) - 0.1, Inf, 1), ncol = 3, byrow = TRUE))

# Plot final map
plot(final_shared_raster, main = "Shared Environmental Niche (All Variables)", col = c("grey80", "darkgreen"), legend = FALSE)


