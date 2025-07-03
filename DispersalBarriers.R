

library(terra)
library(ggplot2)
library(tidyr)
library(ggplot2)
library(viridis)
library(pheatmap)

#We'll be using the high resolution rasters
# Load input rasters
vars       <- rast("~/Desktop/vars/variables_stack_current.tif")
gpp        <- rast("~/Desktop/vars/GPP_1km_res.tif")
names(gpp)<-"gpp"
elevation  <- rast("~/Desktop/vars/srtm_1km_res.tif")
names(elevation)<-"elevation"
hfp        <- rast("~/Desktop/vars/HFP2009ddeg.tif")
names(hfp)<-"hfp"

# Load and clean fractions
fractions  <- rast("~/Desktop/proyecto_GANSO/variables/fractions_esp_2.5km.tif")
fractions  <- fractions[[-c(2)]]  # Remove second layer
names(fractions) <- c("grass", "crops", "shrub", "water")

# Load Spain raster and project to CRS of vars
spain_raster       <- rast("spain_rasterUTM.tif")
spain_raster_proj  <- project(spain_raster, crs(vars))
plot(spain_raster_proj)
# Crop extent from projected Spain raster
spain_ext_proj     <- ext(spain_raster_proj)

# Crop all environmental rasters to that extent
vars_crop      <- crop(vars, spain_ext_proj)
gpp_crop       <- crop(gpp, spain_ext_proj)
elevation_crop <- crop(elevation, spain_ext_proj)
hfp_crop       <- crop(hfp, spain_ext_proj)
fractions_crop <- crop(fractions, spain_ext_proj)
plot(vars_crop)

# Get the CRS string of spain_raster
target_crs <- crs(spain_raster)

# Reproject while keeping original resolution (do NOT use spain_raster as target)
vars_reproj      <- project(vars_crop, target_crs, method = "bilinear")
gpp_reproj       <- project(gpp_crop, target_crs, method = "bilinear")
elevation_reproj <- project(elevation_crop, target_crs, method = "bilinear")
hfp_reproj       <- project(hfp_crop, target_crs, method = "bilinear")
fractions_reproj <- project(fractions_crop, target_crs, method = "bilinear")

# Resample all layers to match resolution/grid of gpp
ref <- gpp_reproj
vars_rs      <- resample(vars_reproj, ref)
elevation_rs <- resample(elevation_reproj, ref)
hfp_rs       <- resample(hfp_reproj, ref)
fractions_rs <- resample(fractions_reproj, ref)

# Combine all into one raster stack
vars_highres_stack <- c(vars_rs, gpp = ref, elevation_rs, hfp_rs, fractions_rs)

# Optional: mask by Spain
plot(vars_highres_stack[[1]])


sp.rast <- resample(spain_raster, vars_highres_stack, method = "near")
vars_highres_stack <- terra::mask(vars_highres_stack, anyNA(sp.rast), maskvalue=T)
#########################################
#load environmental variables with matching resolution to atlas
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
#########################################

get_friction_map <- function(sp_name, focal = 3){
# Load species binary raster (0/1)
sp1.rast <- get.sp.raster(sp_name)
#sp1.rast <- resample(sp1.rast, vars_res, method = "near")

plot(sp1.rast)

# Define 3x3 neighborhood weights
w <- matrix(1, nrow = focal, ncol = focal)

# Identify edge presence cells: 1s with at least one 0 neighbor
absence_neighbors <- focal(sp1.rast == 0, w = w, fun = sum, na.policy = "omit")
plot(absence_neighbors)
edge_pres <- (sp1.rast == 1) & (absence_neighbors > 0)
plot(edge_pres, main = "Edge Presence Cells")

# Identify edge absence cells: 0s with at least one 1 neighbor
presence_neighbors <- focal(sp1.rast == 1, w = w, fun = sum, na.policy = "omit")
edge_abs <- (sp1.rast == 0) & (presence_neighbors > 0)
plot(edge_abs, main = "Edge Absence Cells")

# Extract environmental values
edge_vals <- extract(vars_res, which(values(edge_pres) == 1))
absence_vals <- extract(vars_res, which(values(edge_abs) == 1))

# Convert to data frames
edge_df <- as.data.frame(edge_vals)
absence_df <- as.data.frame(absence_vals)

# Summary statistics
summary(edge_df)
summary(absence_df)

# Prepare long format for plotting
edge_long <- pivot_longer(edge_df, everything(), names_to = "Variable", values_to = "Value")
absence_long <- pivot_longer(absence_df, everything(), names_to = "Variable", values_to = "Value")

edge_long$Group <- "EdgePresence"
absence_long$Group <- "EdgeAbsence"

combined <- rbind(edge_long, absence_long)

# Plot density comparisons
ggplot(combined, aes(x = Value, fill = Group)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Variable, scales = "free") +
  theme_minimal()

# One test per variable
pvals <- mapply(function(x, y) {
  if (length(unique(x)) <= 2 | length(unique(y)) <= 2) return(NA)  # skip binary or nearly constant vars
  t.test(x, y)$p.value
}, edge_df, absence_df)

# Adjust for multiple testing (e.g., Benjamini-Hochberg FDR)
pvals_adj <- p.adjust(pvals, method = "BH")

# Combine into a result table
comp_results <- data.frame(
  Variable = colnames(edge_df),
  p_value = pvals,
  p_adj = pvals_adj
)

comp_results <- comp_results[order(comp_results$p_adj), ]
print(comp_results)


# Cohen's d function
cohen_d <- function(x, y) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  n1 <- length(x)
  n2 <- length(y)
  if (n1 < 2 || n2 < 2) return(NA)  # not enough values to compute variance
  pooled_sd <- sqrt(((n1 - 1) * var(x) + (n2 - 1) * var(y)) / (n1 + n2 - 2))
  if (pooled_sd == 0 || is.na(pooled_sd)) return(NA)
  d <- (mean(x) - mean(y)) / pooled_sd
  return(d)
}

# Calculate Cohen’s d per variable
effect_sizes <- mapply(cohen_d, edge_df, absence_df)

comp_results$cohen_d <- effect_sizes
print(comp_results)

plot(ggplot(comp_results, aes(x = reorder(Variable, cohen_d), y = cohen_d, fill = p_adj < 0.05)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "gray80")) +
  labs(y = "Cohen's d", x = "Variable", fill = "Significant") +
  theme_minimal())

##############
# Keep only significant variables (adjusted p < 0.05)
sig_vars <- comp_results[comp_results$p_adj < 0.05, ]

# Use signed Cohen's d as weights
barrier_vars <- sig_vars$Variable
weights <- sig_vars$cohen_d
names(weights) <- barrier_vars

# Subset the environmental raster stack to just the significant variables
vars_sub <- vars_res[[barrier_vars]]

# Create an empty list to store scaled rasters
scaled_list <- list()

# Loop through each variable layer and rescale to [0, 1]
for (i in 1:nlyr(vars_sub)) {
  lyr <- vars_sub[[i]]
  lyr_min <- as.numeric(global(lyr, "min", na.rm = TRUE))
  lyr_max <- as.numeric(global(lyr, "max", na.rm = TRUE))
  scaled <- (lyr - lyr_min) / (lyr_max - lyr_min)
  scaled_list[[i]] <- scaled
}

# Combine into a single SpatRaster
vars_scaled <- rast(scaled_list)
names(vars_scaled) <- barrier_vars  # re-assign names just in case

# Apply weights (element-wise multiplication)
for (v in barrier_vars) {
  vars_scaled[[v]] <- vars_scaled[[v]] * weights[v]
}
##########################
### PLOT ###
# Sum the weighted rasters to get friction map
friction_map <- sum(vars_scaled)

friction_norm <- (friction_map - minmax(friction_map)[1]) / 
  (minmax(friction_map)[2] - minmax(friction_map)[1])


par(mfrow = c(1, 2))
plot(sp1.rast,main=paste0(sp_name),col = c("gray", "blue"), legend=F)
plot(friction_norm, main = "Normalized Friction Map",
     col = rev(magma(100)))
par(mfrow = c(1, 1))

return(list(
  friction_norm = friction_norm,
  comp_results = comp_results)
)
}


##HILLSHADE RAST
get_hillshade<-function(sp, focal = 3, overlaysp = T){
a<-get_friction_map(sp, focal = focal) #always uneven number
sp1.rast <- get.sp.raster(sp)
sp1.rast <- resample(sp1.rast, vars_res, method = "near")
sp1.rast[sp1.rast != 1] <- NA

b<-a$friction_norm
c<-a$comp_results
library(terra)

# Step 1: Apply exaggeration to elevation (NOT slope)
exaggeration_factor <- 10
b_exaggerated <- exp(b) * exaggeration_factor

# Step 2: Compute slope and aspect from the exaggerated elevation
slope <- terrain(b_exaggerated, "slope", unit = "radians")
aspect <- terrain(b_exaggerated, "aspect", unit = "radians")

# Step 3: Create hillshade with appropriate sun settings
hillshade <- shade(slope, aspect, angle = 180, direction = 330)

# Correctly extract scalar values from the global() result
b_min <- global(b, "min", na.rm = TRUE)[1, 1]
b_max <- global(b, "max", na.rm = TRUE)[1, 1]
# Now normalize between 0 and 1
b_norm <- (b - b_min) / (b_max - b_min)

plot(hillshade, col = rev(gray.colors(100)), main = "Raster 'b' over Hillshade", legend = FALSE)
plot(b_norm, col = rev(viridis::inferno(100)), alpha = 0.6, add = TRUE)
if(overlaysp ==TRUE){plot(sp1.rast, col ="blue", alpha=0.2, add=T)}
##########
#plot high res map infered from lower res map

# Keep only significant variables (adjusted p < 0.05)
sig_vars <- c[c$p_adj < 0.05, ]

# Use signed Cohen's d as weights
barrier_vars <- sig_vars$Variable
weights <- sig_vars$cohen_d
names(weights) <- barrier_vars

# Subset the environmental raster stack to just the significant variables
vars_sub <- vars_highres_stack[[barrier_vars]]

# Create an empty list to store scaled rasters
scaled_list <- list()

# Loop through each variable layer and rescale to [0, 1]
for (i in 1:nlyr(vars_sub)) {
  lyr <- vars_sub[[i]]
  lyr_min <- as.numeric(global(lyr, "min", na.rm = TRUE))
  lyr_max <- as.numeric(global(lyr, "max", na.rm = TRUE))
  scaled <- (lyr - lyr_min) / (lyr_max - lyr_min)
  scaled_list[[i]] <- scaled
}

# Combine into a single SpatRaster
vars_scaled <- rast(scaled_list)
names(vars_scaled) <- barrier_vars  # re-assign names just in case

# Apply weights (element-wise multiplication)
for (v in barrier_vars) {
  vars_scaled[[v]] <- vars_scaled[[v]] * weights[v]
}
##########################
### PLOT ###
# Sum the weighted rasters to get friction map
friction_map <- sum(vars_scaled)
plot(friction_map)


# Step 1: Apply exaggeration to elevation (NOT slope)
exaggeration_factor <- 10
b_exaggerated <- exp(friction_map) * exaggeration_factor

# Step 2: Compute slope and aspect from the exaggerated elevation
slope <- terrain(b_exaggerated, "slope", unit = "radians")
aspect <- terrain(b_exaggerated, "aspect", unit = "radians")

# Step 3: Create hillshade with appropriate sun settings
hillshade <- shade(slope, aspect, angle = 180, direction = 330)

# Correctly extract scalar values from the global() result
b_min <- global(friction_map, "min", na.rm = TRUE)[1, 1]
b_max <- global(friction_map, "max", na.rm = TRUE)[1, 1]
# Now normalize between 0 and 1
b_norm <- (friction_map - b_min) / (b_max - b_min)

plot(hillshade, col = rev(gray.colors(100)), main = "Raster 'b' over Hillshade", legend = FALSE)
plot(b_norm, col = rev(viridis::inferno(100)), alpha = 0.7, add = TRUE)
if(overlaysp ==TRUE){plot(sp1.rast, col ="blue", alpha=0.2, add=T)}

}
#This function computes the map of the friction map in high resolution

get_hillshade("Mustela erminea", overlaysp = T)
##############
# OVERALL ANALYSIS 
#############
# this creates a single map of connectivity for everyone
# Store results
results_list <- list()
species_done <- c()

# Loop over all species
for (sp in unique_species) {
  cat("Running:", sp, "\n")
  try({
    res <- get_friction_map(sp)
    results_list[[sp]] <- res
    species_done <- c(species_done, sp)
  }, silent = TRUE)
}

# Stack all friction rasters
friction_stack <- rast(lapply(results_list, function(x) x$friction_norm))
names(friction_stack) <- species_done
global_friction<- sum(friction_stack)
plot(global_friction)

# Combine comp_results into one data frame
library(dplyr)
all_comp_results <- bind_rows(
  lapply(results_list, function(x) x$comp_results),
  .id = "Species"
)

sig_vars_summary <- all_comp_results %>%
  filter(p_adj < 0.05 & !is.na(cohen_d)) %>%
  group_by(Variable) %>%
  summarise(
    n_species = n(),
    n_positive = sum(cohen_d > 0),
    n_negative = sum(cohen_d < 0),
    mean_d = mean(cohen_d),
    mean_abs_d = mean(abs(cohen_d))
  ) %>%
  arrange(desc(n_species))

ggplot(sig_vars_summary, aes(x = reorder(Variable, mean_d), y = mean_d, fill = mean_d > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "tomato")) +
  labs(
    title = "Mean Directional Effect Size of Variables",
    x = "Environmental Variable",
    y = "Mean Cohen's d",
    fill = "Effect Direction"
  ) +
  theme_minimal()

library(dplyr)

heatmap_data <- all_comp_results %>%
  dplyr::filter(p_adj < 0.05 & !is.na(cohen_d)) %>%
  dplyr::select(Species, Variable, cohen_d)

library(ggplot2)

#order data
# Reorder variables by mean effect size (signed)
var_order <- heatmap_data %>%
  group_by(Variable) %>%
  summarise(mean_d = mean(cohen_d, na.rm = TRUE)) %>%
  arrange(mean_d) %>%
  pull(Variable)

# Reorder species by their total average effect (optional)
species_order <- heatmap_data %>%
  group_by(Species) %>%
  summarise(mean_d = mean(cohen_d, na.rm = TRUE)) %>%
  arrange(mean_d) %>%
  pull(Species)

# Set factor levels for ordering
heatmap_data$Variable <- factor(heatmap_data$Variable, levels = var_order)
heatmap_data$Species <- factor(heatmap_data$Species, levels = species_order)


heatmap_matrix <- heatmap_data %>%
  tidyr::pivot_wider(names_from = Variable, values_from = cohen_d) %>%
  tibble::column_to_rownames("Species") %>%
  as.matrix()

# Optional: Fill NAs with 0 (or use another imputation method)
heatmap_matrix[is.na(heatmap_matrix)] <- 0

# Compute similarity (correlation) between species' effect profiles
cor_mat <- cor(t(heatmap_matrix), method = "pearson")  # species as rows, variables as columns

# Hierarchical clustering using 1 - correlation distance
hc <- hclust(dist(1 - cor_mat))

# Get ordered species names from the clustering
ordered_species <- rownames(cor_mat)[hc$order]

library(ggplot2)

# Convert to long format
plot_data <- heatmap_data %>%
  filter(Species %in% ordered_species)  # keep only those in clustering

# Set the species factor with clustering order
plot_data$Species <- factor(plot_data$Species, levels = ordered_species)

ggplot(plot_data, aes(x = Variable, y = Species, fill = cohen_d)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Cohen's d (Species Clustered by Effect Profile)",
    x = "Environmental Variable",
    y = "Species"
  )
