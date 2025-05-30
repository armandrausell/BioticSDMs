

library(virtualspecies)
library(raster)
library(terra)


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

# Step 2: Resample hfp_projected to match the resolution and extent of spain_raster
vars_res <- resample(vars_p, spain_raster, method = "bilinear")  # or "near" for categorical

# Optional: plot to check result
plot(vars_res)

vars_res <- mask(vars_res, anyNA(spain_raster), maskvalue=TRUE)
plot(vars_res)

#Function
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

#######

species_correlation_matrix <- function(species_names, method = "jaccard") {
  library(terra)
  
  n <- length(species_names)
  cor_matrix <- matrix(NA, nrow = n, ncol = n)
  rownames(cor_matrix) <- colnames(cor_matrix) <- species_names
  
  # Step 1: Cache all rasters
  message("Caching species rasters...")
  raster_list <- lapply(species_names, get.sp.raster)
  names(raster_list) <- species_names
  
  # Step 2: Ensure all rasters are aligned to the first one
  ref_raster <- raster_list[[1]]
  raster_list <- lapply(raster_list, function(r) resample(r, ref_raster, method = "near"))
  
  # Step 3: Loop through unique pairs and compute similarity
  for (i in 1:n) {
    for (j in i:n) {
      r1 <- values(raster_list[[i]])
      r2 <- values(raster_list[[j]])
      mask <- !is.na(r1) & !is.na(r2)
      x <- r1[mask]
      y <- r2[mask]
      
      if (method == "jaccard") {
        intersection <- sum(x == 1 & y == 1)
        union <- sum(x == 1 | y == 1)
        val <- if (union > 0) intersection / union else NA
      } else if (method == "mcc") {
        TP <- sum(x == 1 & y == 1)
        TN <- sum(x == 0 & y == 0)
        FP <- sum(x == 0 & y == 1)
        FN <- sum(x == 1 & y == 0)
        numerator <- (TP * TN) - (FP * FN)
        denominator <- sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))
        val <- if (denominator > 0) numerator / denominator else NA
      } else {
        stop("Unsupported method. Use 'jaccard' or 'mcc'.")
      }
      
      cor_matrix[i, j] <- cor_matrix[j, i] <- val
    }
  }
  
  return(as.data.frame(cor_matrix))
}
cor_df <- species_correlation_matrix(unique_species, method = "jaccard")
# install.packages("corrplot") if needed
library(corrplot)
library(corrplot)

corrplot(
  as.matrix(cor_df),
  method = "color",
  is.corr = FALSE,
  col = colorRampPalette(c("white", "red"))(200),
  cl.lim = c(0, 1),
  type = "upper",
  order = "hclust",         # <--- groups similar species
  tl.cex = 0.5,             # smaller text
  tl.col = "black",
  diag = F
)

library(reshape2)
library(ggplot2)

library(reshape2)
library(ggplot2)

# Convert correlation matrix to full matrix
cor_mat <- as.matrix(cor_df)

# Perform hierarchical clustering on rows/columns
hc <- hclust(dist(1 - cor_mat))  # 1 - similarity = distance
ordered_species <- rownames(cor_mat)[hc$order]

# Convert to long format for ggplot
cor_long <- melt(cor_mat)
names(cor_long) <- c("Species1", "Species2", "Jaccard")

# Convert Species1 and Species2 to factors with clustering order
cor_long$Species1 <- factor(cor_long$Species1, levels = ordered_species)
cor_long$Species2 <- factor(cor_long$Species2, levels = ordered_species)

# Plot clustered heatmap
g<-ggplot(cor_long, aes(Species1, Species2, fill = Jaccard)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", mid = "blue", high = "red", midpoint = 0.5, limits = c(0,1)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
    axis.text.y = element_text(size = 5)
  )
print(g)
ggsave("Figures/Correlogram_distributions.jpg", plot = g, width = 10, height = 6, dpi = 600)

#######
library(ggrepel)
graph_comp_type <- function(sp) {
  # Get the top competitors
  a <- get_top_competitors(species_name = sp, n = 15)[["Top_Pressure_Table"]]
  a <- a[, c("Pressure_on_Species", "Pressure_Exerting_Species","Overlap_Proportion")]
  
  # Clean species names
  sp <- trimws(as.character(sp))
  a$Pressure_Exerting_Species <- trimws(as.character(a$Pressure_Exerting_Species))
  
  # Load focal species raster
  sp_rast <- get.sp.raster(sp)
  sp_vals <- values(sp_rast)
  sp_presence_cells <- which(!is.na(sp_vals) & sp_vals == 1)
  
  # Total cells where focal species is present
  n_sp_cells <- length(sp_presence_cells)
  
  # Create overlap column for the focal species
  a$FocalOverlap <- sapply(a$Pressure_Exerting_Species, function(comp_sp) {
    comp_rast <- get.sp.raster(comp_sp)
    
    # Extract values for focal species cells
    comp_vals <- values(comp_rast)[sp_presence_cells]
    
    # Count how many of those cells are also 1 in the competitor
    shared_cells <- sum(comp_vals == 1, na.rm = TRUE)
    
    # Proportion of A's range that overlaps with B
    return(shared_cells / n_sp_cells)
  })
  
  #Non focal overlap
  a$NonFocalOverlap <- sapply(a$Pressure_Exerting_Species, function(comp_sp) {
  comp_rast <- get.sp.raster(comp_sp)
  
  # Get indices where B is present
  comp_vals <- values(comp_rast)
  comp_presence_cells <- which(!is.na(comp_vals) & comp_vals == 1)
  n_comp_cells <- length(comp_presence_cells)
  
  if (n_comp_cells == 0) return(NA)  # avoid divide-by-zero
  
  # Resample focal raster to match competitor's if needed
  sp_rast_resampled <- resample(sp_rast, comp_rast, method = "near")
  
  # Extract focal species presence at B’s cells
  sp_vals_at_b <- values(sp_rast_resampled)[comp_presence_cells]
  
  # Count overlapping cells
  shared_cells <- sum(sp_vals_at_b == 1, na.rm = TRUE)
  
  # Proportion of B's range that overlaps with A
  return(shared_cells / n_comp_cells)
})
  
  a$MaxOverlap<-pmax(a$FocalOverlap, a$NonFocalOverlap)
   
  cor_long_clean$Species1 <- trimws(as.character(cor_long_clean$Species1))
  cor_long_clean$Species2 <- trimws(as.character(cor_long_clean$Species2))
  
  # Filter cor_long_clean to only rows involving sp
  sp_rows <- subset(cor_long_clean, Species1 == sp | Species2 == sp)
  
  # Create a normalized column to match competitor names
  sp_rows$Competitor <- ifelse(sp_rows$Species1 == sp, sp_rows$Species2, sp_rows$Species1)
  
  # Match each competitor and pull Jaccard value
  a$DistSimilarity <- sapply(a$Pressure_Exerting_Species, function(comp_sp) {
    jaccard_val <- sp_rows$Jaccard[sp_rows$Competitor == comp_sp]
    if (length(jaccard_val) > 0) return(jaccard_val[1]) else return(NA)
  })
  
  a$Pressure_on_Species_Scaled <- (a$Pressure_on_Species - min(a$Pressure_on_Species, na.rm = TRUE)) /
    (max(a$Pressure_on_Species, na.rm = TRUE) - min(a$Pressure_on_Species, na.rm = TRUE))
  
  
  # Calculate midpoint of y-axis
  y_mid <- (max(a$MaxOverlap, na.rm = TRUE) + min(a$MaxOverlap, na.rm = TRUE)) / 2
  
  ggplot(a, aes(x = Pressure_on_Species_Scaled, y = MaxOverlap, label = Pressure_Exerting_Species)) +
    geom_point(color = "darkred", size = 3) +
    geom_text_repel(size = 3, max.overlaps = 10, box.padding = 0.3) +
    
    # Vertical line at 0.5 (middle of scaled pressure)
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "black") +
    
    # Horizontal line at dynamic midpoint of y-axis
    geom_hline(yintercept = y_mid, linetype = "dashed", color = "black") +
    
    labs(
      title = paste("Competitor Pressure vs. Distribution Similarity for", sp),
      x = "Competitor Pressure (scaled 0–1)",
      y = "Distribution Similarity (Jaccard)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    xlim(0, max(a$Pressure_on_Species_Scaled, na.rm = TRUE) * 1.1)
  
  
}


tg<-graph_comp_type("Lanius collurio")%>%print()

ggsave("Figures/Competition_type_lcollurio.jpg", plot = tg, width = 8, height = 6, dpi = 600)

############

#
get_exclusive_comp <- function(sp_name) {
  # Get top competitors
  a <- get_top_competitors(species_name = sp_name, n = 15)[["Top_Pressure_Table"]]
  a <- a[, c("Pressure_on_Species", "Pressure_Exerting_Species","Overlap_Proportion")]
  
  # Clean and prep
  sp <- trimws(as.character(sp_name))
  a$Pressure_Exerting_Species <- trimws(as.character(a$Pressure_Exerting_Species))
  
  
  #########
  # Load focal species raster
  sp_rast <- get.sp.raster(sp)
  sp_vals <- values(sp_rast)
  sp_presence_cells <- which(!is.na(sp_vals) & sp_vals == 1)
  
  # Total cells where focal species is present
  n_sp_cells <- length(sp_presence_cells)
  
  # Create overlap column for the focal species
  a$FocalOverlap <- sapply(a$Pressure_Exerting_Species, function(comp_sp) {
    comp_rast <- get.sp.raster(comp_sp)
    
    # Extract values for focal species cells
    comp_vals <- values(comp_rast)[sp_presence_cells]
    
    # Count how many of those cells are also 1 in the competitor
    shared_cells <- sum(comp_vals == 1, na.rm = TRUE)
    
    # Proportion of A's range that overlaps with B
    return(shared_cells / n_sp_cells)
  })
  
  #Non focal overlap
  a$NonFocalOverlap <- sapply(a$Pressure_Exerting_Species, function(comp_sp) {
    comp_rast <- get.sp.raster(comp_sp)
    
    # Get indices where B is present
    comp_vals <- values(comp_rast)
    comp_presence_cells <- which(!is.na(comp_vals) & comp_vals == 1)
    n_comp_cells <- length(comp_presence_cells)
    
    if (n_comp_cells == 0) return(NA)  # avoid divide-by-zero
    
    # Resample focal raster to match competitor's if needed
    sp_rast_resampled <- resample(sp_rast, comp_rast, method = "near")
    
    # Extract focal species presence at B’s cells
    sp_vals_at_b <- values(sp_rast_resampled)[comp_presence_cells]
    
    # Count overlapping cells
    shared_cells <- sum(sp_vals_at_b == 1, na.rm = TRUE)
    
    # Proportion of B's range that overlaps with A
    return(shared_cells / n_comp_cells)
  })
  
  #Extract the max overlap, that if high, means that there is no exclusive competition
  a$MaxOverlap<-pmax(a$FocalOverlap, a$NonFocalOverlap)
  
  #########
  cor_long_clean$Species1 <- trimws(as.character(cor_long_clean$Species1))
  cor_long_clean$Species2 <- trimws(as.character(cor_long_clean$Species2))
  
  # Filter for rows involving the species
  sp_rows <- subset(cor_long_clean, Species1 == sp | Species2 == sp)
  sp_rows$Competitor <- ifelse(sp_rows$Species1 == sp, sp_rows$Species2, sp_rows$Species1)
  
  # Append Jaccard similarity
  a$DistSimilarity <- sapply(a$Pressure_Exerting_Species, function(comp_sp) {
    jaccard_val <- sp_rows$Jaccard[sp_rows$Competitor == comp_sp]
    if (length(jaccard_val) > 0) return(jaccard_val[1]) else return(NA)
  })
  
  # Scale pressure to 0–1
  a$Pressure_on_Species_Scaled <- (a$Pressure_on_Species - min(a$Pressure_on_Species, na.rm = TRUE)) /
    (max(a$Pressure_on_Species, na.rm = TRUE) - min(a$Pressure_on_Species, na.rm = TRUE))
  
  # Get midpoint of similarity axis
  y_mid <- (max(a$MaxOverlap, na.rm = TRUE) + min(a$MaxOverlap, na.rm = TRUE)) / 2
  
  # Return names in bottom-right quadrant: high pressure, low similarity
  result <- a$Pressure_Exerting_Species[
    a$Pressure_on_Species_Scaled > 0.5 & a$MaxOverlap < y_mid
  ]
  
  return(result)
}

get_exclusive_comp("Sus scrofa")

#############################
# names_sps <- c("Canis lupus", "Lanius collurio", ...) # define your list of species

# Wrapper to return a data.frame row with collapsed competitors into a string of characters
get_exclusive_comp_row <- function(sp_name) {
  comps <- get_exclusive_comp(sp_name)
  data.frame(
    Focal_Species = sp_name,
    Exclusive_Competitors = paste(comps, collapse = ", "),
    stringsAsFactors = FALSE
  )
}

# Apply over all species names
exclusive_df <- do.call(rbind, lapply(unique_species, get_exclusive_comp_row))


######
get_rast_shared_sp <- function(compA, compB) {
  library(terra)
  library(dplyr)
  
  # Step 1: Get unique prey species
  b <- get_prey_species(compA)
  c <- get_prey_species(compB)
  bc <- c(b, c) %>% unique()
  
  # Step 2: Load rasters
  raster_list <- lapply(bc, function(sp_name) {
    tryCatch({
      get.sp.raster(sp_name)
    }, error = function(e) {
      warning(paste("Failed to get raster for", sp_name, ":", e$message))
      return(NULL)
    })
  })
  
  # Step 3: Clean out failed rasters
  raster_list <- raster_list[!sapply(raster_list, is.null)]
  
  if (length(raster_list) == 0) {
    warning("No valid rasters to sum.")
    return(NULL)
  }
  
  # Step 4: Stack and sum all rasters
  shared_stack <- rast(raster_list)
  summed_raster <- app(shared_stack, sum, na.rm = TRUE)%>%plot()
  
  return(summed_raster)
}


get_rast_shared_sp("Milvus migrans", "Felis silvestris")#create a raster map with shared prey


