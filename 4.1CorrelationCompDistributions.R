

library(virtualspecies)
library(raster)
library(terra)


#######

#This function creates a correlation matrix of the distributions of the species among them,
#so how good each distribution correlated with every other distribution
species_correlation_matrix <- function(species_names, method = "jaccard") {
  library(terra)
  
  n <- length(species_names)
  cor_matrix <- matrix(NA, nrow = n, ncol = n)
  rownames(cor_matrix) <- colnames(cor_matrix) <- species_names
  
  # Cache all rasters
  message("Caching species rasters...")
  raster_list <- lapply(species_names, get.sp.raster)
  names(raster_list) <- species_names
  
  # Ensure all rasters are aligned to the first one
  ref_raster <- raster_list[[1]]
  raster_list <- lapply(raster_list, function(r) resample(r, ref_raster, method = "near"))
  
  # Loop through unique pairs and compute similarity
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
library(reshape2)
library(ggplot2)

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

#Based on the competitor index the species had and the amount of distribution similarity,
# a four tile graph is created, where in the top right quadrant we obtain represents
# non-exclusive competition and the bottom represents exclusive competition. Just
# plug the name of a species of interest and see where their competitors fall in the graph.
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
   
  cor_long$Species1 <- trimws(as.character(cor_long$Species1))
  cor_long$Species2 <- trimws(as.character(cor_long$Species2))
  
  # Filter cor_long_clean to only rows involving sp
  sp_rows <- subset(cor_long, Species1 == sp | Species2 == sp)
  
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


tg<-graph_comp_type("Glis glis")%>%print()

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
  cor_long$Species1 <- trimws(as.character(cor_long$Species1))
  cor_long$Species2 <- trimws(as.character(cor_long$Species2))
  
  # Filter for rows involving the species
  sp_rows <- subset(cor_long, Species1 == sp | Species2 == sp)
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

a<-get_exclusive_comp("Glis glis")

#############################
# names_sps <- c("Canis lupus", "Lanius collurio", ...) # define your list of species

#  Returns a data.frame row with collapsed competitors into a string of characters
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

names_gbif_unf<-unique(exclusive_df$Focal_Species)

# Get taxonomic classification to add phylogenetic signal:
library(taxize)

# Check their names
#resolved_names <- gna_verifier(names = names_gbif_unf, best_match_only = TRUE)

# Get classification
#names_gbif<-unique(resolved_names$currentName)
#classification_list <- classification(names_gbif, db = "gbif")

# Try NCBI
classification_list_ncbi <- classification(names_gbif_unf, db = "ncbi")

library(phylotools)
library(vegan)
library(ggtree)
# 2. Convert to Newick format
tax_tree <- class2tree(classification_list_ncbi)
plot(tax_tree$phylo)  # Visualize
names_analyze<-tax_tree$phylo$tip.label
# Plot taxonomic tree
ggtree(tax_tree$phylo) + 
  geom_tiplab() +
  geom_nodelab(aes(label = node), hjust = -0.3) +
  ggtitle("Taxonomic Relationships")


# 3. Extract the tree
tree <- tax_tree$phylo

tax_dist <-tax_tree$distmat%>%as.matrix()
# Heatmap of distances
heatmap(as.matrix(tax_dist), 
        symm = TRUE,
        margins = c(5, 5),
        main = "Taxonomic Distance Matrix")
#Rename rows and columns to ensure potential synonyms are ignored
rownames(tax_dist) <- names_gbif_unf
colnames(tax_dist) <- names_gbif_unf
heatmap(as.matrix(tax_dist), 
        symm = F,
        margins = c(10, 10),
        main = "Taxonomic Distance Matrix")
print(as.matrix(tax_dist))
#######

###########################################################################
## CREATE FUNCTION TO RELATE TAXONOMIC SIGNAL WITH EXCLUSIVE COMPETITION
###########################################################################

#This function calculates how related focal species and competitors are.
get_phylo_sign<-function(species_name){
  
  # Get top competitors
  a <- get_top_competitors(species_name = species_name, n = 5)[["Top_Pressure_Table"]]
  a <- a[, c("Pressure_on_Species", "Pressure_Exerting_Species","Overlap_Proportion")]
  
  # Clean and prep
  sp <- trimws(as.character(species_name))
  a$Pressure_Exerting_Species <- trimws(as.character(a$Pressure_Exerting_Species))
  
  
  #########
  # Load focal species raster
  sp_rast <- get.sp.raster(sp)
  sp_vals <- terra::values(sp_rast)
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
  
  
  # Scale pressure to 0–1
  a$Pressure_on_Species_Scaled <- (a$Pressure_on_Species - min(a$Pressure_on_Species, na.rm = TRUE)) /
    (max(a$Pressure_on_Species, na.rm = TRUE) - min(a$Pressure_on_Species, na.rm = TRUE))
  
  a$CompEx<-a$Pressure_on_Species_Scaled-a$MaxOverlap
  
  comp<-unique(a$Pressure_Exerting_Species)
  

  # After calculating CompEx in dataframe 'a' and before creating results_df:
  
  # Create a named vector of CompEx values for quick lookup
  comp_ex_values <- setNames(a$CompEx, a$Pressure_Exerting_Species)
  
  # Create an empty data frame to store results with additional CompEx column
  results_df <- data.frame(
    species_name = character(),
    competitor = character(),
    TaxonomicDist = numeric(),
    CompEx = numeric(),
    ClassFocal = character(),
    ClassComp = character(),
    stringsAsFactors = FALSE
  )
  
  # Loop through each species and competitor
  for (competitor in comp) {
    if (species_name != competitor) {  # Skip comparing a species to itself
      # Get rows for both species

        # Calculate percentage match (excluding the species column itself if needed)
        # Get taxonomic distance (ensure names match matrix rownames)
        tax_distance <- tax_dist[species_name, competitor]
        
        # Get the CompEx value for this competitor
        comp_ex_value <- ifelse(competitor %in% names(comp_ex_values), 
                                comp_ex_values[competitor], 
                                NA)
        
        class_data <- classification_list_ncbi[[species_name]]
        class_data_comp <- classification_list_ncbi[[competitor]]
        
        # Add to results
        results_df <- rbind(results_df, data.frame(
          species_name = species_name,
          competitor = competitor,
          TaxonomicDist = tax_distance,
          ClassFocal =  class_data$name[class_data$rank == "class"],
          ClassComp = class_data_comp$name[class_data_comp$rank == "class"],
          CompEx = comp_ex_value,
          stringsAsFactors = FALSE
        ))
      
    }
  }
  
  # Now results_df contains all the information together
  
  return(results_df)

}

re<-get_phylo_sign("Accipiter gentilis")

get_phylo_sign_for_all <- function(species_vector) {
  # Initialize an empty list to store results for each species
  all_results <- list()
  
  # Loop through each species in the input vector
  for (species_name in species_vector) {
    # Get the results for this species using your existing function
    species_results <- get_phylo_sign(species_name)
    
    # Only add if we got results
    if (nrow(species_results) > 0) {
      all_results[[species_name]] <- species_results
    }
  }
  
  # Combine all results into a single data frame
  final_results <- do.call(rbind, all_results)
  rownames(final_results) <- NULL  # Clean up row names
  
  return(final_results)
}

# Usage:
all_comparisons <- get_phylo_sign_for_all(names_gbif_unf)

# Subset where ClassFocal == ClassComp
filtered_df <- all_comparisons[all_comparisons$ClassFocal == all_comparisons$ClassComp, ]

 library(ggplot2)
 
 # Basic scatter plot with regression line
 ggplot(filtered_df, aes(x =TaxonomicDist, y = CompEx )) +
   geom_point(alpha = 0.6, color = "steelblue") +  # Semi-transparent points
   geom_smooth(method = "lm",color = "red", se = TRUE) +  # Add regression line with CI
   labs(title = "Correlation between Competitive Exclusion (CompEx) and Taxonomic Similarity",
        x = "Taxonomic Distance",
        y = "Competitive Exclusion Index (CompEx)",
        caption = paste("n =", nrow(all_comparisons), "species pairs")) +
   theme_minimal() +
   facet_wrap(~ClassFocal)
   theme(plot.title = element_text(hjust = 0.5))
 
 
 library(ggpubr)
   # Basic scatter plot with regression line
   ggplot(all_comparisons, aes(x =TaxonomicDist, y = CompEx )) +
     geom_point(alpha = 0.6, color = "steelblue") +  # Semi-transparent points
     geom_smooth(method = "lm",color = "red", se = TRUE) +  # Add regression line with CI
     labs(title = "Correlation between Competitive Exclusion (CompEx) and Taxonomic Similarity",
          x = "Taxonomic Distance",
          y = "Competitive Exclusion Index (CompEx)",
          caption = paste("n =", nrow(all_comparisons), "species pairs")) +
     theme_minimal() +
   theme(plot.title = element_text(hjust = 0.5))
  
   
   
 ggplot(all_comparisons, aes(x = CompEx, y = TaxonomicDist)) +
   geom_hex(bins = 30) +  # Creates hexagonal heatmap
   scale_fill_gradient(low = "blue", high = "red", name = "Count") +
   geom_smooth(method = "lm", color = "white") +
   labs(title = "Density of Species Pairs by CompEx and Taxonomic Similarity",
        x = "Competitive Exclusion Index",
        y = "Percentage Taxonomic Match") +
   theme_dark()
 
 cor_test <- cor.test(all_comparisons$CompEx, all_comparisons$TaxonomicDist, 
                      method = "pearson")
 
 cat(sprintf("Pearson's r = %.2f (p = %.3f)", 
             cor_test$estimate, 
             cor_test$p.value))
#-----------------------------

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


get_rast_shared_sp("Glis glis", "Elyomis quercinus")#create a raster map with shared prey


