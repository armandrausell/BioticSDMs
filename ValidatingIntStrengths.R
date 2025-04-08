

###

# Get all file paths matching FW_*.csv recursively (all the foodwebs)
fw_files <- list.files(path = "~/Desktop/Research/FoodWebsData", pattern = "^FW_.*\\.csv$", 
                       full.names = TRUE, recursive = TRUE)

# Read each CSV and store in a named list
fw_data_list <- lapply(fw_files, read.csv)

# Name the list elements using the base file names (without extension)
names(fw_data_list) <- tools::file_path_sans_ext(basename(fw_files))

# Example: access FW_005 data
fw_data_list$FW_005

# Apply transformation to all matrices in fw_data_list
fw_data_list <- lapply(fw_data_list, function(mat) {
  # Ensure it's a matrix (in case it was read as a data frame)
  mat <- as.matrix(mat)
  
  # Rename rows and columns generically
  rownames(mat) <- paste0("Prey", seq_len(nrow(mat)))
  colnames(mat) <- paste0("Pred", seq_len(ncol(mat)))
  
  # Transpose so rows are predators and columns are prey
  mat_t <- t(mat)
  
  return(mat_t)
})

library(dplyr)

# Step 1: Create binary version of the webs
fw_binary_list <- lapply(fw_data_list, function(mat) {
  mat_bin <- ifelse(mat > 0, 1, 0)
  return(mat_bin)
})

trophic_position_df <- lapply(fw_binary_list, function(mat) {
  # Get all unique species names from rows and columns
  all_species <- union(rownames(mat), colnames(mat))
  
  # Initialize named vectors
  times_preyed_on <- colSums(mat)
  total_prey <- rowSums(mat)
  
  names(times_preyed_on) <- colnames(mat)
  names(total_prey) <- rownames(mat)
  
  # Create data frame with one row per species
  tp_df <- data.frame(
    Species = all_species,
    Times_Preyed_On = ifelse(all_species %in% names(times_preyed_on),
                             times_preyed_on[all_species], 0),
    Total_Prey = ifelse(all_species %in% names(total_prey),
                        total_prey[all_species], 0)
  )
  
  # Calculate Trophic Position
  tp_df$TrophicPosition <- ((tp_df$Times_Preyed_On + 1) / (tp_df$Total_Prey + 1)) * -1
  
  tp_df$TrophicPosition<-rank(tp_df$TrophicPosition)
  return(tp_df)
})
################
compute_rpp_and_tp <- function(binary_matrix) {
  # Extract all species
  predators <- rownames(binary_matrix)
  preys <- colnames(binary_matrix)
  all_species <- union(predators, preys)
  
  # Calculate preyed on and prey vectors
  times_preyed_on <- colSums(binary_matrix)
  total_prey <- rowSums(binary_matrix)
  names(times_preyed_on) <- colnames(binary_matrix)
  names(total_prey) <- rownames(binary_matrix)
  
  # Trophic position dataframe
  tp_df <- data.frame(
    Species = all_species,
    Times_Preyed_On = ifelse(all_species %in% names(times_preyed_on), times_preyed_on[all_species], 0),
    Total_Prey = ifelse(all_species %in% names(total_prey), total_prey[all_species], 0)
  )
  
  # Trophic position calculation
  tp_df$TrophicPosition <- ((tp_df$Times_Preyed_On + 1) / (tp_df$Total_Prey + 1)) * -1
  tp_df$TrophicPosition <- rank(tp_df$TrophicPosition)  # Optional: rank
  
  # Create lookup vectors
  tp <- setNames(tp_df$TrophicPosition, tp_df$Species)
  total_prey_vec <- setNames(tp_df$Total_Prey, tp_df$Species)
  times_preyed_on_vec <- setNames(tp_df$Times_Preyed_On, tp_df$Species)
  
  # Initialize rPP matrix
  rPP_matrix <- matrix(0, nrow = length(predators), ncol = length(preys),
                       dimnames = list(Predator = predators, Prey = preys))
  
  # Fill rPP values where interaction exists
  for (pred in predators) {
    for (prey in preys) {
      if (binary_matrix[pred, prey] == 1 && pred != prey) {
        rPP_matrix[pred, prey] <- (abs(tp[pred] - tp[prey]) * (1 / (total_prey_vec[pred] + 1))) /
          (times_preyed_on_vec[prey] + 1)
      }
    }
  }
  
  # Normalize rows so each predator's outgoing rPPs sum to 1
  rPP_normalized <- rPP_matrix
  for (pred in rownames(rPP_matrix)) {
    row_sum <- sum(rPP_matrix[pred, ], na.rm = TRUE)
    if (row_sum > 0) {
      rPP_normalized[pred, ] <- rPP_matrix[pred, ] / row_sum
    }
  }
  
  # Return both outputs
  return(list(
    rPP_matrix = rPP_normalized,
    trophic_position_df = tp_df
  ))
}
rPP_outputs <- lapply(fw_binary_list, compute_rpp_and_tp)



# Initialize results table
correlation_results <- data.frame(
  Network = character(),
  Pearson = numeric(),
  Spearman = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each network name in the list
for (name in names(fw_data_list)) {
  original_mat <- fw_data_list[[name]]
  rpp_mat <- rPP_outputs[[name]]$rPP_matrix
  
  # Flatten both matrices to vectors
  original_vec <- as.vector(original_mat)
  rpp_vec <- as.vector(rpp_mat)
  
  # Compare only where binary == 1 (existing links)
  mask <- which(original_vec > 0 & !is.na(rpp_vec))
  
  if (length(mask) > 1) {
    binary_links <- original_vec[mask]
    rpp_links <- rpp_vec[mask]
    
    pearson <- cor(binary_links, rpp_links, method = "pearson")
    spearman <- cor(binary_links, rpp_links, method = "spearman")
  } else {
    pearson <- NA
    spearman <- NA
  }
  
  # Assuming you already have binary_links and rpp_links from the previous code
  
  library(ggplot2)
  
  # Create a data frame for plotting
  plot_data <- data.frame(
    Binary_Links = binary_links,
    rPP_Links = rpp_links
  )
  
  # Plot using ggplot
  g<-ggplot(plot_data, aes(x = Binary_Links, y = rPP_Links)) +
    geom_point(alpha = 0.7, color = "blue") +
    labs(
      title = paste0(name," Original vs predicted interaction strenght"),
      x = "Original strenght",
      y = "Predicted strenght"
    ) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  print(g)
  
  correlation_results <- rbind(correlation_results, data.frame(
    Network = name,
    Pearson = pearson,
    Spearman = spearman
  ))
}


