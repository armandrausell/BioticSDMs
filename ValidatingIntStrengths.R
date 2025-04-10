

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
  predators <- rownames(binary_matrix)
  preys <- colnames(binary_matrix)
  all_species <- union(predators, preys)
  
  # Step 1: Compute Times_Preyed_On and Total_Prey
  times_preyed_on <- colSums(binary_matrix)
  total_prey <- rowSums(binary_matrix)
  names(times_preyed_on) <- colnames(binary_matrix)
  names(total_prey) <- rownames(binary_matrix)
  
  # Step 2: Compute Trophic Position
  tp_df <- data.frame(
    Species = all_species,
    Times_Preyed_On = ifelse(all_species %in% names(times_preyed_on), times_preyed_on[all_species], 0),
    Total_Prey = ifelse(all_species %in% names(total_prey), total_prey[all_species], 0)
  )
  tp_df$TrophicPosition <- ((tp_df$Times_Preyed_On + 1) / (tp_df$Total_Prey + 1)) * -1
  tp_df$TrophicPosition <- rank(tp_df$TrophicPosition)
  
  # Lookup vectors
  tp <- setNames(tp_df$TrophicPosition, tp_df$Species)
  total_prey_vec <- setNames(tp_df$Total_Prey, tp_df$Species)
  times_preyed_on_vec <- setNames(tp_df$Times_Preyed_On, tp_df$Species)
  
  # Step 3: Median TP per predator
  median_tp <- sapply(predators, function(pred) {
    prey_names <- preys[binary_matrix[pred, ] == 1]
    if (length(prey_names) == 0) return(NA)
    median(tp[prey_names], na.rm = TRUE)
  })
  
  # Step 4: Max deviation per predator
  max_dev <- sapply(predators, function(pred) {
    prey_names <- preys[binary_matrix[pred, ] == 1]
    if (length(prey_names) == 0) return(0)
    max(abs(tp[prey_names] - median_tp[pred]), na.rm = TRUE)
  })
  
  # Step 5: Compute raw rPP and TCP matrices
  rPP_matrix_raw <- matrix(0, nrow = length(predators), ncol = length(preys),
                           dimnames = list(Predator = predators, Prey = preys))
  tcp_matrix <- matrix(0, nrow = length(predators), ncol = length(preys),
                       dimnames = list(Predator = predators, Prey = preys))
  
  for (pred in predators) {
    for (prey in preys) {
      if (binary_matrix[pred, prey] == 1 && pred != prey) {
        abs_dev <- abs(tp[prey] - median_tp[pred])
        tcp <- if (max_dev[pred] > 0) 1 - (abs_dev / max_dev[pred]) else 1
        tcp_matrix[pred, prey] <- tcp
        
        rPP_matrix_raw[pred, prey] <- (
          log2(abs(tp[pred] - tp[prey]) + 1) *
            (1 / (1 + total_prey_vec[pred] * (1 - tcp)))
        ) / (times_preyed_on_vec[prey] + 1)
      }
    }
  }
  
  # Step 6: Normalize rPP per predator (row)
  rPP_matrix_norm <- rPP_matrix_raw
  for (pred in rownames(rPP_matrix_raw)) {
    row_sum <- sum(rPP_matrix_raw[pred, ], na.rm = TRUE)
    if (row_sum > 0) {
      rPP_matrix_norm[pred, ] <- rPP_matrix_raw[pred, ] / row_sum
    }
  }
  
  # Step 7: Create interaction table with both raw and normalized rPP
  interaction_table <- data.frame()
  
  for (pred in predators) {
    for (prey in preys) {
      if (binary_matrix[pred, prey] == 1 && pred != prey) {
        interaction_table <- rbind(interaction_table, data.frame(
          Predator = pred,
          Prey = prey,
          TrophicPositionPred = tp[pred],
          TrophicPositionPrey = tp[prey],
          TotalPrey_Pred = total_prey_vec[pred],
          Times_Prey_Preyed = times_preyed_on_vec[prey],
          TrophicCoreProximity = tcp_matrix[pred, prey],
          rPP_raw = rPP_matrix_raw[pred, prey],
          rPP = rPP_matrix_norm[pred, prey]
        ))
      }
    }
  }
  
  return(list(
    rPP_matrix = rPP_matrix_norm,
    TP = tp_df,
    InteractionTable = interaction_table
  ))
}
# Apply the function
rPP_outputs <- lapply(fw_binary_list, compute_rpp_and_tp)

matrix_pairs_list_df <- list()  # NEW list to store original & rPP matrices
matrix_pairs_list <- list()  # NEW list to store original & rPP matrices

correlation_results <- data.frame(
  Network = character(),
  Mean_Spearman = numeric(),
  Median_Spearman = numeric(),
  Links = integer(),
  stringsAsFactors = FALSE
)

for (name in names(fw_data_list)) {
  original_mat <- fw_data_list[[name]]
  rpp_mat <- rPP_outputs[[name]]$rPP_matrix
  
  prey_species <- colnames(original_mat)
  predator_species <- rownames(original_mat)
  
  spearman_values <- c()
  
  for (prey in prey_species) {
    original_col <- original_mat[, prey]
    rpp_col <- rpp_mat[, prey]
    
    # Only keep predators that have interaction in original
    valid_mask <- which(original_col > 0 & !is.na(rpp_col))
    
    if (length(valid_mask) > 1) {
      orig_vals <- original_col[valid_mask]
      rpp_vals <- rpp_col[valid_mask]
      
      corr_val <- suppressWarnings(cor(orig_vals, rpp_vals, method = "spearman"))
      spearman_values <- c(spearman_values, corr_val)
    }
  }
  
  # Summary stats
  mean_corr <- mean(spearman_values, na.rm = TRUE)
  median_corr <- median(spearman_values, na.rm = TRUE)
  
  correlation_results <- rbind(
    correlation_results,
    data.frame(
      Network = name,
      Mean_Spearman = mean_corr,
      Median_Spearman = median_corr,
      Links = nrow(original_mat)
    )
  )
}

###############################
# ASSESS DIFFERENCES IN CALCULUS
###############################
difference_list <- list()  # New list to store the result per network

for (name in names(matrix_pairs_list)) {
  original_df <- matrix_pairs_list[[name]]$Original
  rpp_df <- matrix_pairs_list[[name]]$rPP
  
  # Ensure both have the same row and column names
  if (!all(rownames(original_df) == rownames(rpp_df)) ||
      !all(colnames(original_df) == colnames(rpp_df))) {
    warning(paste("Mismatch in row/col names for:", name))
    next
  }
  
  # Convert to matrices for element-wise subtraction
  original_mat <- as.matrix(original_df)
  rpp_mat <- as.matrix(rpp_df)
  
  # Calculate difference: rPP - Original
  diff_mat <- rpp_mat - original_mat
  
  # Convert to long format
  diff_long <- as.data.frame(as.table(diff_mat)) %>%
    rename(Predator = Var1, Prey = Var2, Difference = Freq) %>%
    mutate(Network = name)
  
  difference_list[[name]] <- diff_long
}

# Combine all into a single dataframe
difference_df <- bind_rows(difference_list)

# Optional: arrange by absolute difference
difference_df <- difference_df %>%
  arrange(desc(abs(Difference)))
