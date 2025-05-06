

###

# Get all file paths matching FW_*.csv recursively (all the foodwebs)
fw_files <- list.files(path = "~/Desktop/Research/FoodWebsData/WebsWithSpNames", pattern = "^FW_.*\\.csv$", 
                       full.names = TRUE, recursive = TRUE)

# Read each CSV and store in a named list
fw_data_list <- lapply(fw_files, read.csv)

fw_data_list <- lapply(fw_data_list, function(mat) {
  # Ensure it's a data frame
  df <- as.data.frame(mat)
  
  # Set first column as rownames (species names)
  rownames(df) <- df[[1]]
  
  # Remove first column
  df <- df[, -1, drop = FALSE]
  
  # Clean row and column names: replace spaces or dots with underscores
  rownames(df) <- gsub("[ .]", "_", rownames(df))
  colnames(df) <- gsub("[ .]", "_", colnames(df))
  
  return(df)
})


fw_data_list <- lapply(fw_data_list, function(df) {
  # Ensure df is a data.frame
  df <- as.data.frame(df)
  
  # Get all unique species across rows and columns
  all_species <- union(rownames(df), colnames(df))
  
  # Add missing rows
  missing_rows <- setdiff(all_species, rownames(df))
  if (length(missing_rows) > 0) {
    empty_rows <- matrix(0, nrow = length(missing_rows), ncol = ncol(df))
    rownames(empty_rows) <- missing_rows
    colnames(empty_rows) <- colnames(df)
    df <- rbind(df, empty_rows)
  }
  
  # Add missing columns
  missing_cols <- setdiff(all_species, colnames(df))
  if (length(missing_cols) > 0) {
    empty_cols <- matrix(0, nrow = nrow(df), ncol = length(missing_cols))
    colnames(empty_cols) <- missing_cols
    rownames(empty_cols) <- rownames(df)
    df <- cbind(df, empty_cols)
  }
  
  # Reorder rows and columns to the same species order
  df <- df[all_species, all_species]
  
  return(df)
})


# Name the list elements using the base file names (without extension)
names(fw_data_list) <- tools::file_path_sans_ext(basename(fw_files))

# Example: access FW_005 data
view(fw_data_list$FW_005)

# Apply transformation to all matrices in fw_data_list
fw_data_list <- lapply(fw_data_list, function(mat) {
  # Ensure it's a matrix (in case it was read as a data frame)
  mat <- as.matrix(mat)

  # Transpose so rows are predators and columns are prey
  mat_t <- t(mat)
  
  return(mat_t)
})

library(dplyr)

# Create binary version of the webs
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
  
  #tp_df$TrophicPosition<-rank(tp_df$TrophicPosition)
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
  
  # Median TP per predator
  median_tp <- sapply(predators, function(pred) {
    prey_names <- preys[binary_matrix[pred, ] == 1]
    if (length(prey_names) == 0) return(NA)
    median(tp[prey_names], na.rm = TRUE)
  })
  
  # Median TP per predator
  mean_tp <- sapply(predators, function(pred) {
    prey_names <- preys[binary_matrix[pred, ] == 1]
    if (length(prey_names) == 0) return(NA)
    mean(tp[prey_names], na.rm = TRUE)
  })
  
  # Max deviation per predator
  max_dev <- sapply(predators, function(pred) {
    prey_names <- preys[binary_matrix[pred, ] == 1]
    if (length(prey_names) == 0) return(0)
    max(abs(tp[prey_names] - median_tp[pred]), na.rm = TRUE)
  })
  
  # Median TP of predators for each prey
  median_ptp <- sapply(preys, function(prey) {
    predator_names <- rownames(binary_matrix)[binary_matrix[, prey] == 1]
    if (length(predator_names) == 0) return(NA)
    median(tp[predator_names], na.rm = TRUE)
  })
  
  # Compute max deviation from predator median per prey
  max_dev_ptp <- sapply(preys, function(prey) {
    predator_names <- rownames(binary_matrix)[binary_matrix[, prey] == 1]
    if (length(predator_names) == 0) return(0)
    max(abs(tp[predator_names] - median_ptp[prey]), na.rm = TRUE)
  })
  
  # Compute standard deviation of prey TP per predator
  sd_tp <- sapply(predators, function(pred) {
    prey_names <- preys[binary_matrix[pred, ] == 1]
    if (length(prey_names) == 0) return(NA)
    sd(tp[prey_names], na.rm = TRUE)
  })
  
  # Compute raw rPP and TCP matrices
  rPP_matrix_raw <- matrix(0, nrow = length(predators), ncol = length(preys),
                           dimnames = list(Predator = predators, Prey = preys))
  tcp_matrix <- matrix(0, nrow = length(predators), ncol = length(preys),
                       dimnames = list(Predator = predators, Prey = preys))
  
  for (pred in predators) {
    prey_names <- preys[binary_matrix[pred, ] == 1]
    tp_sd <- sd(tp[prey_names], na.rm = TRUE)  # SD of prey TP for this predator
     for (prey in preys) {
      if (binary_matrix[pred, prey] == 1 && pred != prey) {
        abs_dev <- abs(tp[prey] - median_tp[pred])
        
        # Smooth TCP scaling
        alpha <- 0.2
        beta <- 10
        tcp_original <- if (max_dev[pred] > 0) 1 - (abs_dev / max_dev[pred]) else 1
        scaling_factor <- 1 / (1 + exp(-alpha * (total_prey_vec[pred] - beta)))
        tcp <- 1 - ((1 - tcp_original) * scaling_factor)
        tcp_matrix[pred, prey] <- tcp
        
        # Modulate by similarity to median predator TP for prey
        dev_from_ptp <- abs(tp[pred] - median_ptp[prey])
        max_dev_val <- max_dev_ptp[prey]
        if (is.na(dev_from_ptp) || is.na(max_dev_val) || max_dev_val == 0) {
          pred_similarity <- 1  # no modulation
        } else {
          pred_similarity <- 1 - (dev_from_ptp / max_dev_val)
        }
        
        # Final rPP formula
        #rPP_matrix_raw[pred, prey] <- (
         # log2(abs(tp[pred] - tp[prey]) + 1) *
          #  (1 / (1 + total_prey_vec[pred] * (1 - tcp)))
        #) / (times_preyed_on_vec[prey] + 1)
        
        # Final rPP formula
        #rPP_matrix_raw[pred, prey] <- (
        #  (exp(-(((abs(tp[pred] - tp[prey]) - (tp[pred]/2))^2)) / (2 * (nrow(binary_matrix)*0.2)^2))) *
        #    (1 / (1 + total_prey_vec[pred] * (1 - tcp)))
      #  ) / (times_preyed_on_vec[prey] + 1)
        
        # Formula with tp/2 and sd
       # rPP_matrix_raw[pred, prey] <- (
       #   exp(-(((abs(tp[pred] - tp[prey]) - (tp[pred]/2))^2)) / (2 * (tp_sd^2))) *
       #     (1 / (1 + total_prey_vec[pred] * (1 - tcp)))
       # ) / (times_preyed_on_vec[prey] + 1)
        
        rPP_matrix_raw[pred, prey] <- (
          (exp(-(((abs(tp[pred] - tp[prey]) - (median_tp[pred]))^2)) / (2 * (nrow(binary_matrix)*0.2)^2))) *
            (1 / (1 + total_prey_vec[pred] * (1 - tcp)))
        ) / (times_preyed_on_vec[prey] + 1)
        
      }
    }
  }
  
  # Normalize rPP per prey (column)
  rPP_matrix_norm <- rPP_matrix_raw
  for (prey in colnames(rPP_matrix_raw)) {
    row_sum <- sum(rPP_matrix_raw[,prey ], na.rm = TRUE)
    if (row_sum > 0) {
      rPP_matrix_norm[,prey ] <- rPP_matrix_raw[,prey ] / row_sum
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
          MedianPreyTP = median_tp[pred],
          MedianPredTP = median_ptp[prey],
          AbsDev = abs(tp[prey] - median_tp[pred]),
          rPP_raw = rPP_matrix_raw[pred, prey],
          rPP = rPP_matrix_norm[pred, prey]
        ))
      }
    }
  }
  
  return(list(
    rPP_matrix = rPP_matrix_raw, #
    TP = tp_df,
    InteractionTable = interaction_table
  ))
}
# Apply the function
rPP_outputs <- lapply(fw_binary_list, compute_rpp_and_tp)

matrix_pairs_list_df <- list()  # NEW list to store original & rPP matrices
matrix_pairs_list <- list()  # NEW list to store original & rPP matrices


# Initialize global summary
correlation_results <- data.frame(
  Network = character(),
  Mean_Spearman = numeric(),
  Median_Spearman = numeric(),
  Links = integer(),
  stringsAsFactors = FALSE
)

# Initialize detailed table
prey_spearman_table <- data.frame(
  Network = character(),
  Prey = character(),
  Spearman = numeric(),
  Num_Predators = integer(),
  stringsAsFactors = FALSE
)

# Loop through all food web matrices
for (name in names(fw_data_list)) {
  original_mat <- fw_data_list[[name]]
  rpp_mat <- rPP_outputs[[name]]$rPP_matrix
  
  # Store as data.frames in the list
  matrix_pairs_list[[name]] <- list(
    Original_df = as.data.frame(original_mat),
    rPP_df = as.data.frame(rpp_mat)
  )
  prey_species <- colnames(original_mat)
  predator_species <- rownames(original_mat)
  
  spearman_values <- c()
  
  for (prey in prey_species) {
    original_col <- original_mat[, prey]
    rpp_col <- rpp_mat[, prey]
    
    # Only compare where there's interaction
    valid_mask <- which(original_col > 0 & !is.na(rpp_col))
    
    if (length(valid_mask) > 1) {
      orig_vals <- original_col[valid_mask]
      rpp_vals <- rpp_col[valid_mask]
      
      spearman_corr <- suppressWarnings(cor(orig_vals, rpp_vals, method = "spearman"))
      spearman_values <- c(spearman_values, spearman_corr)
      
      # Store per-prey result with number of predators
      prey_spearman_table <- rbind(
        prey_spearman_table,
        data.frame(
          Network = name,
          Prey = prey,
          Spearman = spearman_corr,
          Num_Predators = length(valid_mask)
        )
      )
    }
  }
  
  # Global summary per network
  correlation_results <- rbind(
    correlation_results,
    data.frame(
      Network = name,
      Mean_Spearman = mean(spearman_values, na.rm = TRUE),
      Median_Spearman = median(spearman_values, na.rm = TRUE),
      Links = nrow(original_mat)
    )
  )
  corr_sum<-sum(correlation_results$Median_Spearman)
}
print(corr_sum)
###############################
# Function to check case by case
Check_comparisons <- function(list_name, prey_species) {
  original_mat <- fw_data_list[[list_name]]
  rpp_mat <- rPP_outputs[[list_name]]$rPP_matrix
  
  # Check that prey_species exists in the matrix
  if (!(prey_species %in% colnames(original_mat))) {
    stop(paste("Prey species", prey_species, "not found in", list_name))
  }
  
  # Extract the prey column
  original_col <- original_mat[, prey_species]
  rpp_col <- rpp_mat[, prey_species]
  
  # Only keep values where interaction is present
  valid_mask <- which(original_col > 0 & !is.na(rpp_col))
  
  if (length(valid_mask) > 1) {
    orig_vals <- original_col[valid_mask]
    rpp_vals <- rpp_col[valid_mask]
    predators <- rownames(original_mat)[valid_mask]
    
    result_df <- data.frame(
      Predator = predators,
      Original = orig_vals,
      rPP = rpp_vals
    )
    
    view(result_df)
    return(result_df)
  } else {
    cat("No valid predator interactions for", prey_species, "in", list_name, "\n")
    return(NULL)
  }
}

Check_comparisons("FW_017_02", "Cephalopoda")


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
