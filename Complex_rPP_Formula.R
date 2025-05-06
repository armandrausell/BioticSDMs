

compute_rpp_and_tp <- function(binary_matrix, simplicity = FALSE) {
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
  
  # Additional step: μ and σ per predator if needed
  if (!simplicity) {
    mu_tp_diff <- sapply(predators, function(pred) {
      prey_names <- preys[binary_matrix[pred, ] == 1]
      if (length(prey_names) == 0) return(0)
      mean(tp[pred] - tp[prey_names])
    })
    
    sigma_tp_diff <- sapply(predators, function(pred) {
      prey_names <- preys[binary_matrix[pred, ] == 1]
      if (length(prey_names) == 0) return(1) # prevent division by zero
      sd(tp[pred] - tp[prey_names])
    })
  }
  
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
        
        # rPP computation
        if (simplicity) {
          rPP_raw <- (log2(abs(tp[pred] - tp[prey]) + 1) * (1 / (1 + total_prey_vec[pred] * (1 - tcp)))) /
            (times_preyed_on_vec[prey] + 1)
        } else {
          gauss_weight <- exp(-((tp[pred] - tp[prey] - mu_tp_diff[pred])^2) / (2 * sigma_tp_diff[pred]^2))
          rPP_raw <- (gauss_weight * (1 / (1 + total_prey_vec[pred] * (1 - tcp)))) /
            (times_preyed_on_vec[prey] + 1)
        }
        
        rPP_matrix_raw[pred, prey] <- rPP_raw
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
  
  # Step 7: Create interaction table
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
          MedianTP_Pred = median_tp[pred],
          AbsDev = abs(tp[prey] - median_tp[pred]),
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
