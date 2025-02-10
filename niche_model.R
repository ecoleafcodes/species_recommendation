library(dplyr)
library(tibble)


compute_niches <- function(data, file_path, method = "count", abundance_var = NULL, 
                           env_vars_parc, env_vars_plot, species_col, parc_col, plot_col) {

  # Read and clean the data
  data <- data %>%
    filter(if_all(everything(), ~ !is.nan(.)))

  # Ensure abundance_var is provided if method is not "count"
  if (method != "count" && is.null(abundance_var)) {
    stop("Error: 'abundance_var' must be provided when method is not 'count'.")
  }

  # Function to compute relative values (count-based or variable-based)
  compute_relative_values <- function(data, group_var, method, abundance_var) {
    data <- data %>%
      group_by(across(all_of(c(group_var, species_col)))) %>%
      summarise(
        value_sum = if (method == "count") n() else sum(.data[[abundance_var]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      group_by(across(all_of(group_var))) %>%
      mutate(relative_value = value_sum / sum(value_sum, na.rm = TRUE)) %>%
      ungroup()
    
    return(data)
  }

  # ---- Compute relative values for parc ----
  relative_parc <- compute_relative_values(data, parc_col, method, abundance_var)

  # Extract environmental variables for parc
  env_vars_parc_df <- data %>%
    select(all_of(c(parc_col, env_vars_parc))) %>%
    distinct()

  # Merge relative values with environmental variables for parc
  final_parc_df <- relative_parc %>%
    left_join(env_vars_parc_df, by = parc_col)

  # ---- Compute relative values for plot ----
  relative_plot <- compute_relative_values(data, plot_col, method, abundance_var)

  # Compute average environmental variables for each plot
  plot_vars_avg <- data %>%
    group_by(across(all_of(plot_col))) %>%
    summarise(across(all_of(env_vars_plot), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

  # Merge relative values with averaged environmental variables for plots
  final_plot_df <- relative_plot %>%
    left_join(plot_vars_avg, by = plot_col)

  # ---- Step 1: Calculate weighted averages for final_parc_df ----
  final_parc_avg <- final_parc_df %>%
    group_by(across(all_of(species_col))) %>%
    summarise(across(all_of(env_vars_parc), ~ weighted.mean(.x, relative_value, na.rm = TRUE), .names = "{.col}"), .groups = "drop")

  # ---- Step 2: Calculate weighted averages for final_plot_df ----
  final_plot_avg <- final_plot_df %>%
    group_by(across(all_of(species_col))) %>%
    summarise(across(all_of(env_vars_plot), ~ weighted.mean(.x, relative_value, na.rm = TRUE), .names = "{.col}"), .groups = "drop")

  # ---- Step 3: Combine both results by species ----
  niche_data <- final_parc_avg %>%
    inner_join(final_plot_avg, by = species_col)

  return(niche_data)
}


predict_species <- function(niche_data, 
                            env_vars_parc, 
                            env_vars_plot, 
                            input_vector, 
                            N = 10) {
  
  # Combine env_vars_parc and env_vars_plot
  env_vars <- c(env_vars_parc, env_vars_plot)
  
  # Check if all specified environmental variables exist in niche_data
  missing_vars <- setdiff(env_vars, colnames(niche_data))
  if (length(missing_vars) > 0) {
    stop(paste("The following variables are missing from niche_data:", paste(missing_vars, collapse = ", ")))
  }
  
  # Standardize the data (Z-score standardization)
  niche_data_standardized <- niche_data %>%
    select(species, all_of(env_vars)) %>%
    mutate(across(all_of(env_vars), ~ ( . - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE)))

  # Standardize the input vector
  input_vector_standardized <- input_vector
  for (var in env_vars) {
    input_vector_standardized[var] <- (input_vector[var] - mean(niche_data[[var]], na.rm = TRUE)) / sd(niche_data[[var]], na.rm = TRUE)
  }
  
  # Extract species vectors from the standardized niche_data
  species_vectors <- niche_data_standardized %>%
    select(species, all_of(env_vars)) %>%
    column_to_rownames(var = "species")  # Make species as row names for easier reference

  # Compute the cosine similarity between the standardized input vector and each species vector
  cosine_similarity_results <- apply(species_vectors, 1, function(species_vector) {
    # Ensure both species_vector and input_vector are in the same order for calculation
    common_vars <- names(species_vector)
    
    dot_product <- sum(species_vector * input_vector_standardized[common_vars], na.rm = TRUE)
    magnitude_species <- sqrt(sum(species_vector^2, na.rm = TRUE))
    magnitude_input <- sqrt(sum(input_vector_standardized[common_vars]^2, na.rm = TRUE))

    # Avoid division by zero
    if (magnitude_species == 0 || magnitude_input == 0) {
      return(0)
    }

    # Cosine similarity calculation
    dot_product / (magnitude_species * magnitude_input)
  })

  # Convert the cosine similarity results into a data frame
  cosine_similarity_df <- data.frame(
    species = names(cosine_similarity_results),
    cosine_similarity = cosine_similarity_results,
    row.names = NULL
  )

  # Select the top N species with the highest cosine similarity
  top_species <- cosine_similarity_df %>%
    arrange(desc(cosine_similarity)) %>%
    slice_head(n = N)

  # Return the top N species
  return(top_species)
}
