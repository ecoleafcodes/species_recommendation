library(dplyr)
library(ggplot2)


# Function to compute relative abundance for a specific plot
compute_relative_abundance <- function(data, plot_col, species_col, abundance_var, plot_id, method) {
  plot_data <- data %>%
    filter(.data[[plot_col]] == plot_id) %>%
    group_by(across(all_of(c(plot_col, species_col)))) %>%
    summarise(
      value_sum = if (method == "count") n() else sum(.data[[abundance_var]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(across(all_of(plot_col))) %>%
    mutate(relative_value = value_sum / sum(value_sum, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(desc(relative_value)) %>%
    mutate(rank = row_number())

  return(plot_data)
}

# Function to get average environmental variables for a specific plot
get_avg_env_vars <- function(data, plot_col, env_vars_plot, plot_id) {
  data %>%
    filter(.data[[plot_col]] == plot_id) %>%
    summarise(across(all_of(env_vars_plot), \(x) mean(x, na.rm = TRUE)))
}

# Function to compute summed relative abundance of overlapping species
compute_overlap_abundance <- function(observed_data, predicted_species) {
  intersection_species <- intersect(observed_data$species, predicted_species$species)
  summed_abundance <- observed_data %>%
    filter(species %in% intersection_species) %>%
    summarise(total_relative_abundance = sum(relative_value, na.rm = TRUE)) %>%
    pull(total_relative_abundance)

  return(summed_abundance)
}


source("niche_model.R")

# Main analysis
file_path <- "DADOS_MARIANA.csv"
data <- read.csv(file_path, sep = ";")

plots <- unique(data$Plot)

abundance_method <- "count"

abundance_var <- "AGBtree"

# Compute niche_data once
niche_data <- compute_niches(
  data = data,
  method = abundance_method, # or count
  abundance_var = abundance_var,
  env_vars_parc = c(),
  env_vars_plot = c("pH", "K", "P", "Ca", "Mg", "Al", "H_Al", "SB", "t", "T", "V", "m", "MO", "Argila", "Silte", "Areia", "AWD1", "AWD2", "AWD3", "AWD4"),
  species_col = "species",
  parc_col = "chave_plot_parc",
  plot_col = "Plot"
)


# Store summed relative abundances
summed_abundances <- sapply(plots, function(plot_id) {
  relative_abundance_plot <- compute_relative_abundance(data, "Plot", "species", abundance_var, plot_id, abundance_method)
  avg_env_vars <- get_avg_env_vars(data, "Plot", c("pH", "K", "P", "Ca", "Mg", "Al", "H_Al", "SB", "t", "T", "V", "m", "MO", "Argila", "Silte", "Areia", "AWD1", "AWD2", "AWD3", "AWD4"), plot_id)

  predicted_species <- predict_species(
    niche_data = niche_data,  # Use precomputed niche_data
    env_vars_parc = c(),
    env_vars_plot = c("pH", "K", "P", "Ca", "Mg", "Al", "H_Al", "SB", "t", "T", "V", "m", "MO", "Argila", "Silte", "Areia", "AWD1", "AWD2", "AWD3", "AWD4"),
    input_vector = avg_env_vars,
    N = nrow(relative_abundance_plot)
  )

  compute_overlap_abundance(relative_abundance_plot, predicted_species)
})

summed_abundances

# Calculate mean summed relative abundance
mean_abundance <- mean(summed_abundances, na.rm = TRUE)

# Plot histogram
histogram_plot <- ggplot(data.frame(summed_abundances), aes(x = summed_abundances)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  geom_vline(aes(xintercept = mean_abundance), color = "red", linetype = "dashed", size = 1) +
  labs(title = "Histogram of Summed Relative Abundances", x = "Summed Relative Abundance", y = "Frequency") +
  theme_minimal()

# Display mean and histogram
print(mean_abundance)
print(histogram_plot)


# Store summed relative abundances and plot IDs
plot_proportions <- data.frame(
  Plot = plots,
  SummedRelativeAbundance = summed_abundances
)

plot_proportions

# ggsave("summed_relative_abundance_histogram.png", plot = histogram_plot, 
#    width = 10,
#    height = 8,
#    units = "cm")

# save.image("analysis_workspace.RData")
