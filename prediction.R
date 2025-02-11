source("niche_model.R")


file_path <- "DADOS_MARIANA.csv"
data <- read.csv(file_path, sep = ";")

plots <- unique(data$Plot)

abundance_method <- "biomass"

abundance_var <- "AGBtree"

env_vars_parc = c()
env_vars_plot = c("pH", "K", "P", "Ca", "Mg", "Al", "H_Al", "SB", "t", "T", "V", "m", "MO", "Argila", "Silte", "Areia", "AWD1", "AWD2", "AWD3", "AWD4", "AWD5", "AWD6")

env_vars <- c(4, 46.91, 2.69, 0.43, 0.31, 1.5, 9.9, 0.86, 2.36, 10.76, 8, 63.56, 4.88, 47, 44, 9, 1.75, 1.85, 1.66, 1.86, 1.7, 1.52)

n_species <- 100

# Compute niche_data once
niche_data <- compute_niches(
  data = data,
  method = abundance_method, # or count
  abundance_var = abundance_var,
  env_vars_parc = env_vars_parc,
  env_vars_plot = env_vars_plot,
  species_col = "species",
  parc_col = "chave_plot_parc",
  plot_col = "Plot"
)


predicted_species <- predict_species(
niche_data = niche_data,  # Use precomputed niche_data
env_vars_parc = env_vars_parc,
env_vars_plot = env_vars_plot,
input_vector = env_vars,
N = n_species
)


print(head(predicted_species))
