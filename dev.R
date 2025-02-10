source("niche_model.R")


file_path <- "DADOS_MARIANA.csv"

data <- read.csv(file_path, sep = ";")


niche_data <- compute_niches(
  data = data,  
  abundance_var = "AGBtree",
  env_vars_parc = c("pH", "K", "P", "Ca", "Mg", "Al", "H_Al", "SB", "t", "T", "V", "m", "MO", "Argila", "Silte", "Areia", "AWD1", "AWD2", "AWD3", "AWD4"),
  env_vars_plot = c("Lat", "Long"),
  species_col = "species",
  parc_col = "chave_plot_parc",
  plot_col = "Plot"
)


input_vector <- c(
  pH = 6.5,
  K = 0.2,
  P = 0.05,
  Ca = 2.0,
  Mg = 1.2,
  Al = 0.1,
  H_Al = 0.3,
  SB = 15,
  t = 22.5,
  T = 18.3,
  V = 3.2,
  m = 5.4,
  MO = 6.7,
  Argila = 35,
  Silte = 40,
  Areia = 25,
  AWD1 = 0.3,
  AWD2 = 0.2,
  AWD3 = 0.1,
  AWD4 = 0.4,
  Lat = -23.5,    
  Long = -45.0   
)

# Call the function with example parameters
top_species <- predict_species(
  niche_data = niche_data,
  env_vars_parc = c("pH", "K", "P", "Ca", "Mg", "Al", "H_Al", "SB", "t", "T", "V", "m", "MO", "Argila", "Silte", "Areia", "AWD1", "AWD2", "AWD3", "AWD4"),
  env_vars_plot = c("Lat", "Long"),
  input_vector = input_vector,
  N = 10
)

# Display the top species
print(top_species)



