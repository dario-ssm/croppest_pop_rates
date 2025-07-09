library(tidyverse)
library(here)
library(readxl)
source(here::here("R/S2-TPC_equations.R"))
source(here::here("R/S1-Functions.R"))

# 1. Load -----------------------------------------------------------------
int_rate_data <- read_excel(here("data/data_source/int_rate_dataset_new.xlsx")) |> 
  mutate(reference = as_factor(reference),
         id_pop = as_factor(id_pop)) |> 
  separate(order, into = c("order", "family"), sep = ": ", remove = TRUE)

# 2. First TPC fitting -----------------------------------------------------------------
for(i in unique(int_rate_data$id_pop)) {
  int_rate_i <-int_rate_data |> 
    filter(id_pop == i) |> 
    filter(!is.na(int_rate))
  species_i <- int_rate_i |> 
    slice(1) |> 
    pull(species)
  reference_i <- int_rate_i |> 
    slice(1) |> 
    pull(reference)
  
  fitted_tpc_i <- fit_tpcs(temp = int_rate_i$temperature,
                          int_rate = int_rate_i$int_rate,
                          model_name = "all")
  plot_tpcs(temp = int_rate_i$temperature,
            int_rate = int_rate_i$int_rate,
            fitted_parameters = fitted_tpc_i,
            species = species_i,
            life_stage = paste0(reference_i,"; ID: ", i))
  ggsave(paste0(here("data/data_sink/figures/supplementary_figs/fit_tpcs/id"), i, ".png"),
         width = 2100,
         height = 2100,
         units = "px")
  }

plot_tpcs(temp = gao2013_example$temperature,
          int_rate = gao2013_example$int_rate,
          fitted_parameters = fit_example,
          species = "Acyrthosiphon gossypii")

# 3. Simulate TPCs -----------------------------------------------------------------
example <- int_rate_data |> 
  filter(id_pop == 24)
fitted_ex <- fit_tpcs(temp = example$temperature,
                      int_rate = example$int_rate,
                      model_name = "all")
plot_tpcs(temp = example$temperature,
          int_rate = example$int_rate,
          fitted_parameters = fitted_ex,
          species = "Aphis gossypii",
          reference = "xia1999")

sim_tpcs_example <- predict_curves(temp = example$temperature,
                                   int_rate = example$int_rate,
                                   fitted_parameters = fitted_ex,
                                   model_name_2boot = c("analytis_kontodimas",
                                                        "atkin",
                                                        "beta", "briere1",
                                                        "johnk",
                                                        "simp_beta",
                                                        "simp_briere1",
                                                        "taylor_sexton",
                                                        "thomas",
                                                        "weibull"))

plot_uncertainties(bootstrap_uncertainties_tpcs = sim_tpcs_example,
                   temp = example$temperature,
                   int_rate = example$int_rate,
                   species = "Aphis gossypii",
                   reference = "xia1999",
                   pop_id = 24)
