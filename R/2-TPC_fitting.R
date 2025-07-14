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
simulated_tpcs <- tibble()

for(i in 1:67) {
  cat(paste0("Beginning population ", i, " of 313\n"))
  int_rate_i <-int_rate_data |> 
    filter(id_pop == i) |> 
    filter(!is.na(int_rate))
  species_i <- int_rate_i |> 
    slice(1) |> 
    pull(species)
  reference_i <- int_rate_i |> 
    slice(1) |> 
    pull(reference) |> 
    as.character()
  
  fitted_tpc_i <- fit_tpcs(temp = int_rate_i$temperature,
                           int_rate = int_rate_i$int_rate,
                           model_name = "all")
  plot_tpcs(temp = int_rate_i$temperature,
            int_rate = int_rate_i$int_rate,
            fitted_parameters = fitted_tpc_i,
            species = species_i,
            reference = paste0(reference_i,"; ID: ", i))
  ggsave(paste0(here("data/data_sink/figures/supplementary_figs/fit_tpcs/id"), i, ".png"),
         width = 2100,
         height = 2100,
         units = "px")
  
  fitted_tpc_equations_i <- unique(fitted_tpc_i$model_name)
  
  possible_error <- tryCatch(expr = {
    sim_tpcs_i <- predict_curves(temp = int_rate_i$temperature,
                               int_rate = int_rate_i$int_rate,
                               fitted_parameters = fitted_tpc_i,
                               model_name_2boot = fitted_tpc_equations_i,n_boots_samples = 100)
  write_rds(sim_tpcs_i, file = paste0(here("data/data_sink/boots_tpcs/boots_tpc_"), i, ".rds"))
  plot_uncertainties(bootstrap_uncertainties_tpcs = sim_tpcs_i,
                     temp = int_rate_i$temperature,
                     int_rate = int_rate_i$int_rate,
                     species = species_i,
                     reference = reference_i,
                     pop_id = i)
  ggsave(paste0(here("data/data_sink/figures/supplementary_figs/sim_tpcs/id"), i, ".png"),
         width = 2100,
         height = 2100,
         units = "px")
  }, # <- inside tryCatch
  error = function(e) e)
if (inherits(possible_error, "error")) {
  fit_nls <- NULL
  warning(paste0("Reference ID", i, "(",reference_i, ") could not fit bootstrapped TPCs"))
}
if (is.null(fit_nls)) {
  simulated_tpcs <- simulated_tpcs
} else {
  simulated_tpcs <- bind_rows(simulated_tpcs, sim_tpcs_i) |> 
  drop_na()
  }
  cat(paste0("Ending population ", i, " of 313\n"))
 }

TPCfit_example <- MuMIn::AICc()(fitted_tpc_i$model_fit[[1]])
summary(TPCfit_example) 



