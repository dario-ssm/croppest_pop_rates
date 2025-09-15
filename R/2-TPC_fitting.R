library(tidyverse)
library(here)
library(readxl)
library(MuMIn)
source(here::here("R/S2-TPC_equations.R"))
source(here::here("R/S1-Functions.R"))

# 1. Load -----------------------------------------------------------------
int_rate_data <- read_excel(here("data/data_source/int_rate_dataset_new.xlsx")) |> 
  mutate(reference = as_factor(reference),
         id_pop = as_factor(id_pop)) |> 
  separate(order, into = c("order", "family"), sep = ": ", remove = TRUE)

# 2. First TPC fitting -----------------------------------------------------------------
tpcs_AICs_params <- tibble()

pb <- progress::progress_bar$new(
  format = "Fitting  TPCs [:bar] :percent",
  total = 312,
  clear = F)
for(i in unique(int_rate_data$id_pop)[103:312]) {
  int_rate_i <-int_rate_data |> 
    filter(id_pop == i) |> 
    filter(!is.na(int_rate))
  species_i <- int_rate_i |> 
    slice(1) |> 
    pull(species)
  reference_i <- int_rate_i |> 
    slice(1) |> 
    pull(reference)
  fitted_tpc_i <- suppressMessages(fit_tpcs(temp = int_rate_i$temperature,
                                            int_rate = int_rate_i$int_rate,
                                            model_name = "all"))
  aics_params_i <- fitted_tpc_i |> 
    group_by(model_name) |> 
    mutate(n_params = n_distinct(param_name)) |> 
    slice(1) |> 
    select(model_name, model_AIC, n_params) |> 
    mutate(id = i)
  tpcs_AICs_params <- bind_rows(tpcs_AICs_params, aics_params_i)
  pb$tick()
  }
tpcs_AICs_params
writexl::write_xlsx(tpcs_AICs_params, here("data/data_sink/tpcs_aics_params.xlsx"))



plot_tpcs(temp = gao2013_example$temperature,
          int_rate = gao2013_example$int_rate,
          fitted_parameters = fit_example,
          species = "Acyrthosiphon gossypii")

# 3. Simulate TPCs -----------------------------------------------------------------
simulated_tpcs <- tibble()

for(i in 142) {
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
  write_rds(sim_tpcs_i, file = paste0(here("data/data_sink/boots_tpc_"), i, ".rds"))
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

# 4. Fit TPCs for AIC filtering -----------------------------------------------------------------
filtered_tpcs <- readxl::read_excel(here("data/data_sink/filter_for_aics.xlsx")) |> 
  filter(tpc_model != "simp_briere1") |> 
  rename(id_pop = id) |> 
  mutate(id_pop = as_factor(id_pop))
joined_int_rate_selected <- left_join(filtered_tpcs, int_rate_data)

tpcs_AICs_params <- tibble()

pb <- progress::progress_bar$new(
  format = "Fitting  TPCs [:bar] :percent",
  total = 312,
  clear = F)

for(i in unique(filtered_tpcs$id_pop)[271:282]) {
  int_rate_i <-int_rate_data |> 
    filter(id_pop == i) |> 
    filter(!is.na(int_rate))
  species_i <- int_rate_i |> 
    slice(1) |> 
    pull(species)
  reference_i <- int_rate_i |> 
    slice(1) |> 
    pull(reference)
  models_i <- filtered_tpcs |> 
    filter(id_pop == i) |> 
    pull(tpc_model)
  
  fitted_tpc_i <- suppressMessages(fit_tpcs(temp = int_rate_i$temperature,
                                            int_rate = int_rate_i$int_rate,
                                            model_name = models_i))
  aics_params_i <- fitted_tpc_i |> 
    group_by(model_name) |> 
    mutate(n_params = n_distinct(param_name)) |> 
    slice(1) |> 
    select(model_name, model_AIC, n_params) |> 
    mutate(id = i)
  tpcs_AICs_params <- bind_rows(tpcs_AICs_params, aics_params_i)
  pb$tick()
}
writexl::write_xlsx(tpcs_AICs_params, here("data/data_sink/tpcs_aics_params.xlsx"))


# 4. Select deltas AICs -----------------------------------------------------------------
int_rate_data_chrs <- int_rate_data |> 
  select(id_pop, reference, species) |> 
  group_by(id_pop) |> 
  slice(1)
tpc_modelselection_aics <- readxl::read_excel(here("data/data_sink/aics_for_deltas.xlsx")) |>
  rename(id_pop = id) |> 
  mutate(id_pop = as_factor(id_pop)) |> 
  select(-species, -reference) |> 
  left_join(int_rate_data_chrs, by = "id_pop") |> 
  relocate(c(1, 9, 10, 2:8)) |> 
  mutate(AIC = as.numeric(AIC)) |> 
  group_by(id_pop) |> 
  mutate(
    min_AIC = min(AIC, na.rm = TRUE),
    delta_AIC = AIC - min_AIC,
    within_delta2 = delta_AIC < 2
  ) |> 
  ungroup()  
  
writexl::write_xlsx(tpc_modelselection_aics,
                    here("data/data_sink/tpcs_selection_filters.xlsx"))


