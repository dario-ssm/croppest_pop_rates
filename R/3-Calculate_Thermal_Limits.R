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


# 2. Thermal limits for selected TPCs -----------------------------------------------------------------

## a) import selected tpcs ----
tpcs_selected <- readxl::read_excel(here("data/data_sink/tpcs_selection_filters_completed.xlsx")) |> 
  filter(step3_uncertainties == "y")

## b) join bootstraps for each selected tpc ----
thermal_limits_tpcs <- tibble()

pb <- progress::progress_bar$new(
  format = "Calculating Thermal Limits [:bar] :percent",
  total = 312,
  clear = F)

for(population_id in tpcs_selected$id_pop) {
  tpc_selected_i <- tpcs_selected |> 
    filter(id_pop == population_id) |> 
    pull(tpc_model)
  int_rate_i <- int_rate_data |> 
    filter(!is.na(int_rate)) |> 
    filter(id_pop == population_id)
  boots_tpc_id_i <- read_rds(paste0(here("data/data_sink/boots_tpcs/boots_tpc_"), population_id, ".rds")) |>
    filter(model_name_iter == tpc_selected_i)
  thermal_limits_pop_i <- tibble()
  nona_iters <- which(!is.na(unique(boots_tpc_id_i$boots_iter)))
  for(iter_i in nona_iters) {
    cat(paste0(iter_i, "/100\n"))
    boots_tpc_id_iter_i <- boots_tpc_id_i |> 
      filter(boots_iter == iter_i)
    thermal_limits_iter_i <- get_therm_lims(boots_tpc = boots_tpc_id_iter_i,
                                            temp = int_rate_i$temperature,
                                            int_rate = int_rate_i$int_rate,
                                            tpc_model = tpc_selected_i,
                                            epsilon = 1e-04)
    therm_limits_tbl_iter_i <- tibble(tmin = thermal_limits_iter_i[1],
                                      tmax = thermal_limits_iter_i[2],
                                      pop_id = population_id)
    thermal_limits_pop_i <- bind_rows(thermal_limits_pop_i, therm_limits_tbl_iter_i) 
  }
  thermal_limits_tpcs <- bind_rows(thermal_limits_tpcs, thermal_limits_pop_i)
  pb$tick()
}
thermal_limits_tpcs

tmin_ests <- thermal_limits_tpcs |>
  select(tmin, pop_id) |> 
  drop_na() |> 
  group_by(pop_id) |>
  mutate(n_tmin_est = n()) |> 
  reframe(tmin_est = mean(tmin, na.rm = TRUE),
          tmin_se =  sd(tmin, na.rm = TRUE) / sqrt(n_tmin_est),
          n_tmin_est) |> 
  group_by(pop_id) |> 
  slice(1) 


tmax_ests <- thermal_limits_tpcs |>
  select(tmax, pop_id) |> 
  drop_na() |> 
  group_by(pop_id) |>
  mutate(n_tmax_est = n()) |> 
  reframe(tmax_est = mean(tmax, na.rm = TRUE),
          tmax_se =  sd(tmax, na.rm = TRUE) / sqrt(n_tmax_est),
          n_tmax_est) |> 
  group_by(pop_id) |> 
  slice(1)

therm_lims_est <- inner_join(tmin_ests, tmax_ests)
view(therm_lims_est)
write_rds(therm_lims_est, here("data/data_sink/therm_lims_est_42.rds"))
