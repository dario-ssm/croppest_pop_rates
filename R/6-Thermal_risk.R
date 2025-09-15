library(tidyverse)
library(here)
library(readxl)
library(sf)
library(terra)
library(leaflet)
source(here::here("R/S2-TPC_equations.R"))
source(here::here("R/S1-Functions.R"))

# 1. Load -----------------------------------------------------------------
int_rate_data <- read_excel(here("data/data_source/int_rate_dataset_new.xlsx")) |> 
  mutate(reference = as_factor(reference),
         id_pop = as_factor(id_pop)) |> 
  separate(order, into = c("order", "family"), sep = ": ", remove = TRUE)

int_rate_data_info <- int_rate_data |> 
  group_by(id_pop) |> 
  slice(1) |> 
  rename(pop_id = id_pop) |> 
  select(reference, species, order, family, lat, lon, pop_id)

thermal_limits <- read_rds(here("data/data_sink/therm_lims_est_42.rds")) |> 
  mutate(pop_id = as_factor(pop_id))

thermal_limits_intrate <- inner_join(int_rate_data_info, thermal_limits)

# 2. Climate data -----------------------------------------------------------------

source_dataset_points_ids <- int_rate_data |> 
  select(lat, lon) |> 
  unique() |> 
  mutate(id_location = 1:n())

list_locations <- int_rate_data |> 
  inner_join(source_dataset_points_ids) |>
  group_by(id_location) |> 
  slice(1) |> 
  filter(!is.na(lon)) 
  

points_intrates <-st_as_sf(x = list_locations,
                           coords = c("lon", "lat")) |> 
  st_set_crs(4326)

#points_vect_intrates <- terra::vect(points_intrates)
present_tmin_climate_wc <-  geodata::worldclim_global(var = "tmin",
                                                      res = 2.5,
                                                      path = tempdir())

points_tmin_present <- terra::extract(present_tmin_climate_wc,
                                      points_intrates, 
                                      id = points_intrates$id_location) |> 
  as_tibble() |> 
  pivot_longer(cols = -ID,
               names_to = "month",
               values_to = "tmin_avg") |> 
  mutate(month = str_sub(month, -2),
         month = as_factor(month))

#### a) Historical ----
present_tmax_climate_wc <-  geodata::worldclim_global(var = "tmax",
                                                      res = 2.5,
                                                      path = tempdir())

points_tmax_present <- terra::extract(present_tmax_climate_wc,
                                      points_intrates, 
                                      id = points_intrates$id_location) |> 
  as_tibble() |> 
  pivot_longer(cols = -ID,
               names_to = "month",
               values_to = "tmax_avg") |> 
  mutate(month = str_sub(month, -2),
         month = as_factor(month))

points_tavg_present <- inner_join(points_tmin_present,
                                  points_tmax_present) |> 
  mutate(tavg = map2_dbl(.x = tmin_avg,
                         .y = tmax_avg,
                         .f = ~mean(c(.x, .y))),
         model = as_factor("WorldClim_1970-2000_v21")) |> 
  rename(tmin = tmin_avg,
         tmax = tmax_avg) |> 
  pivot_longer(cols = 3:5, 
               names_to = "var",
               values_to = "temp_value") |> 
  mutate(time_scenario = as_factor("Present")) |>
  filter(var == "tavg")

### b) Future CMIP6 RCP 4.5 ----

cmip6_avg_tbl <- read_rds(here("data/data_source/cmip6_tavg_2041-2060_ssp245_res25.rds"))

points_tavg_future <- cmip6_avg_tbl |>
  select(ID, month, model, tavg) |> 
  pivot_longer(tavg, names_to = "var", values_to = "temp_value") |> 
  mutate(time_scenario = "Future")
  
### c) Joined (Pres-Future) ----
points_tavg_worldclim <- points_tavg_present |> 
  bind_rows(points_tavg_future)

points_tavg_synth <- points_tavg_worldclim |> 
  group_by(ID, month, time_scenario) |> 
  summarise(temp_value = mean(temp_value)) |> 
  rename(id_location = ID)

int_rate_locations <- int_rate_data |> 
  select(lat, lon) |> 
  unique() |> 
  mutate(id_location = 1:n())

int_rate_data_locations <- int_rate_data |> 
  group_by(lat, lon) |> 
  inner_join(int_rate_locations) |> 
  filter(!is.na(lat)) |> 
  ungroup() |> 
  group_by(id_location) |> 
  mutate(id_location = cur_group_id())

extracted_tavg_intrate <- int_rate_data_locations |> 
  group_by(id_pop, id_location) |> 
  slice(1) |> 
  inner_join(points_tavg_synth) |> 
  select(reference, species, order, family, lat, lon, id_pop, id_location, month, time_scenario, temp_value)

# 3. Project rates -----------------------------------------------------------------
tpcs_selected <- readxl::read_excel(here("data/data_sink/tpcs_selection_filters_completed.xlsx")) |> 
  filter(step3_uncertainties == "y")

tpcs_boots <- tibble()

pb <- progress::progress_bar$new(
  format = "Loading Bootstrapped TPCs :percent",
  total = nrow(tpcs_selected),
  clear = F)

for(population_id in tpcs_selected$id_pop) {
  tpc_selected_i <- tpcs_selected |> 
    filter(id_pop == population_id) |> 
    pull(tpc_model)
  int_rate_i <- int_rate_data |> 
    filter(!is.na(int_rate)) |> 
    filter(id_pop == population_id)
  boots_tpc_id_i <- read_rds(paste0(here("data/data_sink/boots_tpcs/boots_tpc_"), population_id, ".rds")) |>
    filter(model_name_iter == tpc_selected_i) |> 
    mutate(id_pop = population_id)
  calc_parameters_all <- get_therm_lims(boots_tpc = boots_tpc_id_i,
                              temp = int_rate_i$temperature,
                              int_rate = int_rate_i$int_rate,
                              tpc_model = tpc_selected_i,
                              epsilon = 1e-04)
  calc_parameters <- calc_parameters_all |> 
    drop_na() |> 
    summarise(across(where(is.numeric),
                     ~mean(.x[is.finite(.x)], na.rm = TRUE)))
  
  boots_tpc_id_i_params <- boots_tpc_id_i |> 
    mutate(tmin = calc_parameters$tmin,
           tmax = calc_parameters$tmax,
           topt = calc_parameters$topt,
           t_L50 = calc_parameters$t_L50,
           t_R50 = calc_parameters$t_R50,
           rmax = calc_parameters$rmax,
    )
  tpcs_boots <- bind_rows(tpcs_boots, boots_tpc_id_i_params)
  pb$tick()
}
tpcs_boots

#predict performance based on tpcs_boots and extracted monthly tavg

predict_r <- tibble()
pb <- progress::progress_bar$new(
  format = "Predicting population's performance (historical and future) :percent",
  total = 310,
  clear = F)
tpcs_boots_index <- tpcs_boots |> 
  group_by(id_pop) |> 
  slice(1)
iterations_both <- inner_join(tpcs_boots_index, extracted_tavg_intrate)

for(i in unique(iterations_both$id_pop)) {
  tpcs_boots_i <- tpcs_boots |> 
    filter(id_pop == i)  
  tpc_model_i <- tpcs_boots_i$model_name_iter[1]
  extracted_tavg_i <- extracted_tavg_intrate |> 
    filter(id_pop == i)
  reference_i <- as.character(extracted_tavg_i$reference[1])
  species_i <- as.character(extracted_tavg_i$species[1])
  rmax_i <- tpcs_boots_i$rmax[1]
  tmin_i <- tpcs_boots_i$tmin[1]
  t_L50 <- tpcs_boots_i$t_L50[1]
  topt_i <- tpcs_boots_i$topt[1]
  t_R50 <- tpcs_boots_i$t_R50[1]
  tmax_i <- tpcs_boots_i$tmax[1]
  rmax_i <- tpcs_boots_i$rmax[1]
  if(any(is.na(extracted_tavg_i$temp_value))) {
    warning(paste0("id_pop ", i, " had no climate data and was discarded"))
    next
  }
  if(all(is.na(tpcs_boots_i$tmin))) {
    warning(paste0("id_pop ", i, " yield unreliable estimates of Tmin. Please consider selecting a different TPC"))
    
    next}
  
  if(all(is.na(tpcs_boots_i$tmax))) {
    warning(paste0("id_pop ", i, " yield unreliable estimates of Tmax"))
    
    next}
  #loop for each month
  predict_r_pop <- tibble()
  for(month_i in unique(extracted_tavg_intrate$month)) {
    monthly_ext_tavg_i <- extracted_tavg_i |> 
      filter(month == month_i)
    #loop for each time scenario
    predict_r_month <- tibble()
    for(time_i in c("Present", "Future")){
      monthly_ext_tavg_time_i <- monthly_ext_tavg_i |> 
        filter(time_scenario == time_i)
      monthly_time_tvalue <- monthly_ext_tavg_time_i$temp_value
      predicted_r_time_i <- tpcs_boots_i |>
        filter(abs(preds) < 1) |> 
        group_by(id_pop, temp) |> 
        summarise(preds = mean(preds, na.rm = TRUE), .groups = "drop") |> 
        slice_min(temp <= monthly_time_tvalue) |> 
        slice(1) |> 
        select(id_pop, preds)
      if (monthly_time_tvalue < tmin_i || monthly_time_tvalue > tmax_i) {
        predicted_r_time_i$preds <- 0
      }
      predict_r_month_i <- inner_join(monthly_ext_tavg_time_i, predicted_r_time_i, by = "id_pop")
      predict_r_month <- bind_rows(predict_r_month, predict_r_month_i)
    }
    predict_r_month <- predict_r_month |> 
      mutate(rmax = tpcs_boots_i$rmax[1],
             tmin = tpcs_boots_i$tmin[1],
             t_L50 = tpcs_boots_i$t_L50[1],
             topt = tpcs_boots_i$topt[1],
             t_R50 = tpcs_boots_i$t_R50[1],
             tmax = tpcs_boots_i$tmax[1],
             rmax = tpcs_boots_i$rmax[1])
    
    predict_r_pop <- bind_rows(predict_r_pop, predict_r_month)
  }
  predict_r <- bind_rows(predict_r, predict_r_pop)
  pb$tick()
}

predict_r_shift <- predict_r |>
  group_by(reference, species, order, family, lat, lon, id_pop, id_location, month) |> 
  summarise(temp_present = temp_value[time_scenario == "Present"],
            temp_future = temp_value[time_scenario == "Future"],
            preds_present = preds[time_scenario == "Present"],
            preds_future  = preds[time_scenario == "Future"],
            preds_diff = preds_future - preds_present,
            psi = preds_diff/rmax,
            r_max = rmax,
            curve_zone_present = case_when(temp_present < tmin ~ "CE",
                                           temp_present < t_L50 &
                                             temp_present >= tmin ~ "CLP",
                                           temp_present >= t_L50 &
                                             temp_present < topt ~ "OPS",
                                           temp_present >= topt & 
                                             temp_present < t_R50 ~ "OPD",
                                           temp_present >= t_R50 &
                                             temp_present < tmax ~ "HLP",
                                           temp_present > tmax ~ "HE"),
            curve_zone_future = case_when(temp_future < tmin ~ "CE",
                                          temp_future < t_L50 &
                                            temp_future >= tmin ~ "CLP",
                                          temp_future >= t_L50 &
                                            temp_future < topt ~ "OPS",
                                          temp_future >= topt & 
                                            temp_future < t_R50 ~ "OPD",
                                          temp_future >= t_R50 &
                                            temp_future < tmax ~ "HLP",
                                          temp_future >= tmax ~ "HE")
           ) |> 
  slice(1) |> 
  ungroup()

r_max_prop <- predict_r_shift |> 
  group_by(id_pop) |> 
  summarise(preds_present = mean(preds_present, na.rm = TRUE),
            preds_future = mean(preds_future, na.rm = TRUE),
            r_max = mean(r_max, na.rm = TRUE),
            r_hist_perc = (100*preds_present)/r_max,
            r_future_perc = (100*preds_future)/r_max,)
r_max_prop_overall <- r_max_prop |> 
  summarise(r_hist_perc_est = mean(r_hist_perc),
            r_hist_perc_se = sd(r_hist_perc)/sqrt(nrow(r_max_prop)),
            r_future_perc_est = mean(r_future_perc),
            r_future_perc_se = sd(r_future_perc, na.rm = TRUE)/sqrt(nrow(r_max_prop)))

r_count_zones <- predict_r_shift |> 
  select(curve_zone_present, 
         curve_zone_future) |> 
  pivot_longer(cols = 1:2,
               names_to = "time_scenario",
               values_to = "curve_zone") |> 
  mutate(time_scenario = str_sub(time_scenario, 12, -1L)) |> 
  group_by(time_scenario) |> 
  count(curve_zone) |> 
  mutate(n = 100*n/2772) |> 
  mutate(time_scenario = ifelse(time_scenario == "present",
                                "historical",
                                "RCP 4.5 (2041-2060)")) |> 
  mutate(sort = case_when(curve_zone == "CE" ~ 1,
                          curve_zone == "CLP" ~ 2,
                          curve_zone == "OPS" ~ 3,
                          curve_zone == "OPD" ~ 4,
                          curve_zone == "HLP" ~ 5,
                          curve_zone == "HE" ~ 6),
         colorado = case_when(curve_zone == "CE" ~ "#B2F2FD",
                              curve_zone == "CLP" ~ "#61DCA9",
                              curve_zone == "OPS" ~ "#8DB63B",
                              curve_zone == "OPD" ~ "#9B7424",
                              curve_zone == "HLP" ~ "#944046",
                              curve_zone == "HE" ~ "#8C0172"))

ggplot(r_count_zones,
       aes(x = time_scenario, 
           stratum = fct_reorder(curve_zone, sort),
           alluvium = fct_reorder(curve_zone, sort),
           y = n,
           fill = colorado)) +
  ggalluvial::geom_flow(stat = "alluvium", 
                        lode.guidance = "frontback",
                        curve_type = "cubic",
                        alpha = 0.4) +
  ggalluvial::geom_stratum(width = 0.35, color = "gray76") +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_fill_identity()+
  labs(x = "Scenario",
       y = element_blank())

ggsave(here("data/data_sink/figures/supplementary_figs/flowchart_curvezones.svg"),
       width = 2800,
       height = 280,
       units = "px")



write_rds(predict_r_shift, here("data/data_sink/predictions_psi.rds"))


# 4. Analysis risk -----------------------------------------------------------------
predict_r_shift <- read_rds(here("data/data_sink/predictions_psi.rds")) |> 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, na.fail = FALSE) |> 
  st_transform(crs = "+proj=wintri") |> 
  group_by(reference, species, order, family, id_pop, id_location) |> 
  summarise(psi = mean(psi)) |> 
  mutate(win_or_lose = case_when(psi < 0 ~ as_factor("losers"),
                                 psi >= 0 ~ as_factor("winners")))

predict_r_shift
summary(predict_r_shift$psi)
predict_r_shift |> 
  ungroup() |> 
  count(win_or_lose) |> 
  mutate(n = n*100/231)

predict_r_shift |> 
  ungroup() |> 
  group_by(order) |> 
  summarise(psi = mean(psi))


predict_r_shift_preds <- predict_r_shift |> 
  summarise(r_hist = mean(preds_present),
            r_fut = mean(preds_future)) |> 
  mutate(dt_hist = log(2)/r_hist,
         dt_fut = log(2)/r_fut)



join_psi_lat <- read_rds(here("data/data_sink/predictions_psi.rds")) |> 
  group_by(reference, species, order, family, id_pop, id_location) |> 
  summarise(psi = mean(psi),
            lat = mean(lat))  

lat_psi_lm <- lmerTest::lmer(psi~ abs(lat) + (1|reference) + (1|species),
                             data = join_psi_lat)
summary(lat_psi_lm)
# 5. Visualization maps -----------------------------------------------------------------


## a) Static ---------------------------------------------------------------


predict_r_shift <- read_rds(here("data/data_sink/predictions_psi.rds"))


predict_r_shift_sf <- predict_r_shift |> 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, na.fail = FALSE) |> 
  st_transform(crs = "+proj=wintri") |> 
  group_by(reference, species, order, family, id_pop, id_location) |> 
  summarise(psi = mean(psi)) |> 
  mutate(win_or_lose = case_when(psi < 0 ~ as_factor("losers"),
                                 psi >= 0 ~ as_factor("winners")))

worldmap <- rnaturalearth::ne_countries(scale = 50) |> 
  st_transform(crs = "+proj=robin")


ggplot()+
  geom_sf(data = worldmap, fill = "gray45", color = NA)+
  geom_sf(data = predict_r_shift_sf,
          aes(color = psi))+
  theme_bw()+
  theme(legend.position = "bottom")+
  khroma::scale_color_bam(reverse = T, midpoint = 0)+
  ggdark::dark_mode()+
  facet_wrap(~forcats::fct_rev(win_or_lose), ncol = 1)+
  labs(color = expression(italic(psi)[italic(r)]))



ggsave(here("data/data_sink/figures/supplementary_figs/winners_losers.png"),
       width = 2800,
       height = 3600,
       units = "px")

ggplot()+
  geom_sf(data = worldmap, fill = "gray45", color = NA)+
  geom_sf(data = predict_r_shift_sf,
          aes(color = psi))+
  theme_bw()+
  theme(legend.position = "bottom")+
  khroma::scale_color_bam(reverse = T, midpoint = 0)+
  ggdark::dark_mode()+
  labs(color = expression(italic(psi)[italic(r)]))


ggsave(here("data/data_sink/figures/supplementary_figs/psi.png"),
       width = 2800,
       height = 2800,
       units = "px")


#small sample
ggplot()+
  geom_sf(data = worldmap, fill = "gray45", color = NA)+
  geom_sf(data = predict_r_shift_sf,
          aes(color = psi))+
  theme_bw()+
  theme(legend.position = "bottom")+
  khroma::scale_color_bam(reverse = T, midpoint = 0)+
  ggdark::dark_mode()+
  labs(color = expression(italic(psi)[italic(r)]))


ggsave(here("data/data_sink/figures/supplementary_figs/psi_sample.png"),
       width = 21000,
       height = 21000,
       units = "px", dpi = 600)

## b) Interactive ---------------------------------------------------------------

#map it with leaflet

for(i in unique(predict_r_shift$id_pop)){
  predict_r_shift_i <- predict_r_shift  |> 
    filter(id_pop == i) 
  predict_r_shift_curvezone_future <- predict_r_shift_i |> 
    group_by(reference, species, order, family, lat, lon, id_pop, id_location) |> 
    count(curve_zone_future) |> 
    mutate(sort = case_when(curve_zone_future == "CE" ~ 1,
                             curve_zone_future == "CLP" ~ 2,
                             curve_zone_future == "OPS" ~ 3,
                             curve_zone_future == "OPD" ~ 4,
                             curve_zone_future == "HLP" ~ 5,
                             curve_zone_future == "HE" ~ 6),
           colorado = case_when(curve_zone_future == "CE" ~ "#B2F2FD",
                                  curve_zone_future == "CLP" ~ "#61DCA9",
                                  curve_zone_future == "OPS" ~ "#8DB63B",
                                  curve_zone_future == "OPD" ~ "#9B7424",
                                  curve_zone_future == "HLP" ~ "#944046",
                                  curve_zone_future == "HE" ~ "#8C0172"))
  predict_r_shift_curvezone_present <- predict_r_shift_i |> 
    group_by(reference, species, order, family, lat, lon, id_pop, id_location) |> 
    count(curve_zone_present) |> 
    mutate(sort = case_when(curve_zone_present == "CE" ~ 1,
                             curve_zone_present == "CLP" ~ 2,
                             curve_zone_present == "OPS" ~ 3,
                             curve_zone_present == "OPD" ~ 4,
                             curve_zone_present == "HLP" ~ 5,
                             curve_zone_present == "HE" ~ 6),
           colorado = case_when(curve_zone_present == "CE" ~ "#B2F2FD",
                                curve_zone_present == "CLP" ~ "#61DCA9",
                                curve_zone_present == "OPS" ~ "#8DB63B",
                                curve_zone_present == "OPD" ~ "#9B7424",
                                curve_zone_present == "HLP" ~ "#944046",
                                curve_zone_present == "HE" ~ "#8C0172"))
  predict_r_shift_curvezone <- bind_rows(predict_r_shift_curvezone_present,
                                        predict_r_shift_curvezone_future) |> 
    pivot_longer(cols = c(curve_zone_present, curve_zone_future),
                 values_to = "curve_zone",
                 names_to = "time_scenario") |> 
    mutate(time_scenario = str_sub(time_scenario, 12, -1L)) |> 
    mutate(time_scenario = ifelse(time_scenario == "present",
                                  "a) historical",
                                  "b) future")) |> 
    tidyr::drop_na()
  
  species_name <- parse_character(predict_r_shift_curvezone$species)[1]
  reference_name <- as.character(predict_r_shift_curvezone$reference)[1]
  
  barplot_risk_study <- ggplot(data = predict_r_shift_curvezone)+
    geom_bar(aes(x = fct_reorder(as_factor(curve_zone), sort), 
                 y = n, 
                 fill =  as_factor(colorado)),
             stat = "identity"
    )+
    scale_fill_identity()+
    theme_minimal()+
    theme(legend.position = "bottom")+
    labs(title = species_name,
         x = element_blank(),
         y = "#months",
         caption = reference_name)+
    theme(plot.title = element_text(face = "italic"))+
    facet_wrap(~time_scenario)+
    theme(panel.spacing = unit(3, "lines"))  # default is usually 0.5 lines
  
  
  save(barplot_risk_study, 
       file = here(paste0("data/data_sink/figures/riskpoints/barplot_risk_study",i,".RData")))
  ggsave(here(paste0("data/data_sink/figures/riskpoints","/study",i,".jpg")),
         width = 12,
         height = 12,
         units = "cm")
}

list_ggbarplots <- tibble(id = NULL,
                          ggbarplot = NULL)
for(pop_i in unique(predict_r_shift$id_pop)){
  load(file = here(paste0("data/data_sink/figures/riskpoints/barplot_risk_study",pop_i,".RData")))
  ggbarplot_study <- tibble(id_pop = pop_i,
                            ggbarplot = list(barplot_risk_study))
  list_ggbarplots <- bind_rows(list_ggbarplots, ggbarplot_study)
}

leaflet_risk_shift <- predict_r_shift |>  
  inner_join(list_ggbarplots, by = "id_pop") %>% 
  mutate(icons = case_when(curve_zone_present == "CE" ~ "#B2F2FD",
                           curve_zone_present == "CLP" ~ "#61DCA9",
                           curve_zone_present == "OPS" ~ "#8DB63B",
                           curve_zone_present == "OPD" ~ "#9B7424",
                           curve_zone_present == "HLP" ~ "#944046",
                           curve_zone_present == "HE" ~ "#8C0172")
       ) |> 
  group_by(id_pop) |> 
  mutate(psi = mean(psi))

mytext <- paste(
  "<i>Click to see monthly distribution of performance risk</i>", "<br/>",
  "<i>Species</i>: ", leaflet_risk_shift$species,"<br/>", 
  "<i>Order</i>: ", leaflet_risk_shift$order,"<br/>", 
  "<i>Family</i>: ", leaflet_risk_shift$family,"<br/>", 
  "<i>Reference</i>: ", leaflet_risk_shift$reference, "<br/>", 
  "<i>Performance Shift Index </i>: ", leaflet_risk_shift$psi, "<br/>",
  "<i>Performance raw difference </i>: ", leaflet_risk_shift$preds_diff, "<br/>",
  "<i>Performance historic (%)</i>: ", leaflet_risk_shift$preds_present, "<br/>",
  "<i>Performance shift (%)</i>: ", leaflet_risk_shift$preds_future, "<br/>",
  
  
  sep="") %>%
  lapply(htmltools::HTML)
colors_riskcurve <- c("#B2F2FD", "#61DCA9", "#8DB63B", "#9B7424", "#944046", "#8C0172")

colors_shift <- khroma::color(palette = "bam", reverse = TRUE)(100)[1:100]

palette_riskcurve <- colorFactor(palette = colors_riskcurve,
                                 domain = unique(leaflet_risk_shift$curve_zone_future))
palette_shift <- colorNumeric(palette = colors_shift,
                              domain = unique(leaflet_risk_shift$psi))

intrapests_future_risk_leaflet <- leaflet(data = leaflet_risk_shift) %>% 
  addTiles() %>% 
  addCircleMarkers(lng = ~lon, 
                   lat = ~lat, 
                   stroke = FALSE,
                   fillColor = ~palette_shift(psi), 
                   fillOpacity = 0.8, 
                   color = ~palette_shift(psi), 
                   popup = ~reference,
                   label = mytext,
                   group = "reference",
                   labelOptions = labelOptions( 
                     style = list("font-weight" = "normal", 
                                  padding = "3px 8px"), 
                     textsize = "13px", 
                     direction = "auto")) %>% 
  leafpop::addPopupGraphs(leaflet_risk_shift$ggbarplot , 
                          group = "reference", 
                          width = 400, height = 200) %>% 
  addProviderTiles('CartoDB.DarkMatterNoLabels') #%>%
  #addLegend(pal= palette_shift, 
  #          values= ~psi, 
  #          opacity = 0.9,
  #          title = "Performance shift",
  #          position = "bottomleft")
intrapests_future_risk_leaflet <- intrapests_future_risk_leaflet %>% 
  addLegend(pal= palette_shift, 
            values= ~psi, 
            opacity = 0.9,
            title = "Performance shift",
            position = "bottomleft")

save(intrapests_future_risk_leaflet, file = here("data/data_sink/figures/intrapests_future_risk_leaflet.RData"))
htmlwidgets::saveWidget(intrapests_future_risk_leaflet,
           file = here("data/data_sink/figures/intrapests_future_risk_leaflet.html"))


# 6. Examples Figure insets TPC -------------------------------------------


tpcs_selected_inset <- readxl::read_excel(here("data/data_sink/tpcs_selection_filters_completed.xlsx")) |> 
  filter(step3_uncertainties == "y") |> 
  filter(id_pop %in% c(91, 93, 223))

simulated_tpcs <- tibble()


for(i in c(223, 93, 91)) {
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
  tpcs_selected_inset_i <- tpcs_selected_inset |> 
    filter(id_pop == i)
  model_name_i <- tpcs_selected_inset_i$tpc_model
  fitted_tpc_equations_i <- fit_tpcs(temp = int_rate_i$temperature,
                                     int_rate = int_rate_i$int_rate,
                                     model_name = model_name_i)
  possible_error <- tryCatch(expr = {
    sim_tpcs_i <- predict_curves(temp = int_rate_i$temperature,
                                 int_rate = int_rate_i$int_rate,
                                 fitted_parameters = fitted_tpc_i,
                                 model_name_2boot = fitted_tpc_equations_i,n_boots_samples = 100)
    plot_uncertainties(bootstrap_uncertainties_tpcs = sim_tpcs_i,
                       temp = int_rate_i$temperature,
                       int_rate = int_rate_i$int_rate,
                       species = species_i,
                       reference = reference_i,
                       pop_id = i)
    ggsave(paste0(here("data/data_sink/figures/supplementary_figs/inset"), i, ".svg"),
           width = 2100,
           height = 2100,
           units = "px")
  }, # <- inside tryCatch
  error = function(e) e)
  if (inherits(possible_error, "error")) {
    fit_nls <- NULL
    warning(paste0("Reference ID", i, "(",reference_i, ") could not fit bootstrapped TPCs"))
  }
}

