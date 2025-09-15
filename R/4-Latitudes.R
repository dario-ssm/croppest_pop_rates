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

int_rate_data_info <- int_rate_data |> 
  group_by(id_pop) |> 
  slice(1) |> 
  rename(pop_id = id_pop) |> 
  select(reference, species, order, family, lat, lon, pop_id)

thermal_limits <- read_rds(here("data/data_sink/therm_lims_est_42.rds")) |> 
  mutate(pop_id = as_factor(pop_id))

thermal_limits_intrate <- inner_join(int_rate_data_info, thermal_limits)
tpcs_selected <- readxl::read_excel(here("data/data_sink/tpcs_selection_filters_completed.xlsx")) |> 
  filter(step3_uncertainties == "y")

# 2. Tmin ~ Lat -----------------------------------------------------------------
thermal_limits_lats <- thermal_limits_intrate|> 
  filter(is.finite(tmax_est))

## a) 1|reference -----------------------------------------------------------------
tmin_lat_reference <- lmerTest::lmer(tmin_est ~ abs(lat) + (1|reference),
                                     data = thermal_limits_lats)
summary(tmin_lat_reference)

## b) 1|species -----------------------------------------------------------------
tmin_lat_species <- lmerTest::lmer(tmin_est ~ abs(lat) + (1|species),
                                     data = thermal_limits_lats)
summary(tmin_lat_species)

## c) 1|ref + 1|sp -----------------------------------------------------------------
tmin_lat_sp_ref <- lmerTest::lmer(tmin_est ~ abs(lat) + (1|species) + (1|reference),
                                   data = thermal_limits_lats)
summary(tmin_lat_sp_ref)

## d) weights  -----------------------------------------------------------------
w_tmin_est_intrate <- thermal_limits_intrate |> 
  select(tmin_est, tmin_se, lat, species, reference) |> 
  mutate(weights_tmin = 1/((tmin_se)^2)) |> 
  filter(is.finite(weights_tmin))

w_tmin_lat_sp_ref <- lmerTest::lmer(tmin_est ~ abs(lat) + (1|species) + (1|reference),
                                  data = w_tmin_est_intrate,
                                  weights = weights_tmin) 
summary(w_tmin_lat_sp_ref)

## e) plot  -----------------------------------------------------------------
tmin_lat_plot <- sim_and_plot_linears(model_object = tmin_lat_sp_ref,
                                      var_x = abs(thermal_limits_lats$lat),
                                      var_y = thermal_limits_lats$tmin_est,
                                      n_sims = 1000,
                                      your_title = "Minimum",
                                      your_subtitle = "N = 238",
                                      lab_x = "absolute latitude (º)",
                                      lab_y = "thermal limit (ºC)",
                                      color_points = "#0a9396",
                                      color_central = "#005f73",
                                      color_uncertainty = "#94d2bd")
print(tmin_lat_plot)
ggsave(here("data/data_sink/figures/tmin_lat.png"), height = 2600, width = 2600,
       units = "px")

## f) order +1|ref + 1|sp -----------------------------------------------------------------
tmin_lat_order_sp_ref <- lmerTest::lmer(tmin_est ~ abs(lat) + order + (1|species) + (1|reference),
                                        data = thermal_limits_lats)
summary(tmin_lat_order_sp_ref)
anova(tmin_lat_order_sp_ref)

# 3. tmax ~ Lat -----------------------------------------------------------------

## a) 1|reference -----------------------------------------------------------------
tmax_lat_reference <- lmerTest::lmer(tmax_est ~ abs(lat) + (1|reference),
                                     data = thermal_limits_lats) 

summary(tmax_lat_reference)

## b) 1|species -----------------------------------------------------------------
tmax_lat_species <- lmerTest::lmer(tmax_est ~ abs(lat) + (1|species),
                                   data = thermal_limits_lats)
summary(tmax_lat_species)

## c) 1|ref + 1|sp -----------------------------------------------------------------
tmax_lat_sp_ref <- lmerTest::lmer(tmax_est ~ abs(lat) + (1|species) + (1|reference),
                                  data = thermal_limits_lats)
summary(tmax_lat_sp_ref)

## d) weights  -----------------------------------------------------------------
w_tmax_est_intrate <- thermal_limits_intrate |> 
  select(tmax_est, tmax_se, lat, species, reference) |> 
  mutate(weights_tmax = 1/((tmax_se)^2)) |> 
  filter(is.finite(weights_tmax))

w_tmax_lat_sp_ref <- lmerTest::lmer(tmax_est ~ abs(lat) + (1|species) + (1|reference),
                                    data = w_tmax_est_intrate,
                                    weights = weights_tmax) 
summary(w_tmax_lat_sp_ref)

## e) plot  -----------------------------------------------------------------
tmax_lat_plot <- sim_and_plot_linears(model_object = tmax_lat_sp_ref,
                                      var_x = abs(thermal_limits_lats$lat),
                                      var_y = thermal_limits_lats$tmax_est,
                                      n_sims = 1000,
                                      your_title = "Maximum",
                                      your_subtitle = "N = 236",
                                      lab_x = "absolute latitude (º)",
                                      lab_y = "thermal limit (ºC)",
                                      color_points = "#bb3e03",
                                      color_central = "#9b2226",
                                      color_uncertainty = "#ee9b00")
print(tmax_lat_plot)
ggsave(here("data/data_sink/figures/tmax_lat.png"), height = 2600, width = 2600,
       units = "px")

## f) order +1|ref + 1|sp -----------------------------------------------------------------
tmax_lat_order_sp_ref <- lmerTest::lmer(tmax_est ~ abs(lat) + order + (1|species) + (1|reference),
                                  data = thermal_limits_lats)
summary(tmax_lat_order_sp_ref)
anova(tmax_lat_order_sp_ref)

# 4. Composed plots --------------------------------------------------------

## a) pool -----------------------------------------------------------------
tmin_tmax_lat <- cowplot::plot_grid(tmin_lat_plot, tmax_lat_plot)
tmin_tmax_lat
ggsave(here("data/data_sink/figures/tmin_tmax_lat.png"), height = 2600, width = 3200,
       units = "px")
ggsave(here("data/data_sink/figures/tmin_tmax_lat.svg"), height = 1500, width = 1900,
       units = "px")

therm_lims_vert <- thermal_limits_lats |> 
  rename(tmin = tmin_est,
         tmax = tmax_est) |> 
  pivot_longer(cols = c("tmin", "tmax"),
               values_to = "estimates",
               names_to = "thermal_limit")

tmin_tmax_sameplot <- ggplot(data = therm_lims_vert,
                             aes(x = abs(lat),
                                 y = estimates,
                                 color = thermal_limit,
                                 fill = thermal_limit))+
  geom_point(size = 1.5,
             alpha = .75)+
  geom_smooth(method = "lm")+
  theme_bw()+
  scale_fill_manual(values = c("#ee9b00", "#94d2bd"))+
  scale_color_manual(values = c("#bb3e03", "#0a9396"))+
  labs(fill = NULL,
       color = NULL,
       x = "absolute latitude (º)",
       y = "thermal limit (ºC)")
tmin_tmax_sameplot
ggsave(here("data/data_sink/figures/tmax_tmax_sameplot.png"), height = 1400, width = 1000,
       units = "px")

## b) order -----------------------------------------------------------------
tmin_tmax_order_plot <- ggplot(data = therm_lims_vert,
                             aes(x = abs(lat),
                                 y = estimates,
                                 color = thermal_limit,
                                 fill = thermal_limit))+
  geom_point(size = 1.5,
             alpha = .75)+
  geom_smooth(method = "lm")+
  theme_bw()+
  scale_fill_manual(values = c("#ee9b00", "#94d2bd"))+
  scale_color_manual(values = c("#bb3e03", "#0a9396"))+
  facet_wrap(.~order, scales = "free_y")+
  labs(fill = NULL,
       color = NULL,
       x = "absolute latitude (º)",
       y = "thermal limit (ºC)")
tmin_tmax_order_plot
ggsave(here("data/data_sink/figures/tmax_tmax_ordert.png"), height = 3600, width = 3600,
       units = "px")

## c) pool hemispheres -----------------------------------------------------------------
tmin_tmax_hemis_plot <- ggplot(data = therm_lims_vert,
                               aes(x = lat,
                                   y = estimates,
                                   color = thermal_limit,
                                   fill = thermal_limit))+
  geom_point(size = 1.5,
             alpha = .75)+
  geom_smooth(method = "loess")+
  theme_bw()+
  scale_fill_manual(values = c("#ee9b00", "#94d2bd"))+
  scale_color_manual(values = c("#bb3e03", "#0a9396"))+
  labs(fill = NULL,
       color = NULL,
       x = "absolute latitude (º)",
       y = "thermal limit (ºC)")
tmin_tmax_hemis_plot
ggsave(here("data/data_sink/figures/tmin_tmax_hemis_plot.png"), height = 3600, width = 3600,
       units = "px")

## d) facet hemispheres -----------------------------------------------------------------
therm_lims_vert_hemispheres <- therm_lims_vert |> 
  mutate(hemisphere = case_when(lat <= 0 ~ "Northern Hemisphere",
                                lat > 0 ~ "Southern Hemisphere")) |> 
  drop_na()
tmin_tmax_hemis_facplot <- ggplot(data = therm_lims_vert_hemispheres,
                               aes(x = lat,
                                   y = estimates,
                                   color = thermal_limit,
                                   fill = thermal_limit))+
  geom_point(size = 1.5,
             alpha = .75)+
  geom_smooth(method = "lm")+
  facet_wrap(.~hemisphere,
             scales = "free_x")+
  theme_bw()+
  scale_fill_manual(values = c("#ee9b00", "#94d2bd"))+
  scale_color_manual(values = c("#bb3e03", "#0a9396"))+
  labs(fill = NULL,
       color = NULL,
       x = "absolute latitude (º)",
       y = "thermal limit (ºC)")
tmin_tmax_hemis_facplot
ggsave(here("data/data_sink/figures/tmax_tmax_ordert.png"), height = 3600, width = 3600,
       units = "px")

## e) TPC model equation -----------------------------------------------------------------
tpcs_selected_eqs <- tpcs_selected |> 
  select(id_pop, tpc_model) |> 
  rename(pop_id = id_pop)
therm_lims_vert_mods <- inner_join(therm_lims_vert,
                                   tpcs_selected_eqs) |> 
  drop_na()
tmin_tmax_equation_plot <- ggplot(data = therm_lims_vert_mods,
                               aes(x = abs(lat),
                                   y = estimates,
                                   color = thermal_limit,
                                   fill = thermal_limit))+
  geom_point(size = 1.5,
             alpha = .75)+
  geom_smooth(method = "lm")+
  theme_bw()+
  scale_fill_manual(values = c("#ee9b00", "#94d2bd"))+
  scale_color_manual(values = c("#bb3e03", "#0a9396"))+
  facet_wrap(.~tpc_model, scales = "free_y")+
  labs(fill = NULL,
       color = NULL,
       x = "absolute latitude (º)",
       y = "thermal limit (ºC)")
tmin_tmax_equation_plot
ggsave(here("data/data_sink/figures/tmin_tmax_equation_plot.png"), height = 3600, width = 3600,
       units = "px")



# 5. thermal breadth  -----------------------------------------------------------------
thermal_limits_breadth_lat <- thermal_limits_lats |> 
  mutate(thermal_breadth = tmax_est - tmin_est)

thermal_breadth_lat <- lmerTest::lmer(thermal_breadth ~ abs(lat) + (1|reference) + (1|species),
                                      data = thermal_limits_breadth_lat)
summary(thermal_breadth_lat)
thermal_breadth_lat_plot <- sim_and_plot_linears(model_object = thermal_breadth_lat,
                                      var_x = abs(thermal_limits_breadth_lat$lat),
                                      var_y = thermal_limits_breadth_lat$thermal_breadth,
                                      n_sims = 1000,
                                      your_title = "Thermal Breadth across latitudes",
                                      your_subtitle = "Tmin - Tmax,  N = 236",
                                      lab_x = "absolute latitude (º)",
                                      lab_y = "thermal breadth (ºC)",
                                      color_points = "#b56576",
                                      color_central = "#355070",
                                      color_uncertainty = "#eaac8b")
print(thermal_breadth_lat_plot)
ggsave(here("data/data_sink/figures/thermal_breadth_lat_plot.png"), height = 2600, width = 2600,
       units = "px")

