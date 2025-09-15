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

int_rate_data_info <- int_rate_data |> 
  group_by(id_pop) |> 
  slice(1) |> 
  rename(pop_id = id_pop) |> 
  select(reference, species, order, family, lat, lon, pop_id)

thermal_limits <- read_rds(here("data/data_sink/therm_lims_est_42.rds")) |> 
  mutate(pop_id = as_factor(pop_id)) 

thermal_limits_intrate <- inner_join(int_rate_data_info, thermal_limits) |> 
  filter(is.finite(tmax_est))


# 2. Estimates  ----
## a) tmin  ----
tmin_est_int <- lmerTest::lmer(tmin_est ~ 1 + (1|species) + (1|reference),
                               data = thermal_limits_intrate)
summary(tmin_est_int) # 7.067 +/- 0.384
summary(thermal_limits_intrate$tmin_est) # 6.737, median = 5.472
sd(thermal_limits_intrate$tmin_est)/sqrt(238)# se = 0.336

## b) tmax  ----

tmax_est_int <- lmerTest::lmer(tmax_est ~ 1 + (1|species) + (1|reference),
                               data = thermal_limits_intrate)
summary(tmax_est_int) # 34.876 +/- 0.254
summary(thermal_limits_intrate$tmax_est) # mean = 34.98, median = 35.25

tmax_est_int <- lmerTest::lmer(tmax_est ~ 1 + (1|species) + (1|reference),
                               data = thermal_limits_intrate)
summary(tmax_est_int) # 34.876 +/- 0.254
summary(thermal_limits_intrate$tmax_est) # mean = 34.98, median = 35.25


# 4. Orders  ----
## a) tmin  ----
tmin_est_int_order <- lmerTest::lmer(tmin_est ~ order + (1|species) + (1|reference),
                               data = thermal_limits_intrate)
summary(tmin_est_int_order) # 7.067 +/- 0.384
anova(tmin_est_int_order)
tmax_est_int <- lmerTest::lmer(tmax_est ~ 1 + (1|species) + (1|reference),
                               data = thermal_limits_intrate)
summary(tmax_est_int) # 34.876 +/- 0.254
summary(thermal_limits_intrate$tmax_est) # mean = 34.98, median = 35.25

## b) tmax  ----
tmax_est_int_order <- lmerTest::lmer(tmax_est ~ order + (1|species) + (1|reference),
                                     data = thermal_limits_intrate)
summary(tmax_est_int_order) # 
anova(tmax_est_int_order)
tmax_est_int <- lmerTest::lmer(tmax_est ~ 1 + (1|species) + (1|reference),
                               data = thermal_limits_intrate)
summary(tmax_est_int) # 34.876 +/- 0.254
summary(thermal_limits_intrate$tmax_est) # mean = 34.98, median = 35.25

tmax_est_int <- lmerTest::lmer(tmax_est ~ 1 + (1|species) + (1|reference),
                               data = thermal_limits_intrate)
summary(tmax_est_int) # 34.876 +/- 0.254
summary(thermal_limits_intrate$tmax_est) # mean = 34.98, median = 35.25

# 5. F-test (see Herrando-Perez, 2019) ----
variances_therm_lims <- thermal_limits_intrate |> 
  select(pop_id, tmin_est, tmax_est) |> 
  pivot_longer(cols = c("tmin_est", "tmax_est"),
               values_to = "estimates",
               names_to = "thermal_limit") |> 
  filter(is.finite(estimates))
  

#test variance differences
vartest_estimates <- var.test(estimates ~ thermal_limit,
                              variances_therm_lims, 
                              alternative = "greater") 
print(vartest_estimates) #not significative

# 3.  Levene's test (see Hoffmann, 2013) ----
car::leveneTest(y = estimates ~ thermal_limit,
                data = variances_therm_lims) 


#test variance differences
vartest_estimates <- var.test(estimates ~ thermal_limit,
                              variances_therm_lims, 
                              alternative = "greater") 
print(vartest_estimates) #not significative

# 4.  visual inspection ----

trait_clouds <- ggplot(thermal_limits_intrate |> 
                         filter(is.finite(tmax_est)), aes(x = 1.5)) + 
  ggdist::stat_halfeye(aes(y = tmin_est),
                       adjust = .5,
                       width = .5, ## set slab interval to show IQR and 95% data range
                       color = "#0a9396",
                       fill = "#94d2bd"
  ) + 
  ggdist::stat_dots(aes(y = tmin_est),
                    side = "left",
                    dotsize = 1.2,
                    justification = 1.05,
                    binwidth = .5,
                    color = "#0a9396",
                    fill= "#94d2bd",
                    alpha = 1
  )+
  coord_cartesian(xlim = c(1, NA))+
  ggdist::stat_halfeye(aes(y = tmax_est),
                       adjust = .5,
                       width = .5, ## set slab interval to show IQR and 95% data range
                       color = "#bb3e03",
                       fill = "#ee9b00"
  ) + 
  ggdist::stat_dots(aes(y = tmax_est),
                    side = "left",
                    dotsize = 1.2,
                    justification = 1.05,
                    binwidth = .5,
                    color = "#bb3e03",
                    fill= "#ee9b00",
                    alpha = 1
  )+
  ggthemes::theme_few()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  labs(x = element_blank(),
       y = "Temperature (ºC)")
trait_clouds
ggsave(here("data/data_sink/figures/clouds.svg"), height = 1400, width = 1400,
       units = "px")

trait_clouds
  facet_wrap(~order)
  ggsave(here("data/data_sink/figures/clouds_tlims_order.png"), height = 1400, width = 1400,
         units = "px")
  