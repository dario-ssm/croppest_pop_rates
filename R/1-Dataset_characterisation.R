library(tidyverse)
library(here)
library(readxl)
library(ggrepel)
library(sf)
library(leaflet)
# 1. Load -----------------------------------------------------------------
int_rate_data <- read_excel(here("data/data_source/int_rate_dataset_new.xlsx")) |> 
  mutate(reference = as_factor(reference),
         id_pop = as_factor(id_pop)) |> 
  separate(order, into = c("order", "family"), sep = ": ", remove = TRUE)


# id_pop 72 had been manually calculated since the rate was extremely high for a lepidopteran species, seemd a typo error. 
# the values included were recalculated based on the table as r_m = ln(R_0)/T

# 2. Vars. exploration -----------------------------------------------------------------

## a) temperatures ---------------------------------------------------------

### no. temperature treatments
temp_treats <- int_rate_data |> 
  group_by(id_pop) |> 
  distinct(temperature) |> 
  count(name = "temp_treatments")

tempcount_summary <- temp_treats |> 
  ungroup() |> 
  count(temp_treatments) |>
  mutate(perc = 100*n/312) |> 
  mutate(temp_treatments = as_factor(temp_treatments)) |> 
  print()

ggplot(tempcount_summary, aes(x = temp_treatments,
                              y = n,
                              fill = temp_treatments))+
  geom_bar(stat = "identity")+
  khroma::scale_fill_acton(reverse = T,
                           discrete = T)+
  ggthemes::theme_clean()+
  theme(legend.position = "none")+
  labs(x = "Number of temperature treatments in growth chambers",
       y = "No. populations")
ggsave(here("data/data_sink/figures/supplementary_figs/temperatures_dist.png"),
       width = 1600,
       height = 1600,
       units = "px")

### fluctuating vs. constant vs alternating
temp_treats <- int_rate_data |> 
  group_by(id_pop) |>
  slice(1) |> 
  ungroup() |> 
  count(therm_regime) |> 
  print()

### temperatures distribution
ggplot(int_rate_data, aes(x = temperature))+
  geom_histogram(fill = "lightcoral",
                 color = "gray23", linetype = "dashed")+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x = "Controlled temperature in the growth chamber (ºC)")+
  scale_x_continuous(breaks = c(10, 15, 20, 25, 30, 35, 40))
ggsave(here("data/data_sink/figures/supplementary_figs/temperatures_chambers.png"),
       width = 1600,
       height = 1600,
       units = "px")

## b) intrinsic rates of increase ---------------------------------------------------------
### int_rates distribution
int_rate_title <- expression(italic(r)[m]~(day^-1))

ggplot(int_rate_data, aes(x = int_rate))+
  geom_histogram(fill = "darkcyan",
                 color = "gray73", linetype = "dashed",bins = 90)+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x = int_rate_title)

ggsave(here("data/data_sink/figures/supplementary_figs/int_rate_dist.png"),
       width = 1600,
       height = 1600,
       units = "px")

## c) Taxonomic orders ---------------------------------------------------------
###  distribution
int_rate_taxa <- int_rate_data |>  
  group_by(id_pop) |> 
  slice(1) |>
  ungroup() |> 
  count(order, family)

ggplot(int_rate_taxa, aes(y = fct_reorder(family, -n),
                          x = n,
                          fill = order))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("#3E90B5", "#E7CB65", "#E57742", "#AA9ABE",
                               "#A8B85D", "#EB5B4B", "#DECBD6", "#414F44"))+
  theme_bw()+
  labs(x = "Number of populations in the dataset",
       y = "Family")

ggsave(here("data/data_sink/figures/supplementary_figs/int_rate_taxa_family.png"),
       width = 1600,
       height = 1600,
       units = "px")


ggplot(int_rate_taxa, aes(y = fct_reorder(family, -n),
                          x = n,
                          fill = order))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("#3E90B5", "#E7CB65", "#E57742", "#AA9ABE",
                               "#A8B85D", "#EB5B4B", "#DECBD6", "#414F44"))+
  theme_bw()+
  facet_wrap(.~order, scales = "free_y")+
  labs(x = "Number of populations in the dataset",
       y = "Family")

ggsave(here("data/data_sink/figures/supplementary_figs/int_rate_taxa_family_facets.png"),
       width = 2800,
       height = 2800,
       units = "px")

## c) Source locations ---------------------------------------------------------
int_rate_sf <- int_rate_data |> 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, na.fail = FALSE) |> 
  st_transform(crs = "+proj=wintri")  

worldmap <- rnaturalearth::ne_countries(scale = 50) |> 
  st_transform(crs = "+proj=robin")


ggplot()+
  geom_sf(data = worldmap, fill = "gray74", color = NA)+
  geom_sf(data = int_rate_sf,
          aes(color = order))+
  theme_bw()+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#3E90B5", "#E7CB65", "#E57742", "#AA9ABE",
                               "#A8B85D", "#EB5B4B", "#DECBD6", "#414F44"))


ggsave(here("data/data_sink/figures/supplementary_figs/source_map.png"),
       width = 2800,
       height = 2800,
       units = "px")
#map it with leaflet

int_rate_4leaflet <- int_rate_data |> 
group_by(id_pop) |> 
  slice(1)


mytext <- paste(
  "<b>Species</b>: ", int_rate_4leaflet$species,"<br/>", 
  "<b>Reference</b>: ", int_rate_4leaflet$reference, "<br/>", 
  "<i>Order</i>: ", int_rate_4leaflet$order, "<br/>",
  "<i>Family</i>: ", int_rate_4leaflet$family, "<br/>",
  "<i>Host Plant</i>: ", int_rate_4leaflet$host_plant, "<br/>",
  sep="") %>%
  lapply(htmltools::HTML)

colors_order<-c("#3E90B5", "#E7CB65", "#E57742", "#AA9ABE",
                 "#A8B85D", "#EB5B4B", "#DECBD6", "#414F44")

palette_order <- colorFactor(palette = colors_order,
                              domain = unique(int_rate_4leaflet$order))

set.seed(2025)
intrapests_order<- leaflet(data = int_rate_4leaflet) %>% 
  addTiles() %>% 
  addCircleMarkers(lng = ~jitter(lon, factor = 1), 
                   lat = ~jitter(lat, factor = 1), 
                   stroke = FALSE,
                   fillColor = ~palette_order(order), 
                   fillOpacity = 0.75, 
                   color = ~palette_order(order), 
                   popup = ~reference,
                   label = mytext,
                   group = "reference",
                   labelOptions = labelOptions( 
                     style = list("font-weight" = "normal", 
                                  padding = "3px 8px"), 
                     textsize = "13px", 
                     direction = "auto")) %>% 
  addProviderTiles('CartoDB.DarkMatterNoLabels') %>%
  addLegend(pal= palette_order, 
            values= ~order, 
            opacity = 0.9,
            title = "Distribution of populations in the dataset",
            position = "bottomleft")

intrapests_order
htmlwidgets::saveWidget(intrapests_r_shift,
                        here("data/data_sink/figs/intrapests_r_shift.html"))

