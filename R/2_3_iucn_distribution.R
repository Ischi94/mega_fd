library(here)
library(tidyverse)
library(sf)
library(patchwork)

# get functions
source(here("R", 
            "functions.R"))

# read data ---------------------------------------------------------------

# merged presence-absence matrix
dat_presabs <- read_rds(here("data",
                             "presabs_05_res.rds")) %>% 
  # select "grid" cells with >4 species (the trait space is 4D)
  mutate(srich = rowSums(select(., -c(1, 2)))) %>% 
  filter(srich > 4) %>% 
  select(-c(srich, rita_rita)) %>% 
  # correct species names
  rename_with(.cols = -c(1, 2), 
              ~ .x %>% 
                str_to_sentence() %>% 
                str_replace_all("_", " ")) 


# IUCN categories
dat_iucn <- read_csv(here("data",
                          "megafauna_traits.csv")) %>%
  select(species, higher_classification, IUCN)

# World map
world_map_sf <- st_read(here("data", 
                             "worldmap",
                             "ne_10m_land.shp"))

# merge data --------------------------------------------------------------

# prepare data
dat_merged <- dat_presabs %>% 
  pivot_longer(cols = -c(longitude_x, latitude_y), 
               values_to = "presabs", 
               names_to = "species") %>% 
  filter(presabs == 1) %>% 
  left_join(dat_iucn) %>% 
  mutate(ge = case_when(
    IUCN == "LC" ~ 0, 
    IUCN == "NT" ~ 1, 
    IUCN == "VU" ~ 2,
    IUCN == "EN" ~ 3, 
    IUCN == "CR" ~ 4, 
    .default = NA_integer_
  )) 


# average IUCN category per grid ------------------------------------------

plot_av <- dat_merged %>%
  group_by(longitude_x, latitude_y) %>% 
  summarise(mean_ge = mean(ge, na.rm = TRUE)) %>%
  ggplot() +
  geom_raster(aes(x = longitude_x,
                  y = latitude_y, 
                  fill = log(mean_ge))) +
  geom_sf(data = world_map_sf, col = NA, fill = "white", size = 0.1) +
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10,"RdBu")), 
                       name = "Average IUCN category", 
                       breaks = seq(-1.5, 1, by = 0.5), 
                       labels = round(exp(seq(-1.5, 1, by = 0.5)), 1)) +
  labs(x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(-180, 180), ylim= c(-89, 85), expand=FALSE) +
  scale_x_continuous(breaks = seq(-180, 180, 30), 
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-90, 90, 30), 
                     expand = c(0, 0)) 
  

# per category ------------------------------------------------------------

# plot function
plot_metric <- function(metric_col, 
                        label) {
  
  dat_merged %>% 
    filter(IUCN == metric_col) %>% 
    count(longitude_x, latitude_y) %>%
    ggplot() +
    geom_raster(aes(x = longitude_x,
                    y = latitude_y, 
                    fill = n)) +
    geom_sf(data = world_map_sf, col = NA, fill = "grey", size = 0.1) +
    scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10,"RdBu")), 
                         name = label) +
    labs(x = "Longitude", y = "Latitude") +
    coord_sf(xlim = c(-180, 180), ylim= c(-89, 85), expand=FALSE) +
    scale_x_continuous(breaks = seq(-180, 180, 30), 
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(-90, 90, 30), 
                       expand = c(0, 0)) 
}

plot_metric("LC", "Least Concern")

plot_metric("NT", "Near Threatened")

plot_metric("VU", "Vulnerable")

plot_metric("EN", "Endangered")

plot_metric("CR", "Critically Endangered")

dat_merged %>% 
  filter(IUCN == "LC") %>% 
  count(longitude_x, latitude_y) %>% 
  ggplot() +
  geom_raster(aes(x = longitude_x,
                  y = latitude_y, 
                  fill = n)) +
  geom_sf(data = world_map_sf, col = NA, fill = "white", size = 0.1) +
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10,"RdBu")), 
                       name = "Average IUCN category") +
  labs(x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(-180, 180), ylim= c(-89, 85), expand=FALSE) +
  scale_x_continuous(breaks = seq(-180, 180, 30), 
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-90, 90, 30), 
                     expand = c(0, 0)) 
