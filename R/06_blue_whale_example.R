library(here)
library(raster)
library(tidyverse)
library(sf)
library(patchwork)

# get functions
source(here("R", 
            "functions.R"))


# load data ---------------------------------------------------------------


# global metrics
dat_global <- read_rds(here("data", "global_metrics.rds"))

# local metrics
dat_local <- read_rds(here("data", "local_metrics_5_res.rds"))


# world map
world_map_sf <- st_read(here("data", 
                             "worldmap",
                             "ne_10m_land.shp"))

# functional space
dat_space <- read_rds(here("data", "functional_space.rds")) 


# plot per metric -----------------------------------------------------------------------------

# visualise
plot_metric <- dat_local %>% 
  filter(species == "balaenoptera_musculus") %>% 
  pivot_longer(cols = contains("std_local"), 
               names_to = "metric") %>% 
  mutate(metric = word(metric, sep = "_"), 
         metric = ordered(metric, levels = c("FSp", "FDi", "FUn"))) %>% 
  ggplot() +
  geom_raster(aes(x = longitude, 
                  y = latitude, 
                  fill = value)) +
  geom_sf(data = world_map_sf,
          col = NA,
          fill = "grey90",
          size = 0.1) +
  scale_fill_gradient2(
    low = colour_mint,
    mid = colour_purple,
    high = colour_yellow,
    midpoint = 0.5, 
    name = NULL) +
  scale_y_continuous(breaks = seq(-90, 90, 30),
                     expand = c(0, 0)) +
  labs(x = "Longitude", y = "Latitude") +
  coord_sf(expand = FALSE) +
  facet_wrap(~metric, ncol = 1)





# functional space ----------------------------------------------------------------------------


# create functional space plots
list_spaces <- list(max, min) %>%
  map2(
    .x = .,
    .y = c("Highest FUSE cell", "Lowest FUSE cell"),
    .f = ~ dat_local %>%
      filter(species == "balaenoptera_musculus") %>% 
      filter(FUSE_local == .x(FUSE_local, na.rm = TRUE)) %>%
      left_join(dat_local %>% 
                  select(species_pres = species, longitude, latitude), 
                by = c("longitude", "latitude")) %>% 
      distinct(species = species_pres) %>%
      left_join(dat_space %>%
                  as_tibble(rownames = "species")) %>%
      mutate(
        colour_id = if_else(species == "balaenoptera_musculus", "yes", "no"),
        A1_center = mean(A1),
        A2_center = mean(A2)
      ) %>%
      arrange(colour_id) %>%
      ggplot(aes(A1, A2)) +
      geom_polygon(
        fill = NA,
        colour = "grey20",
        data = . %>%
          slice(chull(A1, A2))
      ) +
      geom_point(
        aes(fill = colour_id),
        shape = 21,
        size = 2,
        stroke = 1,
        colour = "grey20"
      ) +
      geom_point(
        aes(A1_center, A2_center),
        fill = colour_yellow,
        size = 4,
        stroke = 1.5,
        shape = 21,
        colour = "white"
      ) +
      scale_fill_manual(values = c(colour_grey, colour_coral)) +
      labs(
        x = "PCoA1",
        y = "PCoA2",
        title = .y,
        colour = NULL
      ) +
      theme(legend.position = "none")
  )



# patch together
plot_bwhale <- plot_metric +
  free(list_spaces[[1]] / list_spaces[[2]]) +
  plot_annotation(tag_levels = "A")

# save plot
ggsave(plot_bwhale, 
       filename = here("figures",
                       "main", 
                       "fig_3.svg"),
       width = 250, height = 150,
       units = "mm")
