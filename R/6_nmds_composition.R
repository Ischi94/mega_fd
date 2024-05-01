library(here)
library(raster)
library(sf)
library(tidyverse)
library(vegan)
library(patchwork)

# get functions
source(here("R", 
            "functions.R"))

# load data ---------------------------------------------------------------


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

# world map
world_map_sf <- st_read(here("data", 
                             "worldmap",
                             "ne_10m_land.shp"))

# calculate distance matrix -----------------------------------------------

# take presence absence matrix with 0.5x0.5 resolution and 
# upscale to 5x5 resolution
dat_rast_agg <- dat_presabs %>%
  rasterFromXYZ() %>%
  aggregate(fact = 10, 
            fun = max)

# clean up
dat_upscaled <- dat_rast_agg %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  mutate(across(everything(), 
                ~ ifelse(.x > 0, 1, .x))) %>% 
  # add_column(longitude_x = coordinates(dat_rast_agg)[, 1], 
  #            latitude_y = coordinates(dat_rast_agg)[, 2], 
  #            .before = 0) %>% 
  # correct species names
  rename_with(.cols = -c(1, 2), 
              ~ .x %>% 
                str_to_sentence() %>% 
                str_replace_all("\\.", " ")) %>% 
  # remove empty cells
  filter(rowSums(across(where(is.numeric))) != 0)

# based on jaccard indices
dat_dist <- dat_upscaled %>% 
  vegdist(method = "jaccard")



# non-metric multidimensional scaling -------------------------------------

# nmds
dat_nmds <-
  metaMDS(dat_dist,
          distance = "jaccard",
          parallel = 10,
          k = 2,
          maxit = 999, 
          trymax = 200,
          wascores = TRUE)


  


# visualise ---------------------------------------------------------------

# non-empty sites
non_empty <- dat_rast_agg %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  reframe(non_empt = rowSums(across(where(is.numeric))) != 0) %>% 
  replace_na(list(non_empt = FALSE)) %>% 
  pull(non_empt)  

rgb(0, 1, 0.5, 0.3)

# assign coordinates
dat_nmds <- dat_nmds$points %>% 
  as_tibble() %>% 
  add_column(longitude_x = coordinates(dat_rast_agg)[non_empty, 1],
             latitude_y = coordinates(dat_rast_agg)[non_empty, 2],
             .before = 0) %>% 
  mutate(long_scld = (longitude_x - min(longitude_x ))/(max(longitude_x)-min(longitude_x )), 
         lat_scld = (latitude_y - min(latitude_y ))/(max(latitude_y)-min(latitude_y)), 
         colour_col = rgb(long_scld, lat_scld, 0.6, 1)) 

# save data
dat_nmds %>% 
  write_rds(here("data", 
                 "mds_score.rds"))


# create nmds grid
dat_grid <- expand_grid(
  MDS1 = seq(min(dat_nmds$MDS1), max(dat_nmds$MDS1), length.out = 100),
  MDS2 = seq(min(dat_nmds$MDS2), max(dat_nmds$MDS2), length.out = 100)) %>% 
  mutate(MDS1_scld = (MDS1 - min(MDS1 ))/(max(MDS1)-min(MDS1 )), 
         MDS2_scld = (MDS2 - min(MDS2 ))/(max(MDS2)-min(MDS2)), 
         colour_col = rgb(MDS1_scld, MDS2_scld, 0.6, 1))

# visualise nmds space
plot_nmds <- dat_grid %>%
  ggplot(aes(MDS1, MDS2)) +
  geom_raster(aes(fill = colour_col)) +
  scale_fill_identity() +
  geom_point(data = dat_nmds, 
             shape = 21,
             stroke = 0.2,
             colour = "white") 

# save plot
ggsave(plot_nmds,
       filename = here("figures",
                       "main",
                       "compositional_map_nmds.pdf"),
       width = 183, height = 100,
       units = "mm",
       bg = "white")

# worldmap
plot_world <- dat_nmds %>%
  mutate(MDS1_scld = (MDS1 - min(MDS1 ))/(max(MDS1)-min(MDS1 )), 
         MDS2_scld = (MDS2 - min(MDS2 ))/(max(MDS2)-min(MDS2)), 
         colour_col = rgb(MDS1_scld, MDS2_scld, 0.6, 1)) %>% 
  # plot
  ggplot() +
  geom_raster(aes(x = longitude_x,
                  y = latitude_y,
                  fill = colour_col)) +
  scale_fill_identity() +
  geom_sf(data = world_map_sf, col = NA, fill = "white", size = 0.1) +
  labs(x = "Longitude", y = "Latitude") +
  scale_x_continuous(breaks = seq(-180, 180, 30),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-90, 90, 30),
                     expand = c(0, 0)) +
  coord_sf(ylim = c(-87, 87)) +
  theme(legend.position = "top")



# save plot
ggsave(plot_world,
       filename = here("figures",
                       "main",
                       "compositional_map.pdf"),
       width = 183, height = 100,
       units = "mm",
       bg = "white")
