library(here)
library(raster)
library(sf)
library(ade4)
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
             shape = 0, 
             size = 0.5,
             colour = "white", 
             alpha = 0.1) +
  theme_void()


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
  labs(x = "Longitude", y = "Latitude", 
       title = "Taxonomic Similarity") +
  scale_x_continuous(breaks = seq(-180, 180, 30),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-90, 90, 30),
                     expand = c(0, 0)) +
  coord_sf(ylim = c(-87, 87)) +
  theme(legend.position = "top") 

# patch together
plot_tax <- plot_world +
  inset_element(plot_nmds, 
                0.02, 0.1, 0.2, 0.4) 




# same for functional richness --------------------------------------------

# distance-trait matrix
load(here("data", 
          "distance_trait_matrix.RData"))

# extract functional space axes
dat_pcoa <- dudi.pco(quasieuclid(distance_trait_matrix),
                 scannf = FALSE,
                 nf = 4) %>% 
  pluck("li") %>% 
  # prepare pcoa for joining
  as_tibble(rownames = "species")



# list of all species present per grid 
spp_per_grid <- apply(dat_presabs, 
                      1, 
                      function(x){colnames(dat_presabs)[which(x==1
                      )]}) %>% 
  # remove empty grids
  compact()

# calculate centroid of pcoa per cell
dat_centroid <- spp_per_grid %>% 
  map_dfr(~ .x %>%
            enframe(value = "species",
                    name = NULL) %>% 
            distinct(species) %>% 
            left_join(dat_pcoa, 
                      by = join_by(species)) %>% 
            summarise(A1 = mean(A1), 
                      A2 = mean(A2)), 
          .progress = TRUE) %>% 
  # add coordinates
  bind_cols(select(dat_presabs, contains("itude"))) %>% 
  mutate(long_scld = (longitude_x - min(longitude_x ))/(max(longitude_x)-min(longitude_x )), 
         lat_scld = (latitude_y - min(latitude_y ))/(max(latitude_y)-min(latitude_y)), 
         colour_col = rgb(long_scld, lat_scld, 0.6, 1)) 

# save data
dat_centroid %>% 
  write_rds(here("data", 
                 "mds_score_functional.rds"))


# create nmds grid
dat_grid_func <- expand_grid(
  A1 = seq(min(dat_centroid$A1), max(dat_centroid$A1), length.out = 100),
  A2 = seq(min(dat_centroid$A2), max(dat_centroid$A2), length.out = 100)) %>% 
  mutate(A1_scld = (A1 - min(A1 ))/(max(A1)-min(A1)), 
         A2_scld = (A2 - min(A2 ))/(max(A2)-min(A2)), 
         colour_col = rgb(A1_scld, A2_scld, 0.6, 1))

# visualise nmds space
plot_nmds_func <- dat_grid_func %>%
  ggplot(aes(A1, A2)) +
  geom_raster(aes(fill = colour_col)) +
  scale_fill_identity() +
  geom_point(data = dat_centroid, 
             shape = 0, 
             size = 0.5,
             colour = "white", 
             alpha = 0.1) +
  theme_void()

# worldmap
plot_world_func <- dat_centroid %>%
  mutate(A1_scld = (A1 - min(A1 ))/(max(A1)-min(A1 )), 
         A2_scld = (A2 - min(A2 ))/(max(A2)-min(A2)), 
         colour_col = rgb(A1_scld, A2_scld, 0.6, 1)) %>% 
  # plot
  ggplot() +
  geom_raster(aes(x = longitude_x,
                  y = latitude_y,
                  fill = colour_col)) +
  scale_fill_identity() +
  geom_sf(data = world_map_sf, col = NA, fill = "white", size = 0.1) +
  labs(x = "Longitude", y = "Latitude", 
       title = "Functional Similarity") +
  scale_x_continuous(breaks = seq(-180, 180, 30),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-90, 90, 30),
                     expand = c(0, 0)) +
  coord_sf(ylim = c(-87, 87)) +
  theme(legend.position = "top") 

# patch together
plot_func <- plot_world_func +
  inset_element(plot_nmds_func, 
                0.02, 0.1, 0.2, 0.4) 

# final plot
plot_final <- plot_tax / 
  plot_func +
  plot_annotation(tag_levels = "a")

# save plot
ggsave(plot_final,
       filename = here("figures",
                       "main",
                       "compositional_map.pdf"),
       width = 183, height = 200,
       units = "mm",
       bg = "white")
