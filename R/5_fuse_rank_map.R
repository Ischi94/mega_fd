library(here)
library(raster)
library(tidyverse)
library(ade4)
library(sf)
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



# distance-trait matrix
load(here("data", 
          "distance_trait_matrix.RData"))

# extract functional space axes
pcoa <- dudi.pco(quasieuclid(distance_trait_matrix),
                 scannf = FALSE,
                 nf = 4)

# World map
world_map_sf <- st_read(here("data", 
                             "worldmap",
                             "ne_10m_land.shp"))

# per cell proportion values
dat_prop <- read_rds(here("data", 
          "proportional_coverage.rds")) %>% 
  filter(cell_res == "1x1")

# assign coordinates -----------------------------------------------

# take presence absence matrix with 0.5x0.5 resolution and 
# upscale to 1x1 resolution
dat_rast_agg <- dat_presabs %>%
  rasterFromXYZ() %>%
  aggregate(fact = 2, 
            fun = max)

# clean up
dat_1x1 <- dat_rast_agg %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  mutate(across(everything(), 
                ~ ifelse(.x > 0, 1, .x))) %>% 
  add_column(longitude_x = coordinates(dat_rast_agg)[, 1], 
             latitude_y = coordinates(dat_rast_agg)[, 2], 
             .before = 0) %>% 
  # correct species names
  rename_with(.cols = -c(1, 2), 
              ~ .x %>% 
                str_to_sentence() %>% 
                str_replace_all("\\.", " "))

# get empty cell indices
empty_cells <- apply(dat_1x1,
                     1,
                     function(x) {
                       colnames(dat_1x1)[which(x == 1)]
                     })  %>%
  # remove empty grids
  map_lgl(is_empty)

# assign values to cells and NA to empty cells
dat_map <- dat_1x1[!empty_cells, c(1, 2)] %>%
  add_column(fuse_prop = dat_prop$fuse_prop) %>%
  bind_rows(dat_1x1[empty_cells, c(1, 2)] %>%
              add_column(fuse_prop = NA_real_)) %>%
  arrange(longitude_x, desc(latitude_y))

# visualise  --------------------------------------------------------------

# create map
plot_map <-dat_map %>%
  # plot
  ggplot() +
  geom_raster(aes(x = longitude_x,
                  y = latitude_y,
                  fill = fuse_prop)) +
  scale_fill_gradient2(low = colour_mint,
                       mid = colour_purple,
                       midpoint = 0.2,
                       high = colour_yellow,
                       na.value = "white",
                       breaks = c(0, 0.2, 0.4),
                       labels = c("0%", "20%", "40%"),
                       name = "FUSE\nGlobal-local \nagreement ") +
  geom_sf(data = world_map_sf, col = NA, fill = "white", size = 0.1) +
  labs(x = "Longitude", y = "Latitude") +
  scale_x_continuous(breaks = seq(-180, 180, 30),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-90, 90, 30),
                     expand = c(0, 0)) +
  coord_sf(ylim = c(-87, 87)) +
  theme(legend.position = "top")


# save plot
ggsave(plot_map,
       filename = here("figures",
                       "main",
                       "agreement_map.pdf"),
       width = 183, height = 100,
       units = "mm",
       bg = "white")

# metrices per grid ------------------------------------------------

# read in metrices and summarise per latitude
dat_metrics <- read_rds(here("data",
                             "functional_metrics.rds")) %>% 
  group_by(latitude_y) %>% 
  summarise(across(c(spR, FRic, fuse, uniq, special), mean)) %>% 
  mutate(across(c(spR, FRic, fuse, uniq, special), ~scale(.x))) %>% 
  pivot_longer(cols = -latitude_y, 
               names_to = "metric") 


plot_latid <- dat_metrics %>%
  ggplot(aes(value, latitude_y, 
             colour = metric)) +
  geom_path(linewidth = 0.5) +
  geom_path(data = dat_map %>% 
              group_by(latitude_y) %>% 
              summarise(value = mean(fuse_prop, 
                                     na.rm = TRUE)) %>% 
              mutate(value = scale(value)) %>% 
              add_column(metric = "fuse_prop"), 
            linewidth = 1.1, 
            colour = "grey10") +
  scale_color_manual(values = c(colour_coral, 
                                colour_mint, 
                                colour_yellow, 
                                colour_purple, 
                                "brown"), 
                     labels = c("FRic", "FUSE", "FSp", 
                                "SR", "FUn"),
                     name = "Metric") +
  annotate(size = 10/.pt, 
           "text", 
           label = "FUSE\nGlobal-local \nagreement", 
           x = 1.5, 
           y = -35) +
  labs(y = "Latitude", 
       x = "Z-score") +
  theme(legend.position = c(0.9, 0.5))

# patch together
plot_final <- plot_map/ plot_latid

# save plot
ggsave(plot_final,
       filename = here("figures",
                       "main",
                       "5_agreement_map.pdf"),
       width = 183, height = 100*2,
       units = "mm",
       bg = "white")
# 
# # computes species richness and functional richness per grid cell
# fd_metrics <- spp_per_grid %>% 
#   map(~get_FV_Sp(ax = c(1:4),
#                  pcoa = pcoa,
#                  Selected_sp = .x), 
#       .progress = TRUE)
# 
# 
# # species richness
# vec_spR <- fd_metrics %>%
#   map_dbl(~pluck(.x, "data") %>% 
#             pull("RS"))
# 
# # functional richness
# vec_FRic <- fd_metrics %>%
#   map_dbl( ~ pluck(.x, "data") %>%
#              pull("FRic"))
# 
# # assign values to cells and NA to empty cells
# dat_1x1[!empty_cells, c(1, 2)] %>%
#   add_column(fuse_prop = vec_spR) %>%
#   bind_rows(dat_1x1[empty_cells, c(1, 2)] %>%
#               add_column(fuse_prop = NA_real_)) %>%
#   arrange(longitude_x, desc(latitude_y)) %>%
#   # plot
#   ggplot() +
#   geom_raster(aes(x = longitude_x,
#                   y = latitude_y,
#                   fill = fuse_prop)) +
#   scale_fill_gradient(low = colour_mint,
#                       high = colour_purple,
#                       na.value = "white",
#                       name = "Coverage") +
#   geom_sf(data = world_map_sf, col = NA, fill = "grey80", size = 0.1) +
#   labs(x = "Longitude", y = "Latitude") +
#   scale_x_continuous(breaks = seq(-180, 180, 30),
#                      expand = c(0, 0)) +
#   scale_y_continuous(breaks = seq(-90, 90, 30),
#                      expand = c(0, 0))
# 
# lm(dat_prop$fuse_prop)