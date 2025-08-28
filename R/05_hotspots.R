library(here)
library(sf)
library(raster)
library(tidyverse)
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

# load metrics, calculated in /supplement/2_1_biodiversity_metrics.R
dat_metrics <- read_rds(here("data",
                             "functional_metrics.rds"))  


# load local metrics, calculated in /supplement/2_2_local_biodiversity_metrics.R
dat_metrics_local <- read_rds(here("data",
                                   "functional_metrics_local.rds"))

# world map
world_map_sf <- st_read(here("data", 
                             "worldmap",
                             "ne_10m_land.shp"))

# global distinctiveness, calculated in /supplement/2_1_4_distinctiveness.R
dat_dist_glob <- read_rds(here("data", "global_distinctiveness.rds"))

# local distinctiveness, calculated in /supplement/2_1_4_distinctiveness.R
dat_dist_loc <- read_rds(here("data", "local_distinctiveness.rds"))



# merge data --------------------------------------------------------------

# list of all species present per grid 
spp_per_grid <- apply(dat_presabs, 
                      1, 
                      function(x){colnames(dat_presabs)[which(x==1
                      )]}) %>% 
  # remove empty grids
  compact()

# join by distinctiveness and calculate mean
dat_dist_glob_sum <- spp_per_grid %>% 
  map_dbl(~ .x %>% 
        enframe(name = NULL, 
                value = "species") %>% 
        left_join(dat_dist_glob, 
                  by = join_by(species)) %>% 
        summarise(FDi = mean(global_di)) %>% 
        pull(FDi), 
        .progress = TRUE)

# assign
dat_metrics <- dat_metrics %>% 
  add_column(FDi = dat_dist_glob_sum) 

# save
dat_metrics %>%
  write_rds(here("data",
                 "functional_metrics_global_FDi.rds"))

# dat_metrics <- read_rds(here("data", "functional_metrics_global_FDi.rds"))

# upscale to 5x5 resolution
dat_rast_agg <- dat_presabs %>%
  rasterFromXYZ() %>%
  aggregate(fact = 10, #5x5 
            fun = max)

# reformat
dat_upscld <- dat_rast_agg %>% 
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


# list of all species present per grid 
spp_per_grid_loc <- apply(dat_upscld, 
                          1, 
                          function(x){colnames(dat_presabs)[which(x==1
                          )]}) %>% 
  # remove empty grids
  compact()


# join by distinctiveness and calculate mean
dat_dist_loc_sum <- spp_per_grid_loc %>% 
  map_dbl(~ .x %>% 
            enframe(name = NULL, 
                    value = "species") %>% 
            left_join(dat_dist_loc, 
                      by = join_by(species)) %>% 
            summarise(FDi = mean(global_di)) %>% 
            pull(FDi), 
          .progress = TRUE) 
  
dat_metrics_local <- dat_metrics_local %>% 
  add_column(FDi = dat_dist_loc_sum)

# save data
dat_metrics_local %>% 
  write_rds(here("data",
                 "functional_metrics_local_FDi.rds"))


# dat_metrics_local <- read_rds(here("data", "functional_metrics_local_FDi.rds"))



# visualise maps ------------------------------------------------------------------------------

# rearrange data
plot_maps <- dat_metrics %>% 
  select(FDi, FSp = special, FUn = uniq, longitude_x, latitude_y) %>% 
  add_column(reso = "Global") %>% 
  pivot_longer(cols = c(FDi, FSp, FUn), 
               names_to = "metric") %>% 
  bind_rows(dat_metrics_local %>% 
              select(FDi, FSp = FSp_std_local_1, FUn = FUn_std_local_1, longitude_x, latitude_y) %>% 
              add_column(reso = "Local") %>% 
              pivot_longer(cols = c(FDi, FSp, FUn), 
                           names_to = "metric")) %>% 
  # summarise(mean(value, na.rm = T), median(value, na.rm = T))
  ggplot() +
  geom_raster(aes(x = longitude_x, 
                  y = latitude_y, 
                  fill = value)) +
  geom_sf(data = world_map_sf,
          col = NA,
          fill = "grey90",
          size = 0.1) +
  scale_fill_gradient2(
    low = colour_mint,
    mid = colour_purple,
    high = colour_yellow,
    midpoint = 0.4, 
    name = NULL) +
  scale_y_continuous(breaks = seq(-90, 90, 30),
                     expand = c(0, 0)) +
  labs(x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(-180, 180),  ylim = c(-89, 85), expand = FALSE) +
  facet_grid(rows = vars(metric), 
             cols = vars(reso))

# save
ggsave(plot_maps,
       filename = here("figures", "main", "fig_4.png"), 
       width = 183, height = 140,
       device = "png", 
       units = "mm", bg = "white")


# calculate hot-spots ------------------------------------------------------

# use 2.5% as cutoff
dat_hot <- dat_metrics %>%
  mutate(FSp_hot = special > quantile(special,
                                    probs =  0.975),
         FUn_hot = uniq > quantile(uniq,
                                    probs =  0.975),
         FDi_hot = FDi > quantile(FDi,
                                    probs =  0.975)) 

# same for local
# use 2.5% as cutoff
dat_hot_loc <- dat_metrics_local %>%
  mutate(FSp_hot = FSp_std_local_1 > quantile(FSp_std_local_1,
                                              probs =  0.975, 
                                              na.rm = TRUE),
         FUn_hot = FUn_std_local_1 > quantile(FUn_std_local_1,
                                               probs =  0.975, 
                                               na.rm = TRUE),
         FDi_hot = FDi > quantile(FDi,
                                  probs =  0.975, 
                                  na.rm = TRUE)) 


# count overlap -----------------------------------------------------------

# upscale
dat_rast_hot <- dat_hot %>%
  select(longitude_x, latitude_y, contains("hot"))  %>% 
  rasterFromXYZ() %>%
  aggregate(fact = 10, #5x5
            fun = sum)

# set up function
get_overlap <- function(metric){
  
  # extract
  dat_rast_hot %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(met_sel = {{ metric }} >= 1) %>%
    add_column(longitude_x = coordinates(dat_rast_hot)[, 1],
               latitude_y = coordinates(dat_rast_hot)[, 2],
               .before = 0) %>% 
    inner_join(dat_hot_loc %>% 
                 select(longitude_x, latitude_y, met_sel_loc = {{ metric }})) %>% 
    filter((met_sel + met_sel_loc) == 2) %>% 
    {{(nrow(.) / nrow(dat_hot_loc %>% 
                    filter({{ metric }} == TRUE) %>% 
                      bind_rows(dat_rast_hot %>%
                                  as.data.frame() %>%
                                  as_tibble()  %>% 
                                  filter({{ metric }} >= 1))))}} %>% 
    {{.*100}} %>% 
    round() %>% 
    paste0(., "%")
}

# visualise ---------------------------------------------------------------

# create plot
plot_hot <- dat_rast_hot %>%
  as.data.frame() %>%
  as_tibble() %>% 
  add_column(longitude_x = coordinates(dat_rast_hot)[, 1],
             latitude_y = coordinates(dat_rast_hot)[, 2],
             .before = 0) %>% 
  pivot_longer(cols = c(contains("hot")), 
               names_to = "hotspot") %>% 
  filter(value >= 1) %>% 
  mutate(hotspot = str_remove_all(hotspot, "_hot")) %>% 
  ggplot(aes(x = longitude_x,
             y = latitude_y)) +
  geom_sf(data = world_map_sf, col = NA, fill = "grey90", size = 0.1, 
          inherit.aes = FALSE) +
  geom_tile(aes(fill = hotspot),
            colour = "grey20",
            data = dat_hot_loc %>%
              pivot_longer(cols = c(contains("hot")),
                           names_to = "hotspot") %>% 
              filter(value == TRUE) %>% 
              mutate(hotspot = str_remove_all(hotspot, "_hot"))) +
  geom_point(aes(colour = hotspot), 
             shape = 21, 
             size = 1.5,
             stroke = 0.75,
             fill = "white") +
  geom_text(aes(label = perc_ov),
            size = 10/.pt,
            data = tibble(hotspot = c("FDi", "FSp", "FUn"),
                          perc_ov = c(get_overlap(FDi_hot),
                                      get_overlap(FSp_hot),
                                      get_overlap(FUn_hot)),
                          longitude_x = -130, latitude_y = -30)) +
  geom_text(aes(label = perc_ov), 
             data = tibble(longitude_x = c(-110, -110, -107), 
                           latitude_y = -30, 
                           perc_ov = "Overlap", 
                           hotspot = c("FDi", "FSp", "FUn")), 
            colour = "grey20",
             size = 10/.pt) + 
  geom_segment(aes(xend = long_end, yend = lat_end), 
               arrow = arrow(length = unit(0.1, "inches")),
               colour = "grey40",
               data = tibble(longitude_x = c(-140, -40), long_end = c(-125, -60),
                             latitude_y = c(30, 40), lat_end = c(60, 57),
                             hotspot = "FDi")) +
  geom_text(aes(label = perc_ov), 
            data = tibble(longitude_x = c(-148, -40), 
                          latitude_y = c(20, 30), 
                          perc_ov = c("Local hotspot", "Global hotspot"), 
                          hotspot = "FDi"), 
            colour = "grey20",
            size = 10/.pt) + 
  scale_fill_manual(values = c(colour_purple, colour_mint, 
                               colour_yellow), 
                    name = "Hotspots") +
  scale_colour_manual(values = c(colour_purple, colour_mint, 
                               colour_yellow), 
                    name = "Hotspots") +
  labs(x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(-180, 180),  ylim = c(-89, 85), expand = FALSE) +
  scale_x_continuous(breaks = seq(-180, 180, 60), 
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-90, 90, 30), 
                     expand = c(0, 0)) +
  theme(legend.position = "none",
        legend.key.size = unit(2, "mm")) +
  facet_wrap(~hotspot, 
             ncol = 1)


# save
ggsave(plot_hot,
       filename = here("figures", "main", "fig_5.png"), 
       width = 183, height = 280,
       device = "png", 
       units = "mm", bg = "white")



