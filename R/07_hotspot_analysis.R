library(here)
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


# provinces of the world
MEOW_sf <- read_sf(here("data", 
                        "meow", 
                        "meow_ecos.shp"))


# world map
world_map_sf <- st_read(here("data", 
                             "worldmap",
                             "ne_10m_land.shp"))

# merged presence-absence matrix
dat_presabs <- read_rds(here("data",
                             "presabs_05_res.rds")) %>% 
  # select "grid" cells with >5 species (the trait space is 5D)
  mutate(srich = rowSums(select(., -c(1, 2)))) %>% 
  filter(srich > 5) %>% 
  select(-srich) 


# global distribution of metrices -------------------------------------------------------------

# aggregate
dat_glob_agg <- dat_presabs %>%
  pivot_longer(cols = -c(longitude_x, latitude_y), 
               names_to = "species") %>% 
  filter(value == 1) %>% 
  left_join(dat_global %>% 
              select(-FUSE)) %>% 
  mutate(longitude_x = floor(longitude_x / 5) * 5 + 2.5,
         latitude_y = floor(latitude_y / 5) * 5 + 2.5) %>%
  group_by(longitude_x, latitude_y) %>%
  summarise(across(c(FSp, FUn, FDi), 
                   ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") %>% 
  add_column(reso = "Global") 

# rearrange data
dat_comp <- dat_glob_agg %>% 
  bind_rows(dat_local %>% 
              select(longitude_x = longitude, latitude_y = latitude,
                     FDi = FDi_local,  FSp = FSp_local, FUn = FUn_local) %>% 
              group_by(longitude_x, latitude_y) %>% 
              summarise(across(c(FDi, FSp, FUn), mean), 
                        .groups = "drop") %>% 
              add_column(reso = "Local")) %>% 
  mutate(across(c(FDi, FSp, FUn), \(x) (x-min(x))/(max(x)-min(x)))) %>% 
  pivot_longer(cols = -c(longitude_x, latitude_y, reso), 
               names_to = "metric")

# visualise
plot_maps <- dat_comp %>%
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
       filename = here("figures", "main", "fig_4.pdf"), 
       width = 183, height = 140,
       units = "mm", bg = "white")


# calculate hot-spots ------------------------------------------------------

# set up function to alternate hot-spot cutoff

get_hotspots <- function(cutoff) {
  
  dat_hot <- dat_local %>%
    rename(longitude_x = longitude, latitude_y = latitude) %>% 
    left_join(dat_global %>% 
                select(species, 
                       FDi_global = FDi, FSp_global = FSp, FUn_global = FUn)) %>% 
    group_by(longitude_x, latitude_y) %>% 
    summarise(across(c(FUn_local:FUn_global), mean), 
              .groups = "drop") %>% 
    pivot_longer(cols = -c(longitude_x, latitude_y), 
                 names_to = "metric") %>% 
    separate_wider_delim(cols = metric, delim = "_", names = c("metric", "reso")) %>% 
    group_by(reso, metric) %>% 
    mutate(is_hot = value > quantile(value, probs =  cutoff)) %>% 
    ungroup() %>% 
    select(-value)
  
  
  # count overlap
  dat_prop <- dat_hot %>%
    filter(is_hot) %>% 
    count(metric, longitude_x, latitude_y) %>% 
    filter(n == 2) %>% 
    mutate(metric = factor(metric, levels = c("FUn", "FSp", "FDi"))) %>% 
    count(metric, name = "shared", 
          .drop = FALSE) %>% 
    left_join(dat_hot %>% 
                filter(is_hot) %>% 
                count(metric)) %>% 
    mutate(prop_ov = (shared / n) * 100, 
           prop_ov = round(prop_ov, 0), 
           prop_ov = paste0(prop_ov, "%")) %>% 
    add_column(cutoff = round((1-cutoff)*100, 1))
  
  list(dat_hot, dat_prop)
  
}

# use 2.5% as cutoff
dat_2.5 <- get_hotspots(0.975)
  
# use 5% as cutoff
dat_5 <- get_hotspots(0.95)

# use 10% as cutoff
dat_10 <- get_hotspots(0.9)

# use 10% as cutoff
dat_20 <- get_hotspots(0.8)


# merge and save
dat_2.5[[2]] %>% 
  bind_rows(dat_5[[2]]) %>% 
  bind_rows(dat_10[[2]]) %>% 
  bind_rows(dat_20[[2]]) %>% 
  write_csv(here("data", "varying_hotspot_thresholds.csv"))


# create plot
plot_hot <- dat_2.5[[1]] %>%
  filter(is_hot) %>% 
  ggplot(aes(x = longitude_x,
             y = latitude_y)) +
  geom_sf(data = world_map_sf, col = NA, fill = "grey90", size = 0.1, 
          inherit.aes = FALSE) +
  geom_tile(aes(fill = metric),
            data = . %>% filter(reso == "local"),
            colour = "grey20") +
  geom_point(aes(colour = metric),
             data = . %>% filter(reso == "global"),
             shape = 21,
             size = 1.5,
             stroke = 0.75,
             fill = "white") +
  geom_text(aes(label = prop_ov),
            size = 10/.pt,
            data = dat_2.5[[2]] %>% 
              mutate(longitude_x = -130, 
                     latitude_y = -30)) +
  geom_text(aes(label = perc_ov),
            data = tibble(longitude_x = -110,
                          latitude_y = -30,
                          perc_ov = "Overlap",
                          metric = c("FDi", "FSp", "FUn")),
            colour = "grey20",
            size = 10/.pt) +
  geom_segment(aes(xend = long_end, yend = lat_end), 
               arrow = arrow(length = unit(0.1, "inches")),
               colour = "grey40",
               data = tibble(longitude_x = c(-148, -40), long_end = c(-127, -63),
                             latitude_y = c(34, 40), lat_end = c(48, 55),
                             metric = "FDi")) +
  geom_text(aes(label = perc_ov), 
            data = tibble(longitude_x = c(-148, -40), 
                          latitude_y = c(25, 31), 
                          perc_ov = c("Local hotspot", "Global hotspot"), 
                          metric = "FDi"), 
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
  facet_wrap(~metric, 
             ncol = 1)


# save
ggsave(plot_hot,
       filename = here("figures", "main", "fig_5.pdf"), 
       width = 183, height = 280,
       units = "mm", bg = "white")



