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

# load metrics calculated in 1
dat_metrics <- read_rds(here("data",
                             "functional_metrics.rds"))  


# load local metrics calculated in 1
dat_metrics_local <- read_rds(here("data",
                                   "functional_metrics_local.rds"))

# world map
world_map_sf <- st_read(here("data", 
                             "worldmap",
                             "ne_10m_land.shp"))

# global distinctiveness
dat_dist_glob <- read_rds(here("data", "global_distinctiveness.rds"))

# local distinctiveness
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

# dat_metrics <- read_rds(here("data",
#                              "functional_metrics_global_FDi.rds"))

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

# set up function
get_overlap <- function(metric){
  
  # upscale
  dat_rast_hot <- dat_hot %>%
    select(longitude_x, latitude_y, {{ metric }}) %>% 
    rasterFromXYZ() %>%
    aggregate(fact = 10, #5x5
              fun = sum)
  
  # extract
  dat_rast_hot %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(met_sel = {{ metric }} >= 1) %>%
    add_column(longitude_x = coordinates(dat_rast_hot)[, 1],
               latitude_y = coordinates(dat_rast_hot)[, 2],
               .before = 0) %>% 
    inner_join(dat_hot_loc %>% 
                 select(longitude_x, latitude_y, met_sel_loc = FSp_hot)) %>% 
    filter((met_sel + met_sel_loc) == 2) %>% 
    {{(nrow(.) / nrow(dat_hot_loc %>% 
                    filter({{ metric }} == TRUE)))}} %>% 
    {{.*100}} %>% 
    round() %>% 
    paste0(., "%")
}


# visualise ---------------------------------------------------------------

# create plot
plot_hot <- dat_hot %>%
  pivot_longer(cols = c(contains("hot")), 
               names_to = "hotspot") %>% 
  filter(value == TRUE) %>% 
  mutate(hotspot = str_remove_all(hotspot, "_hot")) %>% 
  ggplot() +
  geom_sf(data = world_map_sf, col = NA, fill = "grey90", size = 0.1) +
  geom_tile(aes(x = longitude_x,
                y = latitude_y),
            colour = "grey95",
            fill = "grey10",
            data = dat_hot_loc %>%
              pivot_longer(cols = c(contains("hot")),
                           names_to = "hotspot") %>% 
              filter(value == TRUE) %>% 
              mutate(hotspot = str_remove_all(hotspot, "_hot"))) +
  geom_tile(aes(x = longitude_x,
                y = latitude_y,
                fill = hotspot), 
            alpha = 0.8) +
  geom_text(aes(longitude_x, 
                latitude_y, 
                label = perc_ov), 
            size = 10/.pt,
            data = tibble(hotspot = c("FDi", "FSp", "FUn"), 
                          perc_ov = c(get_overlap(FDi_hot),
                                      get_overlap(FSp_hot),
                                      get_overlap(FUn_hot)), 
                          longitude_x = -130, latitude_y = -30)) +
  scale_fill_manual(values = c(colour_purple, colour_mint, 
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
       filename = here("figures", "main", "hotspot_scaling.pdf"), 
       width = 183, height = 300,
       device = "pdf", 
       units = "mm", bg = "white")



# fuse hotspots -----------------------------------------------------------

# calculate quantile threshold
dat_fuse <- dat_metrics %>% 
  mutate(fuse_top_hot = fuse_top > quantile(fuse_top,
                                  probs =  0.975), 
         fuse_hot = fuse > quantile(fuse,
                                    probs =  0.975), 
         overlap_hot = if_else(fuse_top_hot + fuse_hot == 2, 
                               TRUE, FALSE)) %>% 
  pivot_longer(cols = c(contains("hot")), 
               names_to = "hotspot") %>% 
  mutate(hotspot = as.factor(hotspot))

plot_fuse_hot <- dat_fuse %>%
  filter(value == TRUE) %>% 
  ggplot() +
  geom_tile(aes(x = longitude_x,
                  y = latitude_y, 
                  fill = hotspot)) +
  geom_sf(data = world_map_sf, col = NA, fill = "grey", size = 0.1) +
  scale_fill_manual(values = c("#c65831", "#169199"),
                    name = "Hotspots",
                    labels = c("FUSE", "TOP 25% FUSE SR")) +
  labs(x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(-180, 180),  ylim = c(-89, 85), expand = FALSE) +
  scale_x_continuous(breaks = seq(-180, 180, 30), 
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-90, 90, 30), 
                     expand = c(0, 0)) +
  theme(legend.position = "top", 
        legend.key.size = unit(2, "mm"))


# venn diagram  -----------------------------------------------------------

plot_venn <- dat_hot %>%
  filter(hotspot != "overlap_hot") %>% 
  filter(value == TRUE) %>% 
  select(longitude_x, latitude_y, hotspot) %>% 
  full_join(dat_fuse %>% 
              filter(hotspot != "overlap_hot") %>% 
              filter(value == TRUE) %>% 
              select(longitude_x, latitude_y, hotspot)) %>% 
  mutate(across(c(longitude_x, latitude_y), abs)) %>% 
  ggplot(aes(longitude_x, latitude_y, 
             colour = hotspot, 
             fill = hotspot)) +
  stat_ellipse(geom = "polygon", 
               alpha = 0.5) +
  labs(x = "Absolute Longitude", 
       y = "Absolute Latitude", 
       colour = "Hotspot", 
       fill = "Hotspot") +
  scale_color_manual(values = c(colour_mint, "#c65831", "#169199",
                                colour_yellow, colour_purple),
                     labels = c("FRic", "FUSE", "TOP 25% FUSE SR",  
                                "SR", "FUn")) +
  scale_fill_manual(values = c(colour_mint, "#c65831", "#169199",
                               colour_yellow, colour_purple),
                    labels = c("FRic", "FUSE", "TOP 25% FUSE SR",  
                               "SR", "FUn")) +
  theme(legend.position = c(0.9, 0.6), 
        legend.key.size = unit(5, "mm"))

# combine figures ---------------------------------------------------------

# patch together
plot_final <- plot_hot/ plot_fuse_hot/ plot_venn

# save
ggsave(plot_final, 
       filename = here("figures",
                       "main",
                       "3_hotspots.pdf"),
       width = 183, height = 100*3,
       units = "mm",
       bg = "white")
