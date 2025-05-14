library(here)
library(sf)
library(raster)
library(ade4)
library(patchwork)
library(tidyverse)


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



# distance-trait matrix
load(here("data", 
          "distance_trait_matrix.RData"))


# IUCN categories
dat_iucn <- read_csv(here("data",
                          "megafauna_traits.csv")) %>%
  select(species, higher_classification, IUCN)

# World map
world_map_sf <- st_read(here("data", 
                             "worldmap",
                             "ne_10m_land.shp"))

# provinces of the world
MEOW_sf <- read_sf(here("data", 
                        "meow", 
                        "meow_ecos.shp"))

# uniqueness and specialisation comparison
dat_realm <- read_rds(here("data",
                           "global_regional_trends.rds"))

# upscaling ---------------------------------------------------------------

# use 5x5 resolution
dat_rast_agg <- dat_presabs %>%
  rasterFromXYZ() %>%
  aggregate(fact = 10, #5x5 
            fun = max)

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


# calculate metrices ------------------------------------------------------


# extract functional space axes
pcoa <- dudi.pco(quasieuclid(distance_trait_matrix),
                 scannf = FALSE,
                 nf = 4)



# list of all species present per grid 
spp_per_grid <- apply(dat_upscld, 
                      1, 
                      function(x){colnames(dat_upscld)[which(x==1
                      )]}) %>% 
  # remove empty grids
  compact()


# FUSE

# prepare iucn categories
ge <- dat_iucn %>%
  mutate(ge = case_when(
    IUCN == "LC" ~ 0, 
    IUCN == "NT" ~ 1, 
    IUCN == "VU" ~ 2,
    IUCN == "EN" ~ 3, 
    IUCN == "CR" ~ 4, 
    .default = NA_integer_
  )) %>% 
  filter(species %in% rownames(pcoa$li)) 

# arrange correctly
ge <- ge[match(rownames(pcoa$li), ge$species), ] %>% 
  pull(ge) 

# Compute FUSE metrics 
dat_fuse <- get_FUSE(Mat_dist = distance_trait_matrix,
                     Coords = pcoa$li,
                     GE = ge) %>% 
  as_tibble(rownames = "species") 



# compute mean FUSE per grid cell
dat_mean_fuse <- spp_per_grid %>% 
  map(~ dat_fuse %>% 
        select(species, FUSE,
               FUn_std, FSp_std) %>% 
        filter(species %in% .x) %>% 
        summarise(mean_fuse = mean(FUSE, na.rm = TRUE), 
                  mean_uniq = mean(FUn_std, na.rm = TRUE), 
                  mean_spec = mean(FSp_std, na.rm = TRUE)), 
      .progress = TRUE)

# extract

# fuse
vec_fuse <- map_dbl(dat_mean_fuse, 
                    ~ pull(.x, "mean_fuse"), 
                    .progress = TRUE)

# uniqueness
vec_uniq <- map_dbl(dat_mean_fuse,
                    ~ pull(.x, "mean_uniq"), 
                    .progress = TRUE)

# specialisation
vec_spec <- map_dbl(dat_mean_fuse,
                    ~ pull(.x, "mean_spec"), 
                    .progress = TRUE)


# combine ----------------------------------------------------------------

# combine in dataframe with coordinates
dat_metrics_global <- tibble(fuse = vec_fuse, 
                             uniq = vec_uniq, 
                             special = vec_spec) %>% 
  add_column(longitude_x = dat_upscld %>%
               filter(rowSums(select(., -contains("itude"))) > 0) %>%
               pull(longitude_x), 
             latitude_y = dat_upscld %>%
               filter(rowSums(select(., -contains("itude"))) > 0) %>%
               pull(latitude_y), .before = 1) 

# read in local mean metrics
dat_metric_local <- read_rds(here("data",
                                  "functional_metrics_local.rds")) %>% 
  select(contains("_1")) %>% 
  rename_with(~str_remove(.x, "_1"))

# merge together
plot_metrics <- bind_cols(dat_metrics_global, dat_metric_local) %>%
  # calculate contrast
  mutate(FUSE = fuse - FUSE_local, 
         FUn = uniq - FUn_std_local, 
         FSp = special - FSp_std_local,
         .keep = "unused") %>% 
  pivot_longer(-contains("itude"), 
               names_to = "metric") %>% 
  ggplot() +
  geom_raster(aes(x = longitude_x,
                  y = latitude_y, 
                  fill = value)) +
  geom_sf(data = world_map_sf, col = NA, fill = "white", size = 0.1) +
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10,"RdBu")), 
                       name = "Global - Local") +
  labs(x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(-180, 180), ylim= c(-89, 85), expand=FALSE) +
  scale_x_continuous(breaks = seq(-180, 180, 60), 
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-90, 90, 30), 
                     expand = c(0, 0)) +
  facet_wrap(~ metric, 
             ncol = 1) +
  theme(legend.position = "bottom")


# save plot
ggsave(plot_metrics,
       filename = here("figures",
                       "main",
                       "global_minus_local.pdf"),
       width = 140, height = 200,
       units = "mm",
       bg = "white")
