library(here)
library(sf)
library(raster)
library(ade4)
library(patchwork)
library(tidyverse)
library(future.apply)

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

# extract functional space axes
pcoa <- dudi.pco(quasieuclid(distance_trait_matrix),
                 scannf = FALSE,
                 nf = 4) 

# World map
world_map_sf <- st_read(here("data", 
                             "worldmap",
                             "ne_10m_land.shp"))

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



# list of all species present per grid 
spp_per_grid <- apply(dat_upscld, 
                      1, 
                      function(x){colnames(dat_presabs)[which(x==1
                      )]}) %>% 
  # remove empty grids
  compact()

# calculate functional richness
list_FRic <- spp_per_grid %>%
  map( ~ get_FV_Sp(
    ax = c(1:4),
    pcoa = pcoa,
    Selected_sp = .x
  ), .progress = TRUE) %>%
  map_dbl( ~ pluck(.x, "data") %>%
             pull("FRic"))



# calculate null distribution sensu Villéger, Grenouillet and Bros (2013)

# define a function that performs the operations and handles errors
safe_operation <- function() {
  tryCatch({
    spp_per_grid %>%
      map_dbl(length) %>%
      map(~ sample(colnames(dat_presabs)[-c(1, 2)], .x)) %>%
      map(~ get_FV_Sp(ax = c(1:4), pcoa = pcoa, Selected_sp = .x)) %>%
      map_dbl(~ pluck(.x, "data") %>%
                pull("FRic"))
  }, error = function(e) {
    NA_integer_
  })
}

# use 15 cores
plan(multisession, workers = 15)

# use future_replicate with the safe_operation function
list_null <- future_replicate(n = 1000, safe_operation(), simplify = FALSE)

# set up vector
null_larger <- vector("double", 
                    length = length(list_FRic))

null_smaller <- vector("double", 
                      length = length(list_FRic))

# iterate over cells
for (i in seq_along(list_FRic)) {
  
  # extract empirical results
  emp_res <- list_FRic[[i]]
  
  # null distribution
  null_dist <- map_dbl(list_null, i, .default = NA)

  # how many null values are larger than empirical
  null_larger[[i]] <- sum(null_dist > emp_res, na.rm = TRUE)
  
  # same for smaller
  null_smaller[[i]] <- sum(null_dist < emp_res, na.rm = TRUE)
  

  }


# get coordinates
dat_final <- dat_upscld %>% 
  rowwise() %>% 
  filter(sum(c_across(-c(longitude_x, latitude_y))) > 0) %>% 
  ungroup() %>% 
  select(c(longitude_x, latitude_y)) %>% 
  add_column(smaller = null_smaller, 
             larger = null_larger)

# save data
dat_final %>% 
  write_rds(here("data", 
                 "null_comparison.rds"))

plot_null <- dat_final %>% 
  ggplot() +
  geom_raster(aes(x = longitude_x,
                  y = latitude_y, 
                  fill = smaller)) +
  geom_sf(data = world_map_sf, col = NA, fill = "white", size = 0.1) +
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10,"RdBu")), 
                       breaks = seq(0, 1000, by = 250), 
                       limits = c(0, 1000),
                       name = NULL, 
                       labels = c("Clustered", "", "", "", "Overdispersed")) +
  labs(x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(-180, 180), ylim= c(-89, 85), expand=FALSE) +
  scale_x_continuous(breaks = seq(-180, 180, 30), 
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-90, 90, 30), 
                     expand = c(0, 0)) 


# save plot
ggsave(plot_null, 
       filename = here("figures",
                       "main",
                       "overdispersion.pdf"),
       width = 183, height = 100,
       units = "mm",
       bg = "white")
