library(here)
library(sf)
library(raster)
library(ade4)
library(funrar)
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

# World map
world_map_sf <- st_read(here("data", 
                             "worldmap",
                             "ne_10m_land.shp"))

# provinces of the world
MEOW_sf <- read_sf(here("data", 
                        "meow", 
                        "meow_ecos.shp"))


# global distinctiveness -----------------------------------------------

# on a global level
dat_glob <- distinctiveness_global(distance_trait_matrix)

# save data
dat_glob %>% 
  write_rds(here("data", 
                 "global_distinctiveness.rds"))

# provincial distinctiveness ----------------------------------------------


# list of all species present per 0.5*0.5 grid 
spp_per_grid <- apply(dat_presabs, 
                      1, 
                      function(x){colnames(dat_presabs)[which(x==1
                      )]}) 

# per province matrix
matrix_realm <- MEOW_sf %>%
  group_by(REALM) %>% 
  summarize(geometry = st_union(geometry)) %>% 
  st_intersects(st_as_sf(dat_presabs, coords = c("longitude_x", 
                                                 "latitude_y"), 
                         crs = "WGS84"))
# get REALM names 
realm_char <- MEOW_sf %>%
  group_by(REALM) %>% 
  summarize(geometry = st_union(geometry)) %>% 
  pull(REALM)


# set up empty list
list_dist_realm <- vector(mode = "list",
                         length = length(realm_char))

# loop through realms
for (i in 1:length(realm_char)) {
  
  # species present per grid 
  realm_spp <- spp_per_grid[matrix_realm[[i]]] %>%
    flatten_chr() %>% 
    enframe(value = "species",
            name = NULL) %>% 
    distinct(species) %>% 
    pull(species)
  
  
  # subset distance trait matrix to realm species pool
  spec_to_keep <- distance_trait_matrix %>%
    as.matrix() %>%
    rownames() %in% realm_spp
  dist_mat_realm <- as.matrix(distance_trait_matrix)[spec_to_keep, spec_to_keep]
  
  
  # compute dist metric
  list_dist_realm[[i]] <- distinctiveness_global(dist_mat_realm) %>% 
    as_tibble() %>% 
    add_column(realm = realm_char[[i]])

}

# combine in one dataframe
dat_province <- bind_rows(list_dist_realm) 
  
# save data
dat_province %>% 
  write_rds(here("data", 
                 "provincial_distinctiveness.rds"))

# local distinctiveness ----------------------------------------------


# upscale to 5x5 resolution
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


# list of all species present per grid 
spp_per_grid_loc <- apply(dat_upscld, 
                      1, 
                      function(x){colnames(dat_presabs)[which(x==1
                      )]}) %>% 
  # remove empty grids
  compact()

# set up empty list
list_dist_local <- vector(mode = "list", length = length(spp_per_grid_loc))

# loop through cells
for (i in 1:length(spp_per_grid_loc)) {
  
  # species present per grid 
  spp_per_cell <- spp_per_grid_loc[[i]] %>%
    enframe(value = "species",
            name = NULL) %>% 
    distinct(species) %>% 
    pull(species)
  
  
  # subset distance trait matrix to realm species pool
  spec_to_keep <- distance_trait_matrix %>%
    as.matrix() %>%
    rownames() %in% spp_per_cell
  dist_mat_local <- as.matrix(distance_trait_matrix)[spec_to_keep, spec_to_keep]
  
  
  # compute dist metric
  list_dist_local[[i]] <- distinctiveness_global(dist_mat_local) %>% 
    as_tibble() %>% 
    add_column(latitude = dat_upscld[[i, 2]], 
               longitude = dat_upscld[[i, 1]])
  
}

# combine in one dataframe
dat_local <- bind_rows(list_dist_local, 
                       .id = "id") 

# save data
dat_local %>% 
  write_rds(here("data", 
                 "local_distinctiveness.rds"))
            