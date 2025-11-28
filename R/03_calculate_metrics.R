library(here)
library(sf)
library(raster)
library(funrar)
library(tidyverse)

# get functions
source(here("R", 
            "functions.R"))

# read data ---------------------------------------------------------------

# merged presence-absence matrix
dat_presabs <- read_rds(here("data",
                             "presabs_05_res.rds")) %>% 
  # select "grid" cells with >5 species (the trait space is 4D)
  mutate(srich = rowSums(select(., -c(1, 2)))) %>% 
  filter(srich > 5) %>% 
  select(-c(srich, rita_rita)) 

# IUCN categories
dat_iucn <- read_csv(here("data",
                          "megafauna_traits.csv")) %>% 
  select(species, higher_classification, IUCN) %>% 
  mutate(species = str_to_lower(species), 
         species = str_replace_all(species, " ", "_"))

# distance trait matrix
dat_dist <- read_rds(here("data", 
                          "distance_trait_matrix.rds"))

# functional space
dat_space <- read_rds(here("data", "functional_space.rds")) 
  
# provinces of the world
MEOW_sf <- read_sf(here("data", 
                        "meow", 
                        "meow_ecos.shp"))

# global FUSE --------------------------------------------------------------------

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
  filter(species %in% rownames(dat_space)) 

# arrange correctly
ge <- ge[match(rownames(dat_space), ge$species), ] %>% 
  pull(ge) 

# compute metrics 
dat_global <- get_metrics(Mat_dist = dat_dist,
                     Coords = dat_space,
                     GE = ge) %>% 
  as_tibble(rownames = "species") %>% 
  # compute distinctiveness
  left_join(as_tibble(distinctiveness_global(dat_dist, di_name = "FDi_std")), 
            by = join_by(species))

# save 
dat_global %>% 
  write_rds(here("data", 
                 "global_metrics.rds"))

# calculate local metrices ------------------------------------------------------


# list of all species present per grid 
# list of all species present per grid 
spp_per_grid <- apply(dat_presabs, 
                      1, 
                      function(x){colnames(dat_presabs)[which(x==1
                      )]}) %>% 
  # remove empty grids
  compact()

# set up empty list
list_fuse_local <- vector(mode = "list", length = length(spp_per_grid))

# loop through realms
for (i in 1:length(spp_per_grid)) {
  
  # get species pool per realm 
  spp_per_cell <- spp_per_grid[[i]] %>%
    enframe(value = "species",
            name = NULL) %>% 
    distinct(species)  
  
  # subset distance trait matrix to realm species pool
  spec_to_keep <- dat_dist %>%
    as.matrix() %>%
    rownames() %in% spp_per_cell$species
  dist_mat_realm <- as.matrix(dat_dist)[spec_to_keep, spec_to_keep]
  
  # subset trait space to realm species pool
  pcoa_realm <- na.omit(dat_space[rownames(dat_space) %in% spp_per_cell$species, ]) 
  
  
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
    filter(species %in% rownames(pcoa_realm)) 
  
  # arrange correctly
  pcoa_realm <- pcoa_realm[match(ge$species, rownames(pcoa_realm)), ]
  dist_mat_realm <- dist_mat_realm[match(ge$species, rownames(dist_mat_realm)), ]
  ge <- ge[match(rownames(pcoa_realm), ge$species), ] %>% 
    pull(ge) 
  
  
  identical(row.names(dist_mat_realm), row.names(pcoa_realm))
  
  # Compute FUSE metrics 
  list_fuse_local[[i]] <- get_metrics(Mat_dist = dist_mat_realm,
                                      Coords = pcoa_realm,
                                      GE = ge) %>% 
    as_tibble(rownames = "species") %>% 
    arrange(species) %>%
    select(species, FUSE, FUn_std, FSp_std) %>% 
    rename_with(.cols = -species, 
                .fn = ~ paste0(.x, "_local")) %>% 
    # compute distinctiveness
    left_join(as_tibble(distinctiveness_global(dist_mat_realm, di_name = "FDi_std_local")), 
              by = join_by(species))
  
  print(i)
}

# get grids without any species
vec_no_spec <- map_dbl(apply(dat_presabs, 
                             1, 
                             function(x){colnames(dat_presabs)[which(x==1
                             )]}), length) == 0

# unlist and save
list_fuse_local %>%
  bind_rows(.id = "id") %>%
  left_join(tibble(id = unique(.$id),
                   latitude = dat_presabs$latitude_y[!vec_no_spec],
                   longitude = dat_presabs$longitude_x[!vec_no_spec])) %>%
  write_rds(here("data",
                 "local_metrics_1_res.rds"),
            compress = "gz")



# local metrics at lower resolution -----------------------------------------------------------


# use 5x5 resolution
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
             .before = 0) 



# list of all species present per grid 
spp_per_grid <- apply(dat_upscld, 
                      1, 
                      function(x){colnames(dat_presabs)[which(x==1
                      )]}) %>% 
  # remove empty grids
  compact()

# set up empty list
list_fuse_local <- vector(mode = "list", length = nrow(dat_upscld))

# loop through realms
for (i in 1:length(spp_per_grid)) {
  
  # get species pool per realm 
  spp_per_cell <- spp_per_grid[[i]] %>%
    enframe(value = "species",
            name = NULL) %>% 
    distinct(species)  
  
  # subset distance trait matrix to realm species pool
  spec_to_keep <- dat_dist %>%
    as.matrix() %>%
    rownames() %in% spp_per_cell$species
  dist_mat_realm <- as.matrix(dat_dist)[spec_to_keep, spec_to_keep]
  
  # subset trait space to realm species pool
  pcoa_realm <- na.omit(dat_space[rownames(dat_space) %in% spp_per_cell$species, ]) 
  
  
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
    filter(species %in% rownames(pcoa_realm)) 
  
  # arrange correctly
  pcoa_realm <- pcoa_realm[match(ge$species, rownames(pcoa_realm)), ]
  dist_mat_realm <- dist_mat_realm[match(ge$species, rownames(dist_mat_realm)), ]
  ge <- ge[match(rownames(pcoa_realm), ge$species), ] %>% 
    pull(ge) 
  
  
  identical(row.names(dist_mat_realm), row.names(pcoa_realm))
  
  # Compute FUSE metrics 
  list_fuse_local[[i]] <- get_metrics(Mat_dist = dist_mat_realm,
                                   Coords = pcoa_realm,
                                   GE = ge) %>% 
    as_tibble(rownames = "species") %>% 
    arrange(species) %>%
    select(species, FUSE, FUn_std, FSp_std) %>% 
    rename_with(.cols = -species, 
                .fn = ~ paste0(.x, "_local")) %>% 
    # compute distinctiveness
    left_join(as_tibble(distinctiveness_global(dist_mat_realm, di_name = "FDi_std_local")), 
              by = join_by(species))
  
  print(i)
}

# get grids without any species
vec_no_spec <- map_dbl(apply(dat_upscld, 
                             1, 
                             function(x){colnames(dat_presabs)[which(x==1
                             )]}), length) == 0

# unlist and save
list_fuse_local %>%
  bind_rows(.id = "id") %>%
  left_join(tibble(id = unique(.$id),
                   latitude = dat_upscld$latitude_y[!vec_no_spec],
                   longitude = dat_upscld$longitude_x[!vec_no_spec])) %>% 
  write_rds(here("data",
                 "local_metrics_5_res.rds"),
            compress = "gz")



# metrics per realm ----------------------------------------------------------

# list of all species present per 0.5x0.5 grid 
spp_per_grid <- apply(dat_presabs, 
                      1, 
                      function(x){colnames(dat_presabs)[which(x==1
                      )]}) %>% 
  # remove empty grids
  compact()

# assign grid cells to realms 
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
list_fuse_realm <- vector(mode = "list", 
                          length = length(realm_char))

# loop through realms
for (i in 1:length(realm_char)) {
  
  # get species pool per realm 
  spp_per_realm <- spp_per_grid[matrix_realm[[i]]] %>%
    flatten_chr() %>% 
    enframe(value = "species",
            name = NULL) %>% 
    distinct(species) %>% 
    add_column(realm = realm_char[i]) 
  
  # subset distance trait matrix to realm species pool
  spec_to_keep <- dat_dist %>%
    as.matrix() %>%
    rownames() %in% spp_per_realm$species
  dist_mat_realm <- as.matrix(dat_dist)[spec_to_keep, spec_to_keep]
  
  # subset trait space to realm species pool
  pcoa_realm <- na.omit(dat_space[rownames(dat_space) %in% spp_per_realm$species, ]) 
  
  
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
    filter(species %in% rownames(pcoa_realm)) 
  
  # arrange correctly
  pcoa_realm <- pcoa_realm[match(ge$species, rownames(pcoa_realm)), ]
  dist_mat_realm <- dist_mat_realm[match(ge$species, rownames(dist_mat_realm)), ]
  ge <- ge[match(rownames(pcoa_realm), ge$species), ] %>% 
    pull(ge) 
  
  
  identical(row.names(dist_mat_realm), row.names(pcoa_realm))
  
  # Compute FUSE metrics 
  list_fuse_realm[[i]] <- get_metrics(Mat_dist = dist_mat_realm,
                                   Coords = pcoa_realm,
                                   GE = ge) %>% 
    as_tibble(rownames = "species") %>% 
    arrange(species) %>%
    select(species, FUSE, FUn_std, FSp_std) %>% 
    rename_with(.cols = -species, 
                .fn = ~ paste0(.x, "_local")) %>% 
    add_column(realm = realm_char[i]) %>% 
    # compute distinctiveness
    left_join(as_tibble(distinctiveness_global(dist_mat_realm, di_name = "FDi_std_local")), 
              by = join_by(species))
  
}

# combine in one dataframe
dat_realm <- bind_rows(list_fuse_realm) %>% 
  rename_with(.cols = contains("std"), 
              .fn = ~ str_remove_all(.x, "std_"))

# save 
dat_realm %>% 
  write_rds(here("data", 
                 "realm_metrics.rds"))

