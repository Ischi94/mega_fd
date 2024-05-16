library(here)
library(tidyverse)
library(sf)
library(ade4)
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
dat_pcoa <- dudi.pco(quasieuclid(distance_trait_matrix),
                     scannf = FALSE,
                     nf = 4) %>% 
  pluck("li") %>% 
  # prepare pcoa for joining
  as_tibble(rownames = "species")


# provinces of the world
MEOW_sf <- read_sf(here("data", 
                        "meow", 
                        "meow_ecos.shp"))

# uniqueness and specialisation comparison
dat_realm <- read_rds(here("data",
                           "global_regional_trends.rds"))

# species pool per realm -------------------------------------------------------

# list of all species present per grid 
spp_per_grid <- apply(dat_presabs, 
                      1, 
                      function(x){colnames(dat_presabs)[which(x==1
                      )]}) 

# group by realm
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
list_spp_realm <- vector(mode = "list",
                         length = length(realm_char))

# loop through realms
for (i in 1:length(realm_char)) {
  
  list_spp_realm[[i]] <- spp_per_grid[matrix_realm[[i]]] %>%
    flatten_chr() %>% 
    enframe(value = "species",
            name = NULL) %>% 
    distinct(species) %>% 
    add_column(realm = realm_char[i]) 
}


# assign pcoa to realms ---------------------------------------------------

# get all four axis per realm
dat_spp_realm <- list_spp_realm %>% 
  map(~ left_join(.x, dat_pcoa)) %>% 
  bind_rows()

# calculate overall convex hull for axis 1 and 2
dat_hull_1_ov <- dat_spp_realm %>% 
  slice(chull(A1, A2)) %>% 
  select(-realm)

# calculate convex hulls for axis 1 and 2
dat_hull_1 <- dat_spp_realm %>% 
  group_by(realm) %>% 
  slice(chull(A1, A2)) %>% 
  ungroup()

# calculate overall convex hull for axis 1 and 3
dat_hull_2_ov <- dat_spp_realm %>% 
  slice(chull(A1, A3)) %>% 
  select(-realm)

# calculate convex hulls for axis 1 and 3
dat_hull_2 <- dat_spp_realm %>% 
  group_by(realm) %>% 
  slice(chull(A1, A3)) %>% 
  ungroup()

# calculate overall convex hull for axis 2 and 3
dat_hull_3_ov <- dat_spp_realm %>% 
  slice(chull(A2, A3)) %>% 
  select(-realm)

# calculate convex hulls for axis 2 and 3
dat_hull_3 <- dat_spp_realm %>% 
  group_by(realm) %>% 
  slice(chull(A2, A3)) %>% 
  ungroup()


# specialisation ----------------------------------------------------------

dat_spp_realm %>%
  # filter(realm %in% (dat_spp_realm %>% 
  #          filter(species == "Chelonia mydas") %>% 
  #          distinct(realm) %>% 
  #          pull(realm)))
  group_by(realm) %>% 
  # center of traitspace
  summarise(A1 = mean(A1),
            A2 = mean(A2),
            A3 = mean(A3)) %>% 
  full_join(dat_spp_realm %>% 
              summarise(A1 = mean(A1), 
                        A2 = mean(A2), 
                        A3 = mean(A3)) %>% 
              add_column(realm = "global")) %>% 
  ggplot(aes(A2, A3)) +
  geom_text(aes(label = realm)) +
  geom_point(data = dat_spp_realm %>% 
               filter(species == "Chelonia mydas")) 

  
dat_spp_realm