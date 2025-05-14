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

# global/ local fuse ranking
dat_rank <- read_rds(here("data", 
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


# assign pcoa to realms
dat_spp_realm <- list_spp_realm %>% 
  map(~ left_join(.x, dat_pcoa)) %>% 
  bind_rows()

# calculate overall convex hull for axis 1 and 2
dat_hull_1_ov <- dat_spp_realm %>% 
  slice(chull(A1, A2)) %>% 
  select(-realm)

# compare to higher ranking species ---------------------------------------

# get the species that rank higher locally
dat_higher_rank <- dat_rank %>% 
  right_join(dat_rank %>% 
              filter(species == "Chelonia mydas") %>% 
              rename(FUSE_chelonia = FUSE_local) %>% 
              select(realm, FUSE_chelonia)) %>% 
  group_by(realm) %>% 
  mutate(is_higher_local = if_else(FUSE_local > FUSE_chelonia, 
                                   1, 0)) %>% 
  ungroup() %>% 
  distinct(realm, species, is_higher_local)

# visualise 
plot_higher_rank <- dat_spp_realm %>%
  left_join(dat_higher_rank) %>% 
  mutate(is_higher_local = as.character(is_higher_local)) %>% 
  ggplot(aes(A1, A2)) +
  geom_polygon(fill = NA,
               colour = "grey20",
               data = dat_hull_1_ov) +
  geom_point(aes(colour = is_higher_local)) +
  geom_point(data = dat_spp_realm %>% 
               filter(species == "Chelonia mydas"), 
             colour = colour_coral) +
  scale_color_manual(values = c(colour_grey, colour_yellow, colour_grey)) +
  facet_wrap(~realm) +
  theme(legend.position = "none")

# save plot
ggsave(plot_higher_rank, 
       filename = here("figures",
                       "main", 
                       "green_sea_turtle",
                       "higher_ranking_sp.pdf"),
       width = 183*1.2, height = 100*1.5,
       units = "mm",
       bg = "white")



