library(here)
library(tidyverse)
library(sf)
library(letsR)


# IUCN --------------------------------------------------------------------


# IUCN data - species distributions 
dat_iucn <- read_sf(here("data", 
                         "iucn", 
                         "IUCN_megafauna_spp_distributions.shp"))

# Set coordinates reference system
st_transform(dat_iucn, st_crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) 


# Generate presence-absence matrix at 0.5 resolution
dat_presabs <- lets.presab(dat_iucn,
                           show.matrix = TRUE,
                           xmn = -180, xmx = 180,
                           ymn = -90, ymx = 90,
                           resol = 0.5,
                           cover = 0.01,
                           count = TRUE)

# save pres-abs matrix
write.table(dat_presabs, 
            here("data",
                 "IUCN",
                 "megafauna_IUCN_presabs_0.5res.txt"),
            sep = ", ",
            row.names = F)


# dat_presabs <- read_delim(here("data",
#                                "IUCN",
#                                "megafauna_IUCN_presabs_0.5res.txt"))


# aquamaps ----------------------------------------------------------------

# global pool of marine megafauna species
megafauna_spp <- read_csv(here("data",
                               "megafauna_species.csv"))


# get the species that are missing from IUCN
missing_spp <- dat_presabs %>% 
  colnames() %>% 
  .[-c(1, 2)] %>% 
  setdiff(megafauna_spp$species, .)

# species ID occurrence probabilities per half-degree cells
dat_aqua <- read_csv(here("data", 
                          "Aquamaps",
                          "hcaf_species_native.csv"))

# get list of synonyms
dat_syn <- read_csv2(here("data", 
                         "synonyms_lookup_table.csv"))

# species names synonyms
dat_aqua_id <- read_csv(here("data",
                             "Aquamaps",
                             "speciesoccursum.csv")) %>% 
  mutate(species = paste(Genus, Species)) %>% 
  select(SpeciesID, species) %>% 
  left_join(dat_syn, by = c("species" = "synonym.1")) %>% 
  mutate(species_corrected = coalesce(species.y, species)) %>% 
  select(-c(species, species.y)) %>%
  rename(species = "species_corrected")

# add to datafile
dat_aqua <- dat_aqua %>% 
  left_join(dat_aqua_id) %>% 
  # remove species that are in IUCN
  filter(species %in% missing_spp) 

# which species are still missing 
missing_spp %>% 
  setdiff(., unique(dat_aqua$species)) # 10 species are missing

# transform probabilities in 1-0 (pres-abs) data, 
# based on a 50% probability threshold
dat_aqua <- dat_aqua %>% 
  select(CenterLat, CenterLong , Probability, species) %>%
  mutate(Probability  = replace(Probability, Probability > 0, 1)) %>%
  # mutate(Probability  = replace(Probability, Probability <= 0.5, 0)) %>% 
  # transform into a spatial object
  filter(Probability == 1) %>%
  rename(latitude = "CenterLat" ,
         longitude = "CenterLong")

dat_presabs_aqua <- dat_aqua %>% 
  select(longitude, 
         latitude) %>% 
  lets.presab.points(.,
                     species = dat_aqua$species, 
                     show.matrix = TRUE,
                     xmn = -180,
                     xmx = 180,
                     ymn = -90,
                     ymx = 90,
                     crs = "+proj=longlat +datum=WGS84",
                     resol = 0.5) %>%
  as_tibble()


# combine with IUCN
dat_presabs <- expand_grid(longitude_x = seq(-180, 180, by = 0.25), 
            latitude_y = seq(-90, 90, by = 0.25)) %>% 
  left_join(janitor::clean_names(dat_presabs)) %>% 
  left_join(janitor::clean_names(dat_presabs_aqua)) %>% 
  replace(is.na(.), 0) 


# save as rds due to size
dat_presabs %>% 
  write_rds(here("data", 
                 "presabs_05_res.rds"), 
            compress = "gz")
