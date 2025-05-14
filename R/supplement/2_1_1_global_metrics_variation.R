library(here)
library(tidyverse)
library(sf)
library(ade4)
library(patchwork)
library(furrr)

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
                str_replace_all("_", " ")) %>% 
  # unique cell id
  unite(grid, longitude_x, latitude_y) %>% 
  column_to_rownames("grid")



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



# read in fuse metrices 
dat_fuse <- read_rds(here("data", 
                          "fuse_metrics_per_spp.rds"))

# extract functional space axes
pcoa <- dudi.pco(quasieuclid(distance_trait_matrix),
                 scannf = FALSE,
                 nf = 4)



# list of all species present per grid 
spp_per_grid <- apply(dat_presabs, 
                      1, 
                      function(x){colnames(dat_presabs)[which(x==1
                      )]}) %>% 
  # remove empty grids
  compact()


# calculate metrices ------------------------------------------------------

# set up multicores
plan(multisession, workers = 15)

# compute max FUSE per grid cell
dat_max_fuse <- spp_per_grid %>% 
  future_map(~ dat_fuse %>% 
        select(species, FUSE,
               FUn_std, FSp_std) %>% 
        filter(species %in% .x) %>% 
        summarise(max_fuse = max(FUSE, na.rm = TRUE), 
                  max_uniq = max(FUn_std, na.rm = TRUE), 
                  max_spec = max(FSp_std, na.rm = TRUE), 
                  median_fuse = median(FUSE, na.rm = TRUE), 
                  median_uniq = median(FUn_std, na.rm = TRUE), 
                  median_spec = median(FSp_std, na.rm = TRUE),
                  sd_fuse = sd(FUSE, na.rm = TRUE), 
                  sd_uniq = sd(FUn_std, na.rm = TRUE), 
                  sd_spec = sd(FSp_std, na.rm = TRUE)), 
      .progress = TRUE)






# visualise ---------------------------------------------------------------

# plot function
plot_metric <- function(dataset, 
                        label) {
  
  dataset %>%
    as_tibble() %>% 
    add_column(grid = rownames(dat_presabs)) %>% 
    separate_wider_delim(grid, 
                         names = c("longitude_x",
                                   "latitude_y"),
                         delim = "_") %>%
    mutate(across(c(longitude_x, latitude_y), as.double)) %>% 
    ggplot() +
    geom_raster(aes(x = longitude_x,
                    y = latitude_y, 
                    fill = value)) +
    geom_sf(data = world_map_sf, col = NA, fill = "white", size = 0.1) +
    scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10,"RdBu")), 
                         name = label) +
    labs(x = "Longitude", y = "Latitude") +
    coord_sf(xlim = c(-180, 180), ylim= c(-89, 85), expand=FALSE) +
    scale_x_continuous(breaks = seq(-180, 180, 30), 
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(-90, 90, 30), 
                       expand = c(0, 0)) 
}

# get list of names
list_names <- list("FUSE [max]", "FUn [max]", "FSp [max]", 
                   "FUSE [median]", "FUn [median]", "FSp [median]", 
                   "FUSE [sd]", "FUn [sd]", "FSp [sd]")

# set up list to save
list_plots <- vector("list", length = length(list_names))

# iterate through
for (i in seq_along(list_names)) {
  
  list_plots[[i]] <- plot_metric(map_dbl(dat_max_fuse, i), list_names[[i]])

}

# save figures
walk2(list_plots, 
      list_names %>% map(~str_replace_all(.x, " ", "_") %>%
                           str_c(., ".pdf")),
      ~ggsave(plot = .x, filename = .y, 
              path = here("figures", "global_metrics"), 
              width = 183, height = 100, 
              device = "pdf", 
              units = "mm", bg = "white"))
