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

# extract functional space axes
pcoa <- dudi.pco(quasieuclid(distance_trait_matrix),
                 scannf = FALSE,
                 nf = 4) %>% 
  pluck("li")

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



# list of all species present per grid 
spp_per_grid <- apply(dat_upscld, 
                      1, 
                      function(x){colnames(dat_presabs)[which(x==1
                      )]}) %>% 
  # remove empty grids
  compact()

# set up empty list
list_fuse_local <- vector(mode = "list", length = nrow(dat_upscld))

# loop through cells
for (i in 1:length(spp_per_grid)) {
  
  # get species pool per cell 
  spp_per_cell <- spp_per_grid[[i]] %>%
    enframe(value = "species",
            name = NULL) %>% 
    distinct(species)  
  
  # subset distance trait matrix to cell species pool
  spec_to_keep <- distance_trait_matrix %>%
    as.matrix() %>%
    rownames() %in% spp_per_cell$species
  dist_mat_realm <- as.matrix(distance_trait_matrix)[spec_to_keep, spec_to_keep]
  
  # subset trait space to cell species pool
  pcoa_realm <- na.omit(pcoa[spp_per_cell$species, ]) 
  
  
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
  list_fuse_local[[i]] <- get_FUSE(Mat_dist = dist_mat_realm,
                                   Coords = pcoa_realm,
                                   GE = ge) %>% 
    as_tibble(rownames = "species") %>% 
    arrange(species) %>%
    select(species, FUSE, FUn_std, FSp_std) %>% 
    rename_with(.cols = -species, 
                .fn = ~ paste0(.x, "_local")) %>% 
    summarise(across(-species, list(~mean(.x, na.rm = TRUE), 
                                    ~max(.x, na.rm = TRUE), 
                                    ~median(.x, na.rm = TRUE), 
                                    ~sd(.x, na.rm = TRUE))))
  
  print(i)
}

# combine in one dataframe
dat_local_metric <- bind_rows(list_fuse_local) %>% 
  bind_cols(dat_upscld %>% 
              mutate(sumVar = rowSums(select(., -contains("itude")))) %>% 
              filter(sumVar > 0) %>% 
              select(contains("itude")))

# save
dat_local_metric %>% 
  write_rds(here("data", 
                 "functional_metrics_local.rds"), 
            compress = "gz")

# dat_local_metric <- read_rds(here("data", 
#                                   "functional_metrics_local.rds"))


# visualise ---------------------------------------------------------------

# plot function
plot_metric <- function(metric_col, 
                        label) {
  
  dat_local_metric %>%
    ggplot() +
    geom_raster(aes(x = longitude_x,
                    y = latitude_y, 
                    fill = !!(metric_col))) +
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



list_names <- paste(rep(c("FUSE", "FUn", "FSp"), each = 4), 
                    c("[mean]", "[max]", "[median]", "[sd]"))

# set up list to save
list_plots <- vector("list", length = length(list_names))

# iterate through
for (i in seq_along(list_names)) {
  
  list_plots[[i]] <- plot_metric(sym(colnames(dat_local_metric)[[i]]), list_names[[i]])
  
}

# save figures
walk2(list_plots, 
      list_names %>% map(~str_replace_all(.x, " ", "_") %>%
                           str_c(., ".pdf")),
      ~ggsave(plot = .x, filename = .y, 
              path = here("figures", "local_metrics"), 
              width = 183, height = 100, 
              device = "pdf", 
              units = "mm", bg = "white"))
