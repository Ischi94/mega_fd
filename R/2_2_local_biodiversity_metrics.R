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

# loop through realms
for (i in 1:length(spp_per_grid)) {
  
  # get species pool per realm 
  spp_per_cell <- spp_per_grid[[i]] %>%
    enframe(value = "species",
            name = NULL) %>% 
    distinct(species)  
  
  # subset distance trait matrix to realm species pool
  spec_to_keep <- distance_trait_matrix %>%
    as.matrix() %>%
    rownames() %in% spp_per_cell$species
  dist_mat_realm <- as.matrix(distance_trait_matrix)[spec_to_keep, spec_to_keep]
  
  # subset trait space to realm species pool
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
    summarise(across(-species, ~mean(.x, na.rm = TRUE)))
  
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
                        label, 
                        breaks, 
                        limits) {
  
  dat_local_metric %>%
    ggplot() +
    geom_raster(aes(x = longitude_x,
                    y = latitude_y, 
                    fill = !!enquo(metric_col))) +
    geom_sf(data = world_map_sf, col = NA, fill = "white", size = 0.1) +
    scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10,"RdBu")), 
                         name = label, 
                         breaks = breaks, 
                         limits = limits) +
    labs(x = "Longitude", y = "Latitude") +
    coord_sf(xlim = c(-180, 180), ylim= c(-89, 85), expand=FALSE) +
    scale_x_continuous(breaks = seq(-180, 180, 30), 
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(-90, 90, 30), 
                       expand = c(0, 0)) 
}


# FUSE
plot_fuse <- plot_metric(FUSE_local, "FUSE", 
                         breaks = c(0, 0.5, 1),
                         limits = c(0, 1.25)) 

# uniqueness
plot_uniq <- plot_metric(FUn_std_local , "FUn", 
                         breaks = seq(0.1, 0.5, by = 0.1), 
                         limits = c(0.09, 0.52)) 

# specialisation
plot_spec <- plot_metric(FSp_std_local , "FSp", 
                         breaks = seq(0.3, 0.6, by = 0.1), 
                         limits = c(0.19, 0.71))


# patch together
plot_final <- plot_fuse / plot_uniq / plot_spec


# save
ggsave(plot_final, 
       filename = here("figures",
                       "main",
                       "1_1_metrics_local.pdf"),
       width = 183, height = 100*2,
       units = "mm",
       bg = "white")



# visualise fuse metrics --------------------------------------------------


# mean fuse per cell
plot_fuse <- plot_metric(fuse, "FUSE") +
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10,"RdBu")), 
                       name = "FUSE", 
                       breaks = c(0, 0.5, 1),
                       limits = c(0, 1.25)) +
  theme(legend.position = "top")

# latitudinal pattern
plot_lat_fuse <- dat_metrics %>% 
  group_by(latitude_y) %>% 
  summarise(fuse = mean(fuse)) %>% 
  ggplot(aes(fuse, latitude_y, 
             colour = fuse)) +
  scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(10,"RdBu")), 
                         limits = c(0, 1.25)) +  
  geom_path(linewidth = 0.5) +
  labs(y = NULL, 
       x = NULL) +
  theme(legend.position = "none", 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

plot_fuse_top <- plot_metric(fuse_top, "Top 25% FUSE\nSR") +
  theme(legend.position = "bottom")

# latitudinal pattern
plot_lat_fuse_top <- dat_metrics %>%
  group_by(latitude_y) %>% 
  summarise(fuse_top = mean(fuse_top)) %>% 
  ggplot(aes(fuse_top, latitude_y, 
             colour = fuse_top)) +
  scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(10,"RdBu")), 
                         limits = c(0, 16)) +  
  geom_path(linewidth = 0.5) +
  labs(y = NULL, 
       x = NULL) +
  theme(legend.position = "none", 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

# patch together
plot_final_fuse <- (plot_fuse + plot_lat_fuse) /
  (plot_fuse_top + plot_lat_fuse_top) 


# save
ggsave(plot_final_fuse, 
       filename = here("figures",
                       "main",
                       "2_fuse.svg"),
       width = 183*2, height = 100*2,
       units = "mm",
       bg = "white")

