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
                .fn = ~ paste0(.x, "_local"))
  
  print(i)
}

# combine in one dataframe
dat_local_metric <- bind_rows(list_fuse_local) %>%
  group_by(species) %>% 
  summarise(across(contains("local"), list("min" = ~min(.x, na.rm = TRUE), 
                                           "max" = ~ max(.x, na.rm = TRUE), 
                                           "mean" = ~ mean(.x, na.rm = TRUE), 
                                           "median" = ~ median(.x, na.rm = TRUE), 
                                           "sd" = ~ sd(.x, na.rm = TRUE)))) %>% 
  rename_with(~str_remove(., '_local')) %>% 
  rename_with(~str_remove(., '_std'))


# save
dat_local_metric %>% 
  write_csv(here("data", 
                 "metric_variation_per_species.csv"))




# create plots ------------------------------------------------------------

# set up function to create visuals for selected species
create_plot <- function(sel_species, 
                        lat_buffer) {
  
  # select species and coordinates
  dat_sp <- bind_rows(list_fuse_local) %>%
    filter(species == sel_species) %>% 
    bind_cols(dat_upscld %>% 
                rename(sp = sel_species) %>% 
                filter(sp == 1) %>% 
                select(contains("itude")))
  
  # create maps of variation
  list_maps <- dat_sp %>%
    colnames() %>%
    .[2:4] %>%
    map2(.x = .,
         .y = c("FUSE", "FUn", "FSp"),
         .f = ~ dat_sp %>%
           ggplot() +
           geom_raster(aes(x = longitude_x,
                           y = latitude_y,
                           fill = !!sym(.x))) +
           geom_sf(data = world_map_sf, col = NA, fill = "white", size = 0.1) +
           scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10,"RdBu")),
                                name = .y) +
           labs(x = "Longitude", y = "Latitude") +
           coord_sf(xlim = c(min(dat_sp$longitude_x)-5, max(dat_sp$longitude_x)+5),
                    ylim = c(min(dat_sp$latitude_y)-5, max(dat_sp$latitude_y)+lat_buffer), 
                    expand = FALSE))

  # create functional space plots
  list_spaces <- list(max, min) %>%
    map2(.x = ., 
         .y = c(RColorBrewer::brewer.pal(10, "RdBu")[1], 
                RColorBrewer::brewer.pal(10,"RdBu")[10]), 
         .f = ~ dat_sp %>%
          filter(FUSE_local == .x(FUSE_local, na.rm = TRUE)) %>%
          slice_head(n = 1) %>%
          { filter(dat_upscld,
                   longitude_x == .$longitude_x,
                   latitude_y == .$latitude_y) } %>%
          apply(.,
                1,
                function(x){colnames(.)[which(x==1
                )]}) %>%
          .[,] %>%
          enframe(value = "species",
                  name = NULL) %>%
          distinct(species) %>%
          left_join(pcoa %>%
                      as_tibble(rownames = "species")) %>%
          mutate(colour_id = if_else(species == sel_species,
                                     "yes", "no"),
                 A1_center = mean(A1),
                 A2_center = mean(A2)) %>%
          ggplot(aes(A1, A2)) +
          geom_polygon(fill = NA,
                       colour = "grey20",
                       data = . %>%
                         slice(chull(A1, A2))) +
          geom_point(aes(fill = colour_id),
                     shape = 21, size = 2,
                     stroke = 1, colour = "grey20") +
          geom_point(aes(A1_center, A2_center),
                     fill = colour_yellow,
                     size = 4,
                     stroke = 1.5,
                     shape = 21,
                     colour = "white") +
          scale_fill_manual(values = c(colour_grey, colour_coral)) +
          labs(x = "PCoA1",
               y = "PCoA2",
               colour = NULL) +
          theme(legend.position = "none", 
                plot.background = element_rect(colour = .y))
    )

  # combine plots
  list_full <- c(list_maps, list_spaces)

  return(list_full)
}

# identify most varying species
dat_local_metric %>% 
  arrange(desc(FUSE_sd)) %>% 
  select(species)

# sturgeon
plot_list_stella <- create_plot("Acipenser stellatus", 
                                lat_buffer = 5)

# patch together
plot_stella <- plot_list_stella[[1]] +
  inset_element(plot_list_stella[[4]], 
                0.1, 0.6, 0.5, 1) +
  inset_element(plot_list_stella[[5]], 
                0.6, 0.6, 1, 1) +
  plot_annotation(tag_levels = "a")

# save
ggsave(plot_stella, 
       filename = here("figures",
                       "main",
                       "8_sturgeon_variation.pdf"),
       width = 183, height = 100*2,
       units = "mm",
       bg = "white")

# green sea turtle
plot_list_mydas <- create_plot("Chelonia mydas", 
                               lat_buffer = 80) 

# patch together
plot_mydas <- plot_list_mydas[[1]] +
  inset_element(plot_list_mydas[[4]], 
                0.1, 0.6, 0.5, 1) +
  inset_element(plot_list_mydas[[5]], 
                0.6, 0.6, 1, 1) +
  plot_annotation(tag_levels = "a")

# save
ggsave(plot_mydas, 
       filename = here("figures",
                       "main",
                       "9_gs_turtle_variation.pdf"),
       width = 183, height = 100*2,
       units = "mm",
       bg = "white")
