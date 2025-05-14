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
    # compute distinctiveness
    left_join(as_tibble(distinctiveness_global(dist_mat_realm)), 
              by = join_by(species))
  
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

# 
# # save rankings
# list_fuse_local %>%
#   compact() %>%
#   map_dfr(~ .x %>%
#             drop_na(FUSE_local) %>%
#             arrange(desc(FUSE_local)) %>%
#             slice_head(n = 5) %>%
#             mutate(fuse_rank = row_number()) %>%
#             select(species, fuse_rank)) %>%
#   group_by(species) %>%
#   summarise(local_rank = which.max(tabulate(fuse_rank)), 
#             local_rank_min = min(fuse_rank), 
#             local_rank_max = max(fuse_rank)) %>% 
#   write_rds(here("data",
#                  "ranking_variation_per_species.rds"))

 
# # same for un-summarized results
# list_fuse_local %>%
#   bind_rows(.id = "id") %>% 
#   left_join(tibble(id = unique(.$id), 
#                    latidude = dat_upscld$latitude_y[!map_lgl(list_fuse_local, is.null)], 
#                    longitude = dat_upscld$longitude_x[!map_lgl(list_fuse_local, is.null)])) %>% 
#   write_rds(here("data",
#                  "ranking_variation_per_species_umsummarized.rds"),
#             compress = "gz")




# create plots ------------------------------------------------------------

# select species and coordinates
dat_sp <- bind_rows(list_fuse_local) %>%
  filter(species == "Balaenoptera musculus") %>%
  bind_cols(dat_upscld %>%
              rename(sp = "Balaenoptera musculus") %>%
              filter(sp == 1) %>%
              select(contains("itude")))

plot_metric <- dat_sp %>% 
  rename("FDi_std_local" = "global_di") %>% 
  pivot_longer(cols = contains("std_local"), 
               names_to = "metric") %>% 
  mutate(metric = word(metric, sep = "_"), 
         metric = ordered(metric, levels = c("FSp", "FDi", "FUn"))) %>% 
  ggplot() +
  geom_raster(aes(x = longitude_x, 
                  y = latitude_y, 
                  fill = value)) +
  geom_sf(data = world_map_sf,
          col = NA,
          fill = "grey90",
          size = 0.1) +
  scale_fill_gradient2(
    low = colour_mint,
    mid = colour_purple,
    high = colour_yellow,
    midpoint = 0.5, 
    name = NULL) +
  scale_y_continuous(breaks = seq(-90, 90, 30),
                     expand = c(0, 0)) +
  labs(x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(min(dat_sp$longitude_x), max(dat_sp$longitude_x)),
           ylim = c(min(dat_sp$latitude_y), max(dat_sp$latitude_y)),
           expand = FALSE) +
  facet_wrap(~metric, ncol = 1)


# create functional space plots
list_spaces <- list(max, min) %>%
  map2(
    .x = .,
    .y = c("Highest FUSE cell", "Lowest FUSE cell"),
    .f = ~ dat_sp %>%
      filter(FUSE_local == .x(FUSE_local, na.rm = TRUE)) %>%
      slice_head(n = 1) %>%
      {
        filter(dat_upscld,
               longitude_x == .$longitude_x,
               latitude_y == .$latitude_y)
      } %>%
      apply(., 1, function(x) {
        colnames(.)[which(x == 1)]
      }) %>%
      .[, ] %>%
      enframe(value = "species", name = NULL) %>%
      distinct(species) %>%
      left_join(pcoa %>%
                  as_tibble(rownames = "species")) %>%
      mutate(
        colour_id = if_else(species == sel_species, "yes", "no"),
        A1_center = mean(A1),
        A2_center = mean(A2)
      ) %>%
      arrange(colour_id) %>%
      ggplot(aes(A1, A2)) +
      geom_polygon(
        fill = NA,
        colour = "grey20",
        data = . %>%
          slice(chull(A1, A2))
      ) +
      geom_point(
        aes(fill = colour_id),
        shape = 21,
        size = 2,
        stroke = 1,
        colour = "grey20"
      ) +
      geom_point(
        aes(A1_center, A2_center),
        fill = colour_yellow,
        size = 4,
        stroke = 1.5,
        shape = 21,
        colour = "white"
      ) +
      scale_fill_manual(values = c(colour_grey, colour_coral)) +
      labs(
        x = "PCoA1",
        y = "PCoA2",
        title = .y,
        colour = NULL
      ) +
      theme(legend.position = "none")
  )




# patch together
plot_bwhale <- plot_metric /
  free(list_spaces[[1]] + list_spaces[[2]]) +
  plot_layout(heights = c(2, 1))

# save plot
ggsave(plot_bwhale, 
       filename = here("figures",
                       "main", 
                       "fig_3.pdf"),
       width = 120, height = 200,
       units = "mm",
       bg = "white")
