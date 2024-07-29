library(here)
library(tidyverse)
library(sf)
library(ade4)
library(ggdist)
library(vegan)
library(patchwork)

# get functions
source(here("R", 
            "functions.R"))


# load data ---------------------------------------------------------------


# distance-trait matrix
load(here("data", 
          "distance_trait_matrix.RData"))

# fuse metrics
dat_fuse <- read_rds(here("data",
                          "fuse_metrics_per_spp.rds"))

# local metric variation
dat_var <- read_rds(here("data", 
                         "ranking_variation_per_species_umsummarized.rds"))

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

# extract functional space axes
pcoa <- dudi.pco(quasieuclid(distance_trait_matrix),
                 scannf = FALSE,
                 nf = 4) %>% 
  pluck("li")

# provinces of the world
MEOW_sf <- read_sf(here("data", 
                        "meow", 
                        "meow_ecos.shp"))


# IUCN categories
dat_iucn <- read_csv(here("data",
                          "megafauna_traits.csv")) %>%
  select(species, higher_classification, IUCN)

# global distinctiveness
dist_glob <- read_rds(here("data",
                           "global_distinctiveness.rds"))

# provincial distinctiveness
dist_prov <- read_rds(here("data", 
                           "provincial_distinctiveness.rds"))

# local distinctiveness.
dist_local <- read_rds(here("data", 
                           "local_distinctiveness.rds"))

# metrics per realm ----------------------------------------------------------

# assign grid cells to realms 
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
  spec_to_keep <- distance_trait_matrix %>%
    as.matrix() %>%
    rownames() %in% spp_per_realm$species
  dist_mat_realm <- as.matrix(distance_trait_matrix)[spec_to_keep, spec_to_keep]
  
  # subset trait space to realm species pool
  pcoa_realm <- na.omit(pcoa[spp_per_realm$species, ]) 
  
  
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
  list_fuse_realm[[i]] <- get_FUSE(Mat_dist = dist_mat_realm,
                       Coords = pcoa_realm,
                       GE = ge) %>% 
    as_tibble(rownames = "species") %>% 
    arrange(species) %>%
    select(species, FUSE, FUn_std, FSp_std) %>% 
    rename_with(.cols = -species, 
                .fn = ~ paste0(.x, "_local")) %>% 
    add_column(realm = realm_char[i]) %>% 
    left_join(dat_fuse %>% 
                select(species, FUSE, FUn_std, 
                       FSp_std))
  
}

# combine in one dataframe
dat_realm <- bind_rows(list_fuse_realm)

# save 
dat_realm %>% 
  write_rds(here("data", 
                 "global_regional_trends.rds"))
  

# add rank-rank plot ------------------------------------------------------


# global to province
dat_glob_loc <- dat_fuse %>%
  arrange(desc(FUSE)) %>%
  filter(species != "Tridacna gigas") %>% 
  slice_head(n = 5) %>%
  mutate(global_rank = row_number()) %>%
  select(species, global_rank) %>% 
  left_join(dat_realm %>% 
              group_by(realm) %>% 
              select(species, FUSE_local, realm) %>% 
              arrange(desc(FUSE_local)) %>% 
              mutate("local_rank" = row_number()) %>%  
              select(species, local_rank, realm)) %>% 
  pivot_longer(cols = -c(species, realm),
               names_to = "scale",
               values_to = "rank") %>%
  mutate(scale = as.integer(as.factor(scale))) %>% 
  mutate(rank_log = log(rank))



# add local and plot using mode
plot_mode <- dat_glob_loc %>%
  bind_rows(read_rds(here("data",
                          "ranking_variation_per_species.rds")) %>% 
              mutate(scale = 3, rank = local_rank, rank_log = log(rank)) %>% 
              select(species, scale, rank, rank_log) %>% 
              right_join(dat_glob_loc %>% 
                           distinct(species, realm))) %>% 
  ggplot(aes(x = scale,
             y = rank_log)) +
  geom_line(aes(group = interaction(species, realm), 
                colour = species), 
            position = position_dodge(width = 0.1)) +
  geom_label(aes(x = scale, 
                 y = rank_log, 
                 label = rank),
             position = position_dodge(width = 0.2), 
            label.size = 0,
            label.padding	= unit(0.5, "lines"),
            label.r = unit(1, "lines"),
            size = 8/.pt, 
            data = . %>%  
              filter(rank %in% c(1:5, 7, 9, 11, 15, 18, 28))) +
  geom_text(aes(label = species, 
                colour = species),
            data = dat_glob_loc %>%
              filter(scale == 1) %>% 
              distinct(species, rank_log, scale),
            size = 9/.pt,
            hjust = 1,
            position = position_nudge(x = -0.2)) +
  labs(y = NULL,
       title = "FUSE rank",
       x = NULL) +
  scale_color_manual(values = c(colour_coral, 
                                "#8D4D5A", 
                                colour_mint, 
                                "#CDA25D", 
                                "#FF8043"), 
                     na.translate = F, 
                     name = NULL) +
  scale_y_reverse(breaks = NULL) +
  coord_cartesian(xlim = c(-0.01, 3)) +
  scale_x_continuous(breaks = c(1, 2, 3),
                     labels = c("Global", "Province", "Local")) +
  theme(legend.position = "none",
        legend.text = element_text(colour = "grey20", size = 7), 
        legend.key.size = unit(3, "mm"),
        panel.grid = element_blank(), 
        plot.title = element_text(colour = "grey20", size = 10, 
                                  face = "bold", 
                                  hjust = 0.05))


# use highest rank
plot_high <- dat_glob_loc %>%
  bind_rows(read_rds(here("data",
                          "ranking_variation_per_species.rds")) %>% 
              mutate(scale = 3, rank = local_rank_max, rank_log = log(rank)) %>% 
              select(species, scale, rank, rank_log) %>% 
              right_join(dat_glob_loc %>% 
                           distinct(species, realm))) %>% 
  ggplot(aes(x = scale,
             y = rank_log)) +
  geom_line(aes(group = interaction(species, realm), 
                colour = species), 
            position = position_dodge(width = 0.1)) +
  geom_label(aes(x = scale, 
                 y = rank_log, 
                 label = rank),
             label.size = 0,
             label.padding	= unit(0.45, "lines"),
             label.r = unit(1, "lines"),
             size = 8/.pt, 
             data = . %>%  
               filter(rank %in% c(1:5, 7, 9, 11, 15, 18, 28))) +
  geom_text(aes(label = species, 
                colour = species),
            data = dat_glob_loc %>%
              filter(scale == 1) %>% 
              distinct(species, rank_log, scale),
            size = 9/.pt,
            hjust = 1,
            position = position_nudge(x = -0.2)) +
  labs(y = "FUSE rank",
       x = NULL) +
  scale_color_manual(values = c(colour_coral, 
                                "#8D4D5A", 
                                colour_mint, 
                                "#CDA25D", 
                                "#FF8043"), 
                     na.translate = F, 
                     name = NULL) +
  scale_y_reverse(breaks = NULL) +
  coord_cartesian(xlim = c(-0.01, 3)) +
  scale_x_continuous(breaks = c(1, 2, 3),
                     labels = c("Global", "Province", "Local")) +
  theme(legend.position = "none",
        legend.text = element_text(colour = "grey20", size = 7), 
        legend.key.size = unit(3, "mm"),
        panel.grid = element_blank())
            



# rank-rank correlation ---------------------------------------------------


# per province rank-correlation
plot_cor_prov <- dat_realm %>%
  left_join(dist_glob %>% 
              full_join(dist_prov %>% 
                          rename(prov_di = global_di))) %>% 
  group_by(realm) %>%
  summarise(FUn = cor(FUn_std_local, FUn_std, use = "complete.obs", method = "kendall"), 
            FSp = cor(FSp_std_local, FSp_std, use = "complete.obs", method = "kendall"), 
            FDi = cor(prov_di, global_di, use = "complete.obs", method = "kendall"),
            nSp = n()) %>% 
  bind_rows(dat_realm %>% 
              left_join(dist_glob %>% 
                          full_join(dist_prov %>% 
                                      rename(prov_di = global_di))) %>% 
              summarise(FUn = cor(FUn_std_local, FUn_std, use = "complete.obs", method = "kendall"), 
                        FSp = cor(FSp_std_local, FSp_std, use = "complete.obs", method = "kendall"),
                        FDi = cor(prov_di, global_di, use = "complete.obs", method = "kendall"),
                        nSp = n()) %>% 
              add_column(realm = "Overall")) %>% 
  pivot_longer(cols = -c(realm, nSp), 
               names_to = "metric", 
               values_to = "correl") %>% 
  left_join(tibble(lat = MEOW_sf %>%
                     st_centroid() %>%
                     st_coordinates() %>%
                     .[, 2], 
                   realm = MEOW_sf$REALM) %>% 
              group_by(realm) %>% 
              summarise(mean_lat = mean(lat)) %>% 
              add_row(realm = "Overall", 
                      mean_lat = -80)) %>% 
  mutate(realm = fct_reorder(realm, mean_lat)) %>% 
  mutate(shape_id = if_else(realm == "Overall", 23, 21)) %>% 
  ggplot(aes(correl, realm, 
             fill = metric, 
             size = nSp, 
             shape = shape_id)) +
  geom_point(colour = "grey20") +
  scale_fill_manual(values = c(colour_purple,
                               colour_mint,
                               colour_yellow)) +
  scale_shape_identity() +
  scale_size_continuous(guide = "none", 
                        range = c(1, 6)) +
  labs(y = NULL, 
       x = expression(tau), 
       fill = NULL, 
       title = "Global-Province correlation") +
  theme(legend.position = "none", 
        plot.title = element_text(colour = "grey20", size = 10, 
                                  face = "bold", 
                                  hjust = -3.5))


# per local cell rank-correlation
plot_cor_loc <- dat_var %>%
  left_join(dat_fuse %>% 
              select(species, FUSE, 
                     FUn_std, FSp_std)) %>% 
  left_join(dist_glob %>% 
              full_join(dist_local %>% 
                          rename(local_di = global_di)) %>% 
              select(-contains("itude"))) %>% 
  drop_na(FUn_std_local, FUn_std, FSp_std_local, FSp_std, global_di, local_di) %>% 
  group_by(id) %>% 
  summarise(FUn = cor(FUn_std_local, FUn_std, use = "complete.obs", method = "kendall"), 
            FSp = cor(FSp_std_local, FSp_std, use = "complete.obs", method = "kendall"),
            FDi = cor(local_di, global_di, use = "complete.obs", method = "kendall"),
            nSp = n()) %>% 
  left_join(dat_var %>% 
              distinct(id, latidude)) %>% 
  pivot_longer(cols = -c(id, latidude, nSp), 
               names_to = "metric", 
               values_to = "correl") %>% 
  ggplot(aes(correl, latidude, 
             fill = metric, 
             size = nSp)) +
  geom_point(shape = 21, 
             alpha = 0.1, 
             colour = "grey20") +
  scale_fill_manual(values = c(colour_purple,
                               colour_mint,
                               colour_yellow)) +
  scale_size_continuous(guide = "none", 
                        range = c(0.8, 4)) +
  labs(y = "Latitude", 
       x = expression(tau), 
       fill = NULL, 
       title = "Global-Local correlation") +
  theme(legend.position = c(-0.2, -0.12), 
        legend.direction = "horizontal",
        plot.title = element_text(colour = "grey20", size = 10, 
                                  face = "bold", 
                                  hjust = -0.6)) +
  guides(fill = guide_legend(override.aes = list(alpha = 1, 
                                                 size = 2), 
                             reverse = TRUE))


# merge together
plot_final <- free(plot_mode) /
  (plot_cor_prov + plot_cor_loc) +
  plot_layout(heights = c(1.5, 1))

# save plot
ggsave(plot_final, 
       filename = here("figures",
                       "main", 
                       "rank_correlation.pdf"),
       width = 183, height = 200,
       units = "mm",
       bg = "white")
