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
  

# visualise rank-rank -----------------------------------------------------

# fuse
plot_fuse <- dat_realm %>%
  group_by(realm) %>% 
  arrange(desc(FUSE_local)) %>%
  mutate(local_rank = 1:n()) %>% 
  arrange(desc(FUSE)) %>% 
  mutate(global_rank = 1:n()) %>% 
  ungroup() %>% 
  mutate(local_rank = if_else(between(local_rank, 1, 10), 
                              "local", 
                              "none"), 
         global_rank = if_else(between(global_rank, 1, 10), 
                               "global", 
                               "not"), 
         local_rank = case_when(
           local_rank == "local" ~ "local", 
           local_rank == "none" & global_rank == "global" ~ "global", 
           local_rank == "none" & global_rank == "not" ~ "bnothing"
         )) %>% 
  ggplot(aes(FUSE_local, FUSE)) +
  geom_abline(intercept = 0,
              slope = 1, 
              linetype = "dotted",
              colour = "grey50") +

  geom_point(aes(fill = global_rank,
                 colour = local_rank),
             shape = 21,
             alpha = 0.9,
             stroke = 0.8,
             size = 2) +
  scale_fill_manual(values = c(colour_purple,
                               "grey90")) +
  scale_colour_manual(values = c("grey80",
                                 colour_purple, 
                                 colour_mint)) +
  labs(y = "Global FUSE", 
       x = "Regional FUSE") +
  scale_y_continuous(breaks = c(0, 1, 2, 3), 
                     labels = c("0", "1", "2", "3")) +
  scale_x_continuous(breaks = c(0, 1, 2, 3), 
                     labels = c("0", "1", "2", "3")) +
  facet_wrap(~ realm) +
  theme(legend.position = "none")

# save plot
ggsave(plot_fuse, 
       filename = here("figures",
                       "main", 
                       "6_agreement_per_province.pdf"),
       width = 183*1.2, height = 100,
       units = "mm",
       bg = "white")



# add rank-rank plot ------------------------------------------------------

# first from local to global
dat_loc_glob <- dat_realm %>% 
  group_by(realm) %>% 
  select(species, FUSE_local, realm) %>% 
  arrange(desc(FUSE_local)) %>% 
  slice_head(n = 5) %>% 
  mutate("local_rank" = row_number()) %>%  
  select(species, local_rank, realm) %>% 
  left_join(dat_fuse %>%
              arrange(desc(FUSE)) %>%
              mutate(global_rank = row_number()) %>%
              select(species, global_rank)) %>% 
  mutate(local_rank = as.integer(local_rank)) %>% 
  pivot_longer(cols = -c(species, realm),
               names_to = "scale",
               values_to = "rank") %>%
  mutate(scale = as.integer(as.factor(scale))) %>% 
  ungroup() %>% 
  mutate(rank_log = log(rank))

# visualise
# plot_loc_glob <- 
dat_loc_glob %>% 
  ggplot(aes(x = scale,
             y = rank_log)) +
  geom_line(aes(group = species), 
            colour = colour_grey) +
  geom_text(aes(x = 0.7,
                y = rank_log,
                label = rank),
            data = dat_loc_glob %>% 
              filter(scale == 1),
            size = 8/.pt) +
  geom_text(aes(label = species),
            data = dat_loc_glob %>%
              filter(scale == 2),
            size = 7/.pt,
            hjust = 1,
            position = position_nudge(x = -0.8)) +
  labs(y = "FUSE rank",
       x = NULL) +
  scale_x_continuous(trans = "reverse",
                     breaks = c(1, 4),
                     labels = c("Global", "Local")) +
  scale_y_reverse(breaks = NULL) +
  coord_cartesian(xlim = c(6, 0.6)) +
  facet_wrap(~ realm, 
             scales = "free_y") +
  theme(legend.position = "none",
        panel.grid = element_blank())

# save plot
ggsave(plot_loc_glob, 
       filename = here("figures",
                       "main", 
                       "7_rank_local_global.pdf"),
       width = 183*1.2, height = 100*1.5,
       units = "mm",
       bg = "white")


# same for global to local
dat_glob_loc <- dat_fuse %>%
  arrange(desc(FUSE)) %>%
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


# visualise
plot_glob_loc <- dat_glob_loc %>%
  ggplot(aes(x = scale,
             y = rank_log)) +
  geom_line(aes(group = interaction(species, realm), 
                colour = realm), 
            position = position_dodge(width = 0.1)) +
  geom_text(aes(x = 2.3,
                y = rank_log,
                label = rank),
            data = dat_glob_loc %>% 
              filter(scale == 2) %>% 
              distinct(scale, rank, rank_log),
            size = 8/.pt) +
  geom_text(aes(label = species),
            data = dat_glob_loc %>%
              filter(scale == 1) %>% 
              distinct(species, rank_log, scale),
            size = 7/.pt,
            hjust = 1,
            position = position_nudge(x = -0.5)) +
  labs(y = "FUSE rank",
       x = NULL) +
  scale_colour_discrete(na.translate = F, 
                        name = NULL) +
  scale_y_reverse(breaks = NULL) +
  coord_cartesian(xlim = c(-1, 3)) +
  scale_x_continuous(breaks = c(1, 2),
                     labels = c("Global", "Local")) +
  theme(legend.position = "bottom",
        legend.text = element_text(colour = "grey20", size = 7), 
        legend.key.size = unit(3, "mm"),
        panel.grid = element_blank())

# save plot
ggsave(plot_glob_loc, 
       filename = here("figures",
                       "main", 
                       "7_rank_global_local.pdf"),
       width = 183, height = 100,
       units = "mm",
       bg = "white")
