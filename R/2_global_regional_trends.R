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
  
  # Uniqueness calculation
  uni_res <- get_indicator(Mat_dist = dist_mat_realm, 
                           nb_NN = 5)
  
  uniqu <- uni_res$Average_uniqueness[,"Mean"] %>% 
    as_tibble(rownames = "species") %>% 
    mutate(specialisation = as.vector(decostand(value,"range"))) %>% 
    select(-value)
  
  # Specialization calculation
  O <- apply(pcoa_realm, 2, mean)
  spe <- apply(pcoa_realm, 1, function(x){sum((x-O)^2)^0.5}) %>% 
    as_tibble(rownames = "species") %>% 
    mutate(uniqueness = as.vector(decostand(value,"range"))) %>% 
    select(-value)
  
  list_fuse_realm[[i]] <- full_join(uniqu, spe) %>% 
    add_column(realm = realm_char[i]) %>% 
    left_join(dat_fuse %>% 
                select(species, FUn_std, 
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
plot_uniq <- dat_realm %>%
  group_by(realm) %>% 
  arrange(desc(uniqueness)) %>% 
  mutate(local_rank = 1:n()) %>% 
  arrange(desc(FUn_std)) %>% 
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
  ggplot(aes(uniqueness, FUn_std)) +
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
  annotate("text", 
           colour = colour_purple, 
           size = 10/.pt, 
           label = "Global most unique", 
           x = 0.5, 
           y = 1.1, 
           fontface = "bold") +
  annotate("text", 
           colour = colour_mint, 
           size = 10/.pt, 
           label = "Local most unique", 
           x = 1.1, 
           y = 0.5, 
           fontface = "bold", 
           angle = 270) +
  scale_fill_manual(values = c(colour_purple,
                               "grey90")) +
  scale_colour_manual(values = c("grey80",
                                 colour_purple, 
                                 colour_mint)) +
  labs(y = "Global\nFunctional Uniqueness [std]", 
       x = "Regional\nFunctional Uniqueness [std]") +
  scale_y_continuous(breaks = c(0, 0.5, 1), 
                     labels = c("0", "0.5", "1")) +
  scale_x_continuous(breaks = c(0, 0.5, 1), 
                     labels = c("0", "0.5", "1")) +
  theme(legend.position = "none")

# specialisation
plot_spec <- dat_realm %>%
  group_by(realm) %>% 
  arrange(desc(specialisation)) %>% 
  mutate(local_rank = 1:n()) %>% 
  arrange(desc(FSp_std)) %>% 
  mutate(global_rank = 1:n()) %>% 
  ungroup() %>% 
  mutate(local_rank = if_else(between(local_rank, 1, 10), 
                              "local", 
                              "not"), 
         global_rank = if_else(between(global_rank, 1, 10), 
                               "global", 
                               "non"), 
         local_rank = case_when(
           local_rank == "local" ~ "local", 
           local_rank == "not" & global_rank == "global" ~ "global", 
           local_rank == "not" & global_rank == "non" ~ "bnothing"
         )) %>% 
  ggplot(aes(specialisation , FSp_std, 
             colour = realm)) +
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
  annotate("text", 
           colour = colour_purple, 
           size = 10/.pt, 
           label = "Global most specialised", 
           x = 0.5, 
           y = 1.1, 
           fontface = "bold") +
  annotate("text", 
           colour = colour_mint, 
           size = 10/.pt, 
           label = "Local most specialised", 
           x = 1.1, 
           y = 0.5, 
           fontface = "bold", 
           angle = 270) +
  scale_fill_manual(values = c(colour_purple,
                               "grey90")) +
  scale_colour_manual(values = c("grey80",
                                 colour_purple, 
                                 colour_mint)) +
  labs(y = "Global\nFunctional Specialisation [std]", 
       x = "Regional\nFunctional Specialisation [std]") +
  scale_y_continuous(breaks = c(0, 0.5, 1), 
                     labels = c("0", "0.5", "1")) +
  scale_x_continuous(breaks = c(0, 0.5, 1), 
                     labels = c("0", "0.5", "1")) +
  theme(legend.position = "none")


# patch together
plot_final <- plot_uniq /
  plot_spec +
  plot_annotation(tag_levels = "a")


# save plot
ggsave(plot_final, 
       filename = here("figures",
                       "main", 
                       "local_global_overlap.pdf"),
       width = 183, height = 100*2,
       units = "mm",
       bg = "white")


# supplementary plots -----------------------------------------------------


# per uniqueness
plot_uniq_realm <- dat_realm %>%
  group_by(realm) %>% 
  arrange(desc(uniqueness)) %>% 
  mutate(local_rank = 1:n()) %>% 
  arrange(desc(FUn_std)) %>% 
  mutate(global_rank = 1:n()) %>% 
  ungroup() %>% 
  mutate(local_rank = if_else(between(local_rank, 1, 10), 
                              "local", 
                              "none"), 
         global_rank = if_else(between(global_rank, 1, 10), 
                               "global", 
                               "not")) %>% 
  ggplot(aes(uniqueness, FUn_std)) +
  geom_abline(intercept = 0,
              slope = 1, 
              linetype = "dotted",
              colour = "grey50") +
  # geom_point(shape = 21,
  #            alpha = 0.7,
  #            stroke = 0.8,
  #            size = 2,
  #            colour = "grey20",
  #            fill = "grey80") +
  geom_point(aes(fill = global_rank
                 ,colour = local_rank
                 ),
             shape = 21,
             alpha = 0.7,
             stroke = 0.8,
             # colour = "grey20",
             size = 2) +
  scale_fill_manual(values = c("purple3",
                               "grey80")) +
  scale_colour_manual(values = c("orange",
                                 "grey20")) +
  facet_wrap(~ realm) +
  labs(y = "Global\nFunctional Uniqueness [std]", 
       x = "Regional\nFunctional Uniqueness [std]") +
  scale_y_continuous(breaks = c(0, 0.5, 1), 
                     labels = c("0", "0.5", "1")) +
  scale_x_continuous(breaks = c(0, 0.5, 1), 
                     labels = c("0", "0.5", "1")) +
  theme_minimal() +
  theme(legend.position = "none", 
        strip.text = element_text(size = 8))

# save plot uniqueness
ggsave(plot_uniq_realm, 
       filename = here("figures",
                       "supplement",
                       "uniq_comp.png"),
       width = 1920, height = 1080*2,
       units = "px",
       bg = "white", device = ragg::agg_png)

# residuals of uniqueness
dat_realm %>%
  mutate(resid = FUn_std - uniqueness) %>% 
  ggplot(aes(resid)) +
  stat_slab(aes(fill = after_stat(x >= 0))) +
  geom_vline(xintercept = 0)


# specialisation
plot_spec_realm <- dat_realm %>%
  group_by(realm) %>% 
  arrange(desc(specialisation)) %>% 
  mutate(local_rank = 1:n()) %>% 
  arrange(desc(FSp_std)) %>% 
  mutate(global_rank = 1:n()) %>% 
  ungroup() %>% 
  mutate(local_rank = if_else(between(local_rank, 1, 10), 
                              "local", 
                              "not"), 
         global_rank = if_else(between(global_rank, 1, 10), 
                               "global", 
                               "non")) %>% 
  ggplot(aes(specialisation , FSp_std, 
             colour = realm)) +
  geom_abline(intercept = 0,
              slope = 1, 
              linetype = "dotted",
               colour = "grey50") +
  # geom_point(shape = 21,
  #            alpha = 0.7,
  #            stroke = 0.8,
  #            size = 2,
  #            colour = "grey20",
  #            fill = "grey80") +
  geom_point(aes(fill = global_rank
                 , colour = local_rank
                 ),
             shape = 21,
             alpha = 0.7,
             stroke = 0.8,
             # colour = "grey20",
             size = 2) +
  scale_fill_manual(values = c("purple3",
                               "grey80")) +
  scale_colour_manual(values = c("orange",
                                 "grey20")) +
  facet_wrap(~ realm) +
  labs(y = "Global\nFunctional Specialisation [std]", 
       x = "Regional\nFunctional Specialisation [std]") +
  scale_y_continuous(breaks = c(0, 0.5, 1), 
                     labels = c("0", "0.5", "1")) +
  scale_x_continuous(breaks = c(0, 0.5, 1), 
                     labels = c("0", "0.5", "1")) +
  theme_minimal() +
  theme(legend.position = "none", 
        strip.text = element_text(size = 8))
  
# save plot specialisation
ggsave(plot_spec_realm, 
       filename = here("figures",
                       "supplement", 
                       "speci_comp.png"),
       width = 1920, height = 1080*2,
       units = "px",
       bg = "white", device = ragg::agg_png)

# residuals specialisation
dat_realm %>%
  mutate(resid = FSp_std - specialisation) %>% 
  ggplot(aes(resid)) +
  stat_slab(aes(fill = after_stat(x >= 0))) +
  geom_vline(xintercept = 0) +
  facet_wrap(~ realm)



