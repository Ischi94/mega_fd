library(here)
library(raster)
library(tidyverse)
library(ade4)
library(vegan)

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

# IUCN categories
dat_iucn <- read_csv(here("data",
                          "megafauna_traits.csv")) %>%
  select(species, higher_classification, IUCN) %>%
  mutate(ge = case_when(
    IUCN == "LC" ~ 0, 
    IUCN == "NT" ~ 1, 
    IUCN == "VU" ~ 2,
    IUCN == "EN" ~ 3, 
    IUCN == "CR" ~ 4, 
    .default = NA_integer_
  )) 


# upscale raster ----------------------------------------------------------

# set up function to upscale the presence-absence matrix
# by a specified factor
upscale_grid <- function(upscale_factor) {
  
  dat_rast_agg <- dat_presabs %>%
    rasterFromXYZ() %>%
    aggregate(fact = upscale_factor, 
              fun = max)
  
  dat_rast_agg %>% 
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
  
}


# perform upscaling
list_upscaled_grids <- list(2, # 1x1
                            10, # 5x5
                            20, # 10x10
                            100, # 50x50
                            200) %>% # 100x100
  map(upscale_grid, 
      .progress = TRUE)


# calculate metrics per cell ----------------------------------------------

# list of all species present per grid 
list_spp_per_grid <- list_upscaled_grids %>% 
  map(~ apply(.x, 
              1, 
              function(species){colnames(.x)[which(species==1
              )]}) %>% 
        # remove empty grids
        compact())

# this is a nested list, so we need to apply map two times
list_fuse <- list_spp_per_grid %>% 
  map(~ .x %>%  
        map( function(.x) {
          
          # get cell
          first_cell <- sort(.x)
          
          # subset distance trait matrix to realm species pool
          spec_to_keep <- distance_trait_matrix %>%
            as.matrix() %>%
            rownames() %in% first_cell
          dist_mat_realm <- as.matrix(distance_trait_matrix)[spec_to_keep, spec_to_keep]
          # right order
          dist_mat_realm <- dist_mat_realm[sort(rownames(dist_mat_realm)), 
                                           sort(colnames(dist_mat_realm))]
          
          # subset trait space to realm species pool
          pcoa_realm <- na.omit(pcoa[first_cell, ])
          
          # row.names(as.matrix(dist_mat_realm))
          # row.names(pcoa_realm) 
          
          # prepare iucn categories
          ge <- dat_iucn %>% 
            filter(species %in% first_cell) %>% 
            arrange(species) %>% 
            pull(ge)
          
          
          # get metrices
          get_FUSE(Mat_dist = dist_mat_realm,
                   Coords = pcoa_realm,
                   GE = ge) %>% 
            as_tibble(rownames = "species")
          
          
        }, 
        .progress = TRUE))





# calculate proportions ---------------------------------------------------


# get number of cells per list
cell_nr <- map_dbl(list_spp_per_grid, length)

# proportion of 10 highest local FUSE species
# in the 10 highest global FUSE species
dat_fuse_prop <- list_fuse %>% 
  map(~ .x %>%
        map_dbl( ~ .x %>% 
                   left_join(dat_fuse %>% 
                               select(species, global_FUSE = FUSE) %>% 
                               arrange(desc(global_FUSE)) %>% 
                               rowid_to_column("global_rank"), 
                             by = "species") %>% 
                   arrange(desc(FUSE)) %>% 
                   rowid_to_column("local_rank") %>%
                   filter(between(global_rank, 1, 10)) %>%
                   nrow(.) / 10, 
  .progress = TRUE)) %>% 
  # clean up
  map_dfr(~ .x %>% 
            as_tibble() %>% 
            rename("fuse_prop" = value)) %>% 
  # add cell resolution
  add_column(cell_res = rep(c("1x1", "5x5", 
                              "10x10", "50x50", 
                              "100x100"), 
                            times = cell_nr))
  

# proportion of 10 most unique local species
# in the 10 most unique global species
dat_FUn_prop <- list_fuse %>% 
  map(~ .x %>%
        map_dbl(~ .x %>%
                  left_join(
                    dat_fuse %>%
                      select(species, global_FUn = FUn_std) %>%
                      arrange(desc(global_FUn)) %>%
                      rowid_to_column("global_rank"),
                    by = "species"
                  ) %>%
                  arrange(desc(FUn_std)) %>%
                  rowid_to_column("local_rank") %>%
                              filter(between(global_rank, 1, 10)) %>%
                              nrow(.) / 10, 
                 .progress = TRUE)) %>% 
  # clean up
  map_dfr(~ .x %>% 
            as_tibble() %>% 
            rename("FUn_prop" = value)) 



# proportion of 10 most specialised local species
# in the 10 most specialised global species
dat_FSp_prop <- list_fuse %>% 
  map(~ .x %>%
        map_dbl(~ .x %>%
                  left_join(
                    dat_fuse %>%
                      select(species, global_FSp = FSp_std) %>%
                      arrange(desc(global_FSp)) %>%
                      rowid_to_column("global_rank"),
                    by = "species"
                  ) %>%
                  arrange(desc(FSp_std)) %>%
                  rowid_to_column("local_rank") %>%
                  filter(between(global_rank, 1, 10)) %>%
                  nrow(.) / 10, 
                .progress = TRUE)) %>% 
  # clean up
  map_dfr(~ .x %>% 
            as_tibble() %>% 
            rename("FSp_prop" = value)) 

# combine in one dataframe
dat_prop <- dat_fuse_prop %>% 
  bind_cols(dat_FUn_prop) %>% 
  bind_cols(dat_FSp_prop)

# save as rds file
dat_prop %>% 
  write_rds(here("data", 
                 "proportional_coverage.rds"), 
            compress = "gz")


# dat_prop <- read_rds(here("data",
#                           "proportional_coverage.rds"))



# visualise ---------------------------------------------------------------

plot_prop <- dat_prop %>% 
  pivot_longer(-cell_res, 
               values_to = "prop_cov", 
               names_to = "metric") %>% 
  mutate(cell_res = ordered(cell_res, 
                            levels = c("1x1", "5x5", 
                                       "10x10", "50x50", 
                                       "100x100"))) %>% 
  ggplot(aes(cell_res, prop_cov, 
             fill = metric)) +
  geom_boxplot(colour = "grey30", 
               outlier.shape = 21, 
               outlier.fill = "grey90", 
               outlier.size = 1.5, 
               alpha = 0.6, 
               key_glyph = "rect") +
  scale_fill_manual(values = c(colour_mint, 
                               colour_purple, 
                               colour_yellow), 
                    labels = c("Specialisation",
                               "Uniqueness", 
                               "FUSE")) +
  scale_y_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels = function(x) paste(x * 100, "%")) +
  labs(y = "Global-local agreement", 
       x = "Cell resolution",
       fill = "Functional metric") +
  theme(legend.position = c(0.2, 0.8),
        legend.key.size = unit(3, "mm"))

# save plot
ggsave(plot_prop, 
       filename = here("figures",
                       "main",
                       "scaling_plot.pdf"),
       width = 183, height = 100,
       units = "mm",
       bg = "white")

