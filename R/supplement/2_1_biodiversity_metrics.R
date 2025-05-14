library(here)
library(tidyverse)
library(sf)
library(ade4)
library(patchwork)

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




# calculate metrices ------------------------------------------------------


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


# computes species richness, functional richness, and
# functional diversity per grid cell
fd_metrics <- spp_per_grid %>% 
  map(~get_FV_Sp(ax = c(1:4),
                 pcoa = pcoa,
                 Selected_sp = .x), 
      .progress = TRUE)
  

# species richness
vec_spR <- fd_metrics %>%
  map_dbl(~pluck(.x, "data") %>% 
            pull("RS"))

# functional richness
vec_FRic <- fd_metrics %>%
  map_dbl( ~ pluck(.x, "data") %>%
             pull("FRic"))

# functional diversity
vec_FDiv <- fd_metrics %>%
  map_dbl( ~ pluck(.x, "data") %>%
             pull("FDiv"))



# FUSE --------------------------------------------------------------------

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
  filter(species %in% rownames(pcoa$li)) 

# arrange correctly
ge <- ge[match(rownames(pcoa$li), ge$species), ] %>% 
  pull(ge) 

# Compute FUSE metrics 
dat_fuse <- get_FUSE(Mat_dist = distance_trait_matrix,
                Coords = pcoa$li,
                GE = ge) %>% 
  as_tibble(rownames = "species") 

# save 
dat_fuse %>% 
  write_rds(here("data", 
                 "fuse_metrics_per_spp.rds"))


# compute mean FUSE per grid cell
dat_mean_fuse <- spp_per_grid %>% 
  map(~ dat_fuse %>% 
        select(species, FUSE,
               FUn_std, FSp_std) %>% 
        filter(species %in% .x) %>% 
        summarise(mean_fuse = mean(FUSE, na.rm = TRUE), 
                  mean_uniq = mean(FUn_std, na.rm = TRUE), 
                  mean_spec = mean(FSp_std, na.rm = TRUE)), 
      .progress = TRUE)

# extract

# fuse
vec_fuse <- map_dbl(dat_mean_fuse, 
                    ~ pull(.x, "mean_fuse"), 
                    .progress = TRUE)

# uniqueness
vec_uniq <- map_dbl(dat_mean_fuse,
                    ~ pull(.x, "mean_uniq"), 
                    .progress = TRUE)

# specialisation
vec_spec <- map_dbl(dat_mean_fuse,
                    ~ pull(.x, "mean_spec"), 
                    .progress = TRUE)

# 25% top fuse species
dat_fuse_top <- spp_per_grid %>%
  map( ~ enframe(.x,
                 value = "species")) %>% 
  map(~ .x %>% 
        filter(species %in% (dat_fuse %>% # define top fuse species
                               filter(FUSE  > quantile(FUSE ,
                                                       probs =  0.75,
                                                       na.rm = TRUE)) %>% 
                               pull(species))), 
      .progress = TRUE)


# combine ----------------------------------------------------------------

# combine in dataframe with coordinates
dat_metrics <- tibble(spR = vec_spR,
                      FRic = vec_FRic,
                      fuse = vec_fuse, 
                      uniq = vec_uniq, 
                      special = vec_spec) %>% 
  add_column(grid = rownames(dat_presabs)) %>% 
  separate_wider_delim(grid, 
                       names = c("longitude_x",
                                 "latitude_y"),
                       delim = "_") %>%
  mutate(across(c(longitude_x, latitude_y), as.double), 
         FRic = FRic*100) %>% 
  add_column(fuse_top = map_dbl(dat_fuse_top, nrow)) 

# save
dat_metrics %>% 
  write_rds(here("data", 
                 "functional_metrics.rds"), 
            compress = "gz")

# dat_metrics <- read_rds(here("data",
#                              "functional_metrics.rds"))

# visualise ---------------------------------------------------------------

# plot function
plot_metric <- function(metric_col, 
                        label) {
  
  dat_metrics %>%
    ggplot() +
    geom_raster(aes(x = longitude_x,
                    y = latitude_y, 
                    fill = !!enquo(metric_col))) +
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


# functional richness
plot_FRic <- plot_metric(FRic, "FRic (%)") 

# species richness
plot_SR <- plot_metric(spR, "SR") 

# uniqueness
plot_uniq <- plot_metric(uniq, "FUn") 

# specialisation
plot_spec <- plot_metric(special, "FSp")


# patch together
plot_final <- (plot_SR + plot_uniq) /
  (plot_FRic + plot_spec)


# save
ggsave(plot_final, 
       filename = here("figures",
                       "main",
                       "1_metrics.pdf"),
       width = 183*2, height = 100*2,
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

