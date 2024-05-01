library(here)
library(tidyverse)
library(sf)

# get functions
source(here("R", 
            "functions.R"))


# load data ---------------------------------------------------------------

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


# world map
world_map_sf <- st_read(here("data", 
                             "worldmap",
                             "ne_10m_land.shp"))

# provinces of the world
MEOW_sf <- read_sf(here("data", 
                        "meow", 
                        "meow_ecos.shp"))



# fuse per realm ----------------------------------------------------------

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
  
  list_fuse_realm[[i]] <- spp_per_grid[matrix_realm[[i]]] %>%
    flatten_chr() %>% 
    enframe(value = "species",
            name = NULL) %>% 
    distinct(species) %>% 
    left_join(dat_fuse) %>% 
    add_column(realm = realm_char[i]) 
    
}

# add global rank and plot function
plot_comparison <- function(metric, plot_label) {
  
  dat_local <- list_fuse_realm %>% 
    map(~ .x %>% 
          select(species, !!enquo(metric), realm) %>% 
          arrange(desc(!!enquo(metric))) %>% 
          slice_head(n = 5) %>% 
          rownames_to_column("local_rank") %>% 
          select(species, local_rank, realm)) 
  
  
  dat_merged <- dat_local %>% 
    map(~ left_join(.x, dat_fuse %>%
                      arrange(desc(!!enquo(metric))) %>%
                      mutate(global_rank = row_number()) %>%
                      select(species, global_rank))) %>%
    bind_rows() %>%
    mutate(local_rank = as.integer(local_rank)) %>% 
    pivot_longer(cols = -c(species, realm),
                 names_to = "scale",
                 values_to = "rank") %>%
    mutate(scale = as.integer(as.factor(scale))) 
  
  dat_merged %>% 
    ggplot(aes(x = scale,
               y = rank)) +
    geom_line(aes(group = species), 
              colour = colour_grey) +
    geom_text(aes(x = 0.7,
                  y = rank,
                  label = rank),
              data = dat_merged %>% 
                filter(scale == 1),
              size = 8/.pt) +
    geom_text(aes(label = species),
              data = dat_local %>%
                bind_rows() %>%
                add_column(scale = 1) %>%
                rename(rank = local_rank) %>% 
                mutate(scale = as.integer(scale), 
                       rank = as.integer(rank)),
              size = 7/.pt,
              hjust = 1,
              position = position_nudge(x = -1.1)) +
    labs(y = plot_label,
         x = NULL) +
    scale_x_continuous(trans = "reverse",
                       breaks = c(1, 4),
                       labels = c("Global", "Local")) +
    scale_y_reverse(breaks = NULL) +
    coord_cartesian(xlim = c(6, 0.6)) +
    facet_wrap(~ realm, 
               scales = "free_y") +
    theme_minimal() +
    theme(legend.position = "none",
          panel.grid = element_blank())
  
}

plot_uniq <- plot_comparison(FUn_std, "Functional Uniqueness")

plot_fuse <- plot_comparison(FUSE, "FUSE rank")

plot_speci <- plot_comparison(FSp_std, "Functional Specialisation")


# set up saving function
save_plot <- function(plot_name) {
  
  ggsave(plot_name, 
         filename = here("figures",
                         paste0(plot_name$labels$y, ".png")),
         width = 1920, height = 1080*2,
         units = "px",
         bg = "white", device = ragg::agg_png)
  
}

# save plot
list(plot_uniq, plot_fuse, 
     plot_speci) %>% 
  map(save_plot)

# save plot
ggsave(plot_fuse,
       filename = here("figures",
                       "main",
                       "rank_plot.pdf"),
       width = 183, height = 200,
       units = "mm",
       bg = "white")

