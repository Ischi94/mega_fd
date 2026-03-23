library(here)
library(raster)
library(tidyverse)
library(sf)
library(patchwork)

# get functions
source(here("R", 
            "functions.R"))


# load data ---------------------------------------------------------------


# global metrics
dat_global <- read_rds(here("data", "global_metrics.rds"))

# local metrics
dat_local <- read_rds(here("data", "local_metrics_1_res.rds"))

# per realm metrics
dat_realm <- read_rds(here("data", "realm_metrics.rds")) 

# provinces of the world
MEOW_sf <- read_sf(here("data", 
                        "meow", 
                        "meow_ecos.shp"))

# world map
world_map_sf <- st_read(here("data", 
                             "worldmap",
                             "ne_10m_land.shp"))
            
# per province agreement --------------------------------------------------


get_agreement <- function(local_metr, global_metr,
                          title_metr, legend_pos) {
  
  dat_realm %>%
    left_join(dat_global) %>% 
    group_by(realm) %>% 
    arrange(desc(!!enquo(local_metr))) %>%
    mutate(local_rank = 1:n()) %>% 
    arrange(desc(!!enquo(global_metr))) %>% 
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
           ))  %>%
    ggplot(aes(!!enquo(local_metr), !!enquo(global_metr))) +
    geom_abline(intercept = 0,
                slope = 1,
                linetype = "dotted",
                colour = "grey20") +
    geom_point(aes(fill = interaction(local_rank, global_rank),
                   size = interaction(local_rank, global_rank)),
               colour = "grey60",
               shape = 21,
               stroke = 0.3) +
    scale_fill_manual(values = c("grey20", "#2F899D",
                                 colour_coral, "white"),
                      breaks = c("local.global", "global.global", "local.not", "bnothing.not"),
                      labels = c("Shared", "Global", "Local", " ")) +
    scale_size_manual(values = c(2.2, 1.4, 1.4, 1.4),
                      breaks = c("local.global", "global.global", "local.not", "bnothing.not"),
                      guide = "none") +
    labs(y = "Global",
         x = "Provincial",
         fill = "10 highest species",
         title = title_metr) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
    theme(legend.position = legend_pos) +
    guides(fill = guide_legend(
      override.aes = list(colour = "white",
                          size = 2.5)))
  
}


# specialisation
plot_FSp <- get_agreement(FSp_realm, FSp, "FSp", "none")

# distinctiveness
plot_FDi <- get_agreement(FDi_realm, FDi, "FDi", "top")

# uniqueness
plot_FUn <- get_agreement(FUn_realm, FUn, "FUn", "none")

# patch together
plot_province <- plot_FSp + plot_FDi + plot_FUn


# add rank-rank plot ------------------------------------------------------


# global to province
dat_glob_realm <- dat_global %>%
  arrange(desc(FUSE)) %>%
  filter(species %in% c("chelonia_mydas", "probarbus_jullieni", 
                        "arctocephalus_galapagoensis", "enhydra_lutris", 
                        "ursus_maritimus")) %>% 
  mutate(global_rank = row_number()) %>%
  select(species, global_rank) %>% 
  left_join(dat_realm %>% 
              group_by(realm) %>% 
              select(species, FUSE_realm, realm) %>% 
              arrange(realm, desc(FUSE_realm)) %>% 
              mutate("local_rank" = row_number()) %>%  
              select(species, local_rank, realm)) %>% 
  pivot_longer(cols = -c(species, realm),
               names_to = "scale",
               values_to = "rank") %>%
  mutate(scale = as.integer(as.factor(scale))) %>% 
  mutate(rank_log = log(rank))


# add local and plot using mode
plot_mode <- dat_glob_realm %>%
  bind_rows(dat_local %>% 
              group_by(id) %>%
              drop_na(FUSE_local) %>%
              arrange(desc(FUSE_local)) %>%
              slice_head(n = 5) %>%
              mutate(fuse_rank = row_number()) %>%
              select(species, fuse_rank) %>% 
              group_by(species) %>% 
              summarise(local_rank = which.max(tabulate(fuse_rank))) %>% 
              ungroup() %>% 
              mutate(scale = 3, rank = local_rank, rank_log = log(rank)) %>% 
              select(species, scale, rank, rank_log) %>% 
              right_join(dat_glob_realm %>% 
                           distinct(species, realm))) %>%
  mutate(species = str_replace_all(species, "_", " "), 
         species = str_to_sentence(species)) %>% 
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
               filter(rank %in% c(1:8, 13, 16, 23, 29, 34, 41, 49))) +
  geom_text(aes(label = species, 
                colour = species),
            data = dat_glob_realm %>%
              filter(scale == 1) %>% 
              distinct(species, rank_log, scale) %>% 
              mutate(species = str_replace_all(species, "_", " "), 
                     species = str_to_sentence(species)),
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
        panel.grid = element_blank(), 
        plot.title = element_text(colour = "grey20", size = 10, 
                                  face = "bold"))



# global-local agreement ----------------------------------------------------------------------


# calculate proportion of agreement in FUSE rank between local cells and global rank
dat_prop <- dat_local %>% 
  select(id, species, FUSE_local) %>% 
  drop_na(FUSE_local) %>% 
  left_join(dat_global %>% 
              select(species, FUSE) %>% 
              arrange(desc(FUSE)) %>% 
              rowid_to_column("global_rank")) %>% 
  group_by(id) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(fuse_prop = map_dbl(data, 
                             ~ .x %>% 
                               arrange(desc(FUSE_local)) %>% 
                               rowid_to_column("local_rank") %>% 
                               filter(between(global_rank, 1, 10)) %>% 
                               nrow(.)/10))

# create map
plot_map <- dat_prop %>%
  select(-data) %>% 
  left_join(dat_local %>% 
              select(id, latitude, longitude)) %>% 
  ggplot() +
  geom_raster(aes(x = longitude,
                  y = latitude,
                  fill = fuse_prop)) +
  scale_fill_gradient2(low = colour_mint,
                       mid = colour_purple,
                       midpoint = 0.1,
                       high = colour_yellow,
                       na.value = "white",
                       breaks = c(0, 0.1, 0.2, 0.3),
                       labels = c("0%", "10%", "20%", "30%")) +
  geom_sf(data = world_map_sf, col = NA, fill = "white", size = 0.1) +
  labs(x = "Longitude", y = "Latitude", 
       fill = "Global-Local\nFUSE\nAgreement") +
  scale_x_continuous(breaks = seq(-180, 180, 60),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-90, 90, 30),
                     expand = c(0, 0)) +
  coord_sf(ylim = c(-87, 87)) +
  theme(legend.position = "bottom", 
        plot.title = element_text(colour = "grey20", size = 10, 
                                  face = "bold", 
                                  hjust = -0.13))



# patch together 
plot_final <- free(plot_province) /
  (plot_map | plot_mode) +
  plot_annotation(tag_levels = "A")


# save plot
ggsave(plot_final,
       filename = here("figures",
                       "main",
                       "fig_2.png"),
       width = 300, height = 150,
       units = "mm",
       bg = "white")


