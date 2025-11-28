library(here)
library(sf)
library(tidyverse)
library(patchwork)

# get functions
source(here("R", 
            "functions.R"))

# load data ---------------------------------------------------------------


# global metrics
dat_global <- read_rds(here("data", "global_metrics.rds"))

# local metrics
dat_local <- read_rds(here("data", "local_metrics.rds"))

# per realm metrics
dat_realm <- read_rds(here("data", "realm_metrics.rds")) %>% 
  rename_with(cols = contains("local"), 
              .f = ~ str_replace_all(.x, "local", "realm"))

# provinces of the world
MEOW_sf <- read_sf(here("data", 
                        "meow", 
                        "meow_ecos.shp"))

# rank-rank correlation ---------------------------------------------------

# per province rank-correlation
dat_cor_prov <- dat_realm %>%
  left_join(dat_global) %>% 
  group_by(realm) %>%
  summarise(FUn = cor(FUn_realm, FUn_std, use = "complete.obs", method = "kendall"), 
            FSp = cor(FSp_realm, FSp_std, use = "complete.obs", method = "kendall"), 
            FUSE = cor(FUSE_realm, FUSE, use = "complete.obs", method = "kendall"), 
            FDi = cor(FDi_realm, FDi_std, use = "complete.obs", method = "kendall"),
            nSp = n())  %>% 
  bind_rows(dat_realm %>% 
              left_join(dat_global) %>% 
              summarise(FUn = cor(FUn_realm, FUn_std, use = "complete.obs", method = "kendall"), 
                        FSp = cor(FSp_realm, FSp_std, use = "complete.obs", method = "kendall"), 
                        FUSE = cor(FUSE_realm, FUSE, use = "complete.obs", method = "kendall"), 
                        FDi = cor(FDi_realm, FDi_std, use = "complete.obs", method = "kendall"),
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
  mutate(shape_id = if_else(realm == "Overall", 23, 21), 
         metric = ordered(metric, levels = c("FUn", "FSp", 
                                             "FUSE", "FDi"))) 
# plot
plot_cor_prov <- dat_cor_prov %>%
  ggplot(aes(correl, realm, 
             fill = metric, 
             size = nSp, 
             shape = shape_id)) +
  geom_point(colour = "grey20") +
  scale_fill_manual(values = c(colour_yellow,
                               colour_coral, 
                               colour_mint,
                               colour_purple), 
                    breaks = c("FDi", "FUSE", 
                               "FSp", "FUn")) +
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
                                  hjust = -5))


# per local cell rank-correlation
plot_cor_loc <- dat_local %>%
  left_join(dat_global) %>% 
  group_by(id) %>% 
  summarise(FUn = cor(FUn_std_local, FUn_std, use = "complete.obs", method = "kendall"),
            FSp = cor(FSp_std_local, FSp_std, use = "complete.obs", method = "kendall"),
            FUSE = cor(FUSE_local, FUSE, use = "complete.obs", method = "kendall"),
            FDi = cor(FDi_std_local, FDi_std, use = "complete.obs", method = "kendall"),
            nSp = n()) %>% 
  left_join(dat_local %>% 
              distinct(id, latitude)) %>% 
  pivot_longer(cols = -c(id, latitude, nSp), 
               names_to = "metric", 
               values_to = "correl") %>% 
  mutate(metric = ordered(metric, levels = c("FUn", "FSp", 
                                             "FUSE", "FDi"))) %>% 
  ggplot(aes(correl, latitude, 
             fill = metric, 
             size = nSp)) +
  geom_point(shape = 21, 
             alpha = 0.1, 
             colour = "grey20") +
  scale_fill_manual(values = c(colour_yellow,
                               colour_coral, 
                               colour_mint,
                               colour_purple), 
                    breaks = c("FDi", "FUSE", 
                               "FSp", "FUn")) +
  scale_y_continuous(labels = function(x) str_c(x, '°')) +
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
                                  hjust = -0.7)) +
  guides(fill = guide_legend(override.aes = list(alpha = 1, 
                                                 size = 2), 
                             reverse = TRUE))


# merge together
plot_final <- plot_cor_prov + plot_cor_loc + 
  plot_annotation(tag_levels = "A")

# save plot
ggsave(plot_final, 
       filename = here("figures",
                       "main", 
                       "fig_1.png"),
       width = 183, height = 100,
       units = "mm",
       bg = "white")
