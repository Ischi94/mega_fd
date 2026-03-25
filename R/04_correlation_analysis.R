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
dat_local <- read_rds(here("data", "local_metrics_5_res.rds"))

# per realm metrics
dat_realm <- read_rds(here("data", "realm_metrics.rds")) 

# provinces of the world
MEOW_sf <- read_sf(here("data", 
                        "meow", 
                        "meow_ecos.shp"))

# rank-rank correlation ---------------------------------------------------

# per province rank-correlation
dat_cor_prov <- dat_realm %>%
  left_join(dat_global) %>% 
  group_by(realm) %>%
  summarise(FUn = cor(FUn_realm, FUn, use = "complete.obs", method = "kendall"), 
            FSp = cor(FSp_realm, FSp, use = "complete.obs", method = "kendall"), 
            FUSE = cor(FUSE_realm, FUSE, use = "complete.obs", method = "kendall"), 
            FDi = cor(FDi_realm, FDi, use = "complete.obs", method = "kendall"),
            nSp = n())  %>% 
  bind_rows(dat_realm %>% 
              left_join(dat_global) %>% 
              summarise(FUn = cor(FUn_realm, FUn, use = "complete.obs", method = "kendall"), 
                        FSp = cor(FSp_realm, FSp, use = "complete.obs", method = "kendall"), 
                        FUSE = cor(FUSE_realm, FUSE, use = "complete.obs", method = "kendall"), 
                        FDi = cor(FDi_realm, FDi, use = "complete.obs", method = "kendall"),
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
  summarise(FUn = cor(FUn_local, FUn, use = "complete.obs", method = "kendall"),
            FSp = cor(FSp_local, FSp, use = "complete.obs", method = "kendall"),
            FUSE = cor(FUSE_local, FUSE, use = "complete.obs", method = "kendall"),
            FDi = cor(FDi_local, FDi, use = "complete.obs", method = "kendall"),
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



# intravariable correlation -------------------------------------------------------------------

# local 
dat_intcor_loc <- dat_local %>% 
  summarise(FUSE_FUn = cor(FUSE_local, FUn_local, use = "complete.obs", method = "kendall"), 
            FUSE_FSp = cor(FUSE_local, FSp_local, use = "complete.obs", method = "kendall"), 
            FUSE_FDi = cor(FUSE_local, FDi_local, use = "complete.obs", method = "kendall"), 
            FUn_FSp = cor(FUn_local, FSp_local, use = "complete.obs", method = "kendall"), 
            FUn_FDi = cor(FUn_local, FDi_local, use = "complete.obs", method = "kendall"), 
            FSp_FDi = cor(FSp_local, FDi_local, use = "complete.obs", method = "kendall"))


# provincial 
dat_intcor_prov <- dat_realm %>%
  summarise(FUSE_FUn = cor(FUSE_realm, FUn_realm, use = "complete.obs", method = "kendall"), 
            FUSE_FSp = cor(FUSE_realm, FSp_realm, use = "complete.obs", method = "kendall"), 
            FUSE_FDi = cor(FUSE_realm, FDi_realm, use = "complete.obs", method = "kendall"), 
            FUn_FSp = cor(FUn_realm, FSp_realm, use = "complete.obs", method = "kendall"), 
            FUn_FDi = cor(FUn_realm, FDi_realm, use = "complete.obs", method = "kendall"), 
            FSp_FDi = cor(FSp_realm, FDi_realm, use = "complete.obs", method = "kendall"))

# global
dat_intcor_glob <- dat_global %>% 
  summarise(FUSE_FUn = cor(FUSE, FUn, use = "complete.obs", method = "kendall"), 
            FUSE_FSp = cor(FUSE, FSp, use = "complete.obs", method = "kendall"), 
            FUSE_FDi = cor(FUSE, FDi, use = "complete.obs", method = "kendall"), 
            FUn_FSp = cor(FUn, FSp, use = "complete.obs", method = "kendall"), 
            FUn_FDi = cor(FUn, FDi, use = "complete.obs", method = "kendall"), 
            FSp_FDi = cor(FSp, FDi, use = "complete.obs", method = "kendall"))

# merge 
dat_intcor <- dat_intcor_loc %>% 
  add_column(scale = "local") %>% 
  bind_rows(dat_intcor_prov %>% 
              add_column(scale = "province")) %>% 
  bind_rows(dat_intcor_glob %>% 
              add_column(scale = "global")) 
# save
dat_intcor %>% 
  write_rds(here("data", "intravariable_correlation.rds"))


# dat_intcor <- read_rds(here("data", "intravariable_correlation.rds"))

# plot
plot_intcor <- dat_intcor %>% 
  pivot_longer(cols = -scale) %>% 
  separate_wider_delim(name, delim = "_", names = c("var_a", "var_b")) %>% 
  full_join(tibble(var_a = c("FUSE", "FUn", "FSp", "FDi"), 
                   var_b = c("FUSE", "FUn", "FSp", "FDi")) %>% 
              add_column(value = 1) %>% 
              expand_grid(scale = c("local", "province", "global"))) %>% 
  mutate(scale = str_to_sentence(scale), 
         scale = ordered(scale, levels = c("Local", "Province", "Global"))) %>% 
  ggplot(aes(var_a, var_b)) +
  geom_tile(aes(fill = value), color = 'black') +
  geom_text(aes(label = round(value, 2))) +
  scale_fill_gradient2(low = colour_mint, mid = colour_yellow, high = colour_coral, 
                       midpoint = 0.5) +
  labs(x = element_blank(),
       y = element_blank(),
       fill = "Correlation") +
  facet_wrap(~scale)

# save plot
ggsave(plot_intcor, 
       filename = here("figures",
                       "main", 
                       "suppl_1.pdf"),
       width = 183, height = 100,
       units = "mm",
       bg = "white")

