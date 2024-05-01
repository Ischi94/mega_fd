library(here)
library(tidyverse)
library(sf)
library(patchwork)

# get functions
source(here("R", 
            "functions.R"))

# load data ---------------------------------------------------------------

# load metrics calculated in 1
dat_metrics <- read_rds(here("data",
                             "functional_metrics.rds"))


# world map
world_map_sf <- st_read(here("data", 
                             "worldmap",
                             "ne_10m_land.shp"))



# calculate hot-spots ------------------------------------------------------

dat_hot <- dat_metrics %>% 
  mutate(spR_hot = spR > quantile(spR,
                                  probs =  0.975), 
         FRic_hot = FRic > quantile(FRic,
                                    probs =  0.975), 
         uniq_hot = uniq > quantile(uniq,
                                    probs =  0.975),
         overlap_hot = if_else(spR_hot + FRic_hot == 2, TRUE, FALSE)) %>% 
  pivot_longer(cols = c(contains("hot")), 
               names_to = "hotspot")

plot_hot <- dat_hot %>% 
  filter(value == TRUE) %>% 
  ggplot() +
  geom_sf(data = world_map_sf, col = NA, fill = "grey", size = 0.1) +
  geom_raster(aes(x = longitude_x,
                  y = latitude_y, 
                  fill = hotspot)) +
  scale_fill_manual(values = c(colour_mint, colour_coral, 
                               colour_yellow, colour_purple), 
                    name = "Hotspots", 
                    labels = c("FRic", "Overlap", "SR", "FUn")) +
  labs(x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(-180, 180),  ylim = c(-89, 85), expand = FALSE) +
  scale_x_continuous(breaks = seq(-180, 180, 30), 
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-90, 90, 30), 
                     expand = c(0, 0)) +
  theme(legend.position = "top", 
        legend.key.size = unit(2, "mm"))



# venn diagram  -----------------------------------------------------------

plot_venn <- dat_hot %>% 
  filter(hotspot != "overlap_hot") %>% 
  filter(value == TRUE) %>% 
  mutate(across(c(longitude_x, latitude_y), abs)) %>% 
  ggplot(aes(longitude_x, latitude_y, 
             colour = hotspot, 
             fill = hotspot)) +
  stat_ellipse(geom = "polygon", 
               alpha = 0.5) +
  labs(x = "Absolute Longitude", 
       y = "Absolute Latitude", 
       colour = "Hotspot", 
       fill = "Hotspot") +
  scale_color_manual(values = c(colour_mint, 
                                colour_yellow, colour_purple), 
                     labels = c("FRic", "SR", "FUn")) +
  scale_fill_manual(values = c(colour_mint,
                               colour_yellow, colour_purple), 
                    labels = c("FRic", "SR", "FUn")) +
  theme(legend.position = c(0.85, 0.6))


# patch together
plot_final <- plot_hot/ plot_venn

# save
ggsave(plot_final, 
       filename = here("figures",
                       "main",
                       "2_hotspots.pdf"),
       width = 183, height = 100*2,
       units = "mm",
       bg = "white")
