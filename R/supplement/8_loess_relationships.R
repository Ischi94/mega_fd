library(here)
library(tidyverse)
library(patchwork)

# get functions
source(here("R", 
            "functions.R"))

# read data ---------------------------------------------------------------

dat_metrics <- read_rds(here("data", 
                             "functional_metrics.rds"))

plot_fric <- dat_metrics %>% 
  ggplot(aes(spR, FRic)) +
  geom_point(shape = 21, 
             fill = colour_grey, 
             colour = "grey20", 
             alpha = 0.7) +
  geom_smooth(se = FALSE, 
              colour = colour_coral) +
  labs(x = "Species Richness", 
       y = "Functional Richness") +
  theme_minimal()

plot_uniq <- dat_metrics %>% 
  ggplot(aes(spR, uniq)) +
  geom_point(shape = 21, 
             fill = colour_grey, 
             colour = "grey20", 
             alpha = 0.7) +
  geom_smooth(se = FALSE, 
              colour = colour_coral) +
  labs(x = "Species Richness", 
       y = "Functional Uniqueness") +
  theme_minimal()

plot_special <- dat_metrics %>% 
  ggplot(aes(spR, special )) +
  geom_point(shape = 21, 
             fill = colour_grey, 
             colour = "grey20", 
             alpha = 0.7) +
  geom_smooth(se = FALSE, 
              colour = colour_coral) +
  labs(x = "Species Richness", 
       y = "Functional Specialisation") +
  theme_minimal()

# patch together and save
plot_final <- plot_fric/ plot_uniq / plot_special


# save plot
ggsave(plot_final,
       filename = here("figures",
                       "main",
                       "4_loess_plot.pdf"),
       width = 183, height = 200,
       units = "mm",
       bg = "white")
