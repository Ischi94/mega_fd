library(tidyverse) 
library(here)

# load data ---------------------------------------------------------------

# local metrics at 5x5 degree res
dat_local <- read_rds(here("data",
                           "local_metrics_5_res.rds"))

# summarise metrics
dat_fun <- dat_local %>%
  rename_with(~ str_remove(.x, "_?local$")) %>% 
  group_by(species) %>% 
  summarise(across(c(FUSE, FUn, FSp, FDi), 
                   list(sd = ~sd(.x, na.rm = TRUE),
                        min = ~min(.x, na.rm = TRUE),
                        mean = ~mean(.x, na.rm = TRUE),
                        max = ~max(.x, na.rm = TRUE)))) %>% 
  mutate(Species = str_to_sentence(species), 
         Species = str_replace_all(Species, "_", " "), 
         across(where(is.double), ~round(.x, 2)), 
         .before = 1) %>% 
  select(-species)

# save as csv
dat_fun %>% 
  write_csv(here("data", "summary_table.csv"))
  
# calculate how many species become the functionally most important in a given grid cell 
# for each metric
map_dbl(c("FUn", "FSp", "FDi"), 
        function(m) {
  
          dat_local %>%
            rename_with( ~ str_remove(.x, "_?local$")) %>%
            group_by(id) %>%
            filter(.data[[m]] == max(.data[[m]], na.rm = TRUE)) %>%
            ungroup() %>%
            distinct(species) %>%
            nrow()
          
          }) %>% 
  set_names(c("FUn", "FSp", "FDi"))
  

