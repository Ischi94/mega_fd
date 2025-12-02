library(tidyverse) 
library(here)
library(flextable)
library(officer)

# load data ---------------------------------------------------------------

# local metrics at 5x5 degree res
dat_local <- read_rds(here("data",
                           "local_metrics_5_res.rds"))

# summarise metrics
dat_fun <- dat_local %>%
  rename_with(~ str_remove(.x, "_?local$")) %>% 
  rename_with(~ str_remove(.x, "_?std$")) %>% 
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
  

# generate nice looking flextable
table_fun <- dat_fun %>%
  flextable() %>% 
  theme_vanilla() %>% 
  merge_v(j = "Species") %>% 
  italic(j = "Species")

# create word document 

# open docx-file and add flextable
my_doc <- read_docx() %>% 
  body_add_flextable(table_fun) 

# convert to word file/ add input to empty docx
print(my_doc, target = here("data",
                            "supplementary_table.docx"))
