library(tidyverse) 
library(here)
library(flextable)
library(officer)

# load data ---------------------------------------------------------------

# local metrics at 5x5 degree res
dat_local <- read_rds(here("data",
                           "local_metrics_5_res.rds"))

# only maximun FUn
dat_fun <- dat_local %>% 
  select(Species = species, 
         Latitude = latitude, Longitude = longitude,  
         FUn = FUn_std_local) %>%
  filter(FUn == 1) %>% 
  arrange(Species) %>% 
  mutate(Species = str_to_sentence(Species), 
         Species = str_replace_all(Species, "_", " "), 
         FUn = round(FUn, 2)) 

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
