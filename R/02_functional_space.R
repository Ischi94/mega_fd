library(here)
library(tidyverse)
library(ade4)

# get presence absence data
dat_pres_abs <- read_rds(here("data", 
                              "presabs_05_res.rds"))

# load trait data
dat_traits <- read_csv(here("data", "megafauna_traits.csv"))

  

# distance trait matrix -----------------------------------------------------------------------


# prepare for trait analysis
dat_trait_clean <- dat_traits %>% 
  arrange(species) %>% 
  select(species, terrestriality:group_behaviour) %>% 
  mutate(species = str_to_lower(species), 
         species = str_replace_all(species," ", "_")) %>% 
  filter(species %in% colnames(dat_pres_abs)[3:ncol(dat_pres_abs)]) %>% 
  mutate(across(where(is.character), as.factor))

# nominal traits
tabN <- dat_trait_clean %>% 
  column_to_rownames("species") %>% 
  select(-weight)

# quantitative traits
tabQ <- dat_trait_clean %>% 
  column_to_rownames("species") %>% 
  select(weight)

# combining the traits
list_traits <- ktab.list.df(list(tabN, tabQ))

# calculate distance trait matrix
distance_trait_matrix <- dist.ktab(list_traits,
                                 c("N", "Q"),
                                 c("scaledBYrange"))

# # save distance trait matrix
# distance_trait_matrix %>% 
#   write_rds(here("data", 
#                  "distance_trait_matrix.rds"), 
#             compress = "gz")


# select dimensions ---------------------------------------------------------------------------


# compute functional space
p <-  dudi.pco(quasieuclid(distance_trait_matrix),
                       scannf = FALSE,
                       nf = 9) 

# extract scores
scores <- p$li  

# original distance matrix as matrix
D0 <- as.matrix(distance_trait_matrix)

# function to compute the reconstructed distances for k axes
reconstruct_dist <- function(scores, k) {
  S <- as.matrix(scores[, 1:k, drop = FALSE])
  Dc <- dist(S)
  as.matrix(Dc)
}

# compute msd for each possible number of axes
max_k <- ncol(scores)
msd_vec <- numeric(max_k)

for (k in 1:max_k) {
  Dk <- reconstruct_dist(scores, k)
  msd_vec[k] <- mean((D0 - Dk)^2)
}

# visualise elbow inflextion point
plot(1:max_k, msd_vec, type = "b", pch = 19,
     xlab = "number of pcoa axes",
     ylab = "mean squared deviation",
     main = "elbow selection for pcoa dimensions")

# compute functional space with appropriate dimensions
dat_space <-  dudi.pco(quasieuclid(distance_trait_matrix),
               scannf = FALSE,
               nf = 5) %>% 
  pluck("li")

# # save functional space
# dat_space %>% 
#   write_rds(here("data", "functional_space.rds"))
  