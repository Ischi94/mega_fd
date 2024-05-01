library(here)
library(tidyverse)
library(sf)
library(ade4)
library(patchwork)

# get functions
source(here("R", 
            "functions.R"))

# load data ---------------------------------------------------------------


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

# distance-trait matrix
load(here("data", 
          "distance_trait_matrix.RData"))

# extract functional space axes
dat_pcoa <- dudi.pco(quasieuclid(distance_trait_matrix),
                 scannf = FALSE,
                 nf = 4) %>% 
  pluck("li") %>% 
  # prepare pcoa for joining
  as_tibble(rownames = "species")


# provinces of the world
MEOW_sf <- read_sf(here("data", 
                        "meow", 
                        "meow_ecos.shp"))

# uniqueness and specialisation comparison
dat_realm <- read_rds(here("data",
                           "global_regional_trends.rds"))

# species pool per realm -------------------------------------------------------

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
list_spp_realm <- vector(mode = "list",
                         length = length(realm_char))

# loop through realms
for (i in 1:length(realm_char)) {
  
  list_spp_realm[[i]] <- spp_per_grid[matrix_realm[[i]]] %>%
    flatten_chr() %>% 
    enframe(value = "species",
            name = NULL) %>% 
    distinct(species) %>% 
    add_column(realm = realm_char[i]) 
}


# assign pcoa to realms ---------------------------------------------------

# get all four axis per realm
dat_spp_realm <- list_spp_realm %>% 
  map(~ left_join(.x, dat_pcoa)) %>% 
  bind_rows()

# calculate overall convex hull for axis 1 and 2
dat_hull_1_ov <- dat_spp_realm %>% 
  slice(chull(A1, A2)) %>% 
  select(-realm)

# calculate convex hulls for axis 1 and 2
dat_hull_1 <- dat_spp_realm %>% 
  group_by(realm) %>% 
  slice(chull(A1, A2)) %>% 
  ungroup()

# calculate overall convex hull for axis 1 and 3
dat_hull_2_ov <- dat_spp_realm %>% 
  slice(chull(A1, A3)) %>% 
  select(-realm)

# calculate convex hulls for axis 1 and 3
dat_hull_2 <- dat_spp_realm %>% 
  group_by(realm) %>% 
  slice(chull(A1, A3)) %>% 
  ungroup()

# calculate overall convex hull for axis 2 and 3
dat_hull_3_ov <- dat_spp_realm %>% 
  slice(chull(A2, A3)) %>% 
  select(-realm)

# calculate convex hulls for axis 2 and 3
dat_hull_3 <- dat_spp_realm %>% 
  group_by(realm) %>% 
  slice(chull(A2, A3)) %>% 
  ungroup()

# visualise specialisation ---------------------------------------------------------------

# set up function to iterate over axis and realm
plot_specialisation <- function(selected_realm,
                                inset =  TRUE) {
  
  dat_prep <- dat_spp_realm %>%
    filter(realm == selected_realm) %>% 
    # center of traitspace
    mutate(A1_local = mean(A1), 
           A2_local = mean(A2), 
           A3_local = mean(A3),
           A1_global = mean(dat_spp_realm$A1), 
           A2_global = mean(dat_spp_realm$A2), 
           A3_global = mean(dat_spp_realm$A3)) %>%
    left_join(dat_realm %>% 
                filter(realm == selected_realm) %>% 
                select(species, specialisation, FSp_std)) %>% 
    # identify 10  most spezialised species
    arrange(desc(specialisation)) %>% 
    mutate(local_rank = 1:n()) %>% 
    arrange(desc(FSp_std)) %>% 
    mutate(global_rank = 1:n()) %>% 
    ungroup() %>% 
    mutate(local_rank = if_else(between(local_rank, 1, 10), 
                                "local", 
                                "not"), 
           global_rank = if_else(between(global_rank, 1, 10), 
                                 "global", 
                                 "non"), 
           ov_rank = case_when(
             local_rank == "local" & global_rank != "global" ~ "local", 
             local_rank == "not" & global_rank == "global" ~ "global", 
             local_rank == "not" & global_rank == "non" ~ "nothing", 
             local_rank == "local" & global_rank == "global" ~ "both"
           )) %>% 
    filter(ov_rank != "nothing") %>% 
    mutate(ov_rank = ordered(ov_rank, levels = c("global", 
                                                 "local", 
                                                 "both")))
  
  
    plot_1 <- dat_prep %>% 
    ggplot(aes(A1, A2)) +
      geom_polygon(fill = NA,
                   colour = "grey20",
                   data = dat_hull_1_ov) +
      geom_polygon(fill = colour_grey, 
                   data = dat_hull_1 %>% 
                     filter(realm == selected_realm)) +
      geom_point(aes(colour = ov_rank)) +
      geom_point(aes(A1_local, A2_local),
                 fill = colour_mint,
                 stroke = 1.5,
                 shape = 21,
                 colour = "white") +
      geom_point(aes(A1_global, A2_global),
                 fill = colour_purple,
                 stroke = 1.5,
                 shape = 21,
                 colour = "white") +
      geom_segment(aes(x = A1_global,
                       y = A2_global,
                       xend = A1,
                       yend = A2), 
                   data = dat_prep %>% 
                     filter(ov_rank %in% c("both", "global")), 
                   colour = colour_purple) +
      geom_segment(aes(x = A1_local,
                       y = A2_local,
                       xend = A1,
                       yend = A2), 
                   data = dat_prep %>% 
                     filter(ov_rank %in% c("both", "local")), 
                   colour = colour_mint) +
      scale_y_continuous(breaks = c(-0.3, 0, 0.3)) +
      scale_x_continuous(breaks = c(-0.3, 0, 0.3)) +
      scale_color_manual(values = c(colour_purple,
                                    colour_mint, 
                                    colour_yellow)) +
      coord_cartesian(xlim = c(-0.5, 0.5), 
                      ylim = c(-0.5, 0.5)) +
      labs(x = "PCoA1", 
           y = "PCoA2", 
           colour = NULL) +
      theme(legend.position = "none")
 
    plot_2 <- ggplot() +
      annotate("point", 
               x = 0.5, 
               y = 0.5, 
               shape = 21, 
               fill = "white", 
               colour = "black", 
               size = 35) +
      annotate("rect", 
               xmin = 0.2,
               xmax = 0.8,
               ymin = 0,
               ymax = 1, 
               colour = colour_grey, 
               fill = colour_grey) +
      annotate("point", 
               x = c(0.5, 1.38, 1.38),
               y = c(0.5, 0.8, 0.4), 
               colour = c("white", colour_purple, colour_mint),
               size = c(4, 2, 2)) +
      annotate("text", 
               x = 2.3, 
               y = c(0.8, 0.4), 
               label = c("Globally most specialized", 
                         "Locally most specialized"), 
               size = 7/.pt, 
               colour = c(colour_purple, colour_mint)) +
      annotate("text", 
               x = c(0.5, -0.5, 0.5), 
               y = c(1.2, 1.5, 0.8), 
               label = c("Local space",
                         "Global space", 
                         "Centroid"),
               size = 7/.pt, 
               colour = c(colour_grey, "black", "white")) +
      coord_cartesian(xlim = c(-1, 3.2), 
                      ylim = c(-0.5, 1.5)) +
      theme_void()
    
    plot_3 <- dat_prep %>% 
      filter(ov_rank != "nothing") %>%
      ggplot(aes(A1, A3)) +
      geom_polygon(fill = NA,
                   colour = "grey20",
                   data = dat_hull_2_ov) +
      geom_polygon(fill = colour_grey, 
                   data = dat_hull_2 %>% 
                     filter(realm == selected_realm)) +
      geom_point(aes(colour = ov_rank)) +
      geom_point(aes(A1_local, A3_local),
                 fill = colour_mint,
                 shape = 21,
                 stroke = 1.5,
                 colour = "white") +
      geom_point(aes(A1_global, A3_global),
                 fill = colour_purple,
                 shape = 21,
                 stroke = 1.5, 
                 colour = "white") +
      geom_segment(aes(x = A1_global,
                       y = A3_global,
                       xend = A1,
                       yend = A3), 
                   data = dat_prep %>% 
                     filter(ov_rank %in% c("both", "global")), 
                   colour = colour_purple) +
      geom_segment(aes(x = A1_local,
                       y = A3_local,
                       xend = A1,
                       yend = A3), 
                   data = dat_prep %>% 
                     filter(ov_rank %in% c("both", "local")), 
                   colour = colour_mint) +
      scale_y_continuous(breaks = c(-0.3, 0, 0.3)) +
      scale_x_continuous(breaks = c(-0.3, 0, 0.3)) +
      scale_color_manual(values = c(colour_purple,
                                    colour_mint, 
                                    colour_yellow)) +
      coord_cartesian(xlim = c(-0.5, 0.5), 
                      ylim = c(-0.5, 0.5)) +
      labs(x = "PCoA1", 
           y = "PCoA3", 
           colour = NULL) +
      theme(legend.position = "none")
    
    
    plot_4 <- dat_prep %>% 
      filter(ov_rank != "nothing") %>%
      ggplot(aes(A2, A3)) +
      geom_polygon(fill = NA,
                   colour = "grey20",
                   data = dat_hull_3_ov) +
      geom_polygon(fill = colour_grey, 
                   data = dat_hull_3 %>% 
                     filter(realm == selected_realm)) +
      geom_point(aes(colour = ov_rank)) +
      geom_point(aes(A2_local, A3_local),
                 fill = colour_mint,
                 shape = 21,
                 stroke = 1.5,
                 colour = "white") +
      geom_point(aes(A2_global, A3_global),
                 fill = colour_purple,
                 shape = 21,
                 stroke = 1.5, 
                 colour = "white") +
      geom_segment(aes(x = A2_global,
                       y = A3_global,
                       xend = A2,
                       yend = A3), 
                   data = dat_prep %>% 
                     filter(ov_rank %in% c("both", "global")), 
                   colour = colour_purple) +
      geom_segment(aes(x = A2_local,
                       y = A3_local,
                       xend = A2,
                       yend = A3), 
                   data = dat_prep %>% 
                     filter(ov_rank %in% c("both", "local")), 
                   colour = colour_mint) +
      scale_y_continuous(breaks = c(-0.3, 0, 0.3)) +
      scale_x_continuous(breaks = c(-0.3, 0, 0.3)) +
      scale_color_manual(values = c(colour_purple,
                                    colour_mint, 
                                    colour_yellow)) +
      coord_cartesian(xlim = c(-0.5, 0.5), 
                      ylim = c(-0.5, 0.5)) +
      labs(x = "PCoA2", 
           y = "PCoA3", 
           colour = NULL) +
      theme(legend.position = "none")
    
    if (inset == TRUE) {
      (plot_1 + plot_2) /
        (plot_3 + plot_4) 
    } else if (inset == FALSE) {
      (plot_1 + plot_spacer()) /
        (plot_3 + plot_4) +
        plot_annotation(tag_levels = "a",
                        title = selected_realm)
      }
  

  
}


plot_spec <- plot_specialisation("Eastern Indo-Pacific")


# plot_specialisation("Southern Ocean",
#                     inset = TRUE)

            

# visualise uniqueness ----------------------------------------------------

# set up function to plot five closest species
plot_uniqueness <- function(selected_realm,
                            selected_species) {
  
  # grid of axis space
  axis_grid <- list(c("A1", "A2"),
                    c("A1", "A3"),
                    c("A2",  "A3"))
  
  # empty list 
  closest_local <- vector("list", 
                          3)
  
  for (i in 1:length(axis_grid)) {
    
    # five closest species of local space
    closest_local[[i]] <- dat_spp_realm %>%
      filter(realm == selected_realm) %>%
      # identify five closest species locally
      st_as_sf(coords = axis_grid[[i]]) %>% 
      st_distance() %>% 
      as_tibble() %>% 
      magrittr::set_colnames(dat_spp_realm%>%
                               filter(realm == selected_realm) %>% 
                               pull(species)) %>% 
      add_column(species = dat_spp_realm %>%
                   filter(realm == selected_realm) %>% 
                   pull(species)) %>% 
      pivot_longer(!species, names_to = 'closest', values_to = 'dist') %>% 
      filter(!is.na(dist), 
             dist > 0) %>% 
      group_by(species) %>% 
      arrange(dist) %>% 
      slice(1:5) %>% 
      mutate(dist_rank = 1:5) %>% 
      filter(species == selected_species) %>% 
      pull(closest)
    
    
  }
  
  
  # identify five closest species globally
  
  # empty list 
  closest_global <- vector("list", 
                           3)
  # iterate over axis
  for (i in 1:length(axis_grid)) {
    
    closest_global[[i]] <- dat_pcoa %>%
      st_as_sf(coords = axis_grid[[i]]) %>%
      st_distance() %>%
      as_tibble() %>%
      magrittr::set_colnames(dat_pcoa %>%
                               pull(species)) %>%
      add_column(species = dat_pcoa %>%
                   pull(species)) %>%
      pivot_longer(!species, names_to = 'closest', values_to = 'dist') %>%
      filter(!is.na(dist),
             dist > 0) %>%
      group_by(species) %>%
      arrange(dist) %>%
      slice(1:5) %>%
      mutate(dist_rank = 1:5) %>%
      filter(species == selected_species) %>%
      pull(closest)
    
  }
  
  
  # visualise A1 and A2
  plot_1 <- dat_spp_realm %>%
    filter(realm == selected_realm) %>%
    ggplot(aes(A1, A2)) +
    geom_polygon(fill = NA,
                 colour = "grey20",
                 data = dat_hull_1_ov) +
    geom_polygon(fill = colour_grey,
                 data = dat_hull_1 %>%
                   filter(realm == selected_realm)) +
    geom_point(data = dat_spp_realm %>%
                 filter(species %in% closest_local[[1]]),
               colour = colour_mint) +
    geom_segment(aes(x = dat_spp_realm %>%
                       filter(species == selected_species) %>%
                       pull(A1) %>%
                       unique(),
                     y = dat_spp_realm %>%
                       filter(species == selected_species) %>%
                       pull(A2) %>%
                       unique(),
                     xend = A1,
                     yend = A2),
                 colour = colour_mint,
                 data = dat_spp_realm %>%
                   filter(species %in% closest_local[[1]])) +
    geom_point(data = dat_spp_realm %>%
                 filter(species %in% closest_global[[1]]),
               colour = colour_purple) +
    geom_segment(aes(x = dat_spp_realm %>%
                       filter(species == selected_species) %>%
                       pull(A1) %>%
                       unique(),
                     y = dat_spp_realm %>%
                       filter(species == selected_species) %>%
                       pull(A2) %>%
                       unique(),
                     xend = A1,
                     yend = A2),
                 colour = colour_purple,
                 data = dat_spp_realm %>%
                   filter(species %in% closest_global[[1]])) +
    scale_y_continuous(breaks = c(-0.3, 0, 0.3)) +
    scale_x_continuous(breaks = c(-0.3, 0, 0.3)) +
    coord_cartesian(xlim = c(-0.5, 0.5),
                    ylim = c(-0.5, 0.5)) +
    labs(x = "PCoA1",
         y = "PCoA2",
         colour = NULL) +
    theme(legend.position = "none")
  
  plot_2 <- ggplot() +
    annotate("point", 
             x = 0.5, 
             y = 0.5, 
             shape = 21, 
             fill = "white", 
             colour = "black", 
             size = 35) +
    annotate("rect", 
             xmin = 0.2,
             xmax = 0.8,
             ymin = 0,
             ymax = 1, 
             colour = colour_grey, 
             fill = colour_grey) +
    annotate("point", 
             x = c(1.38, 1.38),
             y = c(1, 0.2), 
             colour = c(colour_purple, colour_mint),
             size = c(2, 2)) +
    annotate("text", 
             x = 2.3, 
             y = c(1, 0.2), 
             label = c("Five nearest neighbours \nin global space", 
                       "Five nearest neighbours \nin local space"), 
             size = 7/.pt, 
             colour = c(colour_purple, colour_mint)) +
    annotate("text", 
             x = c(0.5, -0.5), 
             y = c(1.2, 1.5), 
             label = c("Local space",
                       "Global space"),
             size = 7/.pt, 
             colour = c(colour_grey, "black")) +
    coord_cartesian(xlim = c(-1, 3.2), 
                    ylim = c(-0.5, 1.5)) +
    theme_void()
  
  # visualise A1 and A3
  plot_3 <- dat_spp_realm %>%
    filter(realm == selected_realm) %>%
    ggplot(aes(A1, A3)) +
    geom_polygon(fill = NA,
                 colour = "grey20",
                 data = dat_hull_2_ov) +
    geom_polygon(fill = colour_grey,
                 data = dat_hull_2 %>%
                   filter(realm == selected_realm)) +
    geom_point(data = dat_spp_realm %>%
                 filter(species %in% closest_local[[2]]),
               colour = colour_mint) +
    geom_segment(aes(x = dat_spp_realm %>%
                       filter(species == selected_species) %>%
                       pull(A1) %>%
                       unique(),
                     y = dat_spp_realm %>%
                       filter(species == selected_species) %>%
                       pull(A3) %>%
                       unique(),
                     xend = A1,
                     yend = A3),
                 colour = colour_mint,
                 data = dat_spp_realm %>%
                   filter(species %in% closest_local[[2]])) +
    geom_point(data = dat_spp_realm %>%
                 filter(species %in% closest_global[[2]]),
               colour = colour_purple) +
    geom_segment(aes(x = dat_spp_realm %>%
                       filter(species == selected_species) %>%
                       pull(A1) %>%
                       unique(),
                     y = dat_spp_realm %>%
                       filter(species == selected_species) %>%
                       pull(A3) %>%
                       unique(),
                     xend = A1,
                     yend = A3),
                 colour = colour_purple,
                 data = dat_spp_realm %>%
                   filter(species %in% closest_global[[2]])) +
    scale_y_continuous(breaks = c(-0.3, 0, 0.3)) +
    scale_x_continuous(breaks = c(-0.3, 0, 0.3)) +
    coord_cartesian(xlim = c(-0.5, 0.5),
                    ylim = c(-0.5, 0.5)) +
    labs(x = "PCoA1",
         y = "PCoA3",
         colour = NULL) +
    theme(legend.position = "none")
  
  # visualise A2 and A3
  plot_4 <- dat_spp_realm %>%
    filter(realm == selected_realm) %>%
    ggplot(aes(A2, A3)) +
    geom_polygon(fill = NA,
                 colour = "grey20",
                 data = dat_hull_3_ov) +
    geom_polygon(fill = colour_grey,
                 data = dat_hull_3 %>%
                   filter(realm == selected_realm)) +
    geom_point(data = dat_spp_realm %>%
                 filter(species %in% closest_local[[3]]),
               colour = colour_mint) +
    geom_segment(aes(x = dat_spp_realm %>%
                       filter(species == selected_species) %>%
                       pull(A2) %>%
                       unique(),
                     y = dat_spp_realm %>%
                       filter(species == selected_species) %>%
                       pull(A3) %>%
                       unique(),
                     xend = A2,
                     yend = A3),
                 colour = colour_mint,
                 data = dat_spp_realm %>%
                   filter(species %in% closest_local[[3]])) +
    geom_point(data = dat_spp_realm %>%
                 filter(species %in% closest_global[[3]]),
               colour = colour_purple) +
    geom_segment(aes(x = dat_spp_realm %>%
                       filter(species == selected_species) %>%
                       pull(A2) %>%
                       unique(),
                     y = dat_spp_realm %>%
                       filter(species == selected_species) %>%
                       pull(A3) %>%
                       unique(),
                     xend = A2,
                     yend = A3),
                 colour = colour_purple,
                 data = dat_spp_realm %>%
                   filter(species %in% closest_global[[3]])) +
    scale_y_continuous(breaks = c(-0.3, 0, 0.3)) +
    scale_x_continuous(breaks = c(-0.3, 0, 0.3)) +
    coord_cartesian(xlim = c(-0.5, 0.5),
                    ylim = c(-0.5, 0.5)) +
    labs(x = "PCoA2",
         y = "PCoA3",
         colour = NULL) +
    theme(legend.position = "none")
  
  (plot_1 + plot_2) /
    (plot_3 + plot_4)
}


plot_uniq <- plot_uniqueness("Eastern Indo-Pacific", 
                             "Neomonachus schauinslandi")

# patch together and save -------------------------------------------------


plot_final <- (plot_spec / 
  plot_uniq) +
  plot_layout(heights = c(1, 1, 2.5), 
              tag_level = 'new') +
  plot_annotation(tag_levels = list(c("A", rep("", 3), 
                                      "B", rep("", 3))))

# save plot
ggsave(plot_final,
       filename = here("figures",
                       "main",
                       "space_per_realm.pdf"),
       width = 183, height = 100*2,
       units = "mm",
       bg = "white")
