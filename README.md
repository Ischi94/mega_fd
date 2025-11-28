# mega_fd

This is the code and data repository for the manuscript "Mismatches in functional conservation priorities across scales and locations
worldwide" by Griffin et al.

The code is annotated and follows a chronological order. All figures and output are created from code. 

### 01_presabs_matrix.R

This code takes the geographic distributions from the marine megafauna species data from the IUCN database and the AquaMaps database and converts it into a presence-absence matrix. 
  
**Output**  
- IUCN/megafauna_IUCN_presabs_0.5res.txt  
- presabs_05_res.rds  

### 02_functional_space.R

This code uses the presence-absence matrix and the focal traits of the species to calculate a 5 dimensional functional space. 
  
**Output**   
- distance_trait_matrix.rds  
  
  
### 03_calculate_metrics.R

This code uses the presence-absence matrix and the functional space to calculate functional metrics on a global, provincial, and local level.  
  
**Output**   
- global_metrics.rds  
- local_metrics_1_res.rds  
- local_metrics_5_res.rds  
- realm_metrics.rds  
  
  
### 04_correlation_analysis.R

This code calculates correlations of functional metrics on a global, provincial, and local level.  
  
**Output**   
- fig_1.png    
  
  
### 05_rank_comparison.R

This code calculates how much rankings of functional metrics agree over spatial scales.  
  
**Output**   
- fig_2.png  
  
  
### 06_blue_whale_example.R

This code takes an example species (the blue whale) to showcase the variation in functional metrics over its range.    
  
**Output**   
- fig_3.png
  
  
### 07_hotspot_analysis.R

This code visualises the spatial extent of each functional metric on a global and local level, and calculates hotspots and their overlap.    
  
**Output**   
- fig_4.png  
- fig_5.png  
  
  
### 08_supplementary_table.R

This code creates the supplementary table showing where species reach the maximum functional uniqueness across their range.    
  
**Output**   
- data/supplementary_table.docx
  
  
### functions.R

These are the functions used throughout provided by Pimiento, C. et al. Functional diversity of sharks and rays is highly vulnerable and
supported by unique species and locations worldwide. Nature Communications 14,
7691 (2023). https://doi.org/10.1038/s41467-023-43212-3 
  
### Session Info  
  
R version 4.5.1 (2025-06-13 ucrt)  
Platform: x86_64-w64-mingw32/x64  
Running under: Windows 11 x64 (build 22631)  

Matrix products: default LAPACK version 3.12.1

locale:  
[1] LC_COLLATE=German_Germany.utf8  LC_CTYPE=German_Germany.utf8    LC_MONETARY=German_Germany.utf8
[4] LC_NUMERIC=C                    LC_TIME=German_Germany.utf8    

time zone: Europe/Berlin  
tzcode source: internal  

attached base packages:  
[1] stats     graphics  grDevices utils     datasets  methods   base      

other attached packages:  
 [1] officer_0.7.1    flextable_0.9.10 patchwork_1.3.1  funrar_1.5.0     raster_3.6-32   
 [6] sp_2.2-0         ade4_1.7-23      letsR_5.0        terra_1.8-60     sf_1.0-21       
[11] lubridate_1.9.4  forcats_1.0.0    stringr_1.5.1    dplyr_1.1.4      purrr_1.1.0     
[16] readr_2.1.5      tidyr_1.3.1      tibble_3.3.0     ggplot2_3.5.2    tidyverse_2.0.0 
[21] here_1.0.1      

loaded via a namespace (and not attached):  
 [1] gtable_0.3.6            xfun_0.52               lattice_0.22-7          tzdb_0.5.0             
 [5] vctrs_0.6.5             tools_4.5.1             generics_0.1.4          proxy_0.4-27           
 [9] pkgconfig_2.0.3         Matrix_1.7-3            KernSmooth_2.23-26      data.table_1.17.8      
[13] RColorBrewer_1.1-3      uuid_1.2-1              lifecycle_1.0.4         compiler_4.5.1         
[17] farver_2.1.2            textshaping_1.0.1       codetools_0.2-20        fontquiver_0.2.1       
[21] fontLiberation_0.1.0    htmltools_0.5.8.1       class_7.3-23            pillar_1.11.0          
[25] MASS_7.3-65             openssl_2.3.3           rsconnect_1.5.1         classInt_0.4-11        
[29] fontBitstreamVera_0.1.1 zip_2.3.3               tidyselect_1.2.1        digest_0.6.37          
[33] stringi_1.8.7           rprojroot_2.1.1         fastmap_1.2.0           grid_4.5.1             
[37] cli_3.6.5               magrittr_2.0.3          e1071_1.7-16            withr_3.0.2            
[41] gdtools_0.4.4           scales_1.4.0            timechange_0.3.0        rmarkdown_2.29         
[45] ragg_1.4.0              askpass_1.2.1           hms_1.1.3               evaluate_1.0.4         
[49] knitr_1.50              rlang_1.1.6             Rcpp_1.1.0              glue_1.8.0             
[53] geosphere_1.5-20        DBI_1.2.3               xml2_1.3.8              rstudioapi_0.17.1      
[57] R6_2.6.1                systemfonts_1.3.1       units_0.8-7   
