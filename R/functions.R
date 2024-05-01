# general configurations for plots

#set ggplot output
ggplot2::theme_set(ggplot2::theme_minimal(base_size = 10))

ggplot2::theme_update(text = element_text(colour = "grey20", size = 10), 
                      legend.text = element_text(colour = "grey20", size = 10),
                      legend.title = element_text(colour = "grey20", size = 10),
                      axis.text = element_text(colour = "grey40", size = 10),
                      axis.title = element_text(colour = "grey20", 
                                                size = 10, 
                                                face = "bold"),
                      axis.ticks = element_line(colour = "grey50"),
                      strip.text = element_text(colour = "grey20", size = 10),
                      panel.grid.minor = element_blank(),  
                      panel.grid.major = element_blank(),
                      plot.title = element_text(colour = "grey20", size = 10, 
                                                face = "bold"), 
                      plot.subtitle = element_text(colour = "grey20", size = 10,
                                                   face =  "italic"))

# define output sizes
image_width <- 183
image_height <- 100
image_units <- "mm"

# define pallets

# define common color
colour_yellow = "#ffc96a"
colour_purple = "#c198d6ff"
colour_mint = "#5d7a64ff"
colour_grey = "grey55"
colour_coral = "#F95875"

get_FV_Sp <-
  function(ax = ax,
           Fric_tot = richness_funct_global_ST[[1]][2],
           pcoa = pcoa,
           Selected_sp) {
    coord_d <- pcoa$li[, ax]
    coord_d_k <- na.omit(coord_d[as.character(Selected_sp), ])
    rs <- dim(coord_d_k)[1]
    
    # functional Richness
    FRic <- round(geometry::convhulln(coord_d_k, "FA")$vol, 10)
    
    # identity of vertices
    vert0 <- geometry::convhulln(coord_d_k, "Fx TO 'vert.txt'")
    vert1 <- scan("vert.txt", quiet = T)
    vert2 <- (vert1 + 1)[-1]
    
    FE_vert_k <- row.names(coord_d_k)[vert2]
    FE_vertices_cells <- FE_vert_k[!is.na(FE_vert_k)] # arendre
    
    # number of vertices
    nbFE_vertices <- length(FE_vert_k)
    
    ### Calcul du FDiv
    # coord values of vertices
    trvertices <- coord_d_k[vert2, ]
    
    # coordinates of the center of gravity of the vertices (B)
    B <- apply(trvertices, 2, function(x) {
      mean(x, na.rm = TRUE)
    })
    
    # computing euclidian dstances to B (dB)
    dB <- apply(coord_d_k, 1, function(x) {
      (sum((x - B) ^ 2)) ^ 0.5
    })
    
    # mean of dB values, deviations to mean
    meandB <- mean(dB)
    devdB <- dB - meandB
    
    # computation of FDiv
    FDiv <-
      round((sum(devdB) + meandB) / (sum(abs(devdB)) + meandB) , 6)# arendre
    
    # computation of FRic
    Fric_r <- round(FRic / geometry::convhulln(coord_d, "FA")$vol, 6)
    
    # computation of originality of each species: distance to nearest neighbour among the global pool of species
    #dist_sp<-as.matrix(dist(coord_d_k,method="euclidean")) ; dist_sp[which(dist_sp==0)]<-NA
    #orig_sp<-apply(dist_sp, 1, min, na.rm=T )
    #Orig_vect <- c(mean=mean(orig_sp,na.rm=T), min=min(orig_sp, na.rm=T), max=max(orig_sp, na.rm=T), se=sd(orig_sp)/sqrt(length(orig_sp)))
    
    #computing functional specialization
    O <- apply(coord_d_k, 2, mean)
    speS <- apply(coord_d_k, 1, function(x) {
      (sum((x - O) ^ 2)) ^ 0.5
    })
    Specialization <-
      c(
        mean = mean(speS, na.rm = T),
        min = min(speS, na.rm = T),
        max = max(speS, na.rm = T),
        se = sd(speS) / sqrt(length(speS))
      )
    
    #computing FEs
    #FEs <-unique(FE[Selected_sp,2])
    
    #computing FR
    #FR<- (length(Selected_sp))/length(FEs)
    
    #computing FV
    #r <- table(FE[Selected_sp,2])
    #FV <-(length(FEs)-sum(sapply(r,function(x){min(x-1,1)}))/length(FEs))
    
    list(
      data = data.frame(
        RS = rs,
        FRic = FRic,
        FDiv = FDiv,
        nb_vertices = nbFE_vertices,
        Fric_r = Fric_r
      ),
      specialization = Specialization,
      Sp_vertice = FE_vertices_cells
    )
  }  


get_FUSE <- function(Mat_dist,Coords,nb_NN=5,GE,StandGE=F){
  
  require(vegan)
  require(reshape2)
  
  if(!identical(row.names(as.matrix(Mat_dist)),row.names(Coords))){
    stop("Coords lines do not match with the distance matrix")
  } 
  
  nm <- rownames(Coords)
  
  # Specialization calculation
  O <- apply(Coords,2,mean)
  spe <- apply(Coords, 1,function(x){sum((x-O)^2)^0.5})
  
  # Distinctivness calculation 
  dist_sp <- as.matrix(Mat_dist)
  Fdistinct <- apply(dist_sp,1,mean)
  
  # Uniqueness calculation
  uni_res <- get_indicator(Mat_dist=as.matrix(Mat_dist),nb_NN=nb_NN)
  uniqu <- uni_res$Average_uniqueness[,"Mean"]
  
  if(StandGE==TRUE){
    GE <- as.vector(decostand(GE,"range",na.rm=T))
  } 
  
  # FUSE metrics calculation
  FUn_std <- as.vector(decostand(uniqu,"range"))
  FUGE <- log(1+(FUn_std*GE))
  FSp_std <- as.vector(decostand(spe,"range")) 
  FSGE <- log(1+(FSp_std*GE))
  FUS <- FUn_std+FSp_std
  FDist_std <- as.vector(decostand(Fdistinct,"range"))
  FUSE <- setNames(FUGE+FSGE,nm=nm)
  FUSE_alt<-log(1+FUS*GE)
  FDGE <- log(1+(FDist_std*GE))
  
  data.frame(cbind(FUSE,FUSE_alt,FUGE,FDGE,FSGE, FUS,FUn_std,FSp_std,FDist_std))
  
}

get_indicator <- function(Mat_dist,nb_NN){
  
  w <- reshape2::melt(Mat_dist)
  s <- split(w,f=w[,2])
  
  Res <- lapply(s,function(x){get_dist_func(nb_NN=nb_NN,data=x)})
  Res_mean_sd <- do.call(rbind,lapply(1:length(Res),function(i){Res[[i]][[1]]}))
  NN <- lapply(1:length(Res),function(i){Res[[i]][[2]]})
  
  rownames(Res_mean_sd) <- names(NN) <- names(Res)
  list(Average_uniqueness=Res_mean_sd,Nearest_neighbour=NN)
  
}


get_dist_func <- function(nb_NN,data){
  
  data <- data[order(data[,3], decreasing=F),]
  data <- data[-1,]
  mm <- mean(data[1:nb_NN,3])
  sd <- sd(data[1:nb_NN,3])
  sp <- as.character(data[1:nb_NN,1])
  list(c(Mean=mm,Sd=sd),Species=sp)
  
} 
