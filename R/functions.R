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


get_metrics <- function(Mat_dist,Coords,nb_NN=5,GE,StandGE=F){
  
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
