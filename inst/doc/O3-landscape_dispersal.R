## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  library(sf)
#  mylandscape <- st_read(dsn = "myshapefile.shp")
#  library(landsepi)
#  simul_params <- createSimulParams(outputDir = getwd())
#  simul_params <- setLandscape(simul_params, mylandscape)
#  simul_params@Landscape

## ----eval=FALSE---------------------------------------------------------------
#  mylandscape$year_1 <- c(13,2,4,1,1) # croptypes ID allocated to the different polygons
#  mylandscape$year_2 <- c(2,2,13,1,1)

## ----eval=FALSE---------------------------------------------------------------
#  simul_params <- setLandscape(simul_params, mylandscape)
#  simul_params@Landscape

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("RCALI")
#  library(RCALI)
#  library(landsepi)
#  landscape <- landscapeTEST1
#  Npoly <- length(landscape)
#  Npoly
#  plot(landscape)

## ----eval=FALSE---------------------------------------------------------------
#  file_land <- "land_rcali.txt"  ## input for califlopp
#  file_disp <- "disp_rcali.txt"  ## output for califlopp (DO NOT WRITE A PATH)
#  
#  ## Formatting the polygons-file for califlopp
#  cat(Npoly, file=file_land)
#  for (k in 1:Npoly) {
#    ## extract coordinates of polygon vertices
#    coords <- landscape@polygons[[k]]@Polygons[[1]]@coords
#    ## alternatively:
#    # coords <- as.data.frame(landscape$geometry[[k]][[1]])
#    n <- nrow(coords)
#    cat(NULL, file=file_land, append=T, sep="\n")
#    cat(c(k,k,n), file=file_land, append=T, sep="\t")
#    cat(NULL, file=file_land, append=T, sep="\n")
#    cat(coords[1:n,1], file=file_land, append=T, sep="\t")
#    cat(NULL,file=file_land,append=T,sep="\n")
#    cat(coords[1:n,2], file=file_land, append=T, sep="\t")
#  }
#  cat(NULL, file=file_land, append=T, sep="\n")

## ----eval=FALSE---------------------------------------------------------------
#  param <- list(input=2, output=0, method="cub", dp=6000, dz=6000
#                , warn.poly=FALSE, warn.conv=FALSE, verbose=FALSE)
#  califlopp(file=file_land, dispf=1, param=param, resfile=file_disp)

## ----eval=FALSE---------------------------------------------------------------
#  my_df <-function(x, a=40, b=7) ((b-2)*(b-1)/(2*a^2*pi)*(1+(abs(x)/a))^(-b))
#  
#  param <- list(input=2, output=0, method="cub", dp=6000, dz=6000, warn.poly=FALSE,
#                warn.conv=FALSE, verbose=FALSE)
#  califlopp(file=file_land, dispf=my_df, param=param, resfile=file_disp)

## ----eval=FALSE---------------------------------------------------------------
#  ## Import califlopp results
#  disp_df <- getRes(file_disp)
#  ## Double the table because only half of the flows have been calculated
#  emitter <- c(disp_df$poly1, disp_df$poly2)
#  receiver <- c(disp_df$poly2, disp_df$poly1)
#  
#  ## Write a text file containing a vector of areas of all polygons
#  area_e <- c(disp_df$area1, disp_df$area2)
#  area_r <- c(disp_df$area2, disp_df$area1)
#  area <- as.vector(by(area_e, emitter, mean))
#  write(area, file="area.txt", sep=",")
#  
#  ## Generation of the dispersal matrix
#  name_f <- "mean.flow"
#  flow_mean <- c(disp_df[,name_f], disp_df[,name_f])
#  flow_f <- cbind(emitter, receiver, flow_mean, area_e, area_r)
#  
#  ## Remove the doublons (i.e. half the lines where emitter == receiver)
#  flow_f[1:nrow(disp_df),][(disp_df$poly2 - disp_df$poly1) == 0,] <- NA
#  flow_f <- flow_f[is.na(apply(flow_f, 1, sum)) == F,]
#  flow_f <- as.data.frame(flow_f)
#  colnames(flow_f) <- c("emitter", "receiver", "flow", "area_e", "area_r")
#  flow_f <- flow_f[order(flow_f$emitter),]
#  
#  ## lines: emitter
#  ## columns: receiver
#  matrix_f <- NULL
#  for(k in 1:Npoly){
#    ## flow divided by the emitter area
#    matrix_f <- cbind(matrix_f, flow_f$flow[flow_f$receiver==k] / area)
#  }
#  
#  ## Normalisation of the matrix (reflecting boundaries)
#  ## (do not normalise for absorbing boundaries)
#  flowtot_f <- apply(matrix_f,1,sum)
#  for(k in 1:Npoly){
#    matrix_f[k,] <- (matrix_f[k,] / flowtot_f[k]) ## In order to have sum == 1
#  }
#  
#  write(as.vector(matrix_f), file="dispersal.txt", sep=",")

## ----eval=FALSE---------------------------------------------------------------
#  disp_patho <- scan("dispersal.txt", sep=",")

## ----eval=FALSE---------------------------------------------------------------
#  landscape <- landscapeTEST1
#  plot(landscape)
#  plotland(landscape)

## ----eval=FALSE---------------------------------------------------------------
#  poly <- 10
#  colFields <- rep("white", length(landscape))
#  colFields[poly] <- "red"
#  plot(landscape, col = colFields)

## ----eval=FALSE---------------------------------------------------------------
#  ## convert dispersal in matrix
#  mat <- matrix(disp_patho, nrow=sqrt(length(disp_patho)))
#  poly <- 1
#  dispToPlot <- log10(mat[poly,] +1E-20)  ## 1E-20 to avoid log(0)
#  
#  ## Colour palette
#  nCol <- 11
#  whiteYellowRed <- colorRampPalette(c("white", "#FFFF99", "#990000"))
#  col_disp <- whiteYellowRed(nCol)
#  intvls <- seq(min(dispToPlot) - 1, max(dispToPlot) + 1, length.out=nCol)
#  intvls_disp <- findInterval(dispToPlot, intvls)
#  
#  ## Plot
#  plot(landscape, col = col_disp[intvls_disp], main=paste("Dispersal from polygon", poly))

## ----eval=FALSE---------------------------------------------------------------
#  library(ggplot2)
#  ggplot(landscape) + ggtitle(paste("Dispersal from polygon", poly)) +
#      geom_sf(colour="black", aes(fill = dispToPlot)) +
#      scale_fill_gradientn(name="Prob. of\ndispersal", colours=rev(heat.colors(10)), breaks=-1:-10, labels=10^(-1:-10)) +
#      # theme_classic() +
#      theme(axis.line=element_blank(),axis.text.x=element_blank(),
#            axis.text.y=element_blank(),axis.ticks=element_blank(),
#            axis.title.x=element_blank(),
#            axis.title.y=element_blank(),
#            panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
#            panel.grid.minor=element_blank(),plot.background=element_blank())

