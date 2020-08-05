
spatialAutocorrelation <- function(occurrenceRecords,rasterLayers,autocorrelationClassDistance,autocorrelationMaxDistance,autocorrelationSignif) {
  
  set.seed(42)
  
  if(nrow(occurrenceRecords) > 1000) { 
    
    cat( paste0("\n"))
    cat( paste0("\n"))
    cat("More than 1000 occurrence records.","\n")  
    cat("Using a maximum of 1000 random records.","\n")
    
    occurrenceRecords <- occurrenceRecords[sample(1:nrow(occurrenceRecords),1000,replace = FALSE),]
    
  }
  
  occurrenceRecords <- data.frame(Lon=as.numeric(as.character(occurrenceRecords[,1])),Lat=as.numeric(as.character(occurrenceRecords[,2])))
  to.keep <- which(!duplicated(occurrenceRecords[,1]))
  occurrenceRecords <- occurrenceRecords[to.keep,]
    
  # --------
  
  presences.environment <- data.frame(raster::extract(rasterLayers,occurrenceRecords,stringsAsFactors=FALSE))
  to.keep <- which(!is.na(presences.environment[,1]))
  presences.environment <- presences.environment[to.keep,]
  occurrenceRecords <- occurrenceRecords[to.keep,]
  
  # --------

  space <- spDists(as.matrix(occurrenceRecords),as.matrix(occurrenceRecords),longlat=TRUE)
  data <- ecodist::distance( presences.environment , method = "euclidean")
  data <- as.matrix(data)
  
  n.class <- round(autocorrelationMaxDistance / autocorrelationClassDistance)
  
  resultsMatrix <- data.frame(classdistanceFrom=seq(0,autocorrelationMaxDistance-autocorrelationClassDistance,autocorrelationClassDistance),
                              classdistanceTo=seq(autocorrelationClassDistance,autocorrelationMaxDistance,autocorrelationClassDistance),
                              R=NA,
                              pVal=NA)
  
  for( i in 1:nrow(resultsMatrix)) {
    
    d1 = resultsMatrix[i,1]
    d2 = resultsMatrix[i,2]
    
    data.d <- as.vector(data)
    space.d <- as.vector(space)
    
    remove <- which(space.d < d1 | space.d > d2)
    data.d <- data.d[-remove]
    space.d <- space.d[-remove]

    modelobject <- lm(space.d~data.d)
    
    f <- summary(modelobject)$fstatistic
    p <-0
    
    tryCatch( p <- pf(f[1],f[2],f[3],lower.tail=FALSE) , error=function(e) { Error <<- TRUE })
    
    #p <- summary.lm(modelobject)$coefficients[2,"Pr(>|t|)"] 
    resultsMatrix[i,3] <- summary(modelobject)$adj.r.squared
    resultsMatrix[i,4] <- p
    
  }
  
  vect.signif <- as.numeric(resultsMatrix[,4] > autocorrelationSignif)
  if(vect.signif[1] == 1) { vect.signif[1] <- 0}
  vect.signif.pch <- c(19,1)
  
  par(mar = c(4.5, 5.5, 4.5, 4.5) , bg = "#FFFFFF")
  plot(resultsMatrix[,2],resultsMatrix[,3],axes=FALSE,xlab="Distance (Km)",ylab="")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#F6F6F6")
  lines(resultsMatrix[,2],resultsMatrix[,3],lty=2,col="#000000",type="l",xlab="Distance (Km)",ylab="")
  
  points(resultsMatrix[,2],resultsMatrix[,3],pch=vect.signif.pch[vect.signif+1],col="#5B5B5B")
  axis(2,las=2,col="White",col.ticks="Black")
  axis(1,las=0,col="White",col.ticks="Black")
  box()
  title(ylab="Correlation (R)",mgp=c(4,1,0)) 
  abline(h = 0,lty=1)
  
  distance <- round( resultsMatrix[ which(resultsMatrix[,4] >= autocorrelationSignif)  , 2 ][1] )
  if( is.na(distance)) { distance <- autocorrelationMaxDistance }
  cat( paste0("\n"))
  cat( paste0("\n"))
  cat( paste0("First non-correlated distance: ",distance," km"))
  
  return(distance)
  
}


spatialThinning <- function(occurrenceRecords,rasterLayers,minDistance,verbose) {
  
  if( missing(verbose) ) { verbose <- TRUE }
  
  occurrenceRecords <- data.frame(Lon=as.numeric(as.character(occurrenceRecords[,1])),Lat=as.numeric(as.character(occurrenceRecords[,2])))
  to.keep <- which(!duplicated(occurrenceRecords[,1]))
  shape <- subset(rasterLayers,1)
  
  occurrenceRecords.i <- rasterize(occurrenceRecords,shape)
  occurrenceRecords.i <- xyFromCell(occurrenceRecords.i, Which(!is.na(occurrenceRecords.i),cells=TRUE) )
  colnames( occurrenceRecords.i ) <- c( "Lon", "Lat" )
  
  dataThinning <- data.frame(Name="Sp",occurrenceRecords.i)
  
  coordinates.t <- thin(dataThinning,
                        lat.col = "Lat",
                        long.col = "Lon",
                        spec.col = "Name",
                        thin.par = minDistance,
                        reps = 1,
                        write.files = FALSE,
                        locs.thinned.list.return = TRUE,
                        verbose = FALSE)[[1]]
  
  if(verbose) {
    
    cat( paste0("\n"))
    cat( paste0("\n"))    
    cat( paste0("Input Records: ",nrow(occurrenceRecords)))
    cat( paste0("\n"))
    cat( paste0("Final Records: ",nrow(coordinates.t)))
    
  }
  
  # Remove from main dataset of occurrences
  colnames( coordinates.t ) <- c( "Lon", "Lat" )
  
  # Remove log file
  if(length(list.files(".", full.names = TRUE, pattern = "spatial_thin_log.txt"))){
    file.remove( list.files(".", full.names = TRUE, pattern = "spatial_thin_log.txt") ) 
  }
  
  return(coordinates.t)
  
}




