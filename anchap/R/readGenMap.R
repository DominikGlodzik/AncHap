readGenMap <-function(filePath=NULL, build=36) {

  mapPosBp <- list()
  mapRecRt <- list()
  mapCm <- list()
  mapPosBpRaw <- list()
  cumEnds <- vector()

  noBlocks <- 50;
  
  for (c in 1:22) {

    if (build==36) {
      r <- read.table(paste(filePath[[1]], c, filePath[[2]], sep=""), header=TRUE)
      noPoints <- dim(r)[1];
      blockLength <- floor(noPoints/noBlocks)
      m <- matrix(r$COMBINED_rate.cM.Mb.[1: (noBlocks* blockLength)], blockLength,noBlocks )
      mapPosBpRaw[[c]] <- r$position;      
      mapCm[[c]] <- r$Genetic_Map.cM.     
    }    
    else {
      # map build 37
      r <- read.table(paste(filePath[[1]], c,filePath[[2]],sep=""), header=TRUE)
      #r <- read.table(paste('../Data/data/genetic_map_HapMapII_GRCh37.tar/genetic_map_GRCh37_chr', c,'.txt',sep=""), header=TRUE)  
      noPoints <- dim(r)[1];
      blockLength <- floor(noPoints/noBlocks)
      m <- matrix(r$Rate.cM.Mb.[1: (noBlocks* blockLength)], blockLength,noBlocks )
      mapPosBpRaw[[c]] <- r$Position.bp.;
      mapCm[[c]] <- r$Map.cM.      
    }
    
    m2 <- colSums(m)/blockLength
    mapPosBp[[c]] <- 1:noBlocks
    mapRecRt[[c]] <- m2
         

    cumEnds[c] <- max(mapPosBpRaw[[c]])


    
  }


  genMap <- list();

  genMap[['mapPosBp']] <- mapPosBp;
  genMap[['mapRecRt']] <- mapRecRt;
  
  # position [bp] of last marker in the chromosome
  genMap[['cumEnds']] <- cumsum(as.numeric(cumEnds));
  # RAW marker position [bp]
  genMap[['mapPosBpRaw']] <- mapPosBpRaw;
  # genetic map
  genMap[['mapCm']] <- mapCm;
  genMap
}
