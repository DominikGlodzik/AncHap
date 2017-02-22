readGenMapPredictions <-function(filePath=NULL) {

  mapPosBp <- list()
  mapRecRt <- list()
  mapCm <- list()
  mapPosBpRaw <- list()
  cumEnds <- vector()
  for (c in 1:22) {

      r <- read.table(paste(filePath[[1]], c,filePath[[2]],sep=""), header=TRUE)

    noPoints <- dim(r)[1]

    noBlocks <- 50;
    blockLength <- floor(noPoints/noBlocks)

    #m <- matrix(r$Rate.cM.Mb.[1: (noBlocks* blockLength)], blockLength,noBlocks )
    m <- matrix(r$Map.cM.[1: (noBlocks* blockLength)], blockLength,noBlocks )
    m2 <- colSums(m)/blockLength
    mapPosBpRaw[[c]] <- r$Position.bp.;
    mapPosBp[[c]] <- 1:noBlocks
    mapRecRt[[c]] <- m2
         sId <- seq(from = 1, to = noPoints, by =1)
    mapCm[[c]] <- r$Map.cM.[sId]

    cumEnds[c] <- max(r$Position.bp.)
    
  }


  genMap <- list();
  genMap[['mapPosBp']] <- mapPosBp;
  genMap[['mapPosBpRaw']] <- mapPosBpRaw;
  genMap[['mapRecRt']] <- mapRecRt;
  genMap[['cumEnds']] <- cumsum(as.numeric(cumEnds));
  genMap[['mapCm']] <- mapCm;
  genMap
}
