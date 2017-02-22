# calcCoverage
#
# calulates the ratio of region covered by at least sequence
#
# sequences - surrogate sequences, with start and end specified
# allSnps - number of consecutive SNPs in the region of interest
calcCoverage <-function (sequences, allSnps) {

    result <- list()

    noSels <- dim(sequences)[1]

    
    if (noSels>0) {
     sequences <- sequences[order(sequences$start) , ]

    stitches <<- data.frame()
    
    stitchStart <- sequences[1, 'start']
    stitchEnd <- sequences[1, 'end']
    i<-1;
    while (i<=noSels) {
      currStart <- sequences[i, 'start']
      currEnd <-  sequences[i, 'end']     
                                        # see if the current stitch expired
      if((stitchEnd-1)<currStart) {        
        stitch <- data.frame(start=stitchStart, end = stitchEnd)
        stitches <- rbind(stitches, stitch)
        stitchStart <- currStart
        stitchEnd <- currEnd
      }
      else {
                                        # extend the stitch
        stitchEnd <- max (stitchEnd, currEnd)    
      }
      i <- i + 1  
    }

    coverage <- -1
    
    stitch <- data.frame(start=stitchStart, end = stitchEnd)
    stitches <<- rbind(stitches, stitch)

    coverage <-  round(sum(stitches$end-stitches$start)/allSnps, digits=5)

    data.frame(coverage=coverage, noStitches=dim(stitches)[1])

    result[['coverage']] <- coverage;
    result[['noStitches']] <- dim(stitches)[1]
    result[['stitches']] <- stitches
    }
    else {
    result[['coverage']] <- 0;
    result[['noStitches']] <- 0
    result[['stitches']] <- 0
    }
    
    result
    
  }



