calcSwitchErrors <- function(stitches, h1, h2, hm, hp,g=FALSE) {

if (length(h1)!=length(h2)) {
  print('Mismatching haplotype lengths')
  browser()
}
  
# an auxiliary function
searchMin <- function(v,minimumValue) {
# returns the leftmost element of a vect equal or larger than a value
  if (length(v)>0) {
    if (minimumValue<=max(v)){
      mi <- min(which(v>=minimumValue))
      r<- v[mi]}
    else r <- Inf
  }
  else {r<-Inf}
  r
}

switchesPerStitch <- vector()

  indSwitchErrors <- 0
  # loop over the phased stitches
  for (s in 1:(dim(stitches)[1])) {

        snpIndx <- stitches$start[s] : stitches$end[s]
        h1Stitch <- h1[snpIndx]
        h2Stitch <- h2[snpIndx]
        hmStitch <- hm[snpIndx]
        hpStitch <- hp[snpIndx] 
        #gStitch <- as.integer(g[snpIndx])
        #gStitchHomo <- (gStitch==1)|(gStitch==3)
       
        # count the inconsistencies between the pairs of haplotypes
        h1hm <-  (((h1Stitch==3) & (hmStitch==1)) | ((h1Stitch==1)&(hmStitch==3))); 
        h1hp <- (((h1Stitch==3)&(hpStitch==1)) | ((h1Stitch==1)&(hpStitch==3))); 

        h2hm <-  (((h2Stitch==3)&(hmStitch==1)) | ((h2Stitch==1)&(hmStitch==3))); 
        h2hp <- (((h2Stitch==3)&(hpStitch==1)) | ((h2Stitch==1)&(hpStitch==3)));

        # find the loci where the phase would need to be switched
        # done to distinguish switch errors from phasing errors
        wh1hm <- which(h1hm & h2hp)
        wh2hm <- which(h2hm & h1hp)


        # haplotype corresponding to maternal haplotype, one in 1,2
                                        # initial choice of phase

        
        
        opp1 <- sum(h1hm[1:100]) + sum(h2hp[1:100])
        opp2 <- sum(h2hm[1:100]) + sum(h1hp[1:100])

        if (opp1<opp2) { phase<-1} else {phase <-2}

        # looking to inconsistencies with respect to the phase
        #if (phase==1) {
        #  phaseInconsistency <- min(searchMin(wh1hm, 0), searchMin(wh2hp,0))
        #} else if (phase==2){
        #  phaseInconsistency <- min(searchMin(wh2hm, 0), searchMin(wh1hp,0))
        #}

        if (phase==1) {
          phaseInconsistency <- searchMin(wh1hm, 0)
        } else if (phase==2){
          phaseInconsistency <- searchMin(wh2hm, 0)
        }

        #incAll <- min(sum(h1hm),sum(h1hp)) +  min(sum(h2hm),sum(h2hp))
        #if ((incAll>0)&(length(h1Stitch)<500)) {browser()}
         
        switchCount <- 0
        while (!is.infinite(phaseInconsistency)) {

          if (phase==1) {phase<-2} else {phase<-1}
          switchCount <- switchCount+1;

          #if (phase==1) {
          #  phaseInconsistency <- min(searchMin(wh1hm, phaseInconsistency+1), searchMin(wh2hp,phaseInconsistency+1))
          #} else if (phase==2){
          #  phaseInconsistency <- min(searchMin(wh2hm, phaseInconsistency+1), searchMin(wh1hp,phaseInconsistency+1))
          #}
          if (phase==1) {
            phaseInconsistency <- searchMin(wh1hm, phaseInconsistency+1)
          } else if (phase==2){
            phaseInconsistency <- searchMin(wh2hm, phaseInconsistency+1)
          }
          
          
        }

        switchesPerStitch[length(switchesPerStitch)+1] <- switchCount
  }
  # end loop over stitches
  
  #print(paste('noStitches', dim(stitches)[1]))

  
  
#finalResult <- sum(totalSwitches)/(0.33*sum(totalStitchLength))
#mDistance <- 0.01*sum(allData$geneticSnpMap[stitches$end] - allData$geneticSnpMap[stitches$start])
#morganResult <- sum(totalSwitches)/(sum(kimmoWithSeqedParents)*mDistance)


switchesPerStitch

}
