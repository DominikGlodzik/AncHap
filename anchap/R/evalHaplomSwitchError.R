evalHaplomSwitchError <- function(haplosM, matHapMat, patHapMat, snpRegion, kimmoWithSeqedParents, chromEnds, geneticSnpMap, chroms=1:22) {

  
      stitches <- data.frame()
      for (c in chroms) {
        if (chromEnds[c]>0) {
          if (c==1) {prevEnd<-0} else {prevEnd<-chromEnds[c-1]}
          stitch <- data.frame(start=(prevEnd+1), end=chromEnds[c])
                  stitches<-rbind(stitches,stitch)
        }

      }
  
  phasedIds <- rownames(matHapMat)[kimmoWithSeqedParents]
  if (sum(allData$idNames %in% phasedIds)>0) {    

  switches <- vector()
  stitchesGenLength <- vector()
  for (i in 1:length(phasedIds)) {

    i1 <- which(rownames(matHapMat)==phasedIds[i])
    i2 <- which(allData$idNames==phasedIds[i])
    if ((2*i2)<=(dim(haplosM)[1])) {


      hm <- matHapMat[i1,]
      hp <- patHapMat[i1,]

      h1 <-  haplosM[2*i2-1,]
      h2 <- haplosM[2*i2,]
      g <- allData$iGene[[i2]][snpRegion]
        


      r <- calcSwitchErrors(stitches, h1, h2, hm, hp)
       
      switches[(length(switches)+1) :(length(switches)+ length(r))]<-r
      stitchesGenLength[(length(stitchesGenLength)+1) :(length(stitchesGenLength)+ length(r))] <- 0.01*(geneticSnpMap[stitches$end] - geneticSnpMap[stitches$start])
      
    }
  }
 
  r <- sum(switches)/sum(stitchesGenLength)
  } else {r <- NA}

}
