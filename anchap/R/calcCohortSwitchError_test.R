calcCohortSwitchErrorTest<- function(phasedIds, matHapMat1, patHapMat1, matHapMat2, patHapMat2, allData ) {


  
totalSwitches <- vector()
totalStitches <- data.frame() 
for (i in 1:length(phasedIds)) {
#for (i in 1:1) {  

  phasedId <- phasedIds[i]
  print(phasedId)
                       
  kimmoInd <- which(rownames(matHapMat1) == phasedId)
  h1 <- matHapMat1[kimmoInd,]
  h2 <- patHapMat1[kimmoInd,]

  # kimmoInd in 1:597
  kimmoInd <- which(rownames(matHapMat2) == phasedId)
  hm <- matHapMat2[kimmoInd,]
  hp <- patHapMat2[kimmoInd,]


  
  #chr20Stitches <- stitches[(allData$map$chromosome[stitches$start]==20),]
  #snpIndx <- unlist(mapply(seq, chr20Stitches$start, chr20Stitches$end, simplify=TRUE))
  #stitches <- chr20Stitches

  # treating the whole chromosome as one stitch
  stitches <- data.frame(start=allData$chromEnds[19]+1, end=allData$chromEnds[20] )

  r <- calcSwitchErrors(stitches, h1, h2, hm, hp,g)

  print(r)
  
  l <- length(totalSwitches)
  totalSwitches[(l+1):(l+length(r))] <- r
  totalStitches <- rbind(totalStitches, cbind(stitches, data.frame(phasedId=phasedId, kimmoInd=kimmoInd)))
  
  
}

# calculate the result - average number of stiches per morgan

stitchMDistance <- 0.01*sum(allData$geneticSnpMap[totalStitches$end] - allData$geneticSnpMap[totalStitches$start])
morganResult <- sum(totalSwitches)/(stitchMDistance)

morganResult

}
