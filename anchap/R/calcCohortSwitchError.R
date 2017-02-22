calcCohortSwitchError <- function(phasedIds, allData, result,  matHapMat, patHapMat) {


  
totalSwitches <- vector()
totalStitches <- data.frame() 
for (i in 1:length(phasedIds)) {
#for (i in 1:1) {  

  phasedId <- phasedIds[i]
  #print(phasedId)
                       
  datasetInd <- which(allData$idNames == phasedId)
  
  #print(i)
  # extracting the haplotypes from the result
  # phasedIds[i] is in 1:597 
  h1 <- result$phasingResults[[datasetInd]]$h1
  h2 <- result$phasingResults[[datasetInd]]$h2
  g <- allData$iGene[[datasetInd]]

  # kimmoInd in 1:597
  kimmoInd <- which(rownames(matHapMat) == phasedId)
  hm <- matHapMat[kimmoInd,]
  hp <- patHapMat[kimmoInd,]

  # sanity checking
  homoLoci <- ((g==1)|(g==3))
  haploHomo <- (hm[homoLoci]==hp[homoLoci])
  consistentHaplo <- (g[homoLoci]==hp[homoLoci])

  print(paste("Kimmo id", kimmoInd))
  print(sum(!haploHomo)/length(homoLoci))
  print(sum(!consistentHaplo)/length(homoLoci))

  selectedSequences <- rbind(result$stitchingResults[[datasetInd]]$selectedSequencesH1, result$stitchingResults[[datasetInd]]$selectedSequencesH2)
  coverageSel <- calcCoverage(selectedSequences, noSnps)

  stitches <- coverageSel$stitches
  # choose stitches longer than 100 markers
  stitches <- stitches[(stitches$end - stitches$start)>500,]
  chr20Stitches <- stitches[(allData$map$chromosome[stitches$start]==20),]
  #snpIndx <- unlist(mapply(seq, chr20Stitches$start, chr20Stitches$end, simplify=TRUE))
  #stitches <- chr20Stitches


  
  # treating the whole chromosome as one stitch
  #stitches <- data.frame(start=allData$chromEnds[19]+1, end=allData$chromEnds[20] )

  r <- calcSwitchErrors(stitches, h1, h2, hm, hp,g)

  l <- length(totalSwitches)
  totalSwitches[(l+1):(l+length(r))] <- r
  totalStitches <- rbind(totalStitches, cbind(stitches, data.frame(phasedId=phasedId, kimmoInd=kimmoInd)))

}


print(totalStitches)
# calculate the result - average number of stiches per morgan

stitchMDistance <- 0.01*sum(allData$geneticSnpMap[totalStitches$end] - allData$geneticSnpMap[totalStitches$start])
morganResult <- sum(totalSwitches)/(stitchMDistance)

morganResult

}
