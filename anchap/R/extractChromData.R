extractChromData <- function(allData, chromNo) {

  allDataChrom <- list()

  chromSnpFlags <- (allData$map$chromosome==chromNo)
  noChromSnps <- sum(chromSnpFlags)


  allDataChrom$map <- allData$map[chromSnpFlags, ]

  allDataChrom$chromEnds <- integer(22)
  allDataChrom$chromEnds[chromNo] <- noChromSnps

  allDataChrom$cumEnds <- integer(22)
  allDataChrom$cumEnds[chromNo] <- max(allDataChrom$map$position) 
  allDataChrom$cumEnds <- cumsum(allDataChrom$cumEnds)
  
  allDataChrom$chromList <- list()
  for (cn in 1:22) {
     allDataChrom$chromList[[cn]] <- list(NULL)
  }
  allDataChrom$chromList[[chromNo]] <- allData$chromList[[chromNo]]

  allDataChrom$iGene <- list()
  for (i in 1:length(allData$iGene)) {
    allDataChrom$iGene[[i]] <-allData$iGene[[i]][chromSnpFlags]
  }

  allDataChrom$genotypeMatrix <- allData$genotypeMatrix[,chromSnpFlags]
  allDataChrom$geneticSnpMap <- allData$geneticSnpMap[chromSnpFlags]
  
  allDataChrom$probands <- allData$probands
  allDataChrom$idNames <- allData$idNames
  allDataChrom$pedigree <- allData$pedigree
  allDataChrom$orcaIds <- allData$orcaIds

  allDataChrom
 
}
