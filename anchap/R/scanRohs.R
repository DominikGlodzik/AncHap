scanRohs <- function(iGene,chromEnds, destinationFile, rohThresh, probands) {

noSubjects <- max(probands)
subjects <- probands

rohList<- array(data.frame(), c(noSubjects))

for (p in probands) {
  print(paste("Proband ", p))
  
  rohList[[p]] <- extractROHc(iGene[[p]], rohThresh, chromEnds)
}

save(rohList, file=destinationFile);
rohList

}
