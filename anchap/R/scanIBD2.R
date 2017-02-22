scanIBD2 <- function(iGene,chromEnds, destinationFile, lengthThresh, probands, trimMargin) {

noSubjects <- max(probands)
subjects <- probands

surrMat <- array(data.frame(), c(noSubjects,noSubjects))

for (p in probands) {
  print(paste("Proband ", p))
  probGen <- iGene[[p]]
  
  for (i in subjects) {  
    if (p<i) {

      iGenotype <- iGene[[i]]   
      surrs<- extractIBD2(probGen, iGenotype, i, lengthThresh, chromEnds)

      # cut the margins at the borders
      # left
      surrs$start <- surrs$start + trimMargin;
      # right
      surrs$end <- surrs$end - trimMargin;
      
      surrMat[[p,i]]<- surrs    
    }  
  }
}

save(surrMat, file=destinationFile);
surrMat

}
