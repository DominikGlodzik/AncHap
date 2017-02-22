scanSurrs <- function(iGene,chromEnds, destinationFile=NULL, lengthThresh, probands=length(1: length(iGene)), trimMargin, geneticMap=NULL, gDistT=NULL) {


noSubjects <- max(probands)
subjects <- probands

surrMat <- array(data.frame(), c(noSubjects,noSubjects))

for (p in probands) {
  if ((p %% 10)==1) {  print(paste("Finding regions of IBD sharing with proband ", p)) }
  probGen <- iGene[[p]]
  
  for (i in subjects) {  
    if (p<i) {

      iGenotype <- iGene[[i]]
      
      surrs<- extractIBDc(probGen, iGenotype, i, lengthThresh, chromEnds,geneticMap, gDistT)

      
      # cut the margins at the borders
      # left
      surrs$start <- surrs$start + trimMargin;
      # right
      surrs$end <- surrs$end - trimMargin;
      
      # keep only the sequences that after trimming are longer than 0
      surrLength <- surrs$end - surrs$start;
      surrs <- surrs[surrLength>1,];
      
      surrMat[[p,i]]<- surrs    
    }  
  }
}


if (!is.null(destinationFile)) {
  save(surrMat, file=destinationFile);
}
surrMat

}
