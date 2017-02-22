
# finds all ibd sequences for an individual, using the surrogacy matrix
getAllSurrogatesOf <- function (surrMat, proband, margin) {

# recovers ibs sequences between two individuals, from the surrogacy matrix
getIbs <- function(i1, i2) {if (i1<i2) {ibsSeqs <- surrMat[[i1,i2]]}  else if (i1>i2) {ibsSeqs <- surrMat[[i2, i1]]}; ibsSeqs}
  
  # requires surrMat
  # requires
  # requires function filterSurrogates
  noInds <- dim(surrMat)[1]
  allSurrogates <- data.frame();
  noInds <- dim(surrMat)[1]
  for (i in 1:noInds) {
   if(i!=proband) { 
    surrogates<-getIbs(i,proband)

    if (!is.null(surrogates)) {
    
      noSurrogates <- dim(surrogates)[1]
      if (noSurrogates>0) {
        batch <- data.frame(id=i, start=surrogates$start, end=surrogates$end)
                                        #print(batch)
        allSurrogates <- rbind(allSurrogates,batch)
      }
    }

    
   }   
  }

  if ((dim(allSurrogates)[1])>0) {
  allSurrogates<-allSurrogates[allSurrogates$end>0,]
  is.vector(allSurrogates$start)
  allSurrogates.sorted <- allSurrogates[order(allSurrogates$start) , ]
  allSurrogates <-allSurrogates.sorted

  allSurrogates$start <- allSurrogates$start+margin
  allSurrogates$end <- allSurrogates$end-margin
}
  allSurrogates
}




