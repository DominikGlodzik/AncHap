computeSurrDensity<- function(surrMatDf, noSnps=NA) {

if (is.na(noSnps)) {noSnps<-max(surrMatDf$end);}
  
surrCount <- as.vector(matrix(0,noSnps,1))

for (s in 1:(dim(surrMatDf)[1])) {
  surrCount[surrMatDf$start[s]:(surrMatDf$end[s]-1)] <-  surrCount[surrMatDf$start[s]:(surrMatDf$end[s]-1)] + 1
}

surrCount
}
