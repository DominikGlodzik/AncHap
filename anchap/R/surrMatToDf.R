# converts surrogacy matrix - a two dimensional list,
# into a two-dimensional array

surrMatToDf <- function(surrMat) {

noIds <- dim(surrMat)[1]
noHaps <- noIds;



noPairs <- 0;
for (i in 1:noIds) {
  for (j in 1:noIds) {
    noPairsIJ <- dim(surrMat[[i,j]])[1];
    if (length(noPairsIJ)>0) {
      noPairs <- noPairs + noPairsIJ
    }
  }
}

startV <- integer(noPairs)
endV <- integer(noPairs)
id1V <- integer(noPairs)
id2V <- integer(noPairs)

pointerInd <- 0;
for (i in 1:(noHaps)) {
  for (j in 1:(noHaps)) {
    noPairsIJ <- dim(surrMat[[i,j]])[1];
    if ((length(noPairsIJ)>0)&&(noPairsIJ>0)) {
      startV[(pointerInd+1) : (pointerInd+noPairsIJ)] <- surrMat[[i,j]]$start;
      endV[(pointerInd+1) : (pointerInd+noPairsIJ)] <- surrMat[[i,j]]$end;
      #id1V[(pointerInd+1) : (pointerInd+noPairsIJ)] <- surrMatHap[[i,j]]$id1;
      id1V[(pointerInd+1) : (pointerInd+noPairsIJ)] <- i;
      #id2V[(pointerInd+1) : (pointerInd+noPairsIJ)] <- surrMatHap[[i,j]]$id2;      
      id2V[(pointerInd+1) : (pointerInd+noPairsIJ)] <- j;
      # increment the pointer
      pointerInd<- pointerInd +noPairsIJ;      
    }
  }
}


data.frame(start=startV, end = endV, id1=id1V, id2=id2V)
}





