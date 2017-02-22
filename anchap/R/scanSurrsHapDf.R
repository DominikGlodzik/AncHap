#
#
#
#

scanSurrsHapDf <- function(result, lengthThresh, chromEnds, noHaps=2*length(result$phasingResults), dataFraction=-1, geneticMap=NULL, gDistT=NULL) {

  chroms <- 1:22;
  noInds <- length(result$phasingResults)
  
  allSnps <- 1:length(result$phasingResults[[1]]$h1);
  if (dataFraction>0) {
    allSnps <- 1:dataFraction
  }
  
    # make a list of haplotypes
    haplosM <- list()
    for (i in 1:noInds) {
      haplosM[[(i-1)*2+1]] <- as.integer(result$phasingResults[[i]]$h1)[allSnps];
      haplosM[[i*2]] <- as.integer(result$phasingResults[[i]]$h2)[allSnps];
    }# end loop over individuals
  
  surrMat <- array(data.frame(), c(noHaps,noHaps))

    for (p in (1:noHaps)) { # the outer loop over haplotypes
      print(paste("Haplotype ", p))
      h1<- haplosM[[p]]
      for (i in (1:noHaps)) { # the inner loop over haplotypes
        # ensuring haplotypes of a same individual are not compared twice, and that haplotypes of the same individual are not compared twice       

        # making sure they are haplotypes of separate individuals
        if (floor((p+1)/2)<floor((i+1)/2)) {
         h2<- haplosM[[i]]
         surrs <- extractIBDc(h1, h2, ceiling(i/2), lengthThresh, chromEnds, geneticMap, gDistT)
         if ((dim(surrs)[1])>0) {
           colnames(surrs) <- c("start", "end", "length", "id2");
           surrs<- cbind(surrs, data.frame(id1=ceiling(p/2), h1 =  ((p %% 2)==0)+1, h2 =  ((i %% 2)==0)+1))
         }
         
         # testing if the foudn regions indeed have no opposing homozygotes
         #for (s in 1:(dim(surrs)[1])) {
         #  aSurr <- surrs[s,];
         #  region <- aSurr$start:aSurr$end;          
         #  inc <- sum(abs(h1[region]-h2[region])==2)
         #  if (inc>0) browser();
         #}
        
         surrMat[[p,i]] <- surrs;   
        }  
      } # end inner loop over haplotypes
    }# end outer loop over haplotypes


# convert the surrogacy matrix to a data frame
surrMatHapDf <- surrMatToDf(surrMat)
# convert haplotype IDs to indiviudual IDs
surrMatHapDf$hid1 <-  (surrMatHapDf$id1+1)
surrMatHapDf$hid2 <-  (surrMatHapDf$id2+1)
surrMatHapDf$id1 <-  floor((surrMatHapDf$hid1+1)/2)
surrMatHapDf$id2 <-  floor((surrMatHapDf$hid2+1)/2)
 

surrMatHapDf
}# end function
