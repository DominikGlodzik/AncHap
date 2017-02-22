#
#
#
#

scanSurrsHap <- function(result, lengthThresh, chromEnds, noHaps=2*length(result$phasingResults), dataFraction=-1, geneticMap=NULL, gDistT=NULL) {

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
        
        if (floor((p+1)/2)<floor((i+1)/2)) {
         h2<- haplosM[[i]]
         surrs <- extractIBDc(h1, h2, ceiling(i/2), lengthThresh, chromEnds, geneticMap, gDistT)
         if ((dim(surrs)[1])>0) {
           colnames(surrs) <- c("start", "end", "length", "id2");
           surrs<- cbind(surrs, data.frame(id1=ceiling(p/2), h1 =  ((p %% 2)==0)+1, h2 =  ((i %% 2)==0)+1))
         }       
         surrMat[[p,i]] <- surrs;   
        }  
      } # end inner loop over haplotypes
    }# end outer loop over haplotypes



surrMat
}# end function