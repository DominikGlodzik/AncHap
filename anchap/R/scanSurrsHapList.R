scanSurrsHapList <- function(matHapList, patHapList, surrHapThreshold, gDistT, allData) {

noInds <- length(matHapList)
noHaps <- 2*length(matHapList)
noSnps <- length(matHapList[[1]])

  
##########
# create haplotype matrix
# take the result
# 2N x noSnps
    # make a list of haplotypes
    haplosM <- matrix(integer(noHaps*noSnps),(noHaps), noSnps )
    for (i in 1:noInds) {
      haplosM[((i-1)*2+1),] <- matHapList[[i]];
      haplosM[(i*2), ] <- patHapList[[i]];
    }# end loop over individuals
  
##########

##########
# detect sharing from the haplotype matrix
  surrMat <- array(data.frame(), c(noHaps,noHaps))
    for (p in (1:noHaps)) { # the outer loop over haplotypes
      print(paste("Haplotype ", p))
      h1<- haplosM[p,]
      for (i in (1:noHaps)) { # the inner loop over haplotypes
        
        # ensuring haplotypes of a same individual are not compared twice, and that haplotypes of the same individual are not compared twice          
        if (floor((p+1)/2)<floor((i+1)/2)) {
         h2<- haplosM[i,]
          surrs <- extractIBDc(h1, h2, ceiling(i/2), surrHapThreshold, allData$chromEnds, allData$geneticSnpMap, gDistT)

         if ((dim(surrs)[1])>0) {
           colnames(surrs) <- c("start", "end", "length", "id2");
           surrs<- cbind(surrs, data.frame(id1=ceiling(p/2), h1 =  ((p %% 2)==0)+1, h2 =  ((i %% 2)==0)+1))
         }
         surrMat[[p,i]] <- surrs;   
        }  
      } # end inner loop over haplotypes
    }# end outer loop over haplotypes
# convert the surrogacy matrix to a data frame

surrMatHapDf <- surrMatToDf(surrMat)
names(surrMatHapDf) <- c("start","end","hid1","hid2")
# wrt matHapList and patHapList
surrMatHapDf$id1 <- floor((surrMatHapDf$hid1+1)/2) 
surrMatHapDf$id2 <- floor((surrMatHapDf$hid2+1)/2)  




print('The scan for haplotype sharing complete')
##########
surrMatHapDf

}
