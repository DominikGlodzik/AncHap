analyseTrioHaplotyped <- function (iProb, iMo, iFa, phasingsFolder, plotFolder, genomeLength, chromEnds) {

                                        # keeps track of there the sharing between the proband and its parents occurs
  genomeSharing <<- integer(genomeLength)

  # get error rate from two pairs of haplotypes, which sould be the same!
  getErrorRate <- function(id1, id2) {

   allResults <- data.frame()
   blockLength <- 20;

    lastSnp <- 0
    for (c in 1:22) {
      if (chromEnds[c]>0) {

      region <- (lastSnp+1) : chromEnds[c] 

      # how many blocks there are in a region
      noBlocks <- floor(length(region)/blockLength)

      # the region spanning the chromosomes
      chromStart <- lastSnp + 1
      chromEnd <- lastSnp + noBlocks*blockLength  

      
      # reshape the haplotypes
      matrix1h1 <- t(matrix((id1[[1]])[chromStart:chromEnd], ncol=noBlocks ))
      matrix1h2 <- t(matrix((id1[[2]])[chromStart:chromEnd], ncol=noBlocks ))
      matrix2h1 <- t(matrix((id2[[1]])[chromStart:chromEnd], ncol=noBlocks ))
      matrix2h2 <- t(matrix((id2[[2]])[chromStart:chromEnd], ncol=noBlocks ))

      homoLoci1 <- (matrix1h1>0) & ((matrix1h1 - matrix1h2)== 0)
      homoLoci2 <- (matrix2h1>0) & ((matrix2h1 - matrix2h2 )== 0)
     
      heteroLoci1 <- rowSums(matrix1h1 != matrix1h2)
      heteroLoci2 <- rowSums(matrix2h1 != matrix2h2)
      # heterozygous loci in real and phased haplotypes
      importantLoci <- rowSums((matrix1h1 != matrix1h2) & (matrix2h1 != matrix2h2))


      # opposing homozygotes
      oppos1121 <-  (abs(matrix1h1-matrix2h1)==2)
      oppos1122 <- (abs(matrix1h1-matrix2h2)==2)
      oppos1221 <-  (abs(matrix1h2-matrix2h1)==2)
      oppos1222 <- (abs(matrix1h2-matrix2h2)==2)

      opp11 <- rowSums(oppos1121)
      opp12 <- rowSums(oppos1122)
      
      opp21 <- rowSums(oppos1221) # should be the same as 12
      opp22 <- rowSums(oppos1222) # should be the same as 11

      # loop over the blocks
      for (b in 1:noBlocks) {

        # check which of 11-22 12-21 is matching better

        blockOppos11 <- opp11[b] + opp22[b]
        blockOppos12 <- opp12[b] + opp21[b]

#  opp11=opp11[b], opp22=opp22[b], opp12=opp12[b], opp21=opp21[b],
        
        aResult <- data.frame(blockOppos11=blockOppos11, blockOppos12=blockOppos12, heteroLoci1=heteroLoci1[b], heteroLoci2=heteroLoci2[b], importantLocus=importantLoci[b])
        # calculate the error rate best on the better matching
        
        allResults <- rbind(aResult,allResults)
      }   

      # slice the region into 100 SNP bits each
      lastSnp <- chromEnds[c]
    }# end if
  }
   
   # heterozygous loci that had been phased in both haplotypes
   impotantThresh <- 3
   significantResults <- allResults[allResults$importantLocus>impotantThresh,]

   option11 <- ((significantResults$blockOppos11) < (significantResults$blockOppos12))
   option12 <- !option11

   allErrors <- sum(significantResults[option11, 'blockOppos11']) + sum(significantResults[option12, 'blockOppos12'])


   errorResult <- list()
   errorResult[['errorCount']] <- allErrors
   errorResult[['errorProneLoci']] <- sum(significantResults$importantLocus)

   if (errorResult[['errorCount']] > errorResult[['errorProneLoci']]) {browser()}
   
   errorResult
 }


  
  compareHaplos <- function (id1, id2, iParent, height) {
    # compares all pairs of proband-parent hapotypes, and plots the result of the comparison
  
    plotDisagree <- function (pnts, height, c) {
      ys <- integer(length(pnts)) + height
      points(pnts,  ys, lwd=0.05, pch="x", cex=0.2, col=c) 
    }

    
  plotSurrs <- function (surrs, height, c) {
      ys <- c(height, height)
      
      if ((dim(surrs)[1])>0) {
        for (s in 1:(dim(surrs)[1])) {
          
          xs <- c(surrs[s,'start'], surrs[s,'end'])
          lines(xs,ys, col=c)
          iTaken[surrs[s,'start']:surrs[s,'end']] <<- 1

          genomeSharing[surrs[s,'start']:(surrs[s,'end']-1)]<<- genomeSharing[surrs[s,'start']:(surrs[s,'end']-1)] + 1;

        }
      }
    }

    markSurrs <- function(surrs) {
      matching <- integer(genomeLength) 
      for (s in 1:(dim(surrs)[1])) {
        matching[surrs[s,'start']:(surrs[s,'end']-1) ] <- matching[surrs[s,'start']:(surrs[s,'end']-1) ] + 1
      }
      matching      
    }

    
    plotMatching <- function(matching,height) {

      blockSize <- 10
      noBlocks <- floor(length(matching)/blockSize)

      matchingMatrix <- matrix(matching[1:(noBlocks*blockSize)], nrow=blockSize)
      means <- colMeans(matchingMatrix)

      indx <- (0:(noBlocks-1)) * blockSize + 1

      lines(indx  ,0.25*means+height, lwd=0.1, pch=".", cex=0.1)
    }


    r <- getErrorRate(id1, id2)
    
  }
  
  countOppHomoR <- function(s1,s2) {
    if (length(s1)!=length(s2)) {
      stop();
    }
    zeros1 <- (s1==0)
    zeros2 <- (s2==0)
    oppHomos <- (abs(s1-s2)==2) & (!zeros1) & (!zeros2)
    sum(oppHomos)
  }
  

  
  invertAllele <- function(haplotype) {

    haplotypeOut <- integer(length(haplotype))
    
    haplotypeOut[haplotype==1] <- 3
    haplotypeOut[haplotype==3] <- 1
    haplotypeOut
  }
  
  getRealHaplos <- function(iProb, iMo, iFa) {

    realHapMo <- integer(max(allData$chromEnds))
    realHapFa <- integer(max(allData$chromEnds))
    
    gProb <- (allData$iGene[[iProb]])[1 :max(allData$chromEnds)]
    gMo <- (allData$iGene[[iMo]])[1 : (max(allData$chromEnds))]
    gFa <- (allData$iGene[[iFa]])[1 : (max(allData$chromEnds))]

    probHetLoci <- (gProb==2)
    probHomLoci <- (gProb==1)|(gProb==3)
    moHomoLoci <- (gMo==1)|(gMo==3)
    faHomoLoci <- (gFa==1)|(gFa==3)   

    #
    realHapMo[probHomLoci] <- gProb[probHomLoci]
    realHapFa[probHomLoci]  <- realHapMo[probHomLoci] 

    
    # phasing using the mother's genotype
    moPhasingLoci <- probHetLoci & moHomoLoci
    realHapMo[moPhasingLoci] <- gMo[moPhasingLoci]
    realHapFa[moPhasingLoci] <- invertAllele(gMo[moPhasingLoci])
    
    # phasing using the father's genotype
    faPhasingLoci <- probHetLoci & faHomoLoci
    realHapFa[faPhasingLoci] <- gFa[faPhasingLoci]
    realHapMo[faPhasingLoci] <- invertAllele(gFa[faPhasingLoci])

    r <- list()
    r[[1]] <- realHapFa
    r[[2]] <- realHapMo
    
    if (sum(realHapFa==2)>0) {
      print("Warning! heterozygotes in the haplotypes")
      browser()
    }
    
    r
  
  }

  fileName <- paste(plotFolder, "trio","_", iProb ,'.pdf', sep="")
  pdf(fileName)
  
  plot(c(1,max(chromEnds)), c(0 ,11),'n', ylab='', xlab='', bty='o',  yaxt="n")

                                        # load the haplotypes of the proband
  prob <- list()
  filename <- paste(phasingsFolder, "results",iProb, ".Rdata", sep="")
  load(filename)
  prob[[1]] <- h1[1:max(chromEnds)]
  prob[[2]] <- h2[1:max(chromEnds)]
  probPhased <-  1 -(sum(h1==as.raw(0)))/genomeLength
  print("Surrogate-parents phasing rate")
  surrParPhased <- 1 - sum(prob[[1]]==0)/length(prob[[1]])
  print(surrParPhased)

  # get the real haplotypes of the proband
  trioPhasingR <- getRealHaplos(iProb, iMo, iFa)
  
  print("Real-parents phasing rate")
  realParPhased <- 1 - sum(trioPhasingR[[1]]==0)/length(trioPhasingR[[1]])
  print(realParPhased)

  comparisonResult <- compareHaplos(prob,trioPhasingR, iFa, 7)
  
  dev.off()
  print(paste("Generated plot ", fileName))


  result <- list();

  result <- data.frame(errorCount=comparisonResult[['errorCount']], errorProneLoci=comparisonResult[['errorProneLoci']], surrParPhased=surrParPhased, realParPhased =realParPhased )
  

  result
}
