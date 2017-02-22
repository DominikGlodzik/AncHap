# picks the most optimal set of individual to resequence
# accounts for possibly overlapping surrogate regions
#
# Iterates until a number of individuals is found.iTaken[[ids[s]]]
# In each iteration, calculates the coverage of each individual.
# Picks the best individual, and marks is as such.
#
# noPick - number of individuals to resequence
pickReseqIndsFastHaps <- function (stitchingResult=NULL, noPick=noInds, plotPath=NULL, data, noInds=NULL, appFactor=10) {

surrMatDf$start <- ceiling(surrMatDf$start/appFactor)
surrMatDf$end <- ceiling(surrMatDf$end/appFactor)
  
calcSharingSlice <- function(noInds, noSnps, p, surrMarDfSlice, appFactor) {

  auxM <-  matrix(FALSE,noInds,noSnps);
        if ((dim(surrMarDfSlice)[1])>0) {
           for (s in 1:(dim(surrMarDfSlice)[1])) {        
             region <- surrMarDfSlice[s,'start']: (surrMarDfSlice[s,'end']-1)
               auxM[surrMarDfSlice[s,'id'],region] <- TRUE

          }
         }
        auxM[p,]<-1;
  auxM
}

calcImputationPotential <- function(p) {
    surrMarDfSlice <- sharingSliceList[[p]]
    noSharedSeqs<- dim(surrMarDfSlice)[1]

    imputSum <- 0
       toPrune <- integer(noSharedSeqs)
    
        if ( noSharedSeqs>0) {
                                        # whether a sharing sequence should be pruned
          toPrune <- integer(noSharedSeqs)
           for (s in 1:noSharedSeqs) {
             region <- surrMarDfSlice[s,'start']:(surrMarDfSlice[s,'end']-1)            
               newInput <- sum(coverageMatrix[surrMarDfSlice[s,'id'],region]==0)
               imputSum <- imputSum + newInput

             if (newInput==0) {
               toPrune[s]<-1;
               # remove the shared sequence
             }
             
          }
          
          # do the prining

          #print(paste("Pruned", sum(toPrune)))
         }

    newInput <- sum(coverageMatrix[p,]==0)
    imputSum <- imputSum + newInput

    result<-list()
    result$imputSum <-imputSum
    result$toPrune <- toPrune
    result
}



  
print(paste("app factor", appFactor))
  
if((exists('noInds'))&&(is.null(noInds))) {
  noInds <- max(max(surrMatDf$id1),(surrMatDf$id2))
  noPick <- noInds
}
 

  
noSnps <- ceiling(dim(data$genotypeMatrix)[2]/appFactor)
idNames <- rownames(data$genotypeMatrix)

   # ids of the individuals picked so far
  picked <<- vector()
  # useful coverage of the individuals picked so far, takes into account the regions of the genomes that could be imputed with individuals picked before
  pickedCov <<- vector()
  # raw coverage of the individuals picked so far
  pickedRawCov <<- vector()
  # sum of lengths of all genotypes - for normalising
  totalGenome <- noInds*(noSnps)

  # create the data structure iTaken
  # iTaken keeps track of regions that have surrogate coverage and therefore could be imputed

  surrCov <- vector()
  surrCovRaw <- vector() 
  nonSeqCoverage <- vector()
  # looping over number of people to reseqence

  sharingSliceList <- list()
    for (p in 1:noInds) {
      print(paste("Compiling sharing with individual",p))

        psSeqs <-  rbind(stitchingResult[[p]]$selectedSequencesH1, stitchingResult[[p]]$selectedSequencesH2)
        psSeqs <- psSeqs[psSeqs$id<=noInds,]
        psSeqs$start <- ceiling(psSeqs$start/appFactor)
        psSeqs$end <- ceiling(psSeqs$end/appFactor)
        sharingSliceList[[p]] <-psSeqs

        surrCovRaw[p] <- sum(sharingSliceList[[p]]$end - sharingSliceList[[p]]$start)

    }


  # a matrix of size noInds x noSnps
  # 1 if sequenced or has surrogate parents
  coverageMatrix <- matrix(FALSE,noInds,noSnps);


  for (i in 1:noPick) {
    print(paste("picking", i))  
    
    # loop over all individuals
    # check if not picked already
    # otherwise calculate their coverage

    topCoverage <- -1;
    
    for (p in 1:noInds) {
                                        # discard individuals already picked
      if (p %in% picked) {surrCov[p] <- -1}
                                        # calculate the 'coverage' score of the individual
      else {
        # calculate coverage of p
        pSurrParents <- sharingSliceList[[p]]

        #sharingSliceMatrix <- calcSharingSlice(noInds, noSnps, p, pSurrParents, appFactor)    
        #load(paste("../Data/results/pickInds/sharingSliceMatrix",p,'.RData',sep=''))
        #newImputeMatrix<- (sharingSliceMatrix==TRUE) & (coverageMatrix==FALSE) 

        #sum1 <- sum(sharingSliceList[[p]]$end - sharingSliceList[[p]]$start)
        r <- calcImputationPotential(p)        
        psCoverage <- r$imputSum
        #
        pSurrParents<-pSurrParents[!r$toPrune,]
        sharingSliceList[[p]] <- pSurrParents
        
        #sum2 <- sum(sharingSliceList[[p]]$end - sharingSliceList[[p]]$start)

        #print(paste("Pruned", sum2-sum1, "shared seqs"))
         
        if (psCoverage>topCoverage) {topCoverage <- psCoverage;}
        
        # calculate p's coverage, taking into account surrogate parent overalap
        surrCov[p] <- psCoverage;      
      } #end else     
    } # end loop over candidates for picking



    # remove the individuals already picked
    o <- setdiff(order(-surrCov),picked)

    
    
    # pick the top scoring individual in terms of coverage
    pickedNo <- o[1]

    winningMatrix <- calcSharingSlice(noInds, noSnps, pickedNo, sharingSliceList[[pickedNo]], appFactor)
    sharingSliceMatrix <-  (winningMatrix  | coverageMatrix)
    #load(paste("../Data/results/pickInds/sharingSliceMatrix",pickedNo,'.RData', sep=''))
       
    # get statistics
    picked[length(picked)+1] <- pickedNo
    pickedCov[length(pickedCov)+1] <- sum((winningMatrix>0)&(coverageMatrix==FALSE))
    
    #print(pickedCov[length(pickedCov)]/(noInds*noSnps))
    coverageMatrix <- sharingSliceMatrix
    #print(sum(coverageMatrix)/(noSnps*noInds))
    
        #sharingSliceMatrix <- calcSharingSlice(noInds, noSnps, p, pSurrParents, appFactor)    
        #load(paste("../Data/results/pickInds/sharingSliceMatrix",p,'.RData',sep=''))
        #newImputeMatrix<- (sharingSliceMatrix==TRUE) & (coverageMatrix==FALSE) 

        # calculate the raw coverage of p
        pickedRawCov[length(pickedRawCov)+1]  <-surrCovRaw[pickedNo]
    nonSeqCoverage[length(nonSeqCoverage)+1] <- (sum(coverageMatrix) - sum(coverageMatrix[picked,]))/((noInds-length(picked))*noSnps)

  }

  # formatting the output
  result <- data.frame(id=idNames[picked], cohortUsefulCoverage=pickedCov/(noInds*noSnps), cohortCoverage=pickedRawCov/(noInds*noSnps), nonSeqCoverage=nonSeqCoverage)

  # generating the plot summarising the yield as we resequence more individuals

  if (!is.null(plotPath)) {
    pdf(plotPath)
    plot(1:(dim(result)[1]), cumsum(result$cohortUsefulCoverage)/totalGenome, xlab="Number of individuals re-sequenced", ylab="Cumulative cohort coverage", ylim=c(0,1), )
    title(main = "Percentage of genotypes that could be imputed \n versus the number of individuals re-sequenced")
    dev.off();
  }

if (sum(result$cohortUsefulCoverage)!=1) {
  print('Warning! Things dont make sense')
  browser()
}


  result
}
