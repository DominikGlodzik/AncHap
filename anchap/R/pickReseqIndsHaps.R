# picks the most optimal set of individual to resequence
# accounts for possibly overlapping surrogate regions
#
# Iterates until a number of individuals is found.iTaken[[ids[s]]]
# In each iteration, calculates the coverage of each individual.
# Picks the best individual, and marks is as such.
#
# noPick - number of individuals to resequence
pickReseqIndsHaps <- function (hapSurrMatDf=NULL, noPick=noInds, plotPath=NULL, data, noInds=NULL, appFactor=10) {

print(paste("app factor", appFactor))

hapSurrMatDf$start <- ceiling(hapSurrMatDf$start/appFactor)
hapSurrMatDf$end <- ceiling(hapSurrMatDf$end/appFactor)

calcSharingSlice <- function(noInds, noSnps, p, surrMarDfSlice, appFactor) {
  # now the dimension is 2*noInds x noSnps
  auxM <-  matrix(integer(2*noInds*noSnps),2*noInds,noSnps);
        if ((dim(surrMarDfSlice)[1])>0) {
           for (s in 1:(dim(surrMarDfSlice)[1])) {        
             region <- surrMarDfSlice[s,'start']: (surrMarDfSlice[s,'end']-1)

           if (surrMarDfSlice[s,'id1']==p) {
             auxM[surrMarDfSlice[s,'hid2'],region] <- 1
           } else {
             auxM[surrMarDfSlice[s,'hid1'],region] <- 1             
           }

          }
         }
        auxM[2*p-1,]<-1;
        auxM[2*p,]<-1;
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

           if (surrMarDfSlice[s,'id1']==p) {
             newInput <- sum(coverageMatrix[surrMarDfSlice[s,'hid2'],region]==0)
           } else {
             newInput <- sum(coverageMatrix[surrMarDfSlice[s,'hid1'],region]==0)  
           }

             imputSum <- imputSum + newInput

             if (newInput==0) {
               toPrune[s]<-1;
               # remove the shared sequence
             }             
          }
         }

    # add the own, non-sequenced 
    newInput <- sum(coverageMatrix[2*p-1,]==0)+sum(coverageMatrix[2*p,]==0)
    imputSum <- imputSum + newInput

    result<-list()
    result$imputSum <-imputSum
    result$toPrune <- toPrune
    result
}

  
if((exists('noInds'))&&(is.null(noInds))) {
  noInds <- max(max(surrMatDf$id1),(surrMatDf$id2))
  noPick <- noInds
}
 
noSnps <- ceiling(length(data$iGene[[1]])/appFactor)
idNames <- data$idNames

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

        idSeqs1 <- (hapSurrMatDf$id1==p); 
        idSeqs2 <- (hapSurrMatDf$id2==p);

        psSeqs <-  hapSurrMatDf[(idSeqs1 | idSeqs2), ]
        id <- vector()
        id[psSeqs$id1==p] <- psSeqs$id2[psSeqs$id1==p]
        id[psSeqs$id2==p] <- psSeqs$id1[psSeqs$id2==p]
        psSeqs$id <- id;
      
        psSeqs <- psSeqs[psSeqs$id<=noInds,]
        sharingSliceList[[p]] <-psSeqs

        surrCovRaw[p] <- sum(sharingSliceList[[p]]$end - sharingSliceList[[p]]$start)
    }
  cat("compiled the sharing \n")
  browser()


  # a matrix of size noInds x noSnps
  # 1 if sequenced or has surrogate parents
  coverageMatrix <- matrix(integer(2*noInds*noSnps),2*noInds,noSnps);


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

        #sum1 <- sum(sharingSliceList[[p]]$end - sharingSliceList[[p]]$start)
        r <- calcImputationPotential(p)        
        psCoverage <- r$imputSum
        #
        pSurrParents<-pSurrParents[!r$toPrune,]
        sharingSliceList[[p]] <- pSurrParents
        
         
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

if ((sum(result$cohortUsefulCoverage)/2)!=1) {
  print('Warning! Things dont make sense')
  browser()
}


  result
}
