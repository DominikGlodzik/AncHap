# picks the most optimal set of individual to resequence
# accounts for possibly overlapping surrogate regions
#
# Iterates until a number of individuals is found.iTaken[[ids[s]]]
# In each iteration, calculates the coverage of each individual.
# Picks the best individual, and marks is as such.
#
# noPick - number of individuals to resequence
pickReseqInds <- function (surrMatDf, noPick=noInds, plotPath=NULL, data, noInds=NULL) {

if((exists('noInds'))&&(is.null(noInds))) {
  noInds <- max(max(surrMatDf$id1),(surrMatDf$id2))
  noPick <- noInds
}

  

  
noSnps <- dim(data$genotypeMatrix)[2]
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
  # looping over number of people to reseqence

  # a matrix of size noInds x noSnps
  # 1 if sequenced or has surrogate parents
  coverageMatrix <- matrix(0,noInds,noSnps);   

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

        # remember the matrix with individuals picked so far
        coverageMatrixAux <- coverageMatrix;

        # ps sequence will be fully known
        coverageMatrixAux[p,] <- 1;
        
        # find individuals sharing with p
        sp1 <- (surrMatDf$id1==p) & (surrMatDf$id2<=noInds)
        sp2 <- (surrMatDf$id2==p) & (surrMatDf$id1<=noInds)
        
        pSurrParents <- surrMatDf[sp1|sp2,]

        # calculate the raw coverage of p
        surrCovRaw[p] <- noSnps + sum(pSurrParents$end-pSurrParents$start)
          
        # loop over surrogate parents of p
        if ((dim(pSurrParents)[1])>0) {
           for (s in 1:(dim(pSurrParents)[1])) {
          
          region <- pSurrParents[s,'start']: (pSurrParents[s,'end']-1)
          if (pSurrParents[s,'id1']==p) {
            if (pSurrParents[s,'id2']<=noInds) {
              coverageMatrixAux[pSurrParents[s,'id2'],region] <- 1
            }
          }
          else {
            if (pSurrParents[s,'id1']<=noInds) {
              coverageMatrixAux[pSurrParents[s,'id1'],region] <- 1
            }
          }

         
          
        } # end loop over p's shared sequences
        } # end if 

        psCoverage <- sum(sum(coverageMatrixAux-coverageMatrix));  
        if (psCoverage>topCoverage) {topCoverage <- psCoverage; topCoverageMatrix <- coverageMatrixAux}
        
        # calculate p's coverage, taking into account surrogate parent overalap
        surrCov[p] <- psCoverage; 
    

 
      } #end else     
    } # end loop over candidates for picking

    coverageMatrix <- topCoverageMatrix

    # remove the individuals already picked
    o <- setdiff(order(-surrCov),picked)

    # pick the top scoring individual in terms of coverage
    pickedNo <- o[1]

    
    
    # get statistics
    picked[length(picked)+1] <- pickedNo
    pickedCov[length(pickedCov)+1] <- surrCov[pickedNo]
    pickedRawCov[length(pickedRawCov)+1] <- surrCovRaw[pickedNo] 
  }

  # formatting the output
  result <- data.frame(id=idNames[picked], cohortUsefulCoverage=pickedCov/(noInds*noSnps), cohortCoverage=pickedRawCov/(noInds*noSnps))

  # generating the plot summarising the yield as we resequence more individuals

  if (!is.null(plotPath)) {
    pdf(plotPath)
    plot(1:(dim(result)[1]), cumsum(result$cohortUsefulCoverage)/totalGenome, xlab="Number of individuals re-sequenced", ylab="Cumulative cohort coverage", ylim=c(0,1), )
    title(main = "Percentage of genotypes that could be imputed \n versus the number of individuals re-sequenced")
    dev.off();
  }

  result
}
