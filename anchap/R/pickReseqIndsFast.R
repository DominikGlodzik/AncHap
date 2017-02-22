# picks the most optimal set of individual to resequence
# accounts for possibly overlapping surrogate regions
#
# Iterates until a number of individuals is found.iTaken[[ids[s]]]
# In each iteration, calculates the coverage of each individual.
# Picks the best individual, and marks is as such.
#
# noPick - number of individuals to resequence
pickReseqIndsFast <- function (surrMatDf, noPick=NULL, plotPath=NULL, data, noInds=NULL, appFactor=10) {

surrMatDf$start <- ceiling(surrMatDf$start/appFactor)
surrMatDf$end <- ceiling(surrMatDf$end/appFactor)

# calculate a matrix illustrationg imputation potential for an individual
calcSharingSlice <- function(noInds, noSnps, p, surrMarDfSlice, appFactor) {
  # matrix no individuals by no SNPs
  auxM <-  matrix(FALSE,noInds,noSnps);
        if ((dim(surrMarDfSlice)[1])>0) {
           for (s in 1:(dim(surrMarDfSlice)[1])) {        
             region <- surrMarDfSlice[s,'start']: (surrMarDfSlice[s,'end']-1)
             # set the areas that could be imputed to true
             if (surrMarDfSlice[s,'id1']==p) {
               auxM[surrMarDfSlice[s,'id2'],region] <- TRUE
             }
             else {
               auxM[surrMarDfSlice[s,'id1'],region] <- TRUE
             }
          }
         }
        # set the whole of own row to true
        auxM[p,]<-1;
  auxM
}

# calculate imputation potential for an individual
# iterate over all shared sequences
# mark some of them for pruning
calcImputationPotential <- function(p) {
  # sharing for one individual
  surrMarDfSlice <- sharingSliceList[[p]]
  noSharedSeqs<- dim(surrMarDfSlice)[1]
  # stores the imputation potential of one individual
  imputSum <- 0
  # add the own genotype, which could not be imputed so far
  imputSum <- sum(coverageMatrix[p,]==0)
  # a vector telling which sequences should be removed
  toPrune <- integer(noSharedSeqs)
  if ( noSharedSeqs>0) {
          # whether a sharing sequence should be pruned
          toPrune <- integer(noSharedSeqs)
          # loop over all shared sequences with this individual
           for (s in 1:noSharedSeqs) {        
             region <- surrMarDfSlice[s,'start']:(surrMarDfSlice[s,'end']-1)            
             if (surrMarDfSlice[s,'id1']==p) {
               # how much unimputed markers there are so far in the coverage matrix
               newInput <- sum(coverageMatrix[surrMarDfSlice[s,'id2'],region]==0)
               imputSum <- imputSum + newInput
             }
             else {
               newInput <-  sum(coverageMatrix[surrMarDfSlice[s,'id1'],region]==0)
               imputSum <- imputSum + newInput
             }
             if (newInput==0) {
               toPrune[s]<-1;
               # remove the shared sequence
             }            
          }
         }
    result<-list()
    result$imputSum <-imputSum
    result$toPrune <- toPrune
    result
}

print(paste("app factor", appFactor))

# if number of candidates for resequencing hadn't been set
if(is.null(noInds)) {
  noInds <- max(max(surrMatDf$id1),max(surrMatDf$id2))
}

# by default pick all individuals
if(is.null(noPick)) {
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
  # looping over number of people to reseqence

sharingSliceList <- list()
# preparing sharingSliceList
for (p in 1:noInds) {
  print(paste("Compiling sharing with individual",p))    
                                        # take a subset of the surrogacy matrix
  sp1 <- (surrMatDf$id1==p) & (surrMatDf$id2<=noInds)
  sp2 <- (surrMatDf$id2==p) & (surrMatDf$id1<=noInds)
  pSurrParents <- surrMatDf[sp1|sp2,]
  sharingSliceList[[p]] <- pSurrParents
  surrCovRaw[p] <- sum(pSurrParents$end - pSurrParents$start)
}



# a matrix of size noInds x noSnps
# 1 if sequenced or has surrogate parents
coverageMatrix <- matrix(FALSE,noInds,noSnps);

  # loop until a desired number of individuals is picked
  for (i in 1:noPick) {
    print(paste("picking", i))  
    
    # loop over all individuals
    # check if not picked already
    # otherwise calculate their coverage

    topCoverage <- -1;

    # loop over all individuals
    # calculate their coverage
    # surrCov
    for (p in 1:noInds) {
                                        # discard individuals already picked
      if (p %in% picked) {surrCov[p] <- -1}
                                        # calculate the 'coverage' score of the individual
      else {
        # calculate coverage of p
        pSurrParents <- sharingSliceList[[p]]

        r <- calcImputationPotential(p)        
        psCoverage <- r$imputSum
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
    # or between the winning sharing matrix and the coverage matrix
    # what could be imputed with previous + current individual
    sharingSliceMatrix <-  (winningMatrix  | coverageMatrix)
        
    # get statistics
    picked[length(picked)+1] <- pickedNo
    # calculate useful coverage based on the matrix
    pickedCov[length(pickedCov)+1] <- sum((winningMatrix>0)&(coverageMatrix==FALSE))

    # pickedCov[length(pickedCov)] should equal surrCov[pickedNo]
    #browser()
    
    # update the coverage matrix
    coverageMatrix <- sharingSliceMatrix

    # calculate the raw coverage of p
    pickedRawCov[length(pickedRawCov)+1]  <-surrCovRaw[pickedNo]

  }

  # formatting the output
  result <- data.frame(id=idNames[picked], cohortUsefulCoverage=pickedCov/(noInds*noSnps), cohortCoverage=pickedRawCov/(noInds*noSnps))

  # generating the plot summarising the yield as we resequence more individuals

  if (!is.null(plotPath)) {
    pdf(plotPath)
    plot(1:(dim(result)[1]), cumsum(result$cohortUsefulCoverage), xlab="Number of individuals re-sequenced", ylab="Cumulative cohort coverage", ylim=c(0,1), )
    title(main = "Percentage of genotypes that could be imputed \n versus the number of individuals re-sequenced")
    dev.off();
  }

if ((noInds==noPick)&(sum(result$cohortUsefulCoverage)!=1)) {
  print('Warning! Things dont make sense')
  browser()
}


  result
}
