phaseGenDist <- function(data, outputPaths, individualsPhased = ( 1:(dim(data$genotypeMatrix)[1])), genThresh=5, onlySurrs=FALSE,
                  excludeChildren = TRUE, excludeParents = FALSE,  trimMargin=100, overlapT=10, matchT=0.02,  lengthThresh= 200, trioAnalysis=FALSE, surrMatFile=NULL) {
                                        # whether the data is to be read,
                                        # set to FALSE if you are re-running the script, and if all the data is in memory already
                                        # whether the iGene structure is to be recreated, or reloaded from the disk



print("finished loading data")
# end loading the data

# extracting the surrogacy matrix
cat("extracting the surrogacy matrix ...\n")

surrMat <- scanSurrs(iGene=data$iGene,
                     chromEnds=data$chromEnds,
                     lengthThresh= lengthThresh,
                     destinationFile= outputPaths$surrmatFile,
                     probands=individualsPhased,
                     trimMargin,
                     geneticMap = data$geneticSnpMap,
                     gDistT=genThresh
                    )

if (!is.null(surrMatFile)) {
  save(surrMat, file=surrMatFile)
}
surrMatDf <- surrMatToDf(surrMat)

# stitching and phasing is done independently for each individual
# loop over individuals

phasingResults <- list()
stitchingResults <- list()
trioResults <- list()

if (!onlySurrs) {
  for (individualPhased in individualsPhased) {

    
    
    cat(paste("Stitching and phasing individual", individualPhased))
    
                                        # stitching
    cat(": stitching ...")
    
                                        # whether children are to be ommited when phasing parents
    
    excludedParents <- getExcludedIds(individualPhased,excludeParents, excludeChildren, data$pedigree, data$orcaIds );
                                        # don't exclude any parents
    
    stitchingResult <- doStitching3(individualPhased=individualPhased,
                                    surrMat=surrMat,
                                    iGene=data$iGene,
                                    chromList=data$chromList,
                                    overlapT = overlapT,
                                    matchT=matchT,
                                    map=data$map,
                                    chromEnds=data$chromEnds,
                                    excludedParents=excludedParents)
    
    stitchingResults[[individualPhased]] <- stitchingResult 
    stitches <<- data.frame()
    
                                        # phasing
    cat("phasing ...\n")
    phasingResult <- doPhasing(
                               probandId =  individualPhased,
                               iGene = data$iGene,
                               chromList = data$chromList,
                               selectedSequencesH1 = stitchingResult$selectedSequencesH1,
                               selectedSequencesH2 = stitchingResult$selectedSequencesH2,
                               alsoPlot = TRUE,
                               plotFolder =  outputPaths$plotFolder,
                               outputFolder = outputPaths$phasingOutputFolder,
                               cumEnds = data$cumEnds,
                               map = data$map,
                               chromEnds = data$chromEnds,
                               plotFileExt=genThresh)

    phasingResults[[individualPhased]] <- phasingResult   
  
                                        # end loop over individuals
  
                                        # evaluation: comparing the haplotypes of trios in the cohort

  orcaIdsInStudy <- data$orcaIds;
  orcaIdsInStudy[setdiff(1:length(data$orcaIds), individualsPhased)] <- NA
  
} # end loop over individualPhased

# trio analysis, when pedigree is available
trioResults <- NULL;
if ((trioAnalysis) && (!is.null(data$pedigree))) {
                                        # re-evaluate the trio phasing rate
                                        # will have to prepare the Hm matrix first
    noSnps <- length(data$iGene[[1]]); startSnp <- 1; endSnp <- noSnps;  noInds <- max(individualsPhased)
    
    Hm <- matrix(integer(2*noInds*noSnps), nrow =2*noInds, ncol=noSnps)
                                        # build the input matrices
    for (i in 1:noInds) {
      Hm[(i-1)*2+1,] <- as.integer(phasingResults[[i]]$h1)[startSnp:endSnp];
      Hm[i*2,] <- as.integer(phasingResults[[i]]$h2)[startSnp:endSnp];
    }
                                        # this only works with data on one chromosome!
    trioResults <- analyseTriosSimple(data, Hm)
}
  
  
# calculate mean number of surrogate parents
noAlignedSurrogates <- vector();
noSurrogates <- vector()
phasingYield <- vector();
phasingYieldHet <- vector();
betweenGameteInconsistencies <- vector()
withinGameteInconsistencies <- vector()
inconsistenciesPerIndividual <- vector()
homoSurrAllelesPerIndividual <- vector()
noSnps <- length(data$iGene[[1]]);

#Loop over individuals
for (i in individualsPhased) {

  
  # calculate number of surrogate parents for all individuals
  # aligned parents
  allISurrs <- rbind(stitchingResults[[i]]$selectedSequencesH1, stitchingResults[[i]]$selectedSequencesH2)
  noAlignedSurrogates[i] <- sum(allISurrs$end - allISurrs$start)/noSnps;
  # all parents
  indParents <- ((surrMatDf$id1==i) | (surrMatDf$id2==i))
  # average number of surrogate parents for the individual, across many SNPs
  noSurrogates[i] <- sum(surrMatDf[indParents,'end'] - surrMatDf[indParents, 'start'])/noSnps;
  
  phasingYield[i] <- sum(phasingResults[[i]]$h1>0)/noSnps;
  phasingYieldHet[i] <-  sum((phasingResults[[i]]$h1>0) & (data$iGene[[i]]==2))/sum(data$iGene[[i]]==2);
  phasedLoci <- (phasingResults[[i]]$h1>0) & ( data$iGene[[i]]==2);

  ########################################################################
  ## INCONSISTENCIES
  ######################################################################## 
  # between - phased and haplotypes homozygous at a heterozygous locus
  betweenGameteInconsistencies[i] <- sum(phasedLoci &  (phasingResults[[i]]$h1==phasingResults[[i]]$h2) )/sum(phasedLoci);
  # within - sum of inconsistencies on each haplotype
  withinGameteInconsistencies[i] <- (sum(phasingResults[[i]]$i1) + sum(phasingResults[[i]]$i2));


  inconsistenciesPerIndividual[i] <- sum(phasingResults[[i]]$inconsistencies)
  homoSurrAllelesPerIndividual[i] <- sum(phasingResults[[i]]$parentHomozygotes)
} # end loop over individuals


# count of inconsistencies on bith haplotypes, normalized by region size and total number of aligned surrogate parents for all individuals

# count all surrogate alleles  - total length of all surrogate parents, multiplied by 2 (sharing goes two ways)
allAlleleSurrs <- 2*sum(surrMatDf$end - surrMatDf$start);


withinGameteInconsistencies <- sum(withinGameteInconsistencies,na.rm = TRUE)/ (allAlleleSurrs)

# count the surrogate parents across the genome
surrDensity <- NA
if ((dim(surrMatDf)[1])>0)
{surrDensity <- computeSurrDensity(surrMatDf, noSnps)}

resultSummary <- list()

resultSummary$stats <- data.frame(id = individualsPhased,
                             noSurrogatesAtLocus = noSurrogates[individualsPhased],
                             noSurrogatesAtLocusAligned = noAlignedSurrogates[individualsPhased],
                             avgPhasingYield = phasingYield[individualsPhased],
                                  phasingYieldHet = phasingYieldHet[individualsPhased],
                             inconsistenciesPerIndividual = inconsistenciesPerIndividual[individualsPhased],
                             homoSurrAllelesPerIndividual = homoSurrAllelesPerIndividual[individualsPhased]  
)
resultSummary$withinGameteInconsistencies <- withinGameteInconsistencies


result <- list(surrMat=surrMat, surrMatDf=surrMatDf, stitchingResults=stitchingResults, phasingResults=phasingResults, trioResults=trioResults, resultSummary = resultSummary,  chromEnds=data$chromEnds, surrDensity=surrDensity)
} else { # only haplotype sharing required

noSnps <- length(data$iGene[[1]]);
surrDensity <- computeSurrDensity(surrMatDf, noSnps)
result <- list(surrMat=surrMat, surrMatDf=surrMatDf,surrDensity=surrDensity)
} # end compiling the results for IBD sharing detection only
                                   

result
} # end phaseGenDist
