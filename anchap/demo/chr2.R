# load the library, and prepare the data set
 library(anchap)
 data(chr2Data)
 myData <- chr2Data

################################################################
## first stage of scan, and phasing
# work on a subset of the data set
individualsPhased <- 1:10

myResult <- phase( individualsPhased = individualsPhased,
                data = myData,
                ibdThresh = 500
                )

avgSurrs <- 2*sum(myResult$surrMatDf$end - myResult$surrMatDf$start)/(max(individualsPhased)*length(myData$iGene[[1]]))
cat(paste('On average there were ',avgSurrs,'surrogate parents at a locus. \n'))

################################################################
## scan dor surrogate parents from comparisons of partially phased haplotypes


lengthThreshHap <- 200
surrMatHapDf <- scanSurrsHapDf(result=myResult, lengthThresh=lengthThreshHap, chromEnds=myData$chromEnds)



# calculate average number of surrogate parents, from surrMatHapDf
avgSurrs <- 2*sum(surrMatHapDf$end - surrMatHapDf$start)/(max(individualsPhased)*length(myData$iGene[[1]]))

cat(paste('On average there were ',avgSurrs,'surrogate parents at a locus.\n'))

################################################################
## pick individuals for re-sequencing

pickedInds <- pickReseqIndsFast(surrMatDf=myResult$surrMatDf, data=myData)
