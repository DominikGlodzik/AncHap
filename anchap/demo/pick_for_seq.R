
# load data
##################################
data(toyData)
genotypeData <- toyData
rm(toyData)
##################################


# extract IBD sharing between individuals
##################################
if ('geneticSnpMapFile' %in% names(genotypeData)) {
# the genetic map IS available
  phasingResult <- phaseGenDist(data=genotypeData, # use the data read in earlier
                              outputPaths= NULL, # do not make any plots 
                              ibdThresh=3, # IBD segments have to be at least 3cM long
                              trimMargin= 100, # 100 markers trimmed at each side of an IBD segment
                              onlySurrs=TRUE # do not attempt phasing and the second round of IBD detection
                              )
} else {
# the genetic map is NOT available
  phasingResult <- phase(data=genotypeData, # use the data read in earlier
                       outputPaths= NULL, # do not make any plots 
                       ibdThresh=500, # IBD segments have to be at least 500 consecutive markers long
                       trimMargin= 100, # 100 markers trimmed at each side of an IBD segment
                       onlySurrs=TRUE # do not attempt phasing and the second round of IBD detection
                       )   
}
# check the coverage of IBD sharing along the genome
 plot(phasingResult$surrDensity, xlab = 'SNP index', ylab="number of pairs that share IBD" ); # plot the density of IBD sharing along the genome
##################################


# select individuals for resequencing
##################################
pickedInds <- pickReseqIndsFast(surrMatDf = phasingResult$surrMatDf, # who shares IBD with whom and where in the genome - IBD matrix
                                   plotPath='reseqencing.pdf', # name of the file file with the plot summarising 
                                   data=genotypeData, # genotypes for individuals under study
                                   appFactor=10 # approximation factor, "resolution", 1 is full resolution, 100 would be poor
                                   ) 
##################################









