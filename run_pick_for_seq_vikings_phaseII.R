library(anchap)

# load data
##################################
genotypesFile <- "data/vikings_anchap_filtered"; # bed, fam, bim files (plink)

# file with the genetic distances for Orkney
geneticSnpMapFile <- NULL; # do not specify the genetic map
#geneticSnpMapFile <- list('data/LewisMap36/genetic_map_chr', '_b36.txt'); # genetic map files (optional)

cat('Loading data ... \n')
orkneyData <- loadData2(genotypesFile = genotypesFile, geneticSnpMapFile=geneticSnpMapFile)

genotypeData <- orkneyData
rm(orkneyData)
cat('All data read. \n')
##################################


# extract IBD sharing between individuals
##################################
load(file="phasingResult_3000.Rdata")



  # plot the density of IBD sharing along the genome
##################################


# select individuals for resequencing
##################################
pickedInds <- pickReseqIndsFast(surrMatDf = phasingResult$surrMatDf, # who shares IBD with whom and where in the genome - IBD matrix
                                   plotPath='reseqencing.pdf', # name of the file file with the plot summarising 
                                   data=genotypeData, # genotypes for individuals under study
                                   appFactor=10,
 # approximation factor, "resolution", 1 is full resolution, 100 would be poor
     
 noInds=2111)                                                                      
) 
save(pickedInds,file="pickedInds.Rdata")
##################################

