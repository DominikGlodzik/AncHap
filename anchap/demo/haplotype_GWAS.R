data("phasingResult")

#####
# outputs

outputPaths <- list()
outputPaths$anchapPlot <- "../Data/results/orca749chr2/anchap.pdf"

#####

# do the anlysis at one locus only
loci <- c(10050)

individualsPhased <- 1:length(phasingResult$stitchingResults)

iHaploSeqsMyH1 <- list()
iHaploSeqsMyH2 <- list()
for (i in individualsPhased) {
  iHaploSeqsMyH1[[i]] <- phasingResult$stitchingResults[[i]]$selectedSequencesH1
  iHaploSeqsMyH2[[i]] <- phasingResult$stitchingResults[[i]]$selectedSequencesH2
}

anchapResult <- extractAncestralHaplotypes(
                                           plotFileName = outputPaths$anchapPlot,
                                           loci = loci,
                                           probands = individualsPhased,
                                           iHaploSeqsMyH1 = iHaploSeqsMyH1,
                                           iHaploSeqsMyH2 = iHaploSeqsMyH2
                                           )

# randomly generate the phenotypes
phenotypes <- rnorm(length(individualsPhased), mean = 100, sd=10);

for (locus in loci) {

  cat(paste("Sharing at locus",locus,"\n"))
  
  #print("Proportion of all haplotypes in the cluster")
  #print(anchapResult[[locus]]$aResult)
  print("Sizes of the largest clusters (bottom line)")
  print(anchapResult[[locus]]$clusterSizes[1:5])
  print("Mean cluster size")
  print(mean(anchapResult[[locus]]$clusterSizes))
  print("Distribution of cluster sizes")
  print(quantile(anchapResult[[locus]]$clusterSizes))

  print("Sharing at the locus")
  print("Mean number of surrogates")
  print(mean(anchapResult[[locus]]$clusterSizes))
  print("Distribution of number of surrogates")
  print(quantile(anchapResult[[locus]]$clusterSizes))
  print("--------")


  
  pValue <- doHapAssocTest(
                          clusterAssignmentH1 = anchapResult[[locus]]$clusterAssignmentH1,
                          clusterAssignmentH2 = anchapResult[[locus]]$clusterAssignmentH2,
                          clusterSizes = anchapResult[[locus]]$clusterSizes,
                          phenotypes,
                          MIN.CLUSTER.SIZE = 4)
  cat(paste("pValue ", pValue, "\n" ))
}
