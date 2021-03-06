extractAncestralHaplotypes <- function(
                                       plotFileName,
                                       loci,
                                       probands,
                                       iHaploSeqsMyH1,
                                       iHaploSeqsMyH2) {


  

anchapResult <- list()

pickContaining <- function (seqs, locus) {
   # Selects a subset of sequences that contain a given locus.
   # Args:
   #   seqs : Sequences to choose from.
   #   locus : The locus of interest.
   chosen <- ((seqs$start<=locus) &  (seqs$end>locus))

   if (((dim(seqs)[1])>0)&&(length(chosen)>0)){
     seqs <- seqs[chosen,]
   }
   else{
     seqs <- NULL
   }
   
   seqs
 }

# prepare the plot
pdf(plotFileName)
plot(1:(1000), (1:1000)/1000, col="green", xlab='cluster', ylab='proportion of haplotypes known' )

# Loop over all loci in the analysis
for (locus in loci) {

  # Prepare data structures
  iHaploSeqsH1At <- list()
  iHaploSeqsH2At <- list()
  
###### find surrogate parents of all individuals at the locus
  for (i in probands) {
    seqs <- iHaploSeqsMyH1[[i]]
    seqs <- pickContaining(seqs, locus)
    iHaploSeqsH1At[[i]] <- seqs

    
    
    seqs <-  iHaploSeqsMyH2[[i]]
    seqs <- pickContaining(seqs, locus)
    iHaploSeqsH2At[[i]] <- seqs
  }

####### end  find surrogate parents of all individuals at the locus

  
  # store ids of clusters, for every haplotype in the cohort
  # for example, cluster id for individual 3, haplotype 2 is at clusterAssignmentH2[3]
  clusterAssignmentH1 <- integer(length(probands))
  clusterAssignmentH2 <- integer(length(probands))

  # keep track of whether a haplotype identity with other haplotypes had been looked at,
  # indexed as previously
  usedFlagH1 <- integer(length(probands))
  usedFlagH2 <- integer(length(probands))

  # keep track of how many surrogates sequences there are to any haplotype,
  # indexed as previously  
  noSurrogatesH1 <- integer(length(probands))
  noSurrogatesH2 <- integer(length(probands))
  
  noInconsH1 <- integer(length(probands))
  noInconsH2 <- integer(length(probands))


                                        # total number of haplotypes
  noHaplos <- 2* length(probands);

  # id of a currently created cluster
  clusterId <-1;

  # while not all haplotypes have been assigned to clusters
  while ((sum(usedFlagH1)+sum(usedFlagH2))<noHaplos)
    {
      # if possible, add to current clusters
      # check if all haplotypes in the cluseters have been used
      unusedInClustersH1 <- which((clusterAssignmentH1>0) & (!usedFlagH1))
      unusedInClustersH2 <- which((clusterAssignmentH2>0) & (!usedFlagH2))
      
      if (length(c(unusedInClustersH1,unusedInClustersH2))>0) {
      # there are some unexplored haplotypes already assigned to clusters
        
        # mark the origin of the relation as used
        # if the target of the relation is already in another cluster, merge clusters
        
        if (length(unusedInClustersH1)>0) {
        ######
          originH1 <-  unusedInClustersH1[1]
          currentCluster <- clusterAssignmentH1[originH1]
          usedFlagH1[originH1] <- 1
          
          # find the targets
          targets <- iHaploSeqsH1At[[originH1]]$id
          noSurrogatesH1[originH1] <- length(targets)
       
          if (length(targets)>0) {          
            for (t in 1:length(targets)){
              target <- targets[t]
              if (originH1 %in% iHaploSeqsH1At[[target]]$id) {
                # this is target's H1
                oldCluster <-clusterAssignmentH1[target]
                clusterAssignmentH1[target] <-  currentCluster                     
            }
              else if (originH1 %in% iHaploSeqsH2At[[target]]$id) {
                # this is target's H2
                oldCluster <- clusterAssignmentH2[target]
                clusterAssignmentH2[target] <-  currentCluster
                
              }
              else {
                noInconsH1[originH1] <- noInconsH1[originH1] + 1;
              }
            }
          ######
          }
        
        }
        else if (length(unusedInClustersH2)>0) {
        ######
          originH2 <-  unusedInClustersH2[1]
          currentCluster <- clusterAssignmentH2[originH2]
          usedFlagH2[originH2] <- 1
          
          # find the targets
          targets <- iHaploSeqsH2At[[originH2]]$id
          noSurrogatesH2[originH2] <- length(targets)
          if (length(targets)>0) {
            for (t in 1:length(targets)){
              target <- targets[t]
              if (originH2 %in% iHaploSeqsH1At[[target]]$id) {
                                        # this is target's H1
                oldCluster <-clusterAssignmentH1[target]
                clusterAssignmentH1[target] <-  currentCluster
                
              }
              else if (originH2 %in% iHaploSeqsH2At[[target]]$id) {
                                        # this is target's H2
                oldCluster <-clusterAssignmentH2[target]
                clusterAssignmentH2[target] <-  currentCluster
              }
              else {
                noInconsH2[originH2] <- noInconsH2[originH2] + 1;              
              }
            }
          }
          ######
       
        }
      }
      else {
         # start a new cluster

        unusedH1 <- which(usedFlagH1==0) 
        if (length(unusedH1)>0) {
          # using first available unused haplotype 1
          originH1 <- unusedH1[1]
          clusterAssignmentH1[originH1] <- clusterId;
        }
        else  {
          unusedH2 <- which(usedFlagH2==0)
           # using first available unused haplotype 2
          originH2 <- unusedH2[1]
          clusterAssignmentH2[originH2] <- clusterId;
        }
        clusterId <- clusterId +1;
        
      }
    
    }

  ###### post-analysis statistics
  clusterSizes <- table(c(clusterAssignmentH1,clusterAssignmentH2))
  clusterSizes <- clusterSizes[order(-clusterSizes)]
  noSurrogates <- c(noSurrogatesH1,noSurrogatesH2)

  clustering <- cumsum(clusterSizes)/(2*length(probands))
  clusterIds <- c(1,5,10,20,30,40,50,75,100, 200, 400, 1000)
  clusterIds <- clusterIds[clusterIds<=(2*length(probands))]
  clusterIds <- clusterIds[clusterIds<=(length(clustering))]
  
  
  aResult <- data.frame(clusterIndx = clusterIds, haploRatio =clustering[clusterIds])

 lines( 1:(length(clustering)), clustering, ps=2)
  ###### end post-analysis statistics



anchapResult[[locus]] <- list();
anchapResult[[locus]]$clusterAssignmentH1 <- clusterAssignmentH1 ;
anchapResult[[locus]]$clusterAssignmentH2  <- clusterAssignmentH2;
anchapResult[[locus]]$clusterSizes <- clusterSizes;
anchapResult[[locus]]$aResult <- aResult
}
# finish the plot
dev.off();



anchapResult
}
