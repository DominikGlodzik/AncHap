comparePairHaps <- function(haps1, haps2) {
  
      noBlocks <- dim(haps1[[1]])[1] 
  
      # find opposing homozygotes between combinations of haplotypes
      oppos11 <-  (abs(haps1[[1]]-haps2[[1]])==2)
      oppos12 <-  (abs(haps1[[1]]-haps2[[2]])==2)
      oppos21 <-   (abs(haps1[[2]]-haps2[[1]])==2)
      oppos22 <-  (abs(haps1[[2]]-haps2[[2]])==2)

      # ratio of mismatches between haplotypes, wrt loci where they could mismatch
      opp <- list()
      opp[[1]] <- rowSums(oppos11); 
      opp[[2]] <- rowSums(oppos12); 
      opp[[3]] <- rowSums(oppos21); 
      opp[[4]] <- rowSums(oppos22); 
      
      # for each combination, how many loci there were doubly phased
      signif <- list()
      signif[[1]] <- rowSums(haps1[[1]]>0 & haps2[[1]]>0)
      signif[[2]]<- rowSums(haps1[[1]]>0 & haps2[[2]]>0)
      signif[[3]] <- rowSums(haps1[[2]]>0 & haps2[[1]]>0)
      signif[[4]] <- rowSums(haps1[[2]]>0 & haps2[[2]]>0)


      # checking how many heterozygous loci there are, 
      hetero1 <-  rowSums(abs(haps1[[1]]-haps1[[2]])==2)
      hetero2 <-  rowSums(abs(haps2[[1]]-haps2[[2]])==2)     
            
      # loop over the blocks
      blockOpposMin <- vector(); signifLoci <- vector()
      for (b in 1:noBlocks) {
        # check which of the 4 combinations is matching best
        incVectr <- c(opp[[1]][b]/signif[[1]][b] , opp[[2]][b]/signif[[2]][b], opp[[3]][b]/signif[[3]][b], opp[[4]][b]/signif[[4]][b])

        minIndx <- which(incVectr==min(incVectr))[1]
        if (length(opp[[minIndx]][b])==0) {blockOpposMin[b]<-0; signifLoci[b]<-0
                                         } else {
                                           blockOpposMin[b] <- opp[[minIndx]][b]
                                           signifLoci[b] <- signif[[minIndx]][b]
                                         }
      }

      result <- list();
      result$inconsCount <- sum(blockOpposMin)
      result$signifLociCount<- sum(signifLoci)
      result
}
